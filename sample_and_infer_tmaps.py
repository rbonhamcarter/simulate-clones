import argparse
import numpy as np
import lineageot
import scanpy as sc
import cospar as cs
import anndata
import os
import scipy
from numpy.random import default_rng
rng = default_rng()

import utils as src


parser = argparse.ArgumentParser(description='This script takes a dataset as input, '
                                 + 'randomly samples at each timepoint, strips clonal ' 
                                 + 'information from a fraction of the cells, computes ' 
                                 + 'estimates of the population trajectory using 4 '
                                 + 'different methods and saves the resulting tmaps to disk.')

parser.add_argument('-i', '--input', metavar='I', type=str, help='Filename of the input '
                    + 'dataset, stored as an h5ad file.')
parser.add_argument('-o', '--output', metavar='O', type=str, help='Basename for the output '
                    + 'directory for each tmap estimate.')
parser.add_argument('-r', '--sample_rate', metavar='R', type=float, help='The rate at which '
                    + 'the dataset will be sampled.')
parser.add_argument('-b', '--barcode_rate', metavar='B', type=float, help='The rate at which '
                    + 'clone info will be available in the dataset.')
parser.add_argument('-t', '--compute_true', action='store_true', help='Whether or not to '
                    + 'compute and save the true tmap (on the full adata object).', default=False)


# Define a FUZZ factor to help with numerical instabilities
FUZZ = 1e-8

def main():
    args = parser.parse_args()
    adata_path = args.input
    tmap_dir = args.output
    r = args.sample_rate
    b = args.barcode_rate

    # Load the dataset
    full_adata = anndata.read_h5ad(adata_path)

    # Extract cell type and time point information
    time_points = list(full_adata.obs['Order'].unique())
    t0 = min(time_points)
    time_points.remove(t0)
    assert len(time_points) == 2
    t1 = time_points[0]
    t2 = time_points[1]

    # Define adata restricted to sampled points
    adata = full_adata[full_adata.obs.Order > t0]

    # Compute the true coupling if required-------------------------------------------------------
    if args.compute_true:
        true_tmap = src.infer_true_backward_coupling(adata, t1, t2)

        # Save tmap as a h5ad file 
        temp_tmap_dir = os.path.join(tmap_dir, "true")
        if not os.path.isdir(temp_tmap_dir):
            os.makedirs(temp_tmap_dir)

        # Save tmap
        true_tmap.write(os.path.join(temp_tmap_dir, "tmaps_{}_{}.h5ad".format(int(t1), int(t2))))

    #----------------------------------------------------------------------------------
    # Sample the dataset
    # Sample number_to_sample cells at each timepoint irrespective of r
    sampled_adata = src.randomly_sample_timepoints(adata, time_points, n=1000)

    # Removing the clonal information from a fraction of the cells
    sampled_adata, b_new = src.remove_clone_info_from_fraction(sampled_adata, fraction_to_remove=1-b)
    print(b_new)

    # Compute embeddings
    sc.tl.pca(sampled_adata)
    sc.pp.neighbors(sampled_adata, n_neighbors=20)
    sc.tl.umap(sampled_adata, min_dist=0.3)
    sampled_adata.obsm['X_emb'] = sampled_adata.obsm['X_umap']

    # Infer the tmaps--------------------------------------------------------------------------
    # with cospar-mt --------------------------
    # Initialize the CoSpar adata using the convenience function from cospar
    temp_adata_for_cospar = anndata.AnnData(X=sampled_adata.X,
                                            obs=sampled_adata.obs.rename(columns={
                                                'Order': 'time_info', 
                                                'Cell Types': 'state_info'}), 
                                            var=sampled_adata.var)
    temp_adata_for_cospar.obsm['X_pca'] = sampled_adata.obsm['X_pca']
    temp_adata_for_cospar.obsm['X_emb'] = sampled_adata.obsm['X_umap']
    temp_adata_for_cospar.obsm['X_clone'] = scipy.sparse.csr_matrix(sampled_adata.obsm['X_clone'])
    temp_adata_for_cospar.obsm['X_clone'].eliminate_zeros()
    temp_adata_for_cospar = cs.pp.initialize_adata_object(adata=temp_adata_for_cospar)
    random_id = rng.integers(0, 9, size=6)
    temp_adata_for_cospar.uns['data_des'] = [''.join([str(i) for i in random_id])]
    
    # store relative growth rates for cells at t1 (only used in CoSpar-ST (below))
    total_children = temp_adata_for_cospar[temp_adata_for_cospar.obs.time_info == t1].obs['growth'].sum()
    temp_adata_for_cospar.obs['cell_growth_rate'] = sampled_adata.obs['growth']/total_children + FUZZ

    # Compute the transition map
    adata_w_tmap=cs.tmap.infer_Tmap_from_multitime_clones(
        temp_adata_for_cospar, extend_Tmap_space=True, compute_new=True)

    # Save tmap as a h5ad file 
    temp_tmap_dir = os.path.join(tmap_dir, "cospar-mt")
    if not os.path.isdir(temp_tmap_dir):
        os.makedirs(temp_tmap_dir)

    # Save tmap
    multitime_coupling = anndata.AnnData(X=adata_w_tmap.uns['transition_map'], 
                                            obs=adata_w_tmap.obs[adata_w_tmap.obs.time_info == t1],
                                         var=adata_w_tmap.obs[adata_w_tmap.obs.time_info == t2])
    multitime_coupling.write(os.path.join(temp_tmap_dir, "tmaps_{}_{}.h5ad".format(int(t1), int(t2))))

    # with cospar-st --------------------------
    adata_w_tmap = cs.tmap.infer_Tmap_from_one_time_clones(
        temp_adata_for_cospar, initialize_method='OT', OT_cost='GED',
        later_time_point=t2, initial_time_points=[t1], compute_new=True)

    # Save tmap as a h5ad file 
    temp_tmap_dir = os.path.join(tmap_dir, "cospar-st")
    if not os.path.isdir(temp_tmap_dir):
        os.makedirs(temp_tmap_dir)
    
    # Save tmap
    singletime_coupling = anndata.AnnData(X = adata_w_tmap.uns['transition_map'], 
                                        obs=adata_w_tmap.obs[adata_w_tmap.obs.time_info == t1],
                                        var = adata_w_tmap.obs[adata_w_tmap.obs.time_info == t2])
    singletime_coupling.write(os.path.join(temp_tmap_dir, "tmaps_{}_{}.h5ad".format(int(t1), int(t2))))

    # with lineageot-st------------------------------
    print("Now fitting lineageot-st coupling..")
    temp_adata_for_lineageot = anndata.AnnData(X=sampled_adata.X, obs=sampled_adata.obs.rename(columns={'Order': 'time'}), var=sampled_adata.var)
    temp_adata_for_lineageot.obsm['X_pca'] = sampled_adata.obsm['X_pca']
    temp_adata_for_lineageot.obsm['X_emb'] = sampled_adata.obsm['X_umap']
    temp_adata_for_lineageot.obsm['X_clone'] = sampled_adata.obsm['X_clone']
    
    # Remove the cells at t2 which have no clonal info
    time_mask = (temp_adata_for_lineageot.obs.time < t2)
    clone_mask = (np.sum(temp_adata_for_lineageot.obsm['X_clone'], axis=1) != 0)
    mask = time_mask | clone_mask
    temp_adata_for_lineageot = temp_adata_for_lineageot[mask]
    
    # compute relative growth rates to use for earlier time marginal in lineageot
    data_t1 = temp_adata_for_lineageot[temp_adata_for_lineageot.obs.time == t1]
    growth_1 = data_t1.obs['growth']
    marginal_1 = growth_1/growth_1.sum() + FUZZ
    marginal_1 = np.array(marginal_1)
    
    # Fit the tree using the clone information only at t2
    n_clones = sampled_adata.obsm['X_clone'].shape[1]
    clone_times = t0*np.zeros(n_clones)
    lineage_tree_t2 = lineageot.fit_tree(temp_adata_for_lineageot[temp_adata_for_lineageot.obs.time == t2], method="clones", clone_times=clone_times)
    
    # possibly will need to tune epsilon
    coupling = lineageot.fit_lineage_coupling(temp_adata_for_lineageot, t1, t2, lineage_tree_t2, 
                                              epsilon=0.05, marginal_1=marginal_1, marginal_2=[], time_key="time")
    
    # Saving the fitted coupling in the format Waddington-OT expects
    temp_tmap_dir = os.path.join(tmap_dir, "lineageot-st")
    if not os.path.isdir(temp_tmap_dir):
        os.makedirs(temp_tmap_dir)
        
    file_name = os.path.join(temp_tmap_dir, "tmaps_{}_{}.h5ad".format(int(t1), int(t2)))
    file_dir = os.path.dirname(file_name)

    if not os.path.exists(file_dir):
        os.mkdir(file_dir)
    coupling.write(file_name)

    # with lineageot-mt ---------------------------
    print("Now fitting lineageot-mt coupling..")
    lineage_tree = lineageot.fit_tree(temp_adata_for_lineageot, clone_times=clone_times, method="clones")

    coupling = lineageot.fit_lineage_coupling(temp_adata_for_lineageot, t1, t2, lineage_tree,
                                                epsilon=0.05,
                                                marginal_1=marginal_1, marginal_2=[], time_key="time")
    
    # Saving the fitted coupling in the format Waddington-OT expects to use downstream analysis tools
    temp_tmap_dir = os.path.join(tmap_dir, "lineageot-mt")
    if not os.path.isdir(temp_tmap_dir):
        os.makedirs(temp_tmap_dir)
        
    file_name = os.path.join(temp_tmap_dir, "tmaps_{}_{}.h5ad".format(int(t1), int(t2)))
    file_dir = os.path.dirname(file_name)

    if not os.path.exists(file_dir):
        os.mkdir(file_dir)
    coupling.write(file_name)


if __name__ == '__main__':
    main()