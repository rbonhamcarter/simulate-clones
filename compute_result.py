import argparse
import numpy as np
import anndata
import pandas as pd
import os
import wot
from numpy.random import default_rng
import ot

rng = default_rng()

import utils as src

parser = argparse.ArgumentParser(description='This script takes an anndata and a number of tmaps, '
                                 + 'which are inferred for a dataset for a specific sampling rates ' 
                                 + 'and computes metrics for understanding the sampled dataset and ' 
                                 + 'comparing the metrics.')

parser.add_argument('-gp', '--growth_progenitor', metavar='GP', type=float, help='Max growth for progenitor '
                    + 'cell type (num of descendants between t1 and t2). Used to find anndata filename.')
parser.add_argument('-a', '--anndata', metavar='A', type=str, help='The anndata filename. '
                    + 'Should follow the naming convention: full_adata_growth_max_0-sample_rate.h5ad')
parser.add_argument('-r', '--sample_rate', metavar='R', type=float, help='The rate at which '
                    + 'the dataset will be sampled.')
parser.add_argument('-s', '--simulation', metavar='S', type=int, help='Simulation number.')
parser.add_argument('-t', '--tmap_dir', metavar='T', type=str, help='Filename of the directory '
                    + 'containing all the tmaps. The expected structure is as follows: '
                    + '/simulation_number/0-sample_rate/method/tmaps_1_2.h5ad')
parser.add_argument('-o', '--output_dir', metavar='O', type=str, help='Basename for the output '
                    + 'directory for the results.')


# Specify the barcoding rate and timepoints used as global variables
t1 = 1.0
t2 = 2.0

time_points = [t1, t2]
cell_types = {str(t1): ['progenitor', 'type A'], str(t2): ['type A', 'type B']}
all_cell_types = ['progenitor', 'type A', 'type B']

# Assign colors by cell type for plotting
celltype_color = {'progenitor': 'tab:green',
                  'type A': 'tab:red',
                  'type B': 'tab:blue'}


def main():
    # Parse parameters
    args = parser.parse_args()
    growth_max = str(args.growth_progenitor).replace('.', '-')
    output_dir= args.output_dir
    adata_dir = args.anndata
    r = args.sample_rate
    s = args.simulation
    tmap_dir = args.tmap_dir

    rate_key = str(r).replace('.', '-')
    print(s)
    print(rate_key)

    assert os.path.isdir(tmap_dir)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # Get the method info
    method_keys = os.listdir(os.path.join(tmap_dir, str(s), rate_key))
    method_keys.sort()
    try:
        method_keys.remove('true')
    except ValueError:
        pass

    # Load the full adata
    filename = 'full_adata_{}_{}.h5ad'.format(growth_max, rate_key)
    adata = anndata.read_h5ad(os.path.join(adata_dir, filename))

    # Load the tmaps
    tmaps = {}
    tmap_dirs = {}
    exact_rates = {}

    # Load the true tmaps (always computed on full_adata, invariant to simulation number - always use 0)
    true_dir = os.path.join(tmap_dir, '0', rate_key, 'true')
    tmap_dirs['true'] = true_dir
    tmaps['true'] = anndata.read_h5ad(os.path.join(true_dir, "tmaps_{}_{}.h5ad".format(int(t1), int(t2))))

    # Load the estimated tmaps
    for method in method_keys:
        method_dir = os.path.join(tmap_dir, str(s), rate_key, method)
        tmap_dirs[method]= method_dir
        tmaps[method] = anndata.read_h5ad(os.path.join(method_dir, "tmaps_{}_{}.h5ad".format(int(t1), int(t2))))
    
    # Extract cell number, rate information
    mask_t1 = (adata.obs['Cell Types'] == 'progenitor') & (adata.obs['Order'] == t1)
    mask_t2 = (adata.obs['Cell Types'] == 'progenitor') & (adata.obs['Order'] == t2)

    statistics_df = pd.DataFrame(
        [[s, r, str(t1), tmaps['true'].n_obs, adata[mask_t1].n_obs], 
         [s, r, str(t2), tmaps['true'].n_vars, adata[mask_t2].n_obs]], 
         columns=['Simulation', 'Sample rate', 'Order', 'N cells', 'N progenitor'])

    exact_rates[str(t1)] = tmaps['cospar-mt'].n_obs/tmaps['true'].n_obs
    exact_rates[str(t2)] = tmaps['cospar-mt'].n_vars/tmaps['true'].n_vars

    # Write cell number and t2 rates info to disk
    fname = os.path.join(output_dir, "statistics_{}_{}.csv".format(str(s), rate_key))
    if os.path.isfile(fname):
        prev_file = pd.read_csv(fname)
        new_file = prev_file.append(statistics_df)
        new_file.to_csv(fname, index=False)

    else:
        statistics_df.to_csv(fname, index=False)

    # Extract growth rate statistics
    growth_mean = {}
    growth_std = {}

    growth_mean['progenitor'] = adata[mask_t1].obs['n_children'].mean()
    growth_std['progenitor'] = adata[mask_t1].obs['n_children'].std()

    mask_t1 = (adata.obs['Cell Types'] == 'type A') & (adata.obs['Order'] == t1)
    growth_mean['type A'] = adata[mask_t1].obs['n_children'].mean()
    growth_std['type A'] = adata[mask_t1].obs['n_children'].std()

    growth_df = pd.DataFrame(
        [[s, r, str(t1), 'progenitor', growth_mean['progenitor'], growth_std['progenitor']], 
         [s, r, str(t1), 'type A', growth_mean['type A'], growth_std['type A']]], 
         columns=['Simulation', 'Sample rate', 'Order', 'Cell type', 'Growth mean', 'Growth std'])


    # Write growth rate statistics to disk
    fname = os.path.join(output_dir, "growth_{}_{}.csv".format(str(s), rate_key))
    if os.path.isfile(fname):
        prev_file = pd.read_csv(fname)
        new_file = prev_file.append(growth_df)
        new_file.to_csv(fname, index=False)

    else:
        growth_df.to_csv(fname, index=False)

    # Extract the 1000x1000 sampled adatas for each simulation, define cell sets
    sampled_adata = {}
    cell_sets = {}

    method='true'
    cell_sets[method] = {}
    for t in all_cell_types:
        cell_sets[method][t] = list(
            adata.obs.index[adata.obs['Cell Types'] == t])

    b_new_dict = {}
    for method in method_keys:
        cells_in_sample = list(tmaps[method].obs_names) + list(tmaps[method].var_names)

        # Extract the adata from the sample
        sampled_adata[method] = adata[cells_in_sample]

        # Remove the clone information from the appropriate cells
        cells_to_remove = list(tmaps[method][tmaps[method].obs['retained lineage?'] == False].obs_names) +\
             list(tmaps[method][:, tmaps[method].var['retained lineage?'] == False].var_names)
        sampled_adata[method], b_new = src.remove_clone_info_from_cells(sampled_adata[method], cells_to_remove)
        b_new_dict[method] = b_new    

        # Define the cell sets in the sample
        cell_sets[method] = {}
        for t in all_cell_types:
            cell_sets[method][t] = list(
                sampled_adata[method].obs.index[sampled_adata[method].obs['Cell Types'] == t])

    # Create a cell census dataframe for cells in MT for the two timepoints
    n_clones = adata.obsm['X_clone'].shape[1]
    clone_idx = range(n_clones)
    method = 'cospar-mt'    # choose any method that uses all sampled cells
    MT_data = src.retain_multitime_only(sampled_adata[method], time_points, clone_idx)

    print(sampled_adata[method].n_obs)
    print(MT_data.n_obs)

    # Compute the predicted proportion information
    print(b_new_dict)
    combined_rate_t1 = exact_rates[str(t1)] * b_new_dict[method][str(t1)]
    combined_rate_t2 = exact_rates[str(t2)] * b_new_dict[method][str(t2)]
    print("combined rate t1", combined_rate_t1)
    print("combined_rate_t2", combined_rate_t2)

    p_type_given_MT_dict = src.probability_cell_of_type_given_multitime(
                adata, time_points, 
                combined_rate_t1, combined_rate_t2, 
                growth_key='n_children')

    MT_census_df = pd.DataFrame({'progenitor': [],'type A': [], 'type B': []})

    for t in time_points:
        if t == t1:
            data = MT_data[MT_data.obs.Order == t1].obs
            if data.shape[0] == 0:
                MT_census_df = MT_census_df.append({
                                        'Simulation': s, 
                                        'Sample rate': r,
                                        'Order': t, 
                                        'progenitor': 0, 
                                        'type A': 0, 
                                        'type B': 0,
                                        'type A predicted': p_type_given_MT_dict['type A'],
                                        'progenitor predicted': p_type_given_MT_dict['progenitor']
                                        },   
                                        ignore_index = True)
            else:
                MT_census_df = MT_census_df.append({
                                            'Simulation': s, 
                                            'Sample rate': r,
                                            'Order': t, 
                                            'progenitor': data[data['Cell Types'] == 'progenitor'].shape[0]/data.shape[0], 
                                            'type A': data[data['Cell Types'] == 'type A'].shape[0]/data.shape[0], 
                                            'type B': data[data['Cell Types'] == 'type B'].shape[0]/data.shape[0],
                                            'type A predicted': p_type_given_MT_dict['type A'],
                                            'progenitor predicted': p_type_given_MT_dict['progenitor']
                                            },   
                                            ignore_index = True)

        else:
            data = MT_data[MT_data.obs.Order == t2].obs

            if data.shape[0] == 0:
                MT_census_df = MT_census_df.append({
                                        'Simulation': s, 
                                        'Sample rate': r,
                                        'Order': t, 
                                        'progenitor': 0, 
                                        'type A': 0, 
                                        'type B': 0},   
                                        ignore_index = True)
            else:
                MT_census_df = MT_census_df.append({
                                            'Simulation': s, 
                                            'Sample rate': r,
                                            'Order': t, 
                                            'progenitor': data[data['Cell Types'] == 'progenitor'].shape[0]/data.shape[0], 
                                            'type A': data[data['Cell Types'] == 'type A'].shape[0]/data.shape[0], 
                                            'type B': data[data['Cell Types'] == 'type B'].shape[0]/data.shape[0]},   
                                            ignore_index = True)

        

    # Write the MT census and predicted proportions to disk
    fname = os.path.join(output_dir, "MT_census_{}_{}.csv".format(str(s), rate_key))
    if os.path.isfile(fname):
        prev_file = pd.read_csv(fname)
        new_file = prev_file.append(MT_census_df)
        new_file.to_csv(fname, index=False)

    else:
        MT_census_df.to_csv(fname, index=False)

    
    # Check the proportions of ancestors of the type A cells in the randomly sampled datasets
    tolerance = 0.04
    counts = {}

    counts['progenitor'] = []
    counts['type A'] = []

    data_t2 = sampled_adata[method_keys[0]][sampled_adata[method_keys[0]].obs.Order == t2]
    
    for c in ['type A']:
        data_t2_parents_of_c = data_t2[data_t2.obs['Cell Types'] == c].obs['parent']
        temp_counts = {'progenitor': 0, 'type A': 0}
        
        for p in data_t2_parents_of_c:
            p_type = adata.obs['Cell Types'][p]
            temp_counts[p_type] += 1
            
        counts['progenitor'].append(temp_counts['progenitor']/sum(temp_counts.values()))
        counts['type A'].append(temp_counts['type A']/sum(temp_counts.values()))
        
    for c in ['progenitor', 'type A']:
        assert np.allclose(np.mean(counts[c]), np.mean(counts[c]), 
                           atol=tolerance, rtol=tolerance)

    # Set up tmap models, cell sets, and compute fates and trajectories for each method
    tmap_model = {}
    trajectory_ds = {}
    fate_ds = {}
    source_populations = {}
    target_populations = {}
    transition_table = {}

    # Initialize the tmap models
    method = 'true'
    tmap_model[method] = wot.tmap.TransportMapModel.from_directory(tmap_dirs[method] + '/')

    # Define source and target populations
    source_populations[method] = tmap_model[method].population_from_cell_sets(
        cell_sets[method], at_time=t1)
    target_populations[method] = tmap_model[method].population_from_cell_sets(
        cell_sets[method], at_time=t2)

    # Calculate trajectories from tmaps
    trajectory_ds[method] = tmap_model[method].trajectories(target_populations[method])
    trajectory_ds[method].obs = trajectory_ds[method].obs.join(adata.obs)

    # Calculate fates from tmaps
    fate_ds[method] = tmap_model[method].fates(target_populations[method])
    fate_ds[method].obs = fate_ds[method].obs.join(adata.obs) 

    # Calculate true transition table
    transition_table[method] = tmap_model[method].transition_table(
                    source_populations[method], target_populations[method])

    for method in method_keys:
        # Initialize the tmap models
        tmap_model[method] = wot.tmap.TransportMapModel.from_directory(tmap_dirs[method] + '/')

        # Define source and target populations
        source_populations[method] = tmap_model[method].population_from_cell_sets(
            cell_sets[method], at_time=t1)
        target_populations[method] = tmap_model[method].population_from_cell_sets(
            cell_sets[method], at_time=t2)

        # Calculate trajectories from tmaps
        trajectory_ds[method] = tmap_model[method].trajectories(target_populations[method])
        trajectory_ds[method].obs = trajectory_ds[method].obs.join(sampled_adata[method].obs)

        # Calculate fates from tmaps
        fate_ds[method] = tmap_model[method].fates(target_populations[method])
        fate_ds[method].obs = fate_ds[method].obs.join(sampled_adata[method].obs) 

    # Clean nans (replace with 0.0) from true model objects
    method = 'true'
    trajectory_ds[method].X = np.nan_to_num(trajectory_ds[method].X)
    fate_ds[method].X = np.nan_to_num(fate_ds[method].X)

    # Clean nans (replace with 0.0) from cospar-mt objects
    method = 'cospar-mt'
    trajectory_ds[method].X = np.nan_to_num(trajectory_ds[method].X)
    fate_ds[method].X = np.nan_to_num(fate_ds[method].X)

    # Compute results
    fate_corr = {}
    fate_wdist = {}
    ancestor_corr = {}
    ancestor_wdist = {}

    for method in method_keys:
        # For each method, generate a transition table
        print("Now computing results for the {} trajectories...".format(method))
        transition_table[method] = tmap_model[method].transition_table(
                    source_populations[method], target_populations[method])

        # For each method, compute the correlation between the true and estimated fate probabilities
        fate_corr[method] = np.zeros(len(cell_types[str(t2)]))

        true_fates = fate_ds['true']
        true_fates = true_fates[true_fates.obs['Order'] < t2].copy()

        fates = fate_ds[method]
        fates = fates[fates.obs['Order'] < t2].copy()
        
        # Reduce true fates vector to just the subset sampled
        true_fates_subset = true_fates[fates.obs.index]
        
        # Assert order of cell types is correct
        assert np.all(fates.var_names == true_fates_subset.var_names)
        
        # Remove cells with no descendants
        true_fates_subset = true_fates_subset[true_fates_subset.obs.n_children > 0]
        fates = fates[fates.obs.n_children > 0]
        
        for j in range(len(cell_types[str(t2)])):
            if np.std(fates.X.T[j]) == 0.0:
                    fate_corr[method][j] = 0.0 # corr with a constant is always zero
            else:
                fate_corr[method][j] = np.corrcoef(true_fates_subset.X.T[j], fates.X.T[j])[0, 1]

        # Compute the Wasserstein-2 distance between the true and estimate fate probs
        fate_wdist[method] = np.zeros(len(cell_types[str(t2)]))

        # Compute the cost matrix
        a_x = np.float64(adata[true_fates.obs.index].X)
        b_x = np.float64(adata[fates.obs.index].X)
        M = ot.dist(a_x, b_x)

        # Compute the Wasserstein distance for each possible fate
        for j in range(len(cell_types[str(t2)])):

            # Normalize fate vectors
            a = np.float64(true_fates.X.T[j])
            a /= a.sum()
            assert np.isclose(a.sum(), 1.0)

            b = np.float64(fates.X.T[j])
            if b.sum() == 0:
                break
            b /= b.sum()
            assert np.isclose(b.sum(), 1.0)

            # Ensure all inputs to emd2 are the same type
            a = np.array(a)
            b = np.array(b)
            M = np.array(M)

            # Compute the Wasserstein distance for fate j
            fate_wdist[method][j] = ot.emd2(a, b, M, numItermax=100000*2)

        # For each method, compute the correlation between the true and estimated ancestor distributions
        ancestor_corr[method] = np.zeros(len(cell_types[str(t2)]))

        true_ancestors = trajectory_ds['true'][trajectory_ds['true'].obs['Order'] < t2].copy()

        ancestors = trajectory_ds[method][trajectory_ds[method].obs['Order'] < t2].copy()
        
        # Reduce true fates vector to just the subset sampled
        true_ancestors_subset = true_ancestors[ancestors.obs.index]
        
        # Assert order of cell types is correct
        assert np.all(true_ancestors_subset.var_names == ancestors.var_names)
        
        for j in range(len(cell_types[str(t2)])):
            if np.std(ancestors.X.T[j]) == 0.0:
                    ancestor_corr[method][j] = 0.0 # corr with a constant is always zero
            else:
                ancestor_corr[method][j] = np.corrcoef(true_ancestors_subset.X.T[j], ancestors.X.T[j])[0, 1]

        # Compute the Wasserstein-2 distance between the true and estimated ancestor distributions
        ancestor_wdist[method] = np.zeros(len(cell_types[str(t2)]))

        # Assert order of cell types is correct
        assert np.all(ancestors.var_names == true_ancestors.var_names)

        # Compute the cost matrix
        a_x = np.float64(adata[true_ancestors.obs.index].X)
        b_x = np.float64(adata[ancestors.obs.index].X)
        M = ot.dist(a_x, b_x)

        # Compute the Wasserstein distance for each possible fate
        for j in range(len(cell_types[str(t2)])):

            # Normalize ancestor vectors
            a = np.float64(true_ancestors.X.T[j])
            a /= a.sum()
            assert np.isclose(a.sum(), 1.0)

            b = np.float64(ancestors.X.T[j])
            if b.sum() == 0:
                break
            
            b /= b.sum()
            assert np.isclose(b.sum(), 1.0)

            # Ensure all inputs to emd2 are the same type
            a = np.array(a)
            b = np.array(b)
            M = np.array(M)

            # Compute the Wasserstein distance for fate j
            ancestor_wdist[method][j] = ot.emd2(a, b, M, numItermax=100000*2)

    # Store results in dataframes
    tt_data = []
    fate_data = []
    ancestor_data = []


    # True Transition table
    method='true'
    for c_t1 in transition_table[method].obs.index:
        for c_t2 in transition_table[method].var.index:
            tt_data.append([method, c_t1, c_t2, transition_table[method][c_t1, c_t2].X[0][0]])

    for method in method_keys:
        # Transition table
        for c_t1 in transition_table[method].obs.index:
            for c_t2 in transition_table[method].var.index:
                tt_data.append([method, c_t1, c_t2, transition_table[method][c_t1, c_t2].X[0][0]])

        for i in range(len(cell_types[str(t2)])):
            # Fate probabilities
            c_t2 = true_fates_subset.var.index[i]
            fate_data.append([method, c_t2, fate_corr[method][i], fate_wdist[method][i]])

            # Ancestor distributions
            c_t2 = true_ancestors_subset.var.index[i]
            ancestor_data.append([method, c_t2, ancestor_corr[method][i], ancestor_wdist[method][i]])
        
    transition_table_df = pd.DataFrame(tt_data,
        columns=['Method', 'Cell type t1', 'Cell type t2', 'Transition mass'])
    fate_df = pd.DataFrame(fate_data,
        columns=['Method', 'Cell type t2', 'Fate probability correlation', 'Fate probability W-distance'])
    ancestor_df = pd.DataFrame(ancestor_data,
        columns=['Method', 'Cell type t2', 'Ancestor distribution correlation', 'Ancestor distribution W-distance'])

    # Add simulation number and sampling rate
    transition_table_df['Simulation'] = s
    transition_table_df['Sample rate'] = r
    fate_df['Simulation'] = s
    fate_df['Sample rate'] = r
    ancestor_df['Simulation'] = s
    ancestor_df['Sample rate'] = r

    # Save results to disk
    fnames_dict = {"transition_tables_{}_{}.csv".format(str(s), rate_key): transition_table_df,
                   "fate_probability_results_{}_{}.csv".format(str(s), rate_key): fate_df,
                   "ancestor_distribution_results_{}_{}.csv".format(str(s), rate_key): ancestor_df}
    
    for basename in fnames_dict.keys():
        fname = os.path.join(output_dir, basename)
        object_to_save = fnames_dict[basename]

        if os.path.isfile(fname):
            prev_file = pd.read_csv(fname)
            new_file = prev_file.append(object_to_save)
            new_file.to_csv(fname, index=False)

        else:
            object_to_save.to_csv(fname, index=False)

if __name__ == '__main__':
    main()