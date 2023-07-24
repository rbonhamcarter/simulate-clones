import numpy as np
from numpy.random import default_rng
import anndata
import pandas as pd
import wot
import warnings
rng = default_rng()


FUZZ = 1e-12


def normalize_to_backward_coupling(coupling):
    # Normalize to ensure a backward coupling
    temp = coupling.X
    try:
        temp = temp.toarray()
    except:
        pass

    # Column normalize
    if not np.all(temp.sum(axis=0) != 0.0):
        temp = np.divide(temp, temp.sum(axis=0) + FUZZ)
    else:
        temp = np.divide(temp, temp.sum(axis=0))
    
    # Total normalization
    temp /= np.sum(temp)

    return anndata.AnnData(X=temp, obs=coupling.obs, var=coupling.var)



def collapse_coupling_into_cell_sets(coupling, cell_sets, type_of_set="Cell set"):
    """Collapses the coupling into only transitions that are happening within and between cell sets."""
    n_sets = len(cell_sets)
    set_labels = list(cell_sets.keys())

    obs = pd.DataFrame({type_of_set: set_labels}, index=set_labels)
    var = pd.DataFrame({type_of_set: set_labels}, index=set_labels)
    clone_coupling = anndata.AnnData(X=np.zeros((n_sets, n_sets)), obs=obs, var=var)

    temp = np.zeros((n_sets, n_sets))
    for i in range(n_sets):
        for j in range(n_sets):
            c_source = set_labels[i]
            source_cell_set = cell_sets[c_source]

            c_target = set_labels[j]
            target_cell_set = cell_sets[c_target]

            obs_idx = coupling.obs.index.intersection(source_cell_set)
            var_idx = coupling.var.index.intersection(target_cell_set)

            coupling_c_subset = coupling[obs_idx, var_idx]
            temp[i, j] = np.sum(coupling_c_subset.X)

    clone_coupling.X = temp
    return clone_coupling


def collapse_coupling_into_clones(coupling, clone_key='Clone'):
    clone_labels = set(coupling.obs[clone_key]).intersection(set(coupling.var[clone_key]))
    clone_labels = list(clone_labels)

    obs = pd.DataFrame({'Clone': clone_labels}, index=clone_labels)
    clone_coupling = anndata.AnnData(X=np.zeros((len(clone_labels), 1)), obs=obs)

    temp = np.zeros((len(clone_labels), 1))
    for i in range(len(clone_labels)):
        c = clone_labels[i]

        obs_idx = coupling[coupling.obs[clone_key] == c].obs.index 
        var_idx = coupling[coupling.obs[clone_key] == c].var.index 

        coupling_c_subset = coupling[obs_idx, var_idx]
        temp[i] = np.sum(coupling_c_subset.X)

    clone_coupling.X = temp
    return clone_coupling


def total_variation(x, y, normalize=False):
    if normalize:
        x_normalized = x/x.sum()
        y_normalized = y/y.sum()

        return (1/2)*np.linalg.norm(x_normalized-y_normalized, ord=1)

    else:
        return (1/2)*np.linalg.norm(x-y, ord=1)


def get_multitime_clones(adata, t1, t2):
    clones_filtered_multitime = []
    cells_in_MT = {str(t1): [], str(t2): []}
    cell_clone_dict = {}

    for b in range(adata.obsm['X_clone'].shape[1]):
        cell_idx_in_b = adata.obsm['X_clone'][:,b].nonzero()[0]
        cells_in_b = adata[cell_idx_in_b].obs.index
        
        times = set(adata[cells_in_b].obs['Order'])
            
        if len(times) > 1:
            clones_filtered_multitime.append(b)
            
            for c in cells_in_b:
                t = adata[c].obs['Order'][0]
                cells_in_MT[str(t)].append(c)
                cell_clone_dict[c] = b

    return clones_filtered_multitime, cells_in_MT


def infer_cospar_intraclone_coupling(adata, t1, t2, clones_filtered_multitime=None):
    if not clones_filtered_multitime:
        n_clones = adata.obsm['X_clone'].shape[1]
        clones_filtered_multitime = list(range(n_clones))

    data_t1 = adata[adata.obs.Order == t1]
    data_t2 = adata[adata.obs.Order == t2]

    T = np.zeros((data_t1.n_obs, data_t2.n_obs))
    cell_ids_t1 = []
    cell_ids_t2 = []

    for k in range(len(clones_filtered_multitime)):
        c = clones_filtered_multitime[k]
        adata_w_c_indicies_T1 = data_t1.obsm['X_clone'][:, c].nonzero()[0]
        adata_w_c_indicies_T2 = data_t2.obsm['X_clone'][:, c].nonzero()[0]

        n_cells_c_T1 = adata_w_c_indicies_T1.shape[0]
        n_cells_c_T2 = adata_w_c_indicies_T2.shape[0]

        if (n_cells_c_T1 == 0) or (n_cells_c_T2 == 0):
            pass
        else:
            for i in adata_w_c_indicies_T1:
                T[i, adata_w_c_indicies_T2] = 1.0/(n_cells_c_T1 + n_cells_c_T2)

            cell_ids_t1.extend(data_t1.obs.index[adata_w_c_indicies_T1])
            cell_ids_t2.extend(data_t2.obs.index[adata_w_c_indicies_T2])

    return anndata.AnnData(X=T/T.sum(), obs=data_t1[cell_ids_t1].obs, var=data_t2[cell_ids_t2].obs)
        

def infer_intraclone_coupling(adata, t1, t2, clones_filtered_multitime=None):
    if not clones_filtered_multitime:
        n_clones = adata.obsm['X_clone'].shape[1]
        clones_filtered_multitime = list(range(n_clones))

    data_t1 = adata[adata.obs.Order == t1]
    data_t2 = adata[adata.obs.Order == t2]

    T = np.zeros((data_t1.n_obs, data_t2.n_obs))
    cell_ids_t1 = []
    cell_ids_t2 = []

    for k in range(len(clones_filtered_multitime)):
        c = clones_filtered_multitime[k]
        adata_w_c_indicies_T1 = data_t1.obsm['X_clone'][:, c].nonzero()[0]
        adata_w_c_indicies_T2 = data_t2.obsm['X_clone'][:, c].nonzero()[0]
        n_cells_c_T1 = adata_w_c_indicies_T1.shape[0]
        n_cells_c_T2 = adata_w_c_indicies_T2.shape[0]

        if (n_cells_c_T1 == 0) or (n_cells_c_T2 == 0):
            pass
        else:
            for i in adata_w_c_indicies_T1:
                T[i, adata_w_c_indicies_T2] = 1.0/n_cells_c_T1 

            cell_ids_t1.extend(data_t1.obs.index[adata_w_c_indicies_T1])
            cell_ids_t2.extend(data_t2.obs.index[adata_w_c_indicies_T2])

    return anndata.AnnData(X=T/T.sum(), obs=data_t1[cell_ids_t1].obs, var=data_t2[cell_ids_t2].obs)


def infer_true_backward_coupling(adata, t1, t2):
    """Use this for full populations, not sampled datasets."""

    data_t1 = adata[adata.obs.Order == t1]
    data_t2 = adata[adata.obs.Order == t2]

    coupling = anndata.AnnData(X=np.zeros((data_t1.n_obs, data_t2.n_obs)), 
                           obs=data_t1.obs, var=data_t2.obs)

    # For each cell at t2, use it's parent to fill in the coupling
    for child_cell in data_t2.obs_names:
        parent_cell = data_t2.obs['parent'][child_cell]
        clone = data_t2.obs['Clone'][child_cell]

        assert parent_cell in data_t1.obs_names, "Parent cell does not exist at t1!"
        assert data_t1.obs['Clone'][parent_cell] == clone, "Not the same clone!"
        
        coupling[parent_cell, child_cell] = 1.0
    
    return normalize_to_backward_coupling(coupling)


# Define a function for setting up tmap models, cell sets, and compute fates and trajectories
def set_up_model(sampled_adatas, cell_sets, tmap_dir, sample_keys, t1, t2, simulation_keys=None):
    tmap_model = {}
    trajectory_ds = {}
    fate_ds = {}
    source_populations = {}
    target_populations = {}

    if simulation_keys == None:
        for rate in sample_keys:
            # Initialize the tmap models
            tmap_model[rate] = wot.tmap.TransportMapModel.from_directory(tmap_dir[rate] + '/')

            # Define source and target populations
            source_populations[rate] = tmap_model[rate].population_from_cell_sets(
                cell_sets[rate], at_time=t1)
            target_populations[rate] = tmap_model[rate].population_from_cell_sets(
                cell_sets[rate], at_time=t2)

            # Calculate trajectories from tmaps
            trajectory_ds[rate] = tmap_model[rate].trajectories(target_populations[rate])
            trajectory_ds[rate].obs = trajectory_ds[rate].obs.join(sampled_adatas[rate].obs)

            # Calculate fates from tmaps
            fate_ds[rate] = tmap_model[rate].fates(target_populations[rate])
            fate_ds[rate].obs = fate_ds[rate].obs.join(sampled_adatas[rate].obs) 
            
        return tmap_model, trajectory_ds, fate_ds, source_populations, target_populations

    else:
        for sim in simulation_keys:
            tmap_model[sim] = {}
            trajectory_ds[sim] = {}
            fate_ds[sim] = {}
            source_populations[sim] = {}
            target_populations[sim] = {}

            for rate in sample_keys:
                # Initialize the tmap models
                tmap_model[sim][rate] = wot.tmap.TransportMapModel.from_directory(tmap_dir[sim][rate] + '/')

                # Define source and target populations
                source_populations[sim][rate] = tmap_model[sim][rate].population_from_cell_sets(
                    cell_sets[sim][rate], at_time=t1)
                target_populations[sim][rate] = tmap_model[sim][rate].population_from_cell_sets(
                    cell_sets[sim][rate], at_time=t2)

                # Calculate trajectories from tmaps
                trajectory_ds[sim][rate] = tmap_model[sim][rate].trajectories(target_populations[sim][rate])
                trajectory_ds[sim][rate].obs = trajectory_ds[sim][rate].obs.join(sampled_adatas[sim][rate].obs)

                # Calculate fates from tmaps
                fate_ds[sim][rate] = tmap_model[sim][rate].fates(target_populations[sim][rate])
                fate_ds[sim][rate].obs = fate_ds[sim][rate].obs.join(sampled_adatas[sim][rate].obs) 
                
        return tmap_model, trajectory_ds, fate_ds, source_populations, target_populations