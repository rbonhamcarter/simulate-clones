import numpy as np
import random
from numpy.random import default_rng
rng = default_rng()
from math import log


# Functions for the simulation process -------------------------------------------------------------------------------
def get_next_cell_type(current_cell_type, transition_matrix):
    transition_probabilities = transition_matrix[current_cell_type].X[0]
    next_states = transition_matrix.var_names
    
    u = rng.uniform()
    lower_bound = 0.0
    higher_bound = 0.0
    
    for i in range(transition_probabilities.shape[0]):
        higher_bound += transition_probabilities[i]
        if lower_bound < u <= higher_bound:
            return next_states[i]
        else:
            lower_bound += transition_probabilities[i]
            
    raise 'get_next_cell_type() failed'


def compute_gaussian_variance_for_concentration_value(
    concentration_value, concentration_prob=0.05, method="upper_dev"):

    if method == "upper_dev":
        return concentration_prob*concentration_value**2

    elif method == "chebyshev":
        numerator = -concentration_value**2
        denominator = 2*log(concentration_prob)
        return numerator/denominator

    else:
        raise("Method not implemented")



# Functions for the sampling proces -------------------------------------------------------------------------------
def remove_clone_info_from_fraction(data, fraction_to_remove, time_key="Order"):
    adata = data.copy()
    b = 1 - fraction_to_remove
    timepoints = list(adata.obs[time_key].unique())

    # Shuffle the adata
    shuffled_index = adata.obs.sample(frac=1).index
    adata = adata[shuffled_index]
    
    # Randomly sample the idx to strip
    idx_to_strip = random.sample(range(adata.n_obs), k=int(adata.n_obs*fraction_to_remove))

    # Strip the clone information
    new_clone_matrix = adata.obsm['X_clone'].copy()
    new_clone_matrix[idx_to_strip] = np.zeros_like(new_clone_matrix[idx_to_strip])
    adata.obsm['X_clone'] = new_clone_matrix

    # Check and record the exact proportion stripped
    b_new_dict = {}
    for t in timepoints:
        adata_t = adata[adata.obs[time_key] == t]

        n_cells_without_clone = np.count_nonzero(np.sum(adata_t.obsm['X_clone'], axis=1) == 0)
        b_new = 1 - n_cells_without_clone/adata_t.n_obs

        atol = 20/adata_t.n_obs
        assert np.allclose(b_new, b, rtol=0.05, atol=atol), "b_new={}, b={}".format(b_new, b)

        b_new_dict[str(t)] = b_new

    # Add obs for recording which cells still have clonal info
    adata.obs['retained lineage?'] = [True if i not in idx_to_strip else False for i in range(adata.n_obs)]

    # Reorder the adata
    adata = adata[data.obs.index]
    
    return adata, b_new_dict


def remove_clone_info_from_cells(data, cells, time_key='Order'):
    adata = data.copy()
    
    idx_to_strip = [adata.obs.index.get_loc(c) for c in cells]
    new_clone_matrix = adata.obsm['X_clone'].copy()
    new_clone_matrix[idx_to_strip] = np.zeros_like(new_clone_matrix[idx_to_strip])
    
    adata.obsm['X_clone'] = new_clone_matrix

    # Check and record the exact proportion stripped
    timepoints = list(adata.obs[time_key].unique())
    b_dict = {}
    for t in timepoints:
        adata_t = adata[adata.obs[time_key] == t]
        n_cells_without_clone = np.count_nonzero(np.sum(adata_t.obsm['X_clone'], axis=1) == 0)
        b_new = 1 - n_cells_without_clone/adata_t.n_obs

        b_dict[str(t)] = b_new

    # Add obs for recording which cells still have clonal info
    adata.obs['retained lineage?'] = [True if i not in idx_to_strip else False for i in range(adata.n_obs)]
    
    return adata, b_dict


def remove_descendants(parent_cell, adata, parent_key='parent'):
    '''Removes parent_cell and its descendants for adata. Returns new adata.'''

    # remove the parent
    adata = adata[adata.obs.index != parent_cell]

    # extract the children
    children = adata[adata.obs[parent_key] == parent_cell].obs.index

    if len(children) == 0:
        return adata

    else:
        # remove the descendants of those children
        for new_parent in children:
            adata = remove_descendants(new_parent, adata, parent_key)
    
        # remove children
        adata = adata[adata.obs[parent_key] != parent_cell]
        return adata


def randomly_sample_timepoints(adata, times, rate=None, n=None, time_key='Order', 
                               parent_key='parent', remove_sample=False):
    # if (rate==None & n==None) | (rate!=None & n!=None):
    #     raise ValueError('Must specify exactly one of rate and n.')

    times = np.sort(times)

    temp_adata = adata.copy()

    if rate:
        if isinstance(rate, list):
            assert len(rate) == len(times)
        else:
            rate = [rate]*len(times)
    else:
        if isinstance(n, list):
            assert len(n) == len(times)
        else:
            n = [n]*len(times)

    ids_to_sample = {}
    for j in range(len(times)):
        t = times[j]
        adata_t = temp_adata[temp_adata.obs[time_key] == t]
        n_cells_t = adata_t.n_obs
        
        # Shuffle the dataset
        shuffled_index = adata_t.obs.sample(frac=1).index
        adata_t = adata_t[shuffled_index]
        
        # Sample the subset
        if rate:
            n_cells_to_sample = int(n_cells_t*rate[j])
        else:
            n_cells_to_sample = n[j]

        assert n_cells_to_sample <= n_cells_t

        idx_sampled_at_t = random.sample(range(n_cells_t), k=n_cells_to_sample)
        ids_sampled_at_t = list(adata_t.obs.index[idx_sampled_at_t])

        ids_to_sample[str(t)] = ids_sampled_at_t

        if remove_sample and t != times[-1]:
            n_before = temp_adata.n_obs

            # Remove the descendants of the cells sampled at t
            for i in ids_sampled_at_t:
                temp_adata = remove_descendants(i, temp_adata, parent_key)

            n_removed = n_before - temp_adata.n_obs
            print("Removed {} descendants of sampled cells at time {}.".format(n_removed, t))

        else:
            pass

    ids_to_sample = [item for sublist in ids_to_sample.values() for item in sublist]
    assert np.unique(ids_to_sample).shape[0] == len(ids_to_sample)
    
    # Order ids to be in the same order as they appear in the orginal adata
    ordered_ids = []
    for i in adata.obs.index:
        if i in ids_to_sample:
            ordered_ids.append(i)
    
    return adata[ordered_ids]


def retain_multitime_only(adata, times, clone_idx, time_key='Order'):
    clone_sets = []
    for t in times:
        adata_t = adata[adata.obs[time_key] == t]
        clones_at_t = list(np.where(np.sum(adata_t.obsm['X_clone'].toarray(), axis=0) != 0)[0])
        clone_sets.append(set(clones_at_t))
    
    multitime_clones = clone_sets[0]
    for i in range(1, len(clone_sets)):
        multitime_clones = multitime_clones.intersection(clone_sets[i])
        
    multitime_mask = [i in multitime_clones for i in clone_idx]
    multitime = adata[np.sum(adata.obsm['X_clone'][:, multitime_mask], axis=1) != 0]
    return multitime
    
    
# Functions for computing the probabilities relating to the growth effect bias ----------------------------------------------------

def probability_cell_multitime_and_of_type(adata, times, sampling_rate_t1, sampling_rate_t2, 
                                           time_key='Order', growth_key='n_children', state_key='Cell Types'):
    assert len(times) == 2, "Times must contain two timepoints"
    
    temp_adata = adata.copy()
    
    # Extract just the cells at t1
    temp_adata_t1 = temp_adata[temp_adata.obs[time_key] == times[0]]

    cell_types = list(temp_adata_t1.obs['Cell Types'].unique())
    
    # Extract just the cells at t2
    temp_adata_t2 = temp_adata[temp_adata.obs[time_key] == times[1]]

    cell_types = list(temp_adata_t1.obs[state_key].unique())


    # Compute parameters required for the probability computation from the adata
    n_cells_at_t1_total = temp_adata_t1.n_obs
    proba = {}

    for i in range(len(cell_types)):
        c = cell_types[i]
        cells_in_c_at_t1 = temp_adata_t1[temp_adata_t1.obs[state_key] == c]

        clones_in_l_at_t1 = cells_in_c_at_t1.obs.Clone.unique()

        proba_in_clone = np.zeros(clones_in_l_at_t1.shape[0])

        # NOTE: This computation assumes no heterogeneity in clones at t1
        for j in range(clones_in_l_at_t1.shape[0]):
            l = clones_in_l_at_t1[j]
            cells_in_l_at_t1 = cells_in_c_at_t1[cells_in_c_at_t1.obs['Clone'] == l]

            # TO DO: update to use a prior rather than just using the non-integer number for growth
            # Compute the growth rate of the clone at t1 (i.e. the total number of cells in the clone at t2)
            growth_t1_t2 = cells_in_l_at_t1.obs[growth_key].pow(times[1]-times[0])
            growth_rate_at_t1_per_clone = growth_t1_t2.sum()
            
            # Check the growth is the same as the number of cells in the clone at t2
            n_cells_in_l_at_t2 = temp_adata_t2[temp_adata_t2.obs['Clone'] == l].n_obs
            assert growth_rate_at_t1_per_clone == n_cells_in_l_at_t2
            
            # Check number of cells in clone is as expected
            count_cells_in_clone = cells_in_l_at_t1.obsm['X_clone'].sum(axis=0)
            assert count_cells_in_clone.nonzero()[0].shape[0] == 1
            assert count_cells_in_clone.sum() == cells_in_l_at_t1.n_obs
            n_cells_in_clone_at_t1 = cells_in_l_at_t1.n_obs
            
            # Compute the probability of being the clone
            proba_sample_ancestor = sampling_rate_t1
            clone_proportion = n_cells_in_clone_at_t1/n_cells_at_t1_total
            proba_no_descendant_sampled = (1-sampling_rate_t2)**growth_rate_at_t1_per_clone
            proba_in_clone[j] = proba_sample_ancestor*(1 - proba_no_descendant_sampled)*clone_proportion
        
        # Sum over clones of type c to get prob of being type c
        proba[c] = proba_in_clone.sum()
    
    return proba


# def probability_cell_multitime_and_of_type(adata, times, sampling_rate_t1, sampling_rate_t2, time_key='Order', growth_key='growth'):
#     assert len(times) == 2, "Times must contain two timepoints"
    
#     temp_adata = adata.copy()
    
#     # Extract just the cells at t1
#     temp_adata_t1 = temp_adata[temp_adata.obs[time_key] == times[0]]
#     cell_types = list(temp_adata_t1.obs['Cell Types'].unique())
    
#     # Extract just the cells at t2
#     temp_adata_t2 = temp_adata[temp_adata.obs[time_key] == times[1]]

#     # Compute parameters required for the probability computation from the adata
#     growth_rate_at_t1_per_type = np.zeros(len(cell_types))
#     n_cells_in_clones_at_t1 = np.zeros(len(cell_types))
#     n_cells_t1_in_types = np.zeros(len(cell_types))
#     for i in range(len(cell_types)):
#         c = cell_types[i]
#         cells_in_c_at_t1 = temp_adata_t1[temp_adata_t1.obs['Cell Types'] == c]
        
#         # TO DO: update to use a prior rather than just using the non-integer number for growth
#         # Compute the growth rate as the mean across the cell type at t1
#         growth_t1_t2 = cells_in_c_at_t1.obs[growth_key].pow(times[1]-times[0])
#         growth_rate_at_t1_per_type[i] = growth_t1_t2.mean()
        
#         # TO DO: update to use a prior rather than just using the non-integer number for cells in clones
#         # Compute the mean number of cells in clones across the cell type at t1
#         n_cells_in_clones_at_t1[i] = temp_adata_t1.obsm['X_clone'].sum(axis=0).mean()
        
#         # Count the number of cell of type c at t1
#         n_cells_t1_in_types[i] = cells_in_c_at_t1.n_obs
    
#     # Compute the probability
#     n_cells_at_t1_total = temp_adata_t1.n_obs
#     proba = {}
#     for i in range(len(cell_types)):
#         c = cell_types[i]
#         proba_sample_ancestor = sampling_rate_t1
#         proba_no_descendant_sampled = (1-sampling_rate_t2)**(growth_rate_at_t1_per_type[i]*n_cells_in_clones_at_t1[i])
#         cell_type_proportion = n_cells_t1_in_types[i]/n_cells_at_t1_total
#         proba[c] = proba_sample_ancestor*(1 - proba_no_descendant_sampled)*cell_type_proportion
    
#     return proba


def probability_cell_of_type_given_multitime(adata, times, sampling_rate_t1, sampling_rate_t2, time_key='Order', growth_key='growth'):
    
    proba_both = probability_cell_multitime_and_of_type(adata, times, sampling_rate_t1, sampling_rate_t2, time_key, growth_key)
    proba_multitime =  np.sum(list(proba_both.values()))
    adata_t1 = adata[adata.obs[time_key] == times[0]]
    cell_types = list(adata_t1.obs['Cell Types'].unique())

    proba = {}
    for c in cell_types:
        proba[c] = proba_both[c]/proba_multitime
    return proba 
    
    
# Functions for estimating growth rate from cell signature scores -----------------------------------------------------
def logistic(x, L, k, x0=0):
    f = L / (1 + np.exp(-k * (x - x0)))
    return f

def gen_logistic(p, beta_max, beta_min, pmax, pmin, center, width):
    return beta_min + logistic(p, L=beta_max - beta_min, k=4 / width, x0=center)

def beta(p, beta_max=1.7, beta_min=0.3, pmax=1.0, pmin=-0.5, center=0.25):
    return gen_logistic(p, beta_max, beta_min, pmax, pmin, center, width=0.5)

def delta(a, delta_max=1.7, delta_min=0.3, amax=0.5, amin=-0.4, center=0.1):
    return gen_logistic(a, delta_max, delta_min, amax, amin, center,
                          width=0.2)
