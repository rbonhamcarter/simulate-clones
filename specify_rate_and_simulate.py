import argparse
import numpy as np
import lineageot
import scanpy as sc
import pandas as pd
import cospar as cs
import anndata
import os
import scipy
from numpy.random import default_rng
rng = default_rng()
import math

import clonevalidation.src as src


parser = argparse.ArgumentParser(description='This script takes a sampling rate and max growth'
                                 + ' for each cell type, and simulates a populations such that'
                                 + ' sampling 1000 cells at t1 will be at the desired rate.')

parser.add_argument('-o', '--output_dir', metavar='O', type=str, help='Output directory.')
parser.add_argument('-r', '--sample_rate', metavar='R', type=float, help='Desired sampling rate at t1'
                    + ' used to set the size of the population at t1')
parser.add_argument('-gp', '--growth_progenitor', metavar='GP', type=float, help='Max growth for progenitor '
                    + 'cell type (num of descendants between t1 and t2).')
parser.add_argument('-ga', '--growth_type_a', metavar='GA', type=float, help='Max growth for type A '
                    + 'cell type (num of descendants between t1 and t2).')
parser.add_argument('-gpm', '--growth_progenitor_min', metavar='GP', type=float, help='Min growth for progenitor '
                    + 'cell type (num of descendants between t1 and t2).')
parser.add_argument('-gam', '--growth_type_a_min', metavar='GA', type=float, help='Min growth for type A '
                    + 'cell type (num of descendants between t1 and t2).')
parser.add_argument('-x', '--x_max', metavar='X', type=float, help='Max x range value.')
parser.add_argument('-y', '--y_max', metavar='Y', type=float, help='Max y range value.')
parser.add_argument('-c', '--clone_het', action='store_true', help='Whether or not to '
                    + 'have heteorgeneity in the clones at t1.', default=False)



def round_up_to_even(x):
    return math.ceil(x / 2.) * 2

# States (2D) --------------------------------------------

# Simulate states according to different probility distributions over the continous state space

# at t0 ------------------------------------------------------------------------------------------
def state_at_t0(cell_type, x_range, y_range):
    'All progenitors with same state distribution as at t1'
    if cell_type == 'progenitor':
        x_dist = progenitor_state_distribution_at_t1(a=x_range[0], b=(3/5)*(x_range[1]-x_range[0]) + x_range[0])
    else:
        raise ValueError("Invalid cell_type")
        
    x = x_dist.rvs()
    y = rng.uniform(low=y_range[0], high=y_range[1])
    return np.array([x, y])


# at t1 ------------------------------------------------------------------
def prob_slope(x_min, x_max):
    """
    Returns the positive slope of a linear probability density 
    function given the desired parameters specified in the arguments.
    """
    return 2.0/(x_max-x_min)**2

class type_A_state_distribution_at_t1(scipy.stats.rv_continuous):
    def _pdf(self, x):
        slope = prob_slope(self.a, self.b)
        return 0.0 + slope*(x - self.a)
    
class progenitor_state_distribution_at_t1(scipy.stats.rv_continuous):
    def _pdf(self, x):
        slope = prob_slope(self.a, self.b)
        prob_max = slope*(self.b - self.a)
        return prob_max - slope*(x - self.a)
    

def state_at_t1(cell_type, x_range, y_range):
    if cell_type == 'type A':
        x_dist = type_A_state_distribution_at_t1(a=(2/5)*(x_range[1]-x_range[0]) + x_range[0], b=x_range[1])
        y = rng.uniform(low=y_range[0], high=(3/6)*(y_range[1]-y_range[0]) + y_range[0])
    elif cell_type == 'progenitor':
        x_dist = progenitor_state_distribution_at_t1(a=x_range[0], b=(3/5)*(x_range[1]-x_range[0]) + x_range[0])
        y = rng.uniform(low=y_range[0], high=y_range[1])
    else:
        raise ValueError("Invalid cell_type")
        
    x = x_dist.rvs()
    
    return np.array([x, y])

# at t2 ---------------------------------------------------------------------------
def state_at_t2(cell_type, x_range, y_range):
    if cell_type == 'type A':
        y = rng.uniform(low=y_range[0], high=(3/6)*(y_range[1]-y_range[0]) + y_range[0])
        
    elif cell_type == 'type B':
        y = rng.uniform(low=(4/6)*(y_range[1]-y_range[0]) + y_range[0], high=y_range[1])
        
    else:
        raise ValueError("Invalid cell_type")
        
    center = (4/5)*(x_range[1]-x_range[0]) + x_range[0]
    sd = (x_range[1] - x_range[0])/10
    x = center + np.random.normal(loc=0.0, scale=sd)
    return np.array([x, y])


# Probability of state transitions (depends on position in state space)
def transition_matrix_t0_t1(y, x_range, y_range, clone_het=False):
#     # Re-use the type A state prob density as the P_progA(y) transition density
#     type_A_dist = type_A_state_distribution_at_t1(a=y_range[0], b=y_range[1])
#     P_progA = type_A_dist.pdf(y_range[1]-y)
    if clone_het:
        P_progA = 0.5

    else:
        if y <= y_range[0] + (1/2)*(y_range[1] - y_range[0]):
            P_progA = 1.0
        else:
            P_progA = 0.0
    
    tmatrix = np.array([[P_progA, 1.0-P_progA]])
    obs = pd.DataFrame(index=['progenitor'])
    var = pd.DataFrame(index=['type A', 'progenitor'])
    return anndata.AnnData(X=tmatrix, obs=obs, var=var)


def transition_matrix_t1_t2(y, x_range, y_range):
    P_AA = 1.0
    
#     # Re-use the type A state prob density as the P_progA(y) transition density
#     type_A_dist = type_A_state_distribution_at_t1(a=y_range[0], b=y_range[1])
#     P_progA = type_A_dist.pdf(y_range[1]-y)
    if y <= y_range[0] + (1/2)*(y_range[1] - y_range[0]):
        P_progA = 1.0
    else:
        P_progA = 0.0
    
    tmatrix = np.array([[P_AA, 1.0-P_AA], 
                        [P_progA, 1.0-P_progA]])
    obs = pd.DataFrame(index=['type A', 'progenitor'])
    var = pd.DataFrame(index=['type A', 'type B'])
    return anndata.AnnData(X=tmatrix, obs=obs, var=var)


# Growth ----------------------------------------------
def n_descendants(growth_rate):
    """
    Specify how many descendants a cell with this growth rate estimate has. 
    Growth rate is considered an estimate of n_descendants but may not be
    an integer. We derive n_descendants an integer by considering growth rate
    the mean of a 2-state random variable where the the states are {floor(growth_rate), floor(growth_rate) + 1}
    """
    prob_ceil = growth_rate - math.floor(growth_rate)
    prob_floor = 1.0 - prob_ceil
    n = np.random.choice(a=[math.floor(growth_rate), math.ceil(growth_rate)], p=[prob_floor, prob_ceil])
    return n

# t0 ---> t1, use growth uniform across cell types => number of cells in clones at t1 is uniform across cell types
growth_t0_t1 = 2.0 

# t1 ---> t2, use growth which increases linearly from left to right along the x-axis
def growth_t1_t2_at_x(x, x_min, x_max, cell_type, 
                      growth_t1_t2_max_A, growth_t1_t2_max_prog, growth_t1_t2_min_prog, growth_t1_t2_min_A):
    if cell_type == 'type A':
        slope = (growth_t1_t2_max_A - growth_t1_t2_min_A)/(x_max-x_min)
        return growth_t1_t2_max_A - slope*(x-x_min)
    
    elif cell_type == 'progenitor':
        slope = (growth_t1_t2_max_prog - growth_t1_t2_min_prog)/(x_max-x_min)  
        return growth_t1_t2_max_prog - slope*(x-x_min)
    
    else:
        raise ValueError



def main():
    # Parse parameters
    args = parser.parse_args()
    output_dir= args.output_dir
    r = args.sample_rate
    clone_het = args.clone_het
    max_growth = {'progenitor': args.growth_progenitor, 
                  'type A': args.growth_type_a}
    min_growth = {'progenitor': args.growth_progenitor_min, 
                  'type A': args.growth_type_a_min}

    # Globally define basic state parameters
    x_range = [0.0, args.x_max]
    y_range = [0.0, args.y_max]

    # Basic fixed parameters of the system
    t0 = 3/5
    t1 = 1
    t2 = 2
    time_points = [t1, t2]   # sample timepoints
    cell_types = {str(t0): ['progenitor'], 
                  str(t1): ['type A', 'progenitor'], 
                  str(t2): ['type A', 'type B']}

    # Compute the number of cells at t1 using r assuming 1000 cells are sampled
    # Rounding up will result in, if anything, a lower r, which is more desirable than a higher one
    N_t1 = math.ceil(1000/r)

    # Compute the number of cells at t0
    m = growth_t0_t1
    N_t0 = round_up_to_even(N_t1/m)

    # Cell type proportions at t0
    n_cells_at_t0 = {'progenitor': N_t0}   # should choose even numbers for convenience
    n_clones = sum(n_cells_at_t0.values())

    # Simulate the population at t0------------------------------------
    cell_labels = []
    clone_labels = []
    order = []
    cell_type = []
    state = []
    growth = []

    cell_counter = 0
    clone_counter = 0
    for t in cell_types[str(t0)]:
        for x in range(n_cells_at_t0[t]):
            # Annotations for cell at t0
            cell_labels.append('cell_{}'.format(cell_counter))
            clone_labels.append('clone_{}'.format(clone_counter))
            order.append(t0)
            cell_type.append(t)
            state.append(state_at_t0(t, x_range, y_range))
            growth.append(growth_t0_t1)
            cell_counter += 1
            clone_counter += 1
            
    # Store the simulated population as an anndata
    state_array = np.stack(state)
    assert state_array.shape[0] == len(cell_labels)
    obs = pd.DataFrame({'Order': order, 'Clone': clone_labels, 'Cell Types': cell_type, 'growth': growth}, 
                    index=cell_labels)

    adata_t0 = anndata.AnnData(X=state_array, obs=obs)

    # Simulate the population at t1-----------------------------------------------
    cell_labels = []
    clone_labels = []
    order = []
    cell_type = []
    state = []
    growth = []
    parent = []
    n_children_t1 = []
        
    cell_counter = adata_t0.n_obs   # start from last cell label
    for cell in adata_t0.obs_names:
        cell_clone = adata_t0.obs['Clone'][cell]
        current_cell_state_y = adata_t0[cell].X[0][1]
        current_cell_type = adata_t0.obs['Cell Types'][cell]
        
        current_cell_growth = adata_t0.obs['growth'][cell]
        n_children = n_descendants(current_cell_growth)
        n_children_t1.append(n_children)
        cell_transition_matrix = transition_matrix_t0_t1(current_cell_state_y, x_range, y_range, clone_het)
        
        for i in range(n_children):
            # Annotations for cell's children at t1
            cell_labels.append('cell_{}'.format(cell_counter))
            clone_labels.append(cell_clone)    # same clone as parent
            parent.append(cell)
            order.append(t1)
            cell_t = src.get_next_cell_type(current_cell_type, cell_transition_matrix)
            cell_type.append(cell_t)
            
            current_cell_state = state_at_t1(cell_t, x_range, y_range)
            
            # Growth at t1
            current_cell_growth = growth_t1_t2_at_x(
                current_cell_state[0], x_max=x_range[1], x_min=x_range[0], cell_type=cell_t, 
                growth_t1_t2_max_A=max_growth['type A'], 
                growth_t1_t2_max_prog=max_growth['progenitor'], 
                growth_t1_t2_min_A=min_growth['type A'], 
                growth_t1_t2_min_prog=min_growth['progenitor'])
            
            state.append(current_cell_state)
            growth.append(current_cell_growth)
            cell_counter += 1


    # Check the number of cells of each type are as expected
    n_cells_at_t1 = dict.fromkeys(cell_types[str(t1)])
    for t in cell_types[str(t1)]:
        n_cells_at_t1[t] = np.count_nonzero(np.array(cell_type) == t)
    print("Number of cells of each type at t1: ", n_cells_at_t1)
            

    # Store the simulated population as an anndata
    state_array = np.stack(state)
    assert state_array.shape[0] == len(cell_labels)
    obs = pd.DataFrame({'Order': order, 'Clone': clone_labels, 
                        'Cell Types': cell_type, 'growth': growth, 'parent': parent}, 
                    index=cell_labels)

    adata_t1 = anndata.AnnData(X=state_array, obs=obs)
    adata_t0.obs['n_children'] = n_children_t1

    # Simulate the population at t2 --------------------------------
    cell_labels = []
    clone_labels = []
    order = []
    cell_type = []
    state = []
    growth = []
    parent = []
    n_children_t2 = []
        
    cell_counter = adata_t0.n_obs + adata_t1.n_obs   # start from last cell label

    for cell in adata_t1.obs_names:
        cell_clone = adata_t1.obs['Clone'][cell]
        current_cell_state_y = adata_t1[cell].X[0][1]
        current_cell_type = adata_t1.obs['Cell Types'][cell]
        
        current_cell_growth = adata_t1.obs['growth'][cell]
        n_children = n_descendants(current_cell_growth)
        n_children_t2.append(n_children)
        cell_transition_matrix = transition_matrix_t1_t2(current_cell_state_y, x_range, y_range)

        for i in range(n_children):
            # Annotations for cell's children at t2
            cell_labels.append('cell_{}'.format(cell_counter))
            clone_labels.append(cell_clone)    # same clone as parent
            parent.append(cell)
            order.append(t2)
            cell_t = src.get_next_cell_type(current_cell_type, cell_transition_matrix)
            cell_type.append(cell_t)
            
            current_cell_state = state_at_t2(cell_t, x_range, y_range)
            
            # Assume growth at t2 has the same distribution as growth at t1 for type A
            current_cell_growth = growth_t1_t2_at_x(
                current_cell_state[0], x_max=x_range[1], x_min=x_range[0], cell_type='type A', 
                growth_t1_t2_max_A=max_growth['type A'], 
                growth_t1_t2_max_prog=max_growth['progenitor'], 
                growth_t1_t2_min_A=min_growth['type A'], 
                growth_t1_t2_min_prog=min_growth['progenitor'])
            
            state.append(current_cell_state)
            growth.append(current_cell_growth)
            cell_counter += 1

    # Check the number of cells of each type are as expected
    n_cells_at_t2 = dict.fromkeys(cell_types[str(t2)])
    for t in cell_types[str(t2)]:
        n_cells_at_t2[t] = np.count_nonzero(np.array(cell_type) == t)
    print("Number of cells of each type at t2: ", n_cells_at_t2)


    # Store the simulated population as an anndata
    state_array = np.stack(state)
    assert state_array.shape[0] == len(cell_labels)
    obs = pd.DataFrame({'Order': order, 'Clone': clone_labels, 
                        'Cell Types': cell_type, 'growth': growth, 'parent': parent}, 
                    index=cell_labels)

    adata_t2 = anndata.AnnData(X=state_array, obs=obs)
    adata_t1.obs['n_children'] = n_children_t2

    #-------------------------------------------------------------

    # Combine the adatas
    adata_t0.obs['parent'] = np.nan
    adata_t2.obs['n_children'] = 1.0    # to try to get cospar-st using it correctly
    full_adata = anndata.concat([adata_t0, adata_t1, adata_t2], )
    adata = anndata.concat([adata_t1, adata_t2])

    # Create an cell by clone matrix (clone membership matrix)  
    X_clone = np.zeros((adata.n_obs, n_clones))
    full_X_clone = np.zeros((full_adata.n_obs, n_clones))

    for i in range(adata.n_obs):
        clone_label = adata.obs['Clone'][i]
        clone_idx = eval(clone_label.split('_')[-1])
        X_clone[i, clone_idx] = 1

    for i in range(full_adata.n_obs):
        clone_label = full_adata.obs['Clone'][i]
        clone_idx = eval(clone_label.split('_')[-1])
        full_X_clone[i, clone_idx] = 1
            
    assert np.all(np.sum(X_clone, axis=1) == 1), "At least one cell is assigned to more than one clone."
    assert np.all(np.sum(X_clone, axis=0) > 0), "At least one clone has no cells."
    adata.obsm['X_clone'] = X_clone
    full_adata.obsm['X_clone'] = full_X_clone

    # Save the full simulated adata to disk
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    r_label = str(r).replace('.', '-')
    growth_label = str(max_growth['progenitor']).replace('.', '-')
    output_filename = 'full_adata_{}_{}.h5ad'.format(growth_label, r_label)
    full_adata.write_h5ad(os.path.join(output_dir, output_filename))

if __name__ == '__main__':
    main()