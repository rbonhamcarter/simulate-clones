# `simulate-clones`

This repository contains the source code to simulate simple lineage barcoded single cell datasets. It was designed to be used to investigate the presence of a sampling bias in the cell type proportions represented in clones observed over multiple time points, and what effect that bias may have on trajectory inference methods.

The trajectory inference methods tested in this repository are CoSpar (Wang et. al) and LineageOT (Forrow and Schiebinger). For each method two variations were tested, one using only single-time clonal data, and one using multi-time clonal data. These variations are denoted using the "-ST" and "-MT" subscripts. The version of LineageOT used in this study is currently available on this branch of the fork of the method's repository: https://github.com/rbonhamcarter/LineageOT/tree/multi-time-clones.

## For the simulations supporting the derivation of the probabilities

The script `simulate_proportions_to_verify_equations.ipynb` can be used to run simple simulations of cell type proportions. The cell state space is not simulated, only cell growth, cell type labels, sampling and barcode detection. Plots comparing simulated and predicted cell type proportions like the figure in the paper can be generated using this script.

## For the simulations in the study of the impact of the bias on trajectory inference methods

To gather a set of results for one simulation scenario, run the following set of scripts in order:

1. `specify_rate_and_simulate.py`: specify max and min growth rates for between $t_1$ and $t_2$, and a desired $t_1$ sampling rate. Simulates a population over three timepoints $t_0, t_1, t_2$ such that sampling 1000 cells at $t_1$ will give the desired $t_1$ sampling rate. Default times: $t_0=3, t_1=5, t_2=10$.

2. `sample_and_infer_tmaps.py`: specify an adata, $t_1$ sampling rate and barcoding rate (for $t_1$ and $t_2$). Will compute the true $t_1-t_2$ developmental trajectory coupling on the full adata, and then sample the adata at $t_1$ and $t_2$, strip lineage barcodes from the proportion of the population specifed by the barcoding rate, and infer $t_1-t_2$ developmental trajectory couplings using 4 methods: CoSpar-ST, CoSpar-MT, LineageOT-ST and LineageOT-MT.

3. `compute_result.py`: given a set of couplings from the methods and a true coupling (such as from `sample_and_infer_tmaps.py`), computes a number of metrics and statistics about the full/sampled dataset. Statistics include the true and predicted $t_1$ cell type proportions in the cells in multi-time clones, growth rate and cell type info at each timepoint. Metrics/summary results for the couplings include: transition tables, correlation and Wasserstein distance between true and estimated fate probabilities and true and estimated ancestor distributions for each method (CoSpar-ST, CoSpar-MT, LineageOT-ST and LineageOT-MT).

4. `generate_results.ipynb`: Given a set of results (i.e. replicates of running `sample_and_infer_tmaps.py` and `compute_result.py` for a range of $t_1$ sampling rates), this notebook allows one to explore the statistics of the sampling process and plot the results over sampling rates.

Plots like those in the paper that compare two sets of results (baseline and effect) can be generated using `compare-results.ipynb`.


