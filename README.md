# npr_multiple_nets
The code implements multi-graphon (and graphon) estimator as proposed in
'Nonparametric regression for multiple heterogeneous networks' by S. Chandna and P.A. Maugis.
Available on: arXiv preprint arXiv:2001.04938

Input: A collection of m binary or weighted networks (adjacencies), each of size n x n , observed on the same set of nodes; network-level covariates (optional)
Output: Intensity (or probability) array (of size n x n x m) of pairwise interaction between node pairs in each network

Example code uses brain network data Templeton255 publicly available at: https://neurodata.io/mri/# (Kiar et al., 2016)
