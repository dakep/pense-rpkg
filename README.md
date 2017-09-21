# PENSE R package
This R package implements the Penalized Elastic Net S-Estimator (PENSE) and MM-estimator (PENSEM)
for linear regression.

## Usage
The main functions in the package are
* `pense()` … to compute a robust elastic net S-estimator for linear regression
* `pensem()` … to compute a robust elastic net MM-estimator either directly from the data matrix or
    from an S-estimator previously computed with `pense()`.

Both of these functions perform k-fold cross-validation to choose the optimal penalty level
`lambda`, but the optimal balance between the L1 and the L2 penalties (the `alpha` parameter) needs
to be pre-specified by the user.

The default breakdown point is set to 25%. If the user needs an estimator with a higher breakdown
point, the `delta` argument in the `pense_options()` and `initest_options()` can be set to the
desired breakdown point (.e.g, `delta = 0.5`).

The package also exports an efficient classical elastic net algorithm available via the functions
`elnet()` and `elnet_cv()` which chooses an optimal penalty parameter based on cross-validation.
The elastic net solution is computed either by the augmented LARS algorithm
(`en_options_aug_lars()`) or via the Dual Augmented Lagrangian algorithm (Tomioka, et al. 2011)
selected with `en_options_dal()` which is much faster in case of a large number of predictors
(> 500) and a small number of observations (< 200).


## References
Tomioka, R., Suzuki, T., and Sugiyama, M. (2011). Super-linear convergence of dual augmented lagrangian
algorithm for sparsity regularized estimation. _The Journal of Machine Learning Research_, **12**:1537–1586.
