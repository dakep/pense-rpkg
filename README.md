
# pense R package

<!-- begin badges -->

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/pense)](https://CRAN.R-project.org/package=pense)
[![Build
Status](https://travis-ci.com/dakep/pense-rpkg.svg?branch=master)](https://travis-ci.com/dakep/pense-rpkg)
<!-- end badges -->

This R package implements the penalized elastic net S-estimator (PENSE)
and the penalized M-step (PENSEM) as proposed in [Cohen Freue, et al.
(2019)](https://projecteuclid.org/euclid.aoas/1574910036) as well as the
adaptive extensions developed in [Kepplinger
(2020)](https://hdl.handle.net/2429/75637).

## Migrating from pense versions 1.x to 2.x

Version 2.x, release September 2020, introduces many new features and
improved computational speed. The changes, however, make the package
incompatible with previous versions. The [migration
guide](https://dakep.github.io/pense-rpkg/articles/migration_guide.html)
helps migrating existing code to new versions of the package.

## Usage

The main function for most users is `adapense_cv()`, which computes an
adaptive PENSE fit and estimates prediction performance for many values
of the penalization level.

### Example

``` r
library(pense)

# Generate dummy data with n=50 observations and p=25 possible predictors
# (of which only the first 3 are truly relevant).
n <- 50
p <- 25
set.seed(123)
x <- matrix(rt(n * p, df = 5), ncol = p)
y <- x[, 1] + 0.5 * x[, 2] + 2 * x[, 3] + rt(n, df = 2)

# Compute an adaptive PENSE fit with hyper-parameters alpha=0.8, exponent=2,
# and 50 different values for the penalization level.
# Prediction performance of the 50 fits for different penalization levels
# is estimated with 5-fold cross-validation,
# repeated 10 times (the more the better!).
# On selected platforms use 2 CPU cores.
set.seed(123) # Setting the seed is suggested for reproducibility of the CV results.
fit <- adapense_cv(x, y, alpha = 0.9, cv_k = 5, cv_repl = 10, ncores = 2)

# Visualize the estimated prediction performance using plot()
plot(fit)

# Summarize the model with best prediction performance using summary()
summary(fit)

# Summarize the model with "almost as-good prediction performance" as the best
# model using summary(lambda = "se")
summary(fit, lambda = "se")
```

The user can also compute a non-adaptive PENSE fit by using `pense_cv()`
in place of `adapense_cv()` or a PENSEM fit by using `pensem_cv()`.

### Using a parallel cluster

The cross-validation functions in the pense package (those ending in
`_cv()`) also accept a parallel cluster instead of the `ncores=`
argument. Setting the seed ensures the results are reproducible and the
same as if not using a cluster\!

Continuing the example above, this could look like

``` r
library(parallel)

# Set up the cluster using 3 CPUs on the local machine:
par_clust <- makeCluster(3)

set.seed(123) # Setting the seed is suggested for reproducibility of the CV results.
fit_with_cluster <- adapense_cv(x, y, alpha = 0.9, cv_k = 5, cv_repl = 10,
                                cl = par_clust)

stopCluster(par_clust)
```

Note that the `cl=` argument cannot be mixed with the `ncores=`
argument. Only one of them can be specified\!

## Installation

To install the latest release from CRAN, run the following R code in the
R console:

``` r
install.packages("pense")
```

The most recent stable version as well as the developing version might
not yet be available on CRAN. These can be directly installed from
github using the [devtools](https://cran.r-project.org/package=devtools)
package:

``` r
# Install the most recent stable version:
devtools::install_github("dakep/pense-rpkg")
# Install the (unstable) develop version:
devtools::install_github("dakep/pense-rpkg", ref = "develop")
```

## References

  - Kepplinger, D. (2020). *Robust estimation and variable selection in
    high-dimensional linear regression models.* University of British
    Columbia. [available online](https://hdl.handle.net/2429/75637).
  - Cohen Freue, GV. Kepplinger D, Salibián-Barrera M, Smucler E.
    (2019). Robust elastic net estimators for variable selection and
    identification of proteomic biomarkers. *Annals of Applied
    Statistics*. 13(4). [available
    online](https://projecteuclid.org/euclid.aoas/1574910036).