# pense R package

[![R CMD
check](https://github.com/dakep/pense-rpkg/actions/workflows/check-package-extended.yml/badge.svg)](https://github.com/dakep/pense-rpkg/actions/workflows/check-package-extended.yml)

This R package implements the penalized adaptive elastic net S-estimator
(adaptive PENSE), the non-adaptive version (PENSE) and the penalized
M-step (PENSEM) as proposed in [Cohen Freue, et
al. (2019)](https://doi.org/10.1214/19-AOAS1269) and [Kepplinger
(2023)](https://doi.org/10.1016/j.csda.2023.107730). Hyper-parameter
selection is done by robust information sharing cross-validation
(RIS-CV; [Kepplinger & Wei
2025](https://doi.org/10.1080/00401706.2025.2540970)).

## Usage

The main function for most users is
[`adapense_cv()`](reference/pense_cv.md), which computes an adaptive
PENSE fit and estimates prediction performance for many values of the
penalization level.

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

The user can also compute a non-adaptive PENSE fit by using
[`pense_cv()`](reference/pense_cv.md) in place of
[`adapense_cv()`](reference/pense_cv.md) or a PENSEM fit by using
`pensem_cv()`.

## Detailed examples

The package vignette [Estimating predictive
models](https://dakep.github.io/pense-rpkg/articles/computing_adapense.html)
demonstrate in more detail how to compute adaptive and non-adaptive
PENSE estimates.

## Installation

To install the latest release from CRAN, run the following R code in the
R console:

``` r
install.packages("pense")
```

The most recent development version can be installed directly from
github using the [devtools](https://cran.r-project.org/package=devtools)
package:

``` r
# Install the most recent development version:
devtools::install_github("dakep/pense-rpkg")
```

## References

- D. Kepplinger and S. Wei, “Information sharing for robust and stable
  cross-validation,” *Technometrics*, 2025, [doi:
  10.1080/00401706.2025.2540970](https://doi.org/10.1080/00401706.2025.2540970).
- D. Kepplinger, “Robust variable selection and estimation via adaptive
  elastic net S-estimators for linear regression,” *Computational
  Statistics & Data Analysis,* vol. 183, 2023, [doi:
  10.1016/j.csda.2023.107730](https://doi.org/10.1016/j.csda.2023.107730).
- G. V. Cohen Freue, D. Kepplinger, M. Salibián-Barrera, and E. Smucler,
  “Robust elastic net estimators for variable selection and
  identification of proteomic biomarkers,” *Annals of Applied
  Statistics,* vol. 13, no. 4, pp. 2065–2090, 2019, [doi:
  10.1214/19-AOAS1269](https://doi.org/10.1214/19-AOAS1269).
