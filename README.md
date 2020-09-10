# PENSE R package

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/pense)](https://CRAN.R-project.org/package=pense)

This R package implements the penalized elastic net S-estimator (PENSE) and the penalized M-step (PENSEM) as proposed in [Cohen Freue, et al. (2019)](https://projecteuclid.org/euclid.aoas/1574910036) as well as the adaptive extensions developed in [Kepplinger (2020)](https://hdl.handle.net/2429/75637).

## Usage

The main function for most users is `adapense_cv()`, which computes an adaptive PENSE fit and estimates prediction performance for many values of the penalization level.

## Installation
To install the latest release from CRAN, run the following R code in the R console:
```r
install.packages("pense")
```

The most recent stable version as well as the developing version might not yet be available on CRAN.
These can be directly installed from github using the
[devtools](https://cran.r-project.org/package=devtools) package:
```r
# Install the most recent stable version:
devtools::install_github("dakep/pense-rpkg")
# Install the (unstable) develop version:
devtools::install_github("dakep/pense-rpkg", ref = "develop")
```

## References
- Kepplinger, D. (2020). *Robust estimation and variable selection in high-dimensional linear regression models.* University of British Columbia. [available online](https://hdl.handle.net/2429/75637).
- Cohen Freue, GV. Kepplinger D, Salibián-Barrera M, Smucler E. (2019). Robust elastic net estimators for variable selection and identification of proteomic biomarkers. _Annals of Applied Statistics_. 13(4). [available online](https://projecteuclid.org/euclid.aoas/1574910036).
