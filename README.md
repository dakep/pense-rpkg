# PENSE R package
This repository hosts an R package implementing the PENSE robust Elastic Net estimators for linear regression.

## Usage
The main functions of the package are
* `pense` … compute a robust Elastic Net S-estimator for linear regression
* `mstep` … take the S-estimator returned by `pense` as initial estimate for a M-step to improve efficiency

## Installation
To install the most recent stable version of the package from github, either
use the [devtools](https://cran.r-project.org/package=devtools) package's function
```
install_github("gcohenfr/PENSE-package", auth_token = "XXX")
```
with your personal access token (more details under [settings > tokens](https://github.com/settings/tokens)),
clone the repository, or download the ZIP archive provided by github.

### Dependencies
The package depends on R packages `Rcpp`, `RcppArmadillo`, and `robustbase`.
