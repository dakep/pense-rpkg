Major refactoring of code and new functionality.

## Test environments
* local OS X 10.15.6, R version 4.0.2 Patched
* win-builder (devel and release)
* Travis
  * Ubuntu 16.04.7 LTS, R-oldrel
  * Ubuntu 16.04.7 LTS, R-release
  * Ubuntu 16.04.7 LTS, R-devel
* Rhub
  * Debian Linux, R-release, GCC
  * Ubuntu Linux 16.04 LTS, R-devel with rchk
  * Debian Linux, R-devel, GCC ASAN/UBSAN
  * Fedora Linux, R-devel, GCC
  * Fedora Linux, R-devel, clang, gfortran
  * Oracle Solaris 10, x86, 32 bit, R-release
  * Windows Server 2008 R2 SP1, R-oldrel, 32/64 bit
  * Windows Server 2008 R2 SP1, R-patched, 32/64 bit
## R CMD check results

0 ERRORs | 0 WARNINGs | 1 NOTEs.

### NOTEs
> There are ::: calls to the package's namespace in its code. A package almost never needs to use ::: for
> its own objects:
> ‘.__cluster_exported_objects__’ ‘.cluster_export_object’ ‘.cluster_unexport_object’ ‘.elnet_args’ ‘.pense_args’
> ‘.regmest_args’ ‘predict.pense_cvfit’ ‘predict.pense_fit

Methods implementing repeated cross-validation allow the user to provide a parallel cluster (set-up by the user).
These methods run R code on the clusters and access package-internal functions using `:::`.

> installed size is > 5MB

The package uses heavily-templated code (both internal and from RcppArmadillo), bloating the size of the library.
