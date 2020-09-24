Fixes a critical bug which prevents the methods from being used if more than 50% of response values are identical (e.g., discrete or 0-inflated responses).

## Test environments
* local OS X 10.15.6, R version 4.0.2 Patched
* win-builder (devel and release)
* Travis
  * Ubuntu 16.04.7 LTS, R-oldrel
  * Ubuntu 16.04.7 LTS, R-release
  * Ubuntu 16.04.7 LTS, R-devel
* Rhub
  * Ubuntu Linux 16.04 LTS, R-devel with rchk
  * Debian Linux, R-devel, GCC ASAN/UBSAN
  * Fedora Linux, R-devel, clang, gfortran (with VALGRIND)
  * Oracle Solaris 10, x86, 32 bit, R-release
  * Windows Server 2008 R2 SP1, R-oldrel, 32/64 bit
  * Windows Server 2008 R2 SP1, R-patched, 32/64 bit
## R CMD check results

0 ERRORs | 0 WARNINGs | 1 NOTEs.

> installed size is > 5MB

The package uses heavily-templated code (both internal and from RcppArmadillo), bloating the size of the library.
