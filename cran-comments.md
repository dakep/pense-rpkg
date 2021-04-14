Fixes several bugs which, e.g., prevented building on RHEL and degraded performance in certain situations.

## Test environments
* local x86_64-apple-darwin17.0 version 11.2.3 (20D91), R version 4.0.5 Patched (2021-03-31 r80136)
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
