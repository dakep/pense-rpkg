Introduces new functionality, improves speed, and fixes an OpenMP bug with Intel compilers.

## Test environments
* local x86_64-apple-darwin17.0 (64-bit), 4.1.0 Patched (2021-06-19 r80532)
* win-builder (devel and release)
* Travis
  * Ubuntu 18.04 LTS, R-oldrel
  * Ubuntu 18.04 LTS, R-release
  * Ubuntu 18.04 LTS, R-devel
* Rhub
  * Ubuntu Linux 20.04.1 LTS, R-devel with rchk
  * Debian Linux, R-devel, GCC ASAN/UBSAN
  * Fedora Linux, R-devel, clang, gfortran (with VALGRIND)
  * Fedora Linux, R-devel, GCC
  * Oracle Solaris 10, x86, 32 bit, R-release
  * Windows Server 2008 R2 SP1, R-oldrel, 32/64 bit
  * Windows Server 2008 R2 SP1, R-patched, 32/64 bit
## R CMD check results

0 ERRORs | 0 WARNINGs | 1 NOTEs.

> installed size is > 5MB

The package uses heavily-templated code (both internal and from RcppArmadillo), increasing the size of the library.
