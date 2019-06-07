Fix issue with incorrect BLAS/LAPACK prototypes to restore compatibility with the upcoming release of RcppArmadillo 0.9.500.

## Test environments

* local OS X 10.14.5, R 3.6.0
* win-builder (devel and release)
* Rhub
  * Debian Linux, R-release, GCC
  * Debian Linux, R-devel, GCC ASAN/UBSAN
  * Fedora Linux, R-devel, GCC
  * Fedora Linux, R-devel, clang, gfortran
  * Windows Server 2008 R2 SP1, R-oldrel, 32/64 bit
  * Windows Server 2008 R2 SP1, R-devel, 32/64 bit
## R CMD check results

0 ERRORs | 0 WARNINGs | 2 NOTEs.

* on some systems the size of the shared library is > 5MB.
  This is due to the compiler and linker on these systems and beyond our control.
* on R-release win-builder, the checks emit a note that the package is using non-staged installation for x64.
  This seems to be a bug in the checks performed on win-builder.
