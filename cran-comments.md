This release includes the following minor changes:
* prepare for the upcoming changes to Rcpp (using `STRICT_R_HEADERS`)
* fix a bug when computing PSCs with the augmented LARS algorithm for EN

## Test environments

* local OS X 10.14.2, R 3.5.2
* win-builder (devel and release)
* Rhub
  * Debian Linux, R-release, GCC
  * Debian Linux, R-devel, GCC ASAN/UBSAN
  * Fedora Linux, R-devel, GCC
  * Fedora Linux, R-devel, clang, gfortran
  * Windows Server 2008 R2 SP1, R-oldrel, 32/64 bit
  * Windows Server 2008 R2 SP1, R-patched, 32/64 bit
## R CMD check results

0 ERRORs | 0 WARNINGs | 1 NOTE.

* on some systems the size of the shared library is > 5MB.
    This is due to the compiler and linker on these systems and beyond our control.
