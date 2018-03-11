This release addresses bugs that have been reported by users and the build error on Solaris (x86).
We also improved the internal algorithm and the way we determine the maximum value of the regularization parameter.

## Test environments

* local OS X 10.13.3, R 3.4.3
* CentOS 7, R 3.4.3, icc, ifort
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
