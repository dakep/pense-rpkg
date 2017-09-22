This is the initial release of the package.

## Test environments

* local OS X 10.11.6, R 3.4.1
* CentOS 7, R 3.4.0, icc, ifort
* Windows 7 SP1 (x64), R 3.4.1
* Rhub
  * Debian Linux, R-release, GCC
  * Debian Linux, R-devel, GCC ASAN/UBSAN
  * Fedora Linux, R-devel, GCC
  * Fedora Linux, R-devel, clang, gfortran
  * Windows Server 2008 R2 SP1, R-oldrel, 32/64 bit
  * Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  * Windows Server 2008 R2 SP1, R-patched, 32/64 bit
## R CMD check results

0 ERRORs | 0 WARNINGs | 1 NOTE.

* on some systems the size of the shared library is > 5MB (13MB on CentOS 7 using icc; 6.9MB on Fedora Linux, R-devel, GCC).
    This is due to the compiler and linker on these systems and beyond our control.
