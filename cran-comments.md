Addressing issues raised by the previous submission of v. 1.0.4:
 * added a reference in the DESCRIPTION file to the manuscript for more details
 * added a COPYRIGHTS file to list copyright holders of some code parts
 * added examples for all functions in the package

Fixed some smaller bugs and made the interface more consistent.

## Test environments

* local OS X 10.11.6, R 3.4.1
* CentOS 7, R 3.4.0, icc, ifort
* Windows 7 SP1 (x64), R 3.4.1
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

* on some systems the size of the shared library is > 5MB (6.9MB on Fedora Linux, R-devel, GCC).
    This is due to the compiler and linker on these systems and beyond our control.
