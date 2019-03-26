This release fixes the warnings on `r-devel-linux-x86_64-fedora-clang` and
`r-devel-linux-x86_64-fedora-gcc` about missing autoconf templates.
The previous submission of pense_1.2.2 triggered 2 NOTEs on the incoming checks.
After deliberations on the r-package-devel list (https://hypatia.math.ethz.ch/pipermail/r-package-devel/2019q1/003680.html) it appears this is a bug on the win-builder service and out of my control.

## Test environments

* local OS X 10.14.3, R 3.5.3
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
* on R-devel win-builder, the checks emit a note that the package is using non-staged installation for x64.
  This seems to be a bug in the checks performed on win-builder.
