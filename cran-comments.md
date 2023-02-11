Fixes the build error on r-oldrel-windows-ix86+x86_64 due to the incompatibility of the OpenMP data sharing clause between GCC versions < 9 and version >= 9 ("OpenMP data sharing" in https://gcc.gnu.org/gcc-9/porting_to.html).
I updated the autoconf.hpp.win header file for Windows to detect if a GCC compiler older than version 9 is used.
If GCC < 9 is detected constant variables are omitted from the data sharing clause.
On systems supporting the configure script the compiler's behavior has been determined automatically and the issue did not arise.
The fix has been tested on win-builder (devel, release, oldrelease) and the test environments listed below.

## Test environments
* local x86_64-apple-darwin17.0, R 4.2.2 Patched (2023-02-03 r83757)
* win-builder (devel, release, oldrelease)
* GitHub actions:
  * macOS-latest; R-release
  * windows-latest; R-release
  * ubuntu-latest; R-devel
  * ubuntu-latest; R-release
  * ubuntu-latest; R-oldrel-1
* Rhub:
  * windows-x86_64-oldrel
  * windows-x86_64-devel

## R CMD check results

0 ERRORs | 0 WARNINGs | 1 NOTEs.

> installed size is > 5MB

The package uses heavily-templated code (both internally and from RcppArmadillo), increasing the size of the library.
