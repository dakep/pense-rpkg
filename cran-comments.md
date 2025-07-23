Fixes the CRAN build warning reported by GCC 14.1.0.

## Test environments
* local x86_64-apple-darwin20, R version 4.4.0 Patched (2024-04-30 r86503)
* win-builder (devel, release, oldrelease)
* GitHub actions:
  * macOS-latest; R-release
  * windows-latest; R-release
  * ubuntu-latest; R-devel
  * ubuntu-latest; R-release
  * ubuntu-latest; R-oldrel-1

## R CMD check results

0 ERRORs | 0 WARNINGs | 1 NOTEs.

> installed size is > 5MB

The package uses heavily-templated code (both internally and from RcppArmadillo), increasing the size of the library.
