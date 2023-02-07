Adds new algorithm for LS-EN, revamps the regularization path and improves CV strategies; addresses anticipated revdep issues with the forthcoming release of Rcpp.

## Test environments
* local x86_64-apple-darwin17.0, R 4.2.2 Patched (2023-02-03 r83757)
* win-builder (devel and release)
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
