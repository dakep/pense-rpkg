Fixes the CRAN build warning reported by GCC 14.1.0.

## Test environments
* local x86_64-apple-darwin20, R version 4.4.0 Patched (2024-04-30 r86503)
* win-builder (devel, release, oldrelease)
* mac-builder (release: r-devel-macosx-arm64|4.6.0|macosx|macOS 26.2 (25C56)|Mac mini|Apple M1||en_US.UTF-8|macOS 14.4|clang-1700.6.3.2|GNU Fortran (GCC) 14.2.0)
* GitHub actions:
  * macOS-latest; R-release: aarch64-apple-darwin20|R version 4.5.2 (2025-10-31)
  * windows-latest, R-release: x86_64-w64-mingw32|4.5.2 (2025-10-31 ucrt)
  * ubuntu-latest; R-devel: x86_64-pc-linux-gnu|R Under development (unstable) (2026-01-10 r89298)
  * ubuntu-latest; R-release: x86_64-pc-linux-gnu|R version 4.5.2 (2025-10-31)
  * ubuntu-latest; R-oldrel-1: x86_64-pc-linux-gnu|R version 4.4.3 (2025-02-28)

## R CMD check results

0 ERRORs | 0 WARNINGs | 1 NOTEs

- Some systems report a 403 status for the DOI link https://www.doi.org/10.1080/00401706.2025.2540970
