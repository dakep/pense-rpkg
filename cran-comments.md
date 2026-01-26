Fixes the CRAN build error reported by the MKL check system and the malformed doi in the DESCRIPTION file.

## Test environments
* local x86_64-apple-darwin20, R version 4.5.2 Patched (2025-12-15 r89178)
* win-builder (devel, release, oldrelease)
* GitHub actions:
  * R-hub-mkl; R-devel: x86_64-pc-linux-gnu|R Under development (unstable) (2026-01-25 r89330)|gcc with MKL (libmkl_gf_lp64.so.2)
  * R-hub-atlas; R-devel: x86_64-pc-linux-gnu|R Under development (unstable) (2026-01-25 r89330)|gcc with ATLAS (libsatlas.so.3.10)
  * macOS-latest; R-release: aarch64-apple-darwin20|R version 4.5.2 (2025-10-31)
  * windows-latest, R-release: x86_64-w64-mingw32|4.5.2 (2025-10-31 ucrt)
  * ubuntu-latest; R-devel: x86_64-pc-linux-gnu|R Under development (unstable) (2026-01-25 r89330)
  * ubuntu-latest; R-release: x86_64-pc-linux-gnu|R version 4.5.2 (2025-10-31)
  * ubuntu-latest; R-oldrel-1: x86_64-pc-linux-gnu|R version 4.4.3 (2025-02-28)

## R CMD check results

0 ERRORs | 0 WARNINGs | 1 NOTEs

- Some systems report a 403 status for the DOI link https://doi.org/10.1080/00401706.2025.2540970
