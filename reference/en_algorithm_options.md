# Control the Algorithm to Compute (Weighted) Least-Squares Elastic Net Estimates

The package supports different algorithms to compute the EN estimate for
weighted LS loss functions. Each algorithm has certain characteristics
that make it useful for some problems. To select a specific algorithm
and adjust the options, use any of the `en_***_options` functions.

## Details

- [`en_lars_options()`](en_lars_options.md): Use the tuning-free LARS
  algorithm. This computes *exact* (up to numerical errors) solutions to
  the EN-LS problem. It is not iterative and therefore can not benefit
  from approximate solutions, but in turn guarantees that a solution
  will be found.

- [`en_cd_options()`](en_cd_options.md): Use an iterative coordinate
  descent algorithm which needs \\O(n p)\\ operations per iteration and
  converges sub-linearly.

- [`en_admm_options()`](en_admm_options.md): Use an iterative ADMM-type
  algorithm which needs \\O(n p)\\ operations per iteration and
  converges sub-linearly.

- [`en_dal_options()`](en_dal_options.md): Use the iterative Dual
  Augmented Lagrangian (DAL) method. DAL needs \\O(n^3 p^2)\\ operations
  per iteration, but converges exponentially.

## See also

Other LS-EN algorithm options:
[`en_admm_options()`](en_admm_options.md),
[`en_cd_options()`](en_cd_options.md),
[`en_dal_options()`](en_dal_options.md),
[`en_lars_options()`](en_lars_options.md)
