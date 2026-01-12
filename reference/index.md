# Package index

## Robust fitting of linear regression models

### Fitting and estimating prediction performance

Fit linear regression models for a set of penalization levels and
estimate the prediction performance via cross-validation.

- [`change_cv_measure()`](change_cv_measure.md) : Change the
  Cross-Validation Measure
- [`pense_cv()`](pense_cv.md) [`adapense_cv()`](pense_cv.md) :
  Cross-validation for (Adaptive) PENSE Estimates
- [`regmest_cv()`](regmest_cv.md) [`adamest_cv()`](regmest_cv.md) :
  Cross-validation for (Adaptive) Elastic Net M-Estimates

### Fitting only

- [`pense()`](pense.md) : Compute (Adaptive) Elastic Net S-Estimates of
  Regression
- [`regmest()`](regmest.md) : Compute (Adaptive) Elastic Net M-Estimates
  of Regression

## Plotting and printing

Methods for plotting and summarizing fits.

- [`plot(`*`<pense_cvfit>`*`)`](plot.pense_cvfit.md) : Plot Method for
  Penalized Estimates With Cross-Validation
- [`plot(`*`<pense_fit>`*`)`](plot.pense_fit.md) : Plot Method for
  Penalized Estimates
- [`prediction_performance()`](prediction_performance.md)
  [`print(`*`<pense_pred_perf>`*`)`](prediction_performance.md) :
  Prediction Performance of Adaptive PENSE Fits
- [`summary(`*`<pense_cvfit>`*`)`](summary.pense_cvfit.md)
  [`print(`*`<pense_cvfit>`*`)`](summary.pense_cvfit.md) : Summarize
  Cross-Validated PENSE Fit

## Extracting information

Methods for extracting coefficient estimates and predicting values.

- [`coef(`*`<pense_cvfit>`*`)`](coef.pense_cvfit.md) : Extract
  Coefficient Estimates
- [`coef(`*`<pense_fit>`*`)`](coef.pense_fit.md) : Extract Coefficient
  Estimates
- [`predict(`*`<pense_cvfit>`*`)`](predict.pense_cvfit.md) : Predict
  Method for PENSE Fits
- [`predict(`*`<pense_fit>`*`)`](predict.pense_fit.md) : Predict Method
  for PENSE Fits
- [`residuals(`*`<pense_cvfit>`*`)`](residuals.pense_cvfit.md) : Extract
  Residuals
- [`residuals(`*`<pense_fit>`*`)`](residuals.pense_fit.md) : Extract
  Residuals

## Robust location and scale

Compute robust location and scale estimates.

- [`mloc()`](mloc.md) : Compute the M-estimate of Location
- [`mlocscale()`](mlocscale.md) : Compute the M-estimate of Location and
  Scale
- [`mscale()`](mscale.md) : Compute the M-Scale of Centered Values
- [`tau_size()`](tau_size.md) : Compute the Tau-Scale of Centered Values

## Non-robust methods

Non-robust methods for fitting linear regression models.

- [`elnet()`](elnet.md) : Compute the Least Squares (Adaptive) Elastic
  Net Regularization Path
- [`elnet_cv()`](elnet_cv.md) : Cross-validation for Least-Squares
  (Adaptive) Elastic Net Estimates

## Advanced functionality

### Controlling the Robust EN algorithm

Options to choose and control the algorithm to optimize robust EN
objective functions

- [`cd_algorithm_options()`](cd_algorithm_options.md) : Coordinate
  Descent (CD) Algorithm to Compute Penalized Elastic Net S-estimates
- [`mm_algorithm_options()`](mm_algorithm_options.md) : MM-Algorithm to
  Compute Penalized Elastic Net S- and M-Estimates

### Controlling the EN algorithm

Options to choose and control the algorithm to optimize least-squares EN
problems.

- [`en_admm_options()`](en_admm_options.md) : Use the ADMM Elastic Net
  Algorithm
- [`en_algorithm_options`](en_algorithm_options.md) : Control the
  Algorithm to Compute (Weighted) Least-Squares Elastic Net Estimates
- [`en_cd_options()`](en_cd_options.md) : Use Coordinate Descent to
  Solve Elastic Net Problems
- [`en_dal_options()`](en_dal_options.md) : Use the DAL Elastic Net
  Algorithm
- [`en_lars_options()`](en_lars_options.md) : Use the LARS Elastic Net
  Algorithm

### Rho/Psi functions

- [`mscale_algorithm_options()`](mscale_algorithm_options.md) : Options
  for the M-scale Estimation Algorithm
- [`consistency_const()`](rho-tuning-constants.md)
  [`efficiency_const()`](rho-tuning-constants.md) : Get the Constant for
  Consistency for the M-Scale and for Efficiency for the M-estimate of
  Location
- [`rho_function()`](rho_function.md) : List Available Rho Functions

### Initial estimates

Manually compute and alter initial estimates.

- [`enpy_initial_estimates()`](enpy_initial_estimates.md) : ENPY Initial
  Estimates for EN S-Estimators
- [`enpy_options()`](enpy_options.md) : Options for the ENPY Algorithm
- [`prinsens()`](prinsens.md) : Principal Sensitivity Components
- [`starting_point()`](starting_point.md)
  [`as_starting_point()`](starting_point.md) : Create Starting Points
  for the PENSE Algorithm

### Miscellaneous
