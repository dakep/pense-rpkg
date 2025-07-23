#' @param cv_k number of folds per cross-validation.
#' @param cv_repl number of cross-validation replications.
#' @param cv_type what kind of cross-validation should be performed:
#'   robust information sharing (`ris`) or standard (`naive`) CV.
#' @param cv_metric only for `cv_type='naive'`.
#'    Either a string specifying the performance metric to use, or a function to
#'    evaluate prediction errors in a single CV replication.
#'    If a function, the number of arguments define the data the function receives.
#'    If the function takes a single argument, it is called with a single numeric vector of
#'    prediction errors.
#'    If the function takes two or more arguments, it is called with the predicted values as
#'    first argument and the true values as second argument.
#'    The function must always return a single numeric value quantifying the prediction performance.
#'    The order of the given values corresponds to the order in the input data.
#' @param fit_all only for `cv_type='naive'`.
#'    If `TRUE`, fit the model for all penalization levels.
#'    Can also be any combination of `"min"` and `"{x}-se"`, in which case only models at the
#'    penalization level with smallest average CV accuracy, or within `{x}` standard errors,
#'    respectively.
#'    Setting `fit_all` to `FALSE` is equivalent to `"min"`.
#'    Applies to all `alpha` value.
#' @param cl a [parallel][parallel::makeCluster] cluster. Can only be used in combination with
#'   `ncores = 1`.
#'
#' @details
#' The built-in CV metrics are
#' \describe{
#'   \item{`"tau_size"`}{\eqn{\tau}-size of the prediction error, computed by
#'                       [tau_size()] (default).}
#'   \item{`"mape"`}{Median absolute prediction error.}
#'   \item{`"rmspe"`}{Root mean squared prediction error.}
#'   \item{`"auroc"`}{Area under the receiver operator characteristic curve (actually 1 - AUROC).
#'                    Only sensible for binary responses.}
#' }
