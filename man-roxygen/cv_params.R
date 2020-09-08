#' @param cv_k number of folds per cross-validation.
#' @param cv_repl number of cross-validation replications.
#' @param cv_metric either a string specifying the performance metric to use, or a function to evaluate prediction
#'    errors in a single CV replication. If a function, it is called with a single numeric vector of
#'    prediction errors and must return a scalar number.
#' @param fit_all If `TRUE`, fit the model for all penalization levels. Otherwise, only at penalization
#'    level with smallest average CV performance.
#' @param cl a [parallel][parallel::makeCluster] cluster. Can only be used if `ncores = 1`, because multi-threading
#'    can not be used in parallel R sessions on the same host.
