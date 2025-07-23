#' Change the Cross-Validation Measure
#'
#' For cross-validated fits using the RIS-CV strategy, the measure of prediction
#' accuracy can be adjusted post-hoc.
#'
#' @param x fitted (adaptive) PENSE or M-estimator
#' @param measure the measure to use for prediction accuracy
#' @param max_solutions consider only this many of the best solutions.
#'   If missing, all solutions are considered.
#' @return a `pense.cvfit` object using the updated measure of prediction accuracy
#' @importFrom rlang abort warn
#' @export
change_cv_measure <- function (x,
                               measure = c('wrmspe', 'wmape', 'tau_size', 'wrmspe_cv', 'wmape_cv'),
                               max_solutions = Inf) {
  if (!inherits(x, "pense_cvfit")) {
    abort("Cross-validated fit required")
  }

  measure <- match.arg(measure)

  if (!identical(x$cv_type, "ris")) {
    warn('Only RIS-CV fits (`cv_type="ris"`) support different CV measures.')
    return(x)
  }

  cvres <- x$cv_ris[, c("alpha", "lambda", "lambda_index", "solution_index",
                        paste0("avg_", measure), paste0("sd_", measure))]

  colnames(cvres)[5:6] <- c("cvavg", "cvse")

  sel <- lapply(seq_along(x$alpha), \(alpha_index) {
    alpha <- x$alpha[[alpha_index]]
    vapply(seq_along(x$lambda[[alpha_index]]), FUN.VALUE = integer(1L), FUN = \(lambda_index) {
      candidates <- which(abs(cvres$alpha - alpha) < .Machine$double.eps &
                            cvres$solution_index <= max_solutions &
                            cvres$lambda_index == lambda_index)

      candidates[[which.min(cvres$cvavg[candidates])]]
    })
  }) |>
    unlist(FALSE, FALSE)
  x$cvres <- cvres[sel, ]
  x
}
