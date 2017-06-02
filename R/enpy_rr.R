#' PY (Pena-Yohai) initial estimates for EN S-Estimators
#'
#' Computes the PY initial estimates for EN S-Estimators with
#' principal sensitivity components approximated by the ridge regression
#' solution.
#'
#' @param X data matrix with predictors -- a leading column of 1's will be
#'      added!
#' @param y response vector.
#' @param lambda,alpha The EN penalty parameters (NOT adjusted for the number
#'      of observations in \code{X}).
#' @param options additional options for the initial estimator. See
#'      \code{\link{initest_options}} for details.
#' @param en_options additional options for the EN algorithm. See
#'      \code{\link{en_options}} for details.
#'
#' @return \item{coeff}{A numeric matrix with one initial coefficient per column}
#'         \item{objF}{A vector of values of the objective function for the respective coefficient}
#'
#' @useDynLib pense C_enpy_rr
#' @importFrom Rcpp evalCpp
enpy.rr <- function(X, y, alpha, lambda, options, en_options) {
    dX <- dim(X)

    Xtr <- .Call(C_augtrans, X)
    dX[2L] <- dX[2L] + 1L

    if (alpha >= 1 - .Machine$double.eps) {
        if (dX[2L] >= dX[1L]) {
            stop("`enpy_rr` can not be used for data with more variables than ",
                 "observations if `alpha` is 1.")
        } else if (dX[2L] >= ceiling(usableProp * dX[1L])) {
            stop("With the specified proportion of observations to remove, ",
                 "the number of observations will be smaller than the number ",
                 "of variables.\n",
                 "In this case `enpy_rr` can not be used when `alpha` is 1")
        }
    }

    ies <- .Call(C_enpy_rr, Xtr, y, options, en_options)

    return(list(
        coeff = matrix(ies[[1L]], nrow = dX[2L]),
        objF = ies[[2L]]
    ))
}

