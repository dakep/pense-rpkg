#' Penalized Elasitc Net S-estimators for Regression
#'
#' @param x the design matrix
#'
#' @return An object of class \code{"pense"} with elements
#' \item{call}{the call that produced this object.}
#' \item{lambda.opt}{the optimal value of the regularization parameter according to CV.}
#' \item{coefficients}{the vector of coefficients for the optimal lambda \code{lambda.opt}.}
#' \item{residuals}{the residuals, that is response minus fitted values for \code{lambda.opt}.}
#' \item{scale}{the estimated scale for \code{lambda.opt}.}
#' \item{lambda.grid}{a 2-column matrix with values of lambda in the first column and the
#'                    tau-scale estimated via CV in the second column.}
#' \item{alpha}{the value of \code{alpha} the estimatr was run with.}
#' @export
pense <- function(x, ...) {
    UseMethod("pense")
}

#' @section Parallelization:
#' With the parameter \code{ncores}, the number of available processor cores can be set.
#' The grid of lambda values is split into \code{ncores} (almost) equally sized chunks.
#' On each core, the cold initial estimate is computed for the lower endpoint of the
#' lambda sequence. The fully iterated PENSE estimate is then used as warm-start
#' for the next lambda value and so on. Therefore, if \code{ncores = 1}, only one
#' cold initial estimate is computed and all subsequent values of lambda take the
#' previous estimate as warm-start. Similarly, if \code{ncores = nlambda}, a
#' cold initial estimate is computed for every lambda on the grid and no warm-starts
#' are necessary.
#'
#' @param y the response vector.
#' @param alpha the elastic net mixing parameter with \eqn{0 \leq \alpha \leq 1}.
#'      \code{alpha = 1} is the lasso penatly, and \code{alpha = 0} the ridge penalty.
#' @param lambda a single value or a grid of values for the regularization parameter lambda.
#'      Assumed to be on the same scale as the data and adjusted for S-estimation.
#'      Defaults to \code{NULL}, which means a grid of lambda values is automatically
#'      generated (see parameter \code{nlambda}).
#'      If given and \code{standardize = TRUE}, the lambda values will be adjusted
#'      accordingly.
#' @param nlambda if \code{lambda} is not given or \code{NULL}, a grid of \code{nlambda} lambda
#'      values is generated based on the data.
#' @param standardize should the data be standardized robustly? Estimates
#'      are returned on the original scale. Defaults to \code{TRUE}.
#' @param cv.k the number of cross-validation segements to use to choose the optimal
#'      lambda from the grid. If only a single value of lambda is given, this argument
#'      is ignored.
#' @param ncores the number of processor cores to use to estimate the optimal value of lambda.
#'      See details for more information.
#' @param control a list of control parameters as returned by \code{\link{pense.control}}.
#'
#' @rdname pense
#' @export
pense.default <- function(x, y, alpha = 0.5,
                          lambda = NULL, nlambda = 100,
                          standardize = TRUE,
                          cv.k = 5,
                          ncores = getOption("mc.cores", 2L),
                          control = pense.control()) {
    ##
    ## First all arguments are being sanity-checked
    ##
    dY <- dim(y)
    yl <- length(y)
    dX <- dim(x)

    if (is.data.frame(x)) {
        x <- data.matrix(x)
    }

    if (!is.matrix(x) || !is.numeric(x)) {
        stop("`x` must be a numeric matrix")
    }

    if (is.null(yl) || (!is.null(dY) && length(dY) != 1L) || !is.numeric(y)) {
        stop("`yl` must be a numeric vector")
    }

    if (dX[1L] != yl) {
        stop("The number of observations in `x` and `y` does not match")
    }

    if (any(!is.finite(x))) {
        stop("`x` must not contain infinite, NA or NaN values")
    }

    if (any(!is.finite(y))) {
        stop("`y` must not contain infinite, NA or NaN values")
    }

    if (length(alpha) != 1L || !is.numeric(alpha) || is.na(alpha) || alpha < 0 || alpha > 1) {
        stop("`alpha` must be single number in the range [0, 1]")
    }

    if (!is.null(lambda)) {
        if (!is.numeric(lambda) || !is.null(dim(lambda)) || anyNA(lambda) || any(lambda < 0)) {
            stop("Supplied `lambda` must be a numeric vector of non-negative values ",
                 "without NAs")
        }
        nlambda <- length(lambda)
    }

    if (length(standardize) != 1L || !is.logical(standardize) || anyNA(standardize)) {
        warning("`standardize` must be a single logical value. Using default (TRUE).")
        standardize <- TRUE
    }

    if (nlambda > 0L) {
        if (length(cv.k) != 1L || !is.numeric(cv.k) || anyNA(cv.k) || any(cv.k < 0)) {
            warning("`cv.k` must be a single numeric value >= 2. Using default (5).")
            cv.k <- 5L
        }
    }

    control <- .check.pense.control(control)


    ## store the call
    cl <- match.call()
    cl[[1L]] <- as.name("pense")


    return(structure(list(
        residuals = NA_real_,
        coefficients = NA_real_,
        lambda.opt = NA_real_,
        lambda.grid = matrix(NA_real_, ncol = 2L, nrow = 0L),
        alpha = NA_real_,
        scale = NA_real_,
        call = cl
    ), class = "pense"))
}

