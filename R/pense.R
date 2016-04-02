#' Penalized Elasitc Net S-estimators for Regression
#'
#' @param X the design matrix
#'
#' @return An object of class \code{"pense"} with elements
#' \item{call}{the call that produced this object.}
#' \item{lambda.opt}{the optimal value of the regularization parameter according to CV.}
#' \item{coefficients}{the vector of coefficients for the optimal lambda \code{lambda.opt}.}
#' \item{residuals}{the residuals, that is response minus fitted values for \code{lambda.opt}.}
#' \item{scale}{the estimated scale for \code{lambda.opt}.}
#' \item{lambda.grid}{a 2-column matrix with values of lambda in the first column and the
#'                    tau-scale estimated via CV in the second column.}
#' \item{alpha,standardize,control}{the given arguments.}
#' \item{standardize}{}
#' @export
pense <- function(X, ...) {
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
#' @param alpha the elastic net mixing parameter with \eqn{0 \le \alpha \le 1}.
#'      \code{alpha = 1} is the lasso penatly, and \code{alpha = 0} the ridge penalty.
#' @param nlambda if \code{lambda} is not given or \code{NULL} (default),
#'      a grid of \code{nlambda} lambda values is generated based on the data.
#' @param lambda a single value or a grid of values for the regularization parameter lambda.
#'      Assumed to be on the same scale as the data and adjusted for S-estimation.
#'      Defaults to \code{NULL}, which means a grid of lambda values is automatically
#'      generated (see parameter \code{nlambda}).
#'      If given and \code{standardize = TRUE}, the lambda values will be adjusted
#'      accordingly.
#' @param standardize should the data be standardized robustly? Estimates
#'      are returned on the original scale. Defaults to \code{TRUE}.
#' @param cv.k number of cross-validation segements to use to choose the optimal
#'      lambda from the grid. If only a single value of lambda is given, cross-validation
#'      can still done to estimate the prediction performance at this particular lambda.
#' @param warm.reset how often should the warm-start be reset to a cold initial estimate?
#'      If \code{NULL}, no warm starts will be used.
#' @param ncores,cl the number of processor cores or an actual parallel cluster to use to
#'      estimate the optimal value of lambda. See details for more information.
#' @param control a list of control parameters as returned by \code{\link{pense.control}}.
#'
#' @rdname pense
#' @importFrom stats mad median
#' @export
pense.default <- function(X, y, alpha = 0.5,
                          nlambda = 100, lambda = NULL,
                          standardize = TRUE,
                          cv.k = 5, warm.reset = 10,
                          ncores = getOption("mc.cores", 2L), cl = NULL,
                          control = pense.control()) {
    ##
    ## First all arguments are being sanity-checked
    ##
    y <- drop(y)

    dY <- dim(y)
    yl <- length(y)
    dX <- dim(X)

    if (is.data.frame(X)) {
        X <- data.matrix(X)
    }

    if (!is.matrix(X) || !is.numeric(X)) {
        stop("`X` must be a numeric matrix")
    }

    if (is.null(yl) || (!is.null(dY) && length(dY) != 1L) || !is.numeric(y)) {
        stop("`yl` must be a numeric vector")
    }

    if (dX[1L] != yl) {
        stop("The number of observations in `X` and `y` does not match")
    }

    if (any(!is.finite(X))) {
        stop("`X` must not contain infinite, NA, or NaN values")
    }

    if (any(!is.finite(y))) {
        stop("`y` must not contain infinite, NA or NaN values")
    }

    if (length(alpha) != 1L || !is.numeric(alpha) || is.na(alpha) || alpha < 0 || alpha > 1) {
        stop("`alpha` must be single number in the range [0, 1]")
    }

    if (alpha == 0) {
        alpha <- 1e-3
    }

    if (!is.null(lambda)) {
        if (!is.numeric(lambda) || !is.null(dim(lambda)) || anyNA(lambda) || any(lambda < 0)) {
            stop("Supplied `lambda` must be a numeric vector of non-negative values ",
                 "without NAs")
        }
        lambda <- sort(lambda)
        nlambda <- length(lambda)
    } else if (length(nlambda) != 1L || !is.numeric(nlambda) || anyNA(nlambda) ||
               any(nlambda < 1)) {
        stop("`nlambda` must be a single positive numeric value.")
    }
    nlambda <- as.integer(nlambda)

    if (length(standardize) != 1L || !is.logical(standardize) || anyNA(standardize)) {
        warning("`standardize` must be a single logical value. Using default (TRUE).")
        standardize <- TRUE
    }

    if (length(cv.k) != 1L || !is.numeric(cv.k) || anyNA(cv.k) || any(cv.k < 0)) {
        warning("`cv.k` must be a single numeric value >= 2. Using default (5).")
        cv.k <- 5L
    }
    cv.k <- as.integer(cv.k)

    if (is.null(warm.reset)) {
        warm.reset <- nlambda
    }

    if (length(warm.reset) != 1L || !is.numeric(warm.reset) || anyNA(warm.reset) ||
        any(warm.reset < 1)) {
        stop("`warm.reset` must be a single numeric value greater than 0.")
    }
    warm.reset <- as.integer(warm.reset)

    control <- .check.pense.control(control)

    ## store the call
    cl <- match.call()
    cl[[1L]] <- as.name("pense")


    ## Automatically select the PSC method
    if (control$init.psc.method == "auto") {
        ntrain <- ifelse(cv.k > 1,
                         dX[1L] - dX[1L] %/% cv.k - as.integer((dX[1L] %% cv.k) > 0),
                         dX[1L])

        if (dX[2L] >= ntrain) {
            control$init.psc.method <- "Mn"
        } else {
            control$init.psc.method <- "rr"
        }
    }

    scale.x <- 1
    mux <- 0
    muy <- 0

    ## Standardize data -- only needed to compute the grid of lambda values
    if (standardize == TRUE) {
        scale.x <- apply(X, 2, mad)
        mux <- apply(X, 2, median)
        muy <- median(y)

        Xs <- scale(X, center = mux, scale = scale.x)
        yc <- y - muy

        if (!is.null(lambda)) {
            lambda <- lambda / max(scale.x)
        }
    } else {
        Xs <- X
        yc <- y
    }

    ## Generate grid of lambda-values
    if (is.null(lambda)) {
        lambda <- lambda.grid(Xs, yc, nlambda, control, standardize = FALSE)
    }

    ## Create CV segments
    cv.segments <- if(cv.k > 1L) {
        split(seq_len(dX[1L]), sample(rep_len(seq_len(cv.k), dX[1L])))
    } else {
        list(integer(0))
    }

    subgrid.lengths <- nlambda %/% warm.reset
    overlength <- nlambda %% warm.reset
    subgrid.lengths <- c(rep.int(subgrid.lengths, warm.reset - overlength),
                         rep.int(subgrid.lengths + 1L, overlength))

    lambda.subgrids <- split(lambda, rep.int(seq_len(warm.reset), subgrid.lengths))

    jobgrid <- lapply(cv.segments, function(cvs) {
        lapply(lambda.subgrids, function(lvs) {
            list(segment = cvs,
                 lambda = lvs)
        })
    })

    jobgrid <- unlist(jobgrid, recursive = FALSE, use.names = FALSE)

    ## Setup cluster
    cluster <- setupCluster(ncores, cl,
                            export = c("X", "y"),
                            eval = {
                                library(penseinit)
                            })

    ## Run all jobs (combination of all CV segments and all lambda-grids)
    dojobcv <- function(job, alpha, standardize, control) {

        if (length(job$segment) == 0L) {
            X.train <- X
            y.train <- y
        } else {
            X.train <- X[-job$segment, , drop = FALSE]
            y.train <- y[-job$segment]

        }


        est.all <- penseinit:::pense.coldwarm(X.train, y.train,
                                              alpha, job$lambda, standardize, control)

        residuals <- NULL

        if (length(job$segment) == 0L) {
            residuals <- vapply(est.all, function(est) {
                est$resid
            }, FUN.VALUE = y, USE.NAMES = FALSE)
        } else {
            X.test <- X[job$segment, , drop = FALSE]
            y.test <- y[job$segment]
            residuals <- vapply(est.all, function(est) {
                y.test - est$intercept - X.test %*% est$beta
            }, FUN.VALUE = numeric(length(job$segment)), USE.NAMES = FALSE)
        }

        return(residuals)
    }

    tryCatch({
        prediction.errors <- cluster$lapply(jobgrid, dojobcv, alpha = alpha,
                                            standardize = standardize, control = control)
    },
    finally = {
        cluster$stopCluster()
    })

    ## Collect all prediction errors for each lambda sub-grid and determine the optimal lambda
    cv.performance <- unlist(lapply(split(prediction.errors, rep.int(seq_len(warm.reset), cv.k)),
                                    function(preds) {
                                        preds <- do.call(rbind, preds)
                                        apply(preds, 2, control$cv.objective)
                                    }), use.names = FALSE)

    lambda.grid <- cbind(lambda, cv.performance)
    lambda.opt <- lambda[which.min(cv.performance)]

    ## Fit PENSE with optimal lambda using a cold start
    opt.est <- pense.coldwarm(Xs, yc, alpha, lambda.opt, standardize = FALSE, control)[[1L]]


    ## Un-standardize the coefficients
    if (standardize == TRUE) {
        opt.est$beta <- opt.est$beta / scale.x
        opt.est$intercept <- opt.est$intercept + muy - drop(opt.est$beta %*% mux)
    }

    return(structure(list(
        residuals = opt.est$resid,
        coefficients = nameCoefVec(c(opt.est$intercept, opt.est$beta), X),
        objective = opt.est$objF,
        lambda.opt = lambda.opt,
        lambda.grid = lambda.grid,
        scale = opt.est$scale,
        alpha = alpha,
        standardize = standardize,
        control = control,
        call = cl
    ), class = "pense"))
}

