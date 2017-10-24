#' Elastic Net Estimator for Regression
#'
#' Estimate the elastic net regression coefficients.
#'
#' This solves the minimization problem
#' \deqn{\frac{1}{2 N} RSS + \lambda \left(
#'         \frac{(1 - \alpha)} {2} \| \beta \|_2^2 + \alpha \| \beta \|_1
#'  \right)}{
#'      (1/2N) RSS + \lambda * (
#'          (1 - \alpha) / 2 * L2(\beta)^2 + \alpha * L1(\beta)
#'  )}
#'
#' If weights are supplied, the minimization problem becomes
#' \deqn{\frac{1}{2 N} \sum_{i = 1}^n w_i r_i^2 + \lambda \left(
#'         \frac{(1 - \alpha)} {2} \| \beta \|_2^2 + \alpha \| \beta \|_1
#'   \right)}{
#'      (1/2N) \sum(w * (y - r^2)) + \lambda * (
#'          (1 - \alpha) / 2 * L2(\beta)^2 + \alpha * L1(\beta)
#'   )}
#'
#' @section Algorithms:
#' Currently this function can compute the elastic net estimator using either
#' augmented LARS or the Dual Augmented Lagrangian (DAL) algorithm
#' (Tomioka 2011).
#' Augmented LARS performs LASSO via the LARS algorithm (or OLS if
#' \code{alpha = 0}) on the data matrix augmented with the L2 penalty term.
#' The time complexity of this algorithm increases fast with an increasing
#' number of predictors. The algorithm currently can not leverage a previous or
#' an approximate solution to speed up computations. However, it is always
#' guaranteed to find the solution.
#'
#' DAL is an iterative algorithm directly minimizing the Elastic Net objective.
#' The algorithm can take an approximate solution to the problem to speed
#' up convergence. In the case of very small lambda values and a bad starting
#' point, DAL may not converge to the solution and hence give wrong
#' results. This would be indicated in the returned status code. Time complexity
#' of this algorithm is dominated by the number of observations.
#'
#' DAL is much faster for a small number of observations (< 200) and a large
#' number of predictors, especially if an approximate solution is available.
#'
#' @param x data matrix with predictors
#' @param y response vector
#' @param alpha controls the balance between the L1 and the L2 penalty.
#'      \code{alpha = 0} is the ridge (L2) penalty, \code{alpha = 1} is
#'      the lasso.
#' @param nlambda size of the lambda grid if \code{lambda} is not specified.
#' @param lambda a grid of decreasing lambda values.
#' @param weights an optional vector of weights to be used in the fitting
#'      process. Should be \code{NULL} or a numeric vector. If non-NULL,
#'      weighted EN is used with weights \code{weights}. See also 'Details'.
#' @param xtest data matrix with predictors used for prediction. This is useful
#'      for testing the prediction performance on an independent test set.
#' @param options additional options for the EN algorithm. See
#'      \code{\link{en_options}} for details.
#' @param intercept should an intercept be estimated?
#' @param lambda_min_ratio if the lambda grid should be automatically defined,
#'      the ratio of the smallest to the largest lambda value in the grid. The
#'      default is \code{1e-6} if $n < p$, and
#'      \code{1e-5 * 10^floor(log10(p / n))} otherwise.
#' @param correction should the "corrected" EN estimator be returned.
#'      If \code{TRUE}, as by default, the corrected estimator as defined in
#'      Zou & Hastie (2005) is returned.
#'
#' @return
#'  \item{lambda}{vector of lambda values.}
#'  \item{status}{integer specifying the exit status of the EN algorithm.}
#'  \item{message}{explanation of the exit status.}
#'  \item{coefficients}{matrix of regression coefficients. Each
#'      column corresponds to the estimate for the lambda value at the same
#'      index.}
#'  \item{residuals}{matrix of residuals. Each column corresponds to the
#'      residuals for the lambda value at the same index.}
#'  \item{predictions}{if \code{xtest} was given, matrix of predicted values.
#'      Each column corresponds to the predictions for the lambda value at the
#'      same index.}
#'
#' @importFrom stats weighted.mean
#'
#' @seealso \code{\link{elnet_cv}} for automatic selection of the penalty
#'      parameter based on the cross-validated prediction error.
#'
#' @references
#' Tomioka, R., Suzuki, T. and Sugiyama, M. (2011).
#'     Super-Linear Convergence of Dual Augmented Lagrangian Algorithm
#'     for Sparse Learning.
#'     \emph{Journal of Machine Learning Research}
#'     \bold{12}(May):1537-1586.
#'
#' Zou, H. and Hastie, T. (2005).
#'     Regularization and variable selection via the elastic net.
#'     \emph{Journal of the Royal Statistical Society}.
#'     Series B (Statistical Methodology), \bold{67}(2):301-320.
#'
#' @example examples/elnet-1.R
#'
#' @export
elnet <- function(x, y, alpha, nlambda = 100, lambda, weights, intercept = TRUE,
                  options = en_options_aug_lars(),
                  lambda_min_ratio, xtest, correction = TRUE) {
    y <- drop(y)

    dx <- dim(x)
    dY <- dim(y)
    yl <- length(y)

    weighted <- FALSE
    have_test_set <- FALSE

    if (is.null(yl) || (!is.null(dY) && length(dY) != 1L) || !is.numeric(y)) {
        stop("`yl` must be a numeric vector")
    }

    if (is.null(dx) || length(dx) != 2L || !is.numeric(x) || dx[1L] != yl) {
        stop("`x` must be a numeric matrix with the same number of ",
             "observations as `y`")
    }

    if (!missing(xtest) && !is.null(xtest)) {
        dxtest <- dim(xtest)
        if (is.null(dxtest) || length(dxtest) != 2L || !is.numeric(xtest) ||
            dxtest[2L] != dx[2L]) {
            stop("`xtest` must be a numeric matrix with the same number of ",
                 "predictors as `x`")
        }
        have_test_set <- TRUE
    } else {
        xtest <- NULL
    }

    if (missing(weights)) {
        weights <- NULL
    }

    if (!is.null(weights)) {
        weighted <- TRUE
        weights <- drop(weights)
        weights <- .check_arg(weights, "numeric", range = 0, length = yl,
                              range_test_lower = ">=")
    }

    if (anyNA(x) || anyNA(y) || anyNA(weights)) {
        stop("Missing values are not supported")
    }

    intercept <- isTRUE(intercept)
    correction <- !identical(correction, FALSE)
    alpha <- .check_arg(alpha, "numeric", range = c(0, 1),
                        range_test_lower = ">=", range_test_upper = "<=")

    if (missing(lambda)) {
        if (missing(lambda_min_ratio)) {
            lambda_min_ratio <- NULL
        }
        nlambda <- .check_arg(nlambda, "integer", range = 0)
        lambda <- .build_lambda_grid_cl(x, y, alpha, nlambda, lambda_min_ratio)
    }

    lambda <- .check_arg(lambda, "numeric", range = 0, length = NULL)
    lambda <- sort(lambda, decreasing = TRUE)

    ret_struct <- structure(list(
        status = 0L,
        alpha = alpha,
        intercept = intercept,
        message = "no observations",
        coefficients = matrix(
            NA_real_,
            nrow = dx[2L] + intercept,
            ncol = length(lambda)
        ),
        residuals = matrix(
            NA_real_,
            nrow = 0L,
            ncol = length(lambda)
        ),
        lambda = lambda,
        call = match.call(expand.dots = FALSE),
        predictions = NULL
    ), class = "elnetfit")

    # Check the size of x and y
    if (identical(dx[1L], 0L)) {
        # no observations
        return(ret_struct)
    } else if (identical(dx[2L], 0L)) {
        # no predictors given
        if (intercept) {
            my <- ifelse(weighted, weighted.mean(y, weights), mean(y))

            ret_struct$coefficients <- matrix(
                my,
                nrow = 1L,
                ncol = length(lambda)
            )
            ret_struct$residuals <- matrix(
                y - my,
                nrow = dx[1L],
                ncol = length(lambda)
            )
            if (have_test_set) {
                ret_struct$predictions <- matrix(
                    my,
                    nrow = dxtest[1L],
                    ncol = length(lambda)
                )
            }
        } else {
            ret_struct$coefficients <- matrix(
                NA_real_,
                nrow = 0L,
                ncol = length(lambda)
            )
            ret_struct$residuals <- matrix(
                y,
                nrow = dx[1L],
                ncol = length(lambda)
            )
        }
        return(ret_struct)
    }

    elnetres <- if (weighted) {
        .elnet.wfit(
            x = x,
            y = y,
            weights = weights,
            alpha = alpha,
            lambda = lambda,
            intercept = intercept,
            options = options,
            xtest = xtest,
            correction = correction
        )
    } else {
        .elnet.fit(
            x = x,
            y = y,
            alpha = alpha,
            lambda = lambda,
            intercept = intercept,
            options = options,
            xtest = xtest,
            correction = correction
        )
    }

    for (ri in names(elnetres)) {
        ret_struct[[ri]] <- elnetres[[ri]]
    }

    lambda_ord <- rev(seq_along(ret_struct$lambda))
    ret_struct$lambda <- ret_struct$lambda[lambda_ord]
    ret_struct$coefficients <- ret_struct$coefficients[, lambda_ord, drop = FALSE]
    ret_struct$residuals <- ret_struct$residuals[, lambda_ord, drop = FALSE]

    if (!is.null(ret_struct$predictions)) {
        ret_struct$predictions <- ret_struct$predictions[, lambda_ord, drop = FALSE]
    }

    return(ret_struct)
}

#' Cross-validate Elastic Net
#'
#' Perform k-fold cross-validation for \code{\link{elnet}}.
#'
#' @param x data matrix with predictors
#' @param y response vector
#' @param alpha controls the balance between the L1 and the L2 penalty.
#'      \code{alpha = 0} is the ridge (L2) penalty, \code{alpha = 1} is
#'      the lasso.
#' @param nlambda size of the lambda grid if \code{lambda} is not specified.
#' @param lambda a grid of decreasing lambda values.
#' @param weights an optional vector of weights to be used in the fitting
#'      process. Should be \code{NULL} or a numeric vector. If non-NULL,
#'      weighted EN is used with weights \code{weights}. See also 'Details'.
#' @param intercept should an intercept be estimated?
#' @param cv_k number of cross validation folds.
#' @param cv_measure function (name) which takes a matrix of residuals and
#'      returns a performance measure for each column.
#' @param ncores,cl the number of processor cores or an actual \code{parallel}
#'      cluster. At most \code{cv_k + 1} nodes will be used.
#' @param lambda_min_ratio if the lambda grid should be automatically defined,
#'      the ratio of the smallest to the largest lambda value in the grid. The
#'      default is \code{1e-6} if $n < p$, and
#'      \code{1e-5 * 10^floor(log10(p / n))} otherwise.
#' @param options additional options for the EN algorithm. See
#'      \code{\link{en_options}} for details.
#' @param ... additional arguments passed on to \code{\link{elnet}}.
#'
#' @return
#'  \item{lambda}{vector of lambda values.}
#'  \item{status}{integer specifying the exit status of the EN algorithm.}
#'  \item{message}{explanation of the exit status.}
#'  \item{coefficients}{matrix of regression coefficients. Each
#'      column corresponds to the estimate for the lambda value at the same
#'      index.}
#'  \item{residuals}{matrix of residuals. Each column corresponds to the
#'      residuals for the lambda value at the same index.}
#'  \item{cvres}{data frame with lambda, average cross-validated performance
#'      and the estimated standard deviation.}
#'
#' @importFrom stats weighted.mean sd
#' @seealso \code{\link{elnet}} to compute only the solution path, without
#'      selecting the optimal penalty parameter using CV.
#'
#' @example examples/elnet_cv-1.R
#'
#' @export
elnet_cv <- function(x, y, alpha, nlambda = 100, lambda, weights,
                     intercept = TRUE, cv_k = 10,
                     cv_measure, ncores = getOption("mc.cores", 1L),
                     cl = NULL, options = en_options_aug_lars(),
                     lambda_min_ratio, ...) {
    y <- drop(y)

    dx <- dim(x)
    dY <- dim(y)
    yl <- length(y)

    if (is.null(yl) || (!is.null(dY) && length(dY) != 1L) || !is.numeric(y)) {
        stop("`yl` must be a numeric vector")
    }

    if (is.null(dx) || length(dx) != 2L || !is.numeric(x) || dx[1L] != yl) {
        stop("`x` must be a numeric matrix with the same number of ",
             "observations as `y`")
    }

    if (anyNA(x) || anyNA(y)) {
        stop("Missing values are not supported")
    }

    call <- match.call(expand.dots = FALSE)

    intercept <- isTRUE(intercept)

    if (missing(weights)) {
        weights <- NULL
    }

    alpha <- .check_arg(alpha, "numeric", range = c(0, 1),
                        range_test_lower = ">=", range_test_upper = "<=")

    if (missing(lambda)) {
        if (missing(lambda_min_ratio)) {
            lambda_min_ratio <- NULL
        }
        nlambda <- .check_arg(nlambda, "integer", range = 0)
        lambda <- .build_lambda_grid_cl(x, y, alpha, nlambda, lambda_min_ratio)
    }

    lambda <- .check_arg(lambda, "numeric", range = 0, length = NULL)
    lambda <- sort(lambda, decreasing = TRUE)

    ## Define return structure
    ret_struct <- structure(list(
        status = 0L,
        alpha = alpha,
        intercept = intercept,
        options = options,
        message = "no observations",
        coefficients = matrix(
            NA_real_,
            nrow = dx[2L] + intercept,
            ncol = length(lambda)
        ),
        residuals = matrix(
            NA_real_,
            nrow = 0L,
            ncol = length(lambda)
        ),
        lambda = lambda,
        cvres = data.frame(lambda = lambda, avg = NA_real_, sd = NA_real_),
        call = call,
        lambda_opt = NA_real_
    ), class = c("cv_elnetfit", "elnetfit"))

    # Check the size of x and y
    if (identical(dx[1L], 0L)) {
        # no observations
        return(ret_struct)
    } else if (identical(dx[2L], 0L)) {
        # no predictors given
        if (intercept) {
            my <- if(!is.null(weights)) {
                weighted.mean(y, weights)
            } else {
                mean(y)
            }

            ret_struct$coefficients <- matrix(
                my,
                nrow = 1L,
                ncol = length(lambda)
            )
            ret_struct$residuals <- matrix(
                y - my,
                nrow = dx[1L],
                ncol = length(lambda)
            )
        } else {
            ret_struct$coefficients <- matrix(
                NA_real_,
                nrow = 0L,
                ncol = length(lambda)
            )
            ret_struct$residuals <- matrix(
                y,
                nrow = dx[1L],
                ncol = length(lambda)
            )
        }
        return(ret_struct)
    }

    # Define CV segments randomly
    cv_k <- .check_arg(cv_k, "integer", range = c(1, yl))
    cv_segments <- c(
        list(full = integer(0L)),
        split(seq_along(y), sample(rep_len(seq_len(cv_k), yl)))
    )

    if (missing(cv_measure)) {
        cv_measure <- function(resids) {
            sqrt(colMeans(resids^2))
        }
    } else {
        cv_measure <- match.fun(cv_measure)
    }

    ## Worker CV function.
    ## Requires the following variables in the environment:
    ## x, y, weights, cv_measure
    elnet_job_cv <- function(cv_ind, alpha, lambda, intercept, options, ...) {
        if (length(cv_ind) > 0L) {
            x_test <- x[cv_ind, , drop = FALSE]
            y_test <- y[cv_ind]
            x_train <- x[-cv_ind, , drop = FALSE]
            y_train <- y[-cv_ind]
            weights_train <- weights[-cv_ind]
        } else {
            x_test <- NULL
            y_test <- NULL
            x_train <- x
            y_train <- y
            weights_train <- weights
        }

        cv_fold_res <- elnet(
            x_train,
            y_train,
            alpha,
            lambda = lambda,
            weights = weights_train,
            intercept = intercept,
            xtest = x_test,
            options = options,
            ...
        )

        if (length(cv_ind) > 0L) {
            # Only return the performance measures
            # `predictions` are ordered from smallest to largest lambda
            return(rev(cv_measure(y_test - cv_fold_res$predictions)))
        } else {
            # Return the full EN result
            return(cv_fold_res)
        }
    }

    # Setup cluster
    cluster <- setupCluster(
        ncores,
        cl,
        export = c("x", "y", "weights", "cv_measure"),
        eval = {
            library(pense)
        }
    )

    # Run all jobs (combination of all CV segments and all lambda-grids)
    tryCatch({
        cv_results <- cluster$lapply(
            cv_segments,
            elnet_job_cv,
            alpha = alpha,
            lambda = lambda,
            intercept = intercept,
            options = options,
            ...
        )
    },
    finally = {
        cluster$stopCluster()
    })

    names(cv_results) <- names(cv_segments)

    # Merge the full result into the return structure
    for (ri in names(cv_results$full)) {
        ret_struct[[ri]] <- cv_results$full[[ri]]
    }

    # Add the CV results
    cv_results_mat <- matrix(
        unlist(cv_results[-1L]),
        ncol = cv_k
    )

    ret_struct$cvres <- data.frame(
        lambda = lambda,
        cvavg = rowMeans(cv_results_mat),
        cvsd = apply(cv_results_mat, 1, sd)
    )

    ret_struct$lambda_opt <- with(ret_struct$cvres, lambda[which.min(cvavg)])
    ret_struct$call <- call
    return(ret_struct)
}


## Internal function to compute a lambda grid for elastic net
#' @importFrom stats cov
.build_lambda_grid_cl <- function(x, y, alpha, nlambda, lambda_min_ratio = NULL) {
    if (is.null(lambda_min_ratio)) {
        dx <- dim(x)
        lambda_min_ratio <- min(1e-6, 1e-5 * 10^floor(log10(dx[2L] / dx[1L])))
    } else {
        lambda_min_ratio <- .check_arg(lambda_min_ratio, "numeric",
                                       range = c(0, 1))
    }
    lambda_max <- max(abs(apply(x, 2, cov, y))) / max(0.01, alpha)

    exp(rev(seq(log(lambda_max), log(lambda_min_ratio * lambda_max),
                length.out = nlambda)))
}


## Internal function to fit an EN linear regression WITHOUT parameter checks!
#' @useDynLib pense, .registration = TRUE
#' @importFrom methods is
#' @importFrom stats sd
#' @importFrom Matrix Matrix Diagonal crossprod colMeans
#' @importClassesFrom Matrix dgCMatrix ddiMatrix
.elnet.fit <- function(x, y, alpha, lambda, intercept = TRUE,
                       options = en_options_aug_lars(),
                       warm_coefs, xtest = NULL, correction = TRUE) {
    y <- drop(y)
    dx <- dim(x)

    ## Add leading column of 1's if necessary
    xtr <- if (isTRUE(sd(x[, 1L]) < sqrt(.Machine$double.eps)) &&
               isTRUE(abs(x[1L, 1L] - 1) < .Machine$double.eps)) {
        # The first column has no variation and the first element is a 1
        # --> we have a column of 1's
        t(x)
    } else {
        .Call(C_augtrans, x)
    }

    if (missing(warm_coefs)) {
        warm_coefs <- NULL
        options$warmStart <- FALSE
    } else {
        if (!is(warm_coefs, "dgCMatrix")) {
            warm_coefs <- Matrix(warm_coefs, sparse = TRUE, ncol = 1L)
        }
        options$warmStart <- TRUE
    }

    if (isTRUE(correction)) {
        options$naive <- FALSE
    } else {
        options$naive <- TRUE
    }

    elnetres <- .Call(
        C_elnet_sp,
        xtr,
        y,
        warm_coefs,
        alpha,
        lambda,
        intercept,
        options,
        xtest
    )

    if (!identical(elnetres[[1L]], 0L)) {
        warning("Elastic Net algorithm had non-zero return status.")
    }

    return(elnetres)
}

## Internal function to fit an EN linear regression WITHOUT parameter checks!
#' @useDynLib pense, .registration = TRUE
#' @importFrom methods is
#' @importFrom stats sd
#' @importFrom Matrix Matrix Diagonal crossprod colSums
#' @importClassesFrom Matrix dgCMatrix ddiMatrix
.elnet.wfit <- function(x, y, weights, alpha, lambda, intercept = TRUE,
                        options = en_options_aug_lars(),
                        warm_coefs, xtest = NULL, correction = TRUE) {
    y <- drop(y)
    dx <- dim(x)

    ## Add leading column of 1's if necessary
    xtr <- if (isTRUE(sd(x[, 1L]) < sqrt(.Machine$double.eps)) &&
               isTRUE(abs(x[1L, 1L] - 1) < .Machine$double.eps)) {
        # The first column has no variation and the first element is a 1
        # --> we have a column of 1's
        t(x)
    } else {
        .Call(C_augtrans, x)
    }

    if (missing(warm_coefs)) {
        warm_coefs <- NULL
        options$warmStart <- FALSE
    } else {
        if (!is(warm_coefs, "dgCMatrix")) {
            warm_coefs <- Matrix(warm_coefs, sparse = TRUE, ncol = 1L)
        }
        options$warmStart <- TRUE
    }

    if (isTRUE(correction)) {
        options$naive <- FALSE
    } else {
        options$naive <- TRUE
    }

    elnetres <- .Call(
        C_elnet_weighted_sp,
        xtr,
        y,
        weights,
        warm_coefs,
        alpha,
        lambda,
        intercept,
        options,
        xtest
    )

    if (!identical(elnetres[[1L]], 0L)) {
        warning("Elastic Net algorithm had non-zero return status.")
    }

    elnetres$weights <- weights

    return(elnetres)
}
