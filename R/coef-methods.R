#' Extract Model Coefficients
#'
#' @param object a PENSE or PENSEM estimate to extract coefficients from.
#' @param lambda the value of the penalty parameter. Default is to use the
#'      optimal lambda (\code{object$lambda_opt}).
#' @param exact if \code{lambda} is not part of the lambda grid, should the
#'      estimates be obtained by linear interpolation between the nearest
#'      lambda values (default) or computed exactly.
#' @param sparse return a sparse vector or a dense (base R) numeric vector
#' @param correction should a correction factor be applied to the EN estimate?
#'       See \code{\link{elnet}} for details on the applied correction.
#' @param ... currently not used.
#' @return if \code{sparse = FALSE} a numeric vector of size \eqn{p + 1}.
#'      Otherwise a sparse matrix with one column and \eqn{p + 1} rows.
#'
#' @importFrom Matrix drop
#' @importFrom stats weighted.mean
#' @importClassesFrom Matrix dgCMatrix
#'
#' @example examples/pense-methods.R
#'
#' @export
coef.pense <- function(object, lambda, exact = FALSE, sparse = FALSE, correction = TRUE, ...) {
    exact <- isTRUE(exact)
    sparse <- isTRUE(sparse)
    correction <- as.integer(correction)
    correction <- correction * ((object$alpha < 1) && !is.null(object$adjusted))

    if (missing(lambda) || is.null(lambda)) {
        lambda <- object$lambda_opt
        exact <- FALSE
    }

    lambda <- .check_arg(lambda, "numeric", range = 0)
    max_lambda_ind <- which.max(object$lambda)
    min_lambda_ind <- which.min(object$lambda)

    if (lambda > object$lambda[max_lambda_ind] &&
        all(abs(object$coefficients[-1L, max_lambda_ind]) < .Machine$double.eps)) {
        # requested lambda is larger than the largest lambda in the grid AND
        # the largest lambda in the grid gave an all-zero vector -->
        # we can return the result of the largest lambda
        return(object$coefficients[, max_lambda_ind, drop = !sparse])
    }

    lambda_diff <- object$lambda - lambda
    lambda_diff_abs <- abs(lambda_diff)

    if (any(lambda_diff_abs < .Machine$double.eps)) {
        # the requested lambda value is alrady in the grid
        sel_lambda_index <- which.min(lambda_diff_abs)
        ret_coef <- object$coefficients[, sel_lambda_index, drop = FALSE]
        if (correction > 0L) {
            ret_coef[1L, ] <- object$adjusted$intercept[sel_lambda_index]
            adj_fact <- object$adjusted$factor[sel_lambda_index]
        }
    } else if (!exact) {
        # the requested lambda value is not part of the grid
        # take the average of the closest two solutions
        lambda_diff_neg_ind <- which(lambda_diff < 0)
        lambda_diff_pos_ind <- which(lambda_diff > 0)
        if (length(lambda_diff_neg_ind) > 0L &&
            length(lambda_diff_pos_ind) > 0L) {

            interp_lambda_left <- lambda_diff_neg_ind[which.max(
                lambda_diff[lambda_diff_neg_ind]
            )]
            interp_lambda_right <- lambda_diff_pos_ind[which.min(
                lambda_diff[lambda_diff_pos_ind]
            )]

            interp_w <- 1 / lambda_diff_abs[c(interp_lambda_left, interp_lambda_right)]
            interp_c <- (interp_w[1L] *
                             object$coefficients[, interp_lambda_left, drop = FALSE] +
                             interp_w[2L] *
                             object$coefficients[, interp_lambda_right, drop = FALSE]) /
                sum(interp_w)

            if (correction > 0L) {
                interp_c[1L, ] <- weighted.mean(
                    object$adjusted$intercept[c(interp_lambda_left, interp_lambda_right)],
                    interp_w
                )
                adj_fact <- weighted.mean(
                    object$adjusted$factor[c(interp_lambda_left, interp_lambda_right)],
                    interp_w
                )
            }

            ret_coef <- interp_c
        }
    } else {
        # else:
        # there is no left or right interpolation point available --> use
        # exact solution
        x <- if (!is.null(object$call$x_train)) {
            data.matrix(eval(object$call$x_train))
        # } else if (!is.null(object$sest$call$x)) {
        #     data.matrix(eval(object$sest$call$x))
        } else {
            data.matrix(eval(object$call$x))
        }

        y <- if (!is.null(object$call$y_train)) {
            data.matrix(eval(object$call$y_train))
        # } else if (!is.null(object$sest$call$y)) {
        #     data.matrix(eval(object$sest$call$y))
        } else {
            data.matrix(eval(object$call$y))
        }

        std_data <- standardize_data(x, y, object$standardize)
        xs <- std_data$xs
        yc <- std_data$yc

        adj_fact <- .en_correction_factor(correction, object$alpha, lambda)

        ## is it the S- or MM-estimate
        if ("pensem" %in% class(object)) {
            ## compute the MM-estimate starting from the "best" S-estimate
            sel_lambda_ind <- with(object$sest, which(lambda == lambda_opt))
            init_int <- object$sest$coefficients[1L, sel_lambda_ind]
            init_beta <- object$sest$coefficients[-1L, sel_lambda_ind, drop = FALSE]

            if (isTRUE(object$standardize)) {
                std_coefs <- std_data$standardize_coefs(list(
                    intercept = init_int,
                    beta = init_beta
                ))
                init_int <- std_coefs$intercept
                init_beta <- std_coefs$beta
            }

            estimate <- pensemstep(
                xs,
                yc,
                init_scale = object$init_scale,
                init_int = init_int,
                init_coef = init_beta,
                alpha = object$alpha,
                lambda = lambda,
                options = object$options,
                en_options = object$en_options
            )
        } else {
            ## compute the S-estimate starting from the closest lambda
            closests_lambda_ind <- which.min(lambda_diff_abs)
            init_int <- object$coefficients[1L, closests_lambda_ind]
            init_beta <- object$coefficients[-1L, closests_lambda_ind, drop = FALSE]

            if (isTRUE(object$standardize)) {
                std_coefs <- std_data$standardize_coefs(list(
                    intercept = init_int,
                    beta = init_beta
                ))
                init_int <- std_coefs$intercept
                init_beta <- std_coefs$beta
            }

            estimate <- pen_s_reg(
                xs,
                yc,
                alpha = object$alpha,
                lambda = lambda,
                init_int = init_int,
                init_coef = init_beta,
                options = object$pense_options,
                en_options = object$en_options
            )
        }

        if (isTRUE(object$standardize)) {
            estimate <- std_data$unstandardize_coefs(estimate)
        }

        ret_coef <- rbind(estimate$intercept, estimate$beta)
        rownames(ret_coef) <- rownames(object$coefficients)

        if (correction > 0L) {
            ret_coef[-1L, ] <- ret_coef[-1L, , drop = FALSE] * adj_fact
            residuals <- drop(y - x %*% ret_coef[-1L, , drop = FALSE])
            ret_coef[1L, ] <- weighted.mean(residuals, estimate$weights)
            correction <- 0L
        }
    }

    if (correction > 0L) {
        # The intercept was already adjusted before!
        ret_coef[-1L, ] <- ret_coef[-1L, , drop = FALSE] * adj_fact
    }

    if (!sparse) {
        ret_coef <- as.numeric(ret_coef)
    }

    return(ret_coef)

}

#' Extract Model Coefficients
#'
#' @param object object of type \code{elnetfit} to extract coefficients from.
#' @param lambda the value of the penalty parameter. Default is to use the
#'      optimal lambda (\code{object$lambda_opt}).
#' @param exact if the lambda is not part of the lambda grid, should the
#'      estimates be obtained by linear interpolation between the nearest
#'      lambda values (default) or computed exactly.
#' @param sparse return a sparse vector or a dense (base R) numeric vector
#' @param ... currently not used.
#' @return if \code{sparse = FALSE} a numeric vector of size \eqn{p + 1}.
#'      Otherwise a sparse matrix with one column and \eqn{p + 1} rows.
#'
#' @example examples/elnet_cv-methods.R
#'
#' @export
#' @importFrom stats setNames
coef.elnetfit <- function(object, lambda, exact = FALSE, sparse = FALSE, ...) {
    exact <- isTRUE(exact)
    sparse <- isTRUE(sparse)

    if (missing(lambda) || is.null(lambda)) {
        lambda <- object$lambda_opt
        exact <- FALSE
    }

    lambda <- .check_arg(lambda, "numeric", range = 0)
    max_lambda_ind <- which.max(object$lambda)
    min_lambda_ind <- which.min(object$lambda)

    if (lambda > object$lambda[max_lambda_ind] &&
        all(abs(object$coefficients[-1L, max_lambda_ind]) < .Machine$double.eps)) {
        # requested lambda is larger than the largest lambda in the grid AND
        # the largest lambda in the grid gave an all-zero vector -->
        # we can return the result of the largest lambda
        return(object$coefficients[, max_lambda_ind, drop = !sparse])
    }

    lambda_diff <- object$lambda - lambda
    lambda_diff_abs <- abs(lambda_diff)

    if (any(lambda_diff_abs < .Machine$double.eps)) {
        # the requested lambda value is alrady in the grid
        return(object$coefficients[, which.min(lambda_diff_abs), drop = !sparse])
    }

    # the requested lambda value is not part of the grid
    if (!exact) {
        # take the average of the closest two solutions
        lambda_diff_neg_ind <- which(lambda_diff < 0)
        lambda_diff_pos_ind <- which(lambda_diff > 0)
        if (length(lambda_diff_neg_ind) > 0L &&
            length(lambda_diff_pos_ind) > 0L) {

            interp_lambda_left <- lambda_diff_neg_ind[which.max(
                lambda_diff[lambda_diff_neg_ind]
            )]
            interp_lambda_right <- lambda_diff_pos_ind[which.min(
                lambda_diff[lambda_diff_pos_ind]
            )]

            interp_w <- 1 / lambda_diff_abs[c(interp_lambda_left, interp_lambda_right)]
            interp_c <- (interp_w[1L] *
                             object$coefficients[, interp_lambda_left, drop = !sparse] +
                             interp_w[2L] *
                             object$coefficients[, interp_lambda_right, drop = !sparse]) /
                sum(interp_w)

            if (!sparse) {
                interp_c <- as.numeric(interp_c)
            }
            return(interp_c)
        }
        # else:
        # there is no left or right interpolation point available --> use
        # exact solution
    }

    closests_lambda_ind <- which.min(lambda_diff_abs)

    x <- data.matrix(eval(object$call$x))
    y <- drop(eval(object$call$y))

    init_int <- object$coefficients[1L, closests_lambda_ind]
    init_beta <- object$coefficients[-1L, closests_lambda_ind, drop = FALSE]

    est <- if (!is.null(object$weights)) {
        .elnet.wfit(
            x = x,
            y = y,
            weights = object$weights,
            alpha = object$alpha,
            lambda = lambda,
            intercept = object$intercept,
            options = object$options,
            warm_coefs = object$coefficients[ , closests_lambda_ind, drop = FALSE]
        )
    } else {
        .elnet.fit(
            x = x,
            y = y,
            alpha = object$alpha,
            lambda = lambda,
            intercept = object$intercept,
            options = object$options,
            warm_coefs = object$coefficients[ , closests_lambda_ind, drop = FALSE]
        )
    }

    if (sparse) {
        coef <- est$coefficients
        rownames(coef) <- rownames(object$coefficients)
    } else {
        coef <- setNames(
            as.numeric(est$coefficients),
            rownames(object$coefficients)
        )
    }
    return(coef)
}

