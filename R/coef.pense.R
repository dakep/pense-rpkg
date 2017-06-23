#' Extract Model Coefficients
#'
#' @param object object of type \code{pense} to extract coefficients from.
#' @param lambda the value of the penalty parameter. Default is to use the
#'      optimal lambda (\code{object$lambda_opt}).
#' @param exact if the lambda is not part of the lambda grid, should the
#'      estimates be obtained by linear interpolation between the nearest
#'      lambda values (default) or computed exactly.
#' @param sparse return a sparse vector or a dense (base R) numeric vector
#' @param ... currently not used.
#' @return A numeric vector of size \eqn{p + 1}.
#' @export
coef.pense <- function(object, lambda, exact = FALSE, sparse = FALSE, ...) {
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

    x <- data.matrix(eval(object$call$X))
    y <- drop(eval(object$call$y))

    std_data <- pense:::standardize_data(x, y, object$standardize)
    xs <- std_data$xs
    yc <- std_data$yc

    init_int <- object$coefficients[1L, closests_lambda_ind]
    init_beta <- object$coefficients[-1L, closests_lambda_ind, drop = FALSE]

    if (isTRUE(object$standardize)) {
        init_int <- init_int - std_data$muy +
            as.numeric(std_data$mux %*% init_beta)
        init_beta <- init_beta * std_data$scale_x
    }

    estimate <- pen_s_reg(
        xs,
        yc,
        alpha = object$alpha,
        lambda = lambda / max(std_data$scale_x),
        init_int = init_int,
        init_coef = init_beta,
        options = object$pense_options,
        en_options = object$en_options
    )

    if (isTRUE(object$standardize)) {
        estimate$beta <- estimate$beta / std_data$scale_x
        estimate$intercept <- estimate$intercept + std_data$muy -
            as.numeric(std_data$mux %*% estimate$beta)
    }

    if (sparse) {
        coef <- rbind(estimate$intercept, estimate$beta)
        rownames(coef) <- rownames(object$coefficients)
    } else {
        coef <- setNames(
            c(estimate$intercept, as.numeric(estimate$beta)),
            rownames(object$coefficients)
        )
    }
    return(coef)

}

#' Extract Model Coefficients
#'
#' @param object object of type \code{pensem} to extract coefficients from.
#' @param sparse return a sparse vector or a dense (base R) numeric vector
#' @param ... currently not used.
#' @return A numeric vector of size \eqn{p + 1}.
#' @export
coef.pensem <- function(object, sparse = FALSE, ...) {
    if (isTRUE(sparse)) {
        return(object$coefficients)
    } else {
        return(as.numeric(object$coefficients))
    }
}

#' Extract Residuals from a Fitted Penalized Elastic-Net S-estimator
#'
#' @param object an object of type \code{pense} to extract the residuals from.
#' @param lambda the value of the penalty parameter. Default is to use the
#'      optimal lambda \code{object$lambda_opt}.
#' @param exact if the lambda is not part of the lambda grid, should the
#'      estimates be obtained by linear interpolation between the nearest
#'      lambda values (default) or computed exactly.
#' @param ... currently ignored.
#' @return a numeric vector of residuals for the given lambda.
#' @export
residuals.pense <- function(object, lambda, exact = FALSE, ...) {
    exact <- isTRUE(exact)

    if (missing(lambda) || is.null(lambda)) {
        lambda <- object$lambda_opt
        exact <- FALSE
    }

    lambda_diff_abs <- abs(object$lambda - lambda)
    lambda_match <- which(lambda_diff_abs < .Machine$double.eps)

    if (length(lambda_match) > 0L) {
        return(object$residuals[ , lambda_match[1L]])
    }

    coefs <- coef.pense(object, lambda = lambda, exact = exact, sparse = TRUE)

    x <- data.matrix(eval(object$call$X))
    y <- drop(eval(object$call$y))

    return(drop(y - coefs[1L] - x %*% coefs[-1L]))
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
#' @return A numeric vector of size \eqn{p + 1}.
#' @export
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

    x <- data.matrix(eval(object$call$X))
    y <- drop(eval(object$call$y))

    init_int <- object$coefficients[1L, closests_lambda_ind]
    init_beta <- object$coefficients[-1L, closests_lambda_ind, drop = FALSE]

    est <- if (is.null(object$weights)) {
        .elnet.wfit(
            X = x,
            y = y,
            weights = object$weights,
            alpha = object$alpha,
            lambda = lambda,
            intercept = object$intercept,
            addLeading1s = TRUE,
            options = object$options,
            warm_coefs = object$coefficients[ , closests_lambda_ind, drop = FALSE]
        )
    } else {
        .elnet.fit(
            X = x,
            y = y,
            alpha = object$alpha,
            lambda = lambda,
            intercept = object$intercept,
            addLeading1s = TRUE,
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

