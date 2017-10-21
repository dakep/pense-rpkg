## Internal function to calculate PENSE for a given initial estimate (init.coef)
## in C++
##
#' @useDynLib pense, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom methods is
#' @importFrom Matrix Matrix
#' @importClassesFrom Matrix dgCMatrix
pen_s_reg <- function(x, y, alpha, lambda, init_int, init_coef, warn = TRUE,
                      options, en_options) {
    dx <- dim(x)

    xtr <- .Call(C_augtrans, x)
    dx[2L] <- dx[2L] + 1L

    if (!is(init_coef, "dgCMatrix")) {
        init_coef <- Matrix(init_coef, ncol = 1L, sparse = TRUE)
    }

    ret <- .Call(
        C_pen_s_reg_sp,
        xtr,
        y,
        init_int,
        init_coef,
        alpha,
        lambda,
        options,
        en_options
    )

    ##
    ## Check if the S-step converged.
    ## Be extra careful with the comparison as rel_change may be NaN or NA
    ##
    if (!isTRUE(ret$rel_change < options$eps) && isTRUE(warn)) {
        warning(sprintf("PENSE did not converge for lambda = %.6f", lambda))
    }

    return(ret)
}
