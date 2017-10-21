#' Additional Options for the Initial Estimator
#'
#' Specify additional options for the initial estimator based on the Pena-Yohai
#' estimator.
#'
#' Two different methods to calculate the sensitivity components are
#' implemented:
#' \describe{
#'      \item{\code{"rr"}}{Approximate the PSCs by using the residuals from the
#'          elastic net fit and the hat matrix from the ridge regression.
#'          This method only works if \code{alpha} < 1 or
#'          \code{ncol(x)} < \code{nrow(x)}.}
#'      \item{\code{"exact"}}{Calculate the PSCs from the difference between the
#'          residuals and leave-one-out residuals from elastic net.}
#' }
#'
#' @param keep_solutions how many initial estimates should be kept to perform
#'      full PENSE iterations?
#' @param psc_method The method to use for computing the principal sensitivity
#'      components. See details for the possible choices.
#' @param maxit maximum number of refinement iterations.
#' @param maxit_pense_refinement maximum number of PENSE iterations to refine
#'      initial estimator.
#' @param eps numeric tolerance for convergence.
#' @param psc_keep proportion of observations to keep based on the PSC scores.
#' @param resid_keep_method How to clean the data based on large residuals.
#'      If \code{"proportion"}, observations with the smallest
#'      \code{resid_keep_prop} residuals will be retained.
#'      If \code{"threshold"}, all observations with scaled residuals smaller
#'      than the threshold \code{resid_keep_thresh} will be retained.
#' @param resid_keep_prop,resid_keep_thresh proportion or threshold for
#'      observations to keep based on their residual.
#' @param delta desired breakdown point of the S-estimator.
#' @param mscale_eps,mscale_maxit maximum number of iterations and numeric
#'      tolerance for the M-scale.
#' @param cc tuning constant for the S-estimator. Default is to chosen based
#'      on the breakdown point \code{delta}. Should never have to be changed.
#' @return a checked options list.
#' @export
#' @family specifying additional options
#'
#' @references Pena, D., & Yohai, V.. (1999). A Fast Procedure for Outlier
#' Diagnostics in Large Regression Problems. \emph{Journal of the American
#' Statistical Association}, 94(446), 434-445.
#' \url{http://doi.org/10.2307/2670164}
initest_options <- function (
    keep_solutions = 5,
    psc_method = c("rr", "exact"),
    maxit = 10,
    maxit_pense_refinement = 5,
    eps = 1e-6,
    psc_keep = 0.5,
    resid_keep_method = c("proportion", "threshold"),
    resid_keep_prop = 0.6,
    resid_keep_thresh = 2,
    delta = 0.25,
    mscale_eps = 1e-8,
    mscale_maxit = 200,
    cc
) {
    psc_method <- match.arg(psc_method)
    resid_keep_method <- match.arg(resid_keep_method)
    resid_keep_method <- as.integer(pmatch(resid_keep_method,
                                           c("proportion", "threshold"))) - 1L

    if (missing(cc)) {
        cc <- consistency.rho(delta, 1L)
    }

    return(list(
        keepSolutions = .check_arg(keep_solutions, "integer", range = 0),
        pscMethod = psc_method,
        maxit = .check_arg(maxit, "integer", range = 0),
        maxitPenseRefinement = .check_arg(maxit_pense_refinement, "integer", range = 0),
        eps = .check_arg(eps, "numeric", range = 0),
        keepResidualsMethod = resid_keep_method,
        keepResidualsThreshold = .check_arg(
            resid_keep_thresh,
            "numeric",
            range = 0
        ),
        keepResidualsProportion = .check_arg(
            resid_keep_prop,
            "numeric",
            range = c(0, 1)
        ),
        keepPSCProportion = .check_arg(psc_keep, "numeric", range = c(0, 1)),
        mscaleDelta = .check_arg(
            delta,
            "numeric",
            range = c(0, 0.5),
            range_test_upper = "<="
        ),
        mscaleCC = .check_arg(cc, "numeric", range = 0),
        mscaleEps = .check_arg(mscale_eps, "numeric", range = 0),
        mscaleMaxit = .check_arg(mscale_maxit, "integer", range = 0)
    ))
}

#' Additional Options for the Penalized EN S-estimator
#'
#'
#' @param delta desired breakdown point of the resulting estimator.
#' @param maxit maximum number of iterations allowed.
#' @param eps numeric tolerance for convergence.
#' @param mscale_eps,mscale_maxit maximum number of iterations and numeric
#'      tolerance for the M-scale.
#' @param cc tuning constant for the S-estimator. Default is to chosen based
#'      on the breakdown point \code{delta}. Should never have to be changed.
#' @param verbosity verbosity of the algorithm.
#' @param en_correction should the corrected EN estimator be used to choose
#'      the optimal lambda with CV.
#'      If \code{TRUE}, as by default, the estimator is "bias corrected".
#' @return a checked options list.
#' @export
#' @family specifying additional options
pense_options <- function (
    delta = 0.25,
    maxit = 1000,
    eps = 1e-6,
    mscale_eps = 1e-8,
    mscale_maxit = 200,
    verbosity = 0,
    cc,
    en_correction = TRUE
) {
    if (missing(cc)) {
        cc <- consistency.rho(delta, 1L)
    }

    return(list(
        bdp = .check_arg(delta, "numeric", range = c(0, 0.5),
                         range_test_upper = "<="),
        maxit = .check_arg(maxit, "integer", range = 0),
        eps = .check_arg(eps, "numeric", range = 0),
        cc = .check_arg(cc, "numeric", range = 0),
        mscaleEps = .check_arg(mscale_eps, "numeric", range = 0),
        mscaleMaxit = .check_arg(mscale_maxit, "integer", range = 0),
        naiveEn = !.check_arg(en_correction, "logical"),
        verbosity = .check_arg(verbosity, "integer", range = 0,
                         range_test_lower = ">=")
    ))
}

#' Additional Options for the Penalized EN MM-estimator
#'
#' @param cc tuning constant for the M-estimator.
#' @param maxit maximum number of iterations allowed.
#' @param eps numeric tolerance for convergence.
#' @param adjust_bdp should the breakdown point be adjusted based on the
#'      effective degrees of freedom?
#' @param en_correction should the corrected EN estimator be used to choose
#'      the optimal lambda with CV.
#'      If \code{TRUE}, as by default, the estimator is "bias corrected".
#' @param verbosity verbosity of the algorithm.
#' @return a checked options list.
#' @export
#' @family specifying additional options
mstep_options <- function (
    cc = 3.44,
    maxit = 1000,
    eps = 1e-6,
    adjust_bdp = FALSE,
    verbosity = 0,
    en_correction = TRUE
) {
    return(list(
        maxit = .check_arg(maxit, "integer", range = 0),
        eps = .check_arg(eps, "numeric", range = 0),
        cc = .check_arg(cc, "numeric", range = 0),
        adjustBdp = .check_arg(adjust_bdp, "logical"),
        naiveEn = !.check_arg(en_correction, "logical"),
        verbosity = .check_arg(verbosity, "integer", range = 0,
                               range_test_lower = ">=")
    ))
}

#' Additional Options for the EN Algorithms
#'
#' Specify additional options for the augmented LARS and the DAL algorithm
#' to solve the EN problem.
#'
#' @param use_gram should the Gram matrix be pre-computed.
#' @param eps numeric tolerance for convergence.
#' @return a checked options list.
#' @export
#' @family specifying additional options
#' @rdname en_options
#' @aliases en_options
en_options_aug_lars <- function (
    use_gram = c("auto", "yes", "no"),
    eps = 1e-12
) {
    use_gram <- match.arg(use_gram)
    use_gram <- as.integer(pmatch(use_gram, c("auto", "yes", "no"))) - 1L
    return(list(
        algorithm = .enalgo2IntEnalgo("augmented-lars"),
        warmStart = FALSE,
        useGram = use_gram,
        eps = .check_arg(eps, "numeric", range = 0)
    ))
}

#' @param maxit maximum number of iterations allowed.
#' @param eta_mult multiplier to increase eta at each iteration.
#' @param eta_start the start value for eta.
#' @param eta_start_numerator if \code{eta_start} is missing, it is defined
#'      by \code{eta_start = eta_start_numerator / lambda}.
#' @param preconditioner preconditioner for the numerical solver. If none,
#'      a standard solver will be used, otherwise the faster preconditioned
#'      conjugate gradient is used.
#' @param verbosity verbosity of the algorithm.
#' @export
#' @rdname en_options
en_options_dal <- function (
    maxit = 100,
    eps = 1e-8,
    eta_mult = 2,
    eta_start_numerator = 1e-2,
    eta_start,
    preconditioner = c("approx", "none", "diagonal"),
    verbosity = 0
) {
    eta_start <- if (!missing(eta_start)) {
        .check_arg(eta_start, "numeric", range = 0)
    } else {
        -1
    }

    preconditioner <- match.arg(preconditioner)
    preconditioner <- -1L + as.integer(
        pmatch(preconditioner, c("none", "diagonal", "approx"))
    )

    return(list(
        algorithm = .enalgo2IntEnalgo("dal"),
        warmStart = FALSE,
        maxit = .check_arg(maxit, "integer", range = 0),
        etaMultiplier = .check_arg(eta_mult, "numeric", range = 1),
        etaStartNumerator = .check_arg(eta_start_numerator, "numeric", range = 0),
        etaStart = eta_start,
        eps = .check_arg(eps, "numeric", range = 0),
        preconditioner = .check_arg(preconditioner, "integer"),
        verbosity = .check_arg(
            verbosity,
            "integer",
            range = 0,
            range_test_lower = ">="
        )
    ))
}


## Utility function to check arguments for the correct type and optionally
## length and value range.
##
##
.check_arg <- function (
    x,
    type = c("integer", "numeric", "logical", "character"),
    range,
    length = 1L,
    range_test_lower = ">",
    range_test_upper = "<"
) {
    type <- match.arg(type)
    ok_type <- switch(
        type,
        integer = is.integer(x) || (
            is.numeric(x) && isTRUE(abs(as.integer(x) - x) < .Machine$double.eps)
        ),
        numeric = is.numeric(x),
        logical = is.logical(x),
        character = is.character(x)
    )

    if (!isTRUE(ok_type)) {
        stop(sprintf("%s is expected to be `%s`", deparse(substitute(x)), type))
    }

    x <- switch(
        type,
        integer = as.integer(x),
        numeric = as.double(x),
        x
    )

    ok_length <- if (!is.null(length)) {
        length(x) == length
    } else {
        TRUE
    }

    if (!isTRUE(ok_length)) {
        stop(sprintf("%s is expected to be of length %d",
                     deparse(substitute(x)), length))
    }

    range_fun_lower <- match.fun(range_test_lower)
    range_fun_upper <- match.fun(range_test_upper)

    ok_range <- if (!missing(range)) {
        if (length(range) == 1L) {
            all(range_fun_lower(x, range[1L]))
        } else {
            all(range_fun_lower(x, range[1L]) & range_fun_upper(x, range[2L]))
        }
    } else {
        TRUE
    }

    if (!isTRUE(ok_range)) {
        message <- if (length(range) > 1L) {
            message <- paste(
                range_test_lower,
                range[1L],
                "and",
                range_test_upper,
                range[2L],
                sep = " "
            )
        } else {
            paste(range_test_lower, range[1L], sep = " ")
        }
        stop(sprintf("%s is expected to be %s",
                     deparse(substitute(x)), message))
    }

    return(x)
}

## Get the integer representation of a supported rho function
##
##
.rho2IntRho <- function(rho_fun = c("huber", "bisquare", "gauss")) {
    rho_fun <- switch(
        match.arg(rho_fun),
        huber = 0L,
        bisquare = 1L,
        gauss = 5L
    )

    rho_fun <- rho_fun[which(!is.na(rho_fun))[1L]]

    if (is.na(rho_fun)) {
        rho_fun <- 1L
        warning("Unknown rho function selected. Using Tukey's bisquare.")
    }

    return(rho_fun)
}

## Get the integer representation of a supported EN algorithm
##
.enalgo2IntEnalgo <- function(en_algo = c("augmented-lars", "dal")) {
    en_algo <- switch(
        match.arg(en_algo),
        `augmented-lars` = 0L,
        dal = 1L
    )

    en_algo <- en_algo[which(!is.na(en_algo))[1L]]

    if (is.na(en_algo)) {
        en_algo <- 0L
        warning("Unknown Elastic Net algorithm selected. Using augmented LARS.")
    }

    return(en_algo)
}

