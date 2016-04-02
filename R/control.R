#' Control Parameters for PENSE
#'
#' Set control parameters for PENSE.
#'
#' @param pense.maxit maximum number of iterations for PENSE.
#' @param pense.tol convergence tolerance for PENSE.
#' @param pense.en.tol convergence tolerance for the elastic net in PENSE
#' @param pense.en.maxit maximum number of iterations for each elastic net fit in PENSE
#'
#' @param cv.objective a user-supplied function taking the cross-validated residuals
#'          to choose lambda. Default is to choose the lambda the minimizes the
#'          \link[robustbase:scaleTau2]{tau-scale}.
#'
#' @param init.tol convergence tolerance for initial estimates.
#' @param init.maxit maximum number of iterations for initial estimates.
#' @param init.resid.clean.method how to clean the data based on large residuals.
#'          If \code{"threshold"}, all observations with scaled residuals larger than
#'          \code{init.resid.threshold} will be removed (\code{init.resid.threshold}
#'          corresponds to the constant \eqn{C_1} from equation (21) in Pena & Yohai (1999).
#'          If \code{"proportion"}, observations with the largest
#'          \code{init.resid.proportion} residuals will be removed.
#' @param init.resid.proportion,init.resid.threshold see details of argument
#'          \code{init.resid.clean.method}
#' @param init.psc.method how to calculate the principal sensitivity components. If \code{"auto"},
#'          \code{rr} is used when there are enough observations, otherwise \code{Mn} is
#'          selected. See \code{\link{enpy}} for details.
#' @param init.psc.proportion the proportion of observations to remove based on the score on
#'          the principal sensitivity components.
#' @param init.csteps the number of PENSE iterations (concentration steps) to perform on
#'          the initial estimates.
#' @param init.nkeep how many initial estimates should be kept to perform full PENSE iterations?
#' @param init.en.tol convergence tolerance for elastic net fits to find initial estimates.
#' @param init.en.maxit maximum number of iterations for each elastic net fit for initial estimates.
#'
#' @param mscale.cc tuning parameter for the M-estimate of scale.
#' @param mscale.delta expected value under the normal model of the rho function with tuning
#'          constant equal to \code{mscale.cc}.
#' @param mscale.maxit maximum number of iterations to calculate the M-estimate of scale.
#' @param mscale.tol convergence tolerance for calculating the M-estimate of scale.
#' @param mscale.rho.fun string specifying the rho function of the S-estimator.
#'
#' @return A list with the given arguments, checked for errors.
#'
#' @references Pena, D., & Yohai, V.. (1999). A Fast Procedure for Outlier Diagnostics in Large
#' Regression Problems. \emph{Journal of the American Statistical Association}, 94(446),
#' 434–445. \url{http://doi.org/10.2307/2670164}
#'
#' @export pense.control
#' @importFrom robustbase scaleTau2
pense.control <- function(
    pense.maxit = 500,
    pense.tol = 1e-6,
    pense.en.tol = 1e-9,
    pense.en.maxit = 5e5,

    cv.objective = scaleTau2,

    init.tol = 1e-6,
    init.maxit = 10,
    init.resid.clean.method = c("proportion", "threshold"),
    init.resid.proportion = 0.4,
    init.resid.threshold = 2,
    init.psc.method = c("auto", "rr", "Mn"),
    init.psc.proportion = 0.8,
    init.csteps = 10,
    init.nkeep = 5,
    init.en.tol = 1e-8,
    init.en.maxit = pense.en.maxit,

    mscale.cc = 1.54764,
    mscale.delta = 0.5,
    mscale.maxit = 200,
    mscale.tol = 1e-8,
    mscale.rho.fun = c("bisquare", "huber", "gauss")
) {
    ret <- as.list(environment())

    ret$cv.objective <- match.fun(ret$cv.objective)
    ret$init.psc.method <- match.arg(init.psc.method)
    ret$init.resid.clean.method <- match.arg(init.resid.clean.method)
    ret$mscale.rho.fun <- match.arg(mscale.rho.fun)

    return(.check.pense.control(ret))
}

#' Control Parameters for ENPY
#'
#' Set control parameters for ENPY.
#'
#' @param en.maxit Maximum number of iterations allowed for the elastic net algorithm.
#' @param en.tol Convergence threshold for elastic net.
#' @param en.centering Should the data be centered for the elastic net algorithm. Only those rows
#'          with leading 1's will be considered for centering.
#' @param mscale.maxit Maximum number of iterations allowed for the m-scale algorithm.
#' @param mscale.tol Convergence threshold for the m-scale
#' @param mscale.rho.fun The rho function to use for the m-scale.
#'
#' @return A list with the given arguments.
#'
#' @export
enpy.control <- function(en.maxit = 50000,
                         en.tol = 1e-8,
                         en.centering = TRUE,
                         mscale.maxit = 200,
                         mscale.tol = 1e-8,
                         mscale.rho.fun = c("bisquare", "huber", "gauss")) {
    ret <- as.list(environment())
    ret$mscale.rho.fun <- match.arg(mscale.rho.fun)
    return(ret)
}


#' Internal function to check PENSE control parameters
.check.pense.control <- function(ctrl) {
    simpleCheck <- function(x) {
        if (length(x) != 1L || !is.numeric(x) || anyNA(x) || x <= 0) {
            stop(sprintf("`%s` must be single positive number", deparse(substitute(x))))
        }
    }

    ##
    ## Check string arguments first
    ##
    ctrl$init.resid.clean.method <- match.arg(ctrl$init.resid.clean.method,
                                              c("proportion", "threshold"))
    ctrl$init.psc.method <- match.arg(ctrl$init.psc.method, c("auto", "rr", "Mn"))
    ctrl$mscale.rho.fun <- match.arg(ctrl$mscale.rho.fun, c("bisquare", "huber", "gauss"))

    with(ctrl, {
        simpleCheck(pense.maxit)
        simpleCheck(pense.tol)
        simpleCheck(pense.en.tol)
        simpleCheck(pense.en.maxit)

        simpleCheck(init.tol)
        simpleCheck(init.maxit)
        simpleCheck(init.psc.proportion)
        simpleCheck(init.csteps)
        simpleCheck(init.nkeep)
        simpleCheck(init.en.tol)
        simpleCheck(init.en.maxit)

        simpleCheck(mscale.cc)
        simpleCheck(mscale.delta)
        simpleCheck(mscale.maxit)
        simpleCheck(mscale.tol)

        if (init.psc.proportion > 1) {
            stop("`init.psc.proportion` must be less than 1")
        }
    })

    if (ctrl$init.resid.clean.method == "proportion") {
        with(ctrl, simpleCheck(init.resid.proportion))

        ctrl$init.resid.threshold <- -1

        if (ctrl$init.resid.proportion > 1) {
            stop("`init.resid.proportion` must be less than 1")
        }
    } else {
        with(ctrl, simpleCheck(init.resid.threshold))
        ctrl$init.resid.proportion <- -1
    }

    if (!is.function(ctrl$cv.objective)) {
        stop("`cv.objective` must be a function")
    }

    tmp <- ctrl$cv.objective(1:20)

    if (!is.numeric(tmp) && length(tmp) != 1L) {
        stop("`cv.objective` must return a single number")
    }

    return(ctrl)
}


#' Creates the internal control list for PY initial estimators
#'
#' Takes care of the correct storage mode for the arguments passed to the C/C++ code
#'
#' @param lambda1 numeric >= 0
#' @param lambda2 numeric >= 0
#' @param numIt integer > 0
#' @param eps numeric > 0
#' @param resid.clean.method character
#' @param resid.threshold If \code{resid.clean.method = "threshold"} numeric > 0, otherwise not
#'          referenced
#' @param resid.proportion If \code{resid.clean.method = "proportion"} numeric > 0, otherwise not
#'          referenced
#' @param psc.proportion numeric > 0
#' @param mscale.delta numeric > 0
#' @param mscale.cc numeric > 0
#'
initest.control <- function(
    lambda1,
    lambda2,
    numIt,
    eps = 1e-6,
    resid.clean.method = c("proportion", "threshold"),
    resid.threshold,
    resid.proportion,
    psc.proportion,
    mscale.delta = 0.5,
    mscale.cc = 1.54764,
    enpy.control
) {
    ret <- as.list(environment())
    ret$enpy.control <- NULL
    ret <- c(ret, enpy.control)


    simpleCheck <- function(x, checkLTE = TRUE) {
        if (length(x) != 1L || !is.numeric(x) || is.na(x) || x < 0 || (checkLTE && (x == 0))) {
            stop(sprintf("`%s` must be single positive number", deparse(substitute(x))))
        }
    }

    ret$resid.clean.method <- match.arg(resid.clean.method)
    ret$mscale.rho.fun <- as.integer(pmatch(ret$mscale.rho.fun, c("bisquare", "huber", "gauss"))) - 1L

    ret$mscale.rho.fun <- ret$mscale.rho.fun[which(!is.na(ret$mscale.rho.fun))[1L]]

    if (ret$resid.clean.method == "proportion") {
        simpleCheck(resid.proportion)

        ret$resid.threshold <- -1

        if (resid.proportion > 1) {
            stop("`resid.proportion` must be less than 1")
        }
    } else {
        simpleCheck(resid.threshold)
        ret$resid.proportion <- -1
    }

    simpleCheck(lambda1, FALSE)
    simpleCheck(lambda2, FALSE)
    simpleCheck(numIt)
    simpleCheck(eps)
    simpleCheck(psc.proportion)
    simpleCheck(ret$en.maxit)
    simpleCheck(ret$en.tol)
    simpleCheck(mscale.delta)
    simpleCheck(mscale.cc)
    simpleCheck(ret$mscale.maxit)
    simpleCheck(ret$mscale.tol)
    simpleCheck(ret$mscale.rho.fun, FALSE)

    if (length(ret$en.centering) != 1L || !is.logical(ret$en.centering) || is.na(ret$en.centering)) {
        stop("`en.centering` must be single logical value")
    }

    if (ret$psc.proportion > 1) {
        stop("`psc.proportion` must be less than 1")
    }

    ret$numIt <- as.integer(ret$numIt)
    ret$en.maxit <- as.integer(ret$en.maxit)
    ret$mscale.maxit <- as.integer(ret$mscale.maxit)

    return(ret)
}

