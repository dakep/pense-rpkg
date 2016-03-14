#' Control Parameters for PENSE
#'
#' Set control parameters for PENSE.
#'
#' @return A list with the given arguments.
#'
#' @export
pense.control <- function(
    pense.maxit = 500,
    pense.tol = 1e-6,
    pense.en.tol = 1e-8,
    pense.en.maxit = 5e5,

    init.tol = 1e-6,
    init.maxit = 10,
    init.resid.clean.method = c("proportion", "threshold"),
    init.resid.proportion = 0.4,
    init.resid.threshold = 2,
    init.psc.method = c("rr", "Mn"),
    init.psc.proportion = 0.8,
    init.csteps = 10,
    init.nkeep = 5,
    init.en.tol = pense.en.tol,
    init.en.maxit = pense.en.maxit,

    mscale.cc = 1.54764,
    mscale.delta = 0.5,
    mscale.maxit = 200,
    mscale.tol = 1e-8,
    mscale.rho.fun = c("bisquare", "huber", "gauss")
) {
    ret <- as.list(environment())

    ret$init.psc.method <- match.arg(init.psc.method)
    ret$init.resid.clean.method <- match.arg(init.resid.clean.method)
    ret$mscale.rho.fun <- match.arg(mscale.rho.fun)

    return(ret)
}

#' Control Parameters for ENPY
#'
#' Set control parameters for ENPY.
#'
#' @param enMaxIt Maximum number of iterations allowed for the elastic net algorithm.
#' @param enEPS Convergence threshold for elastic net.
#' @param enCentering Should the data be centered for the elastic net algorithm. Only those rows
#'          with leading 1's will be considered for centering.
#' @param mscaleMaxIt Maximum number of iterations allowed for the m-scale algorithm.
#' @param mscaleEPS Convergence threshold for the m-scale
#' @param mscaleRhoFun The rho function to use for the m-scale.
#'
#' @return A list with the given arguments.
#'
#' @export
enpy.control <- function(enMaxIt = 50000,
                         enEPS = 1e-8,
                         enCentering = TRUE,
                         mscaleMaxIt = 200,
                         mscaleEPS = 1e-8,
                         mscaleRhoFun = c("bisquare", "huber", "gauss")) {
    ret <- as.list(environment())
    ret$mscaleRhoFun <- match.arg(mscaleRhoFun)
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
    ctrl$init.psc.method <- match.arg(ctrl$init.psc.method, c("rr", "Mn"))
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
#' @param residCleanMethod character
#' @param residThreshold If \code{residCleanMethod = "threshold"} numeric > 0, otherwise not
#'          referenced
#' @param residProportion If \code{residCleanMethod = "proportion"} numeric > 0, otherwise not
#'          referenced
#' @param pscProportion numeric > 0
#' @param mscaleB numeric > 0
#' @param mscaleCC numeric > 0
#'
initest.control <- function(
    lambda1,
    lambda2,
    numIt,
    eps = 1e-6,
    residCleanMethod = c("proportion", "threshold"),
    residThreshold,
    residProportion,
    pscProportion,
    mscaleB = 0.5,
    mscaleCC = 1.54764,
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

    ret$residCleanMethod <- match.arg(residCleanMethod)
    ret$mscaleRhoFun <- as.integer(pmatch(ret$mscaleRhoFun, c("bisquare", "huber", "gauss"))) - 1L

    ret$mscaleRhoFun <- ret$mscaleRhoFun[which(!is.na(ret$mscaleRhoFun))[1L]]

    if (ret$residCleanMethod == "proportion") {
        simpleCheck(residProportion)

        ret$residThreshold <- -1

        if (residProportion > 1) {
            stop("`residProportion` must be less than 1")
        }
    } else {
        simpleCheck(residThreshold)
        ret$residProportion <- -1
    }

    simpleCheck(lambda1, FALSE)
    simpleCheck(lambda2, FALSE)
    simpleCheck(numIt)
    simpleCheck(eps)
    simpleCheck(pscProportion)
    simpleCheck(ret$enMaxIt)
    simpleCheck(ret$enEPS)
    simpleCheck(mscaleB)
    simpleCheck(mscaleCC)
    simpleCheck(ret$mscaleMaxIt)
    simpleCheck(ret$mscaleEPS)
    simpleCheck(ret$mscaleRhoFun, FALSE)

    if (length(ret$enCentering) != 1L || !is.logical(ret$enCentering) || is.na(ret$enCentering)) {
        stop("`enCentering` must be single logical value")
    }

    if (ret$pscProportion > 1) {
        stop("`pscProportion` must be less than 1")
    }

    ret$numIt <- as.integer(ret$numIt)
    ret$enMaxIt <- as.integer(ret$enMaxIt)
    ret$mscaleMaxIt <- as.integer(ret$mscaleMaxIt)

    return(ret)
}
