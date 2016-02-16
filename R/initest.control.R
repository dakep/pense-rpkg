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


#' Creates the internal control list for PY initial estimators
#'
#' Computes the PY initial estimates for EN with approximated principal sensitivity components
#' by the ridge regression solution.
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
#' @param enMaxIt integer > 0
#' @param enEPS numeric > 0
#' @param enCentering logical
#' @param mscaleB numeric > 0
#' @param mscaleCC numeric > 0
#' @param mscaleMaxIt integer > 0
#' @param mscaleEPS numeric > 0
#' @param mscaleRhoFun character
#'
initest.control <- function(lambda1,
                            lambda2,
                            numIt,
                            eps = 1e-6,
                            residCleanMethod = c("proportion", "threshold"),
                            residThreshold,
                            residProportion,
                            pscProportion,
                            mscaleB = 0.5,
                            mscaleCC = 1.54764,

                            enpy.control) {
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
