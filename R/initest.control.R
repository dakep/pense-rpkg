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

                            enMaxIt = 1000,
                            enEPS = 1e-6,

                            mscaleB = 0.5,
                            mscaleCC = 1.54764,
                            mscaleMaxIt = 200,
                            mscaleEPS = 1e-8,
                            mscaleRhoFun = c("bisquare", "huber", "gauss")) {
    ret <- as.list(environment())


    simpleCheck <- function(x, checkLTE = TRUE) {
        if (length(x) != 1L || !is.numeric(x) || is.na(x) || x < 0 || (checkLTE && (x == 0))) {
            stop(sprintf("`%s` must be single positive number", deparse(substitute(x))))
        }
    }

    ret$residCleanMethod <- match.arg(residCleanMethod)
    ret$mscaleRhoFun <- as.integer(pmatch(mscaleRhoFun, c("bisquare", "huber", "gauss"))) - 1L

    ret$mscaleRhoFun <- ret$mscaleRhoFun[which(!is.na(ret$mscaleRhoFun))[1L]]

    if (ret$residCleanMethod == "proportion") {
        simpleCheck(residProportion)
        ret$residThreshold <- -1
    } else {
        simpleCheck(residThreshold)
        ret$residProportion <- -1
    }

    simpleCheck(lambda1, FALSE)
    simpleCheck(lambda2, FALSE)
    simpleCheck(numIt)
    simpleCheck(eps)
    simpleCheck(pscProportion)
    simpleCheck(enMaxIt)
    simpleCheck(enEPS)
    simpleCheck(mscaleB)
    simpleCheck(mscaleCC)
    simpleCheck(mscaleMaxIt)
    simpleCheck(mscaleEPS)
    simpleCheck(ret$mscaleRhoFun, FALSE)


    ret$numIt <- as.integer(ret$numIt)
    ret$enMaxIt <- as.integer(ret$enMaxIt)
    ret$mscaleMaxIt <- as.integer(ret$mscaleMaxIt)

    return(ret)
}
