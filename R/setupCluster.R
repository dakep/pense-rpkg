## Convenience function for parallel computing
##
## This function handles setting up, using, and closing a potential cluster.
## If no cluster of computing nodes is requested, it will create a proxy to the local
## \code{lapply} function.
##
#' @importFrom parallel clusterEvalQ clusterExport clusterApplyLB stopCluster
#' @importFrom parallel makePSOCKcluster clusterSetRNGStream
#' @importFrom methods is
setupCluster <- function(ncores = 1L, cl = NULL, eval, export, envir = parent.frame()) {
    retlocal <- list(
        lapply = lapply,
        ncores = 1L,
        stopCluster = function() {},
        setSeed = set.seed,
        exportedObject = function(obj) {
            return(obj)
        }
    )

    ret <- retlocal

    ##
    ## Set up a potential cluster
    ##
    if (!is.numeric(ncores) || length(ncores) != 1L || ncores < 1) {
        warning("`ncores` must be a positive integer of length one.")
    } else {
        ret$ncores <- ncores
    }

    tryCatch({
        withCallingHandlers({
            if (is(cl, "cluster") || ncores > 1L) {
                if(!is(cl, "cluster")) {
                    cl <- makePSOCKcluster(ncores)
                    ret$stopCluster <- function() {
                        parallel::stopCluster(cl)
                    }
                }

                ret$ncores <- length(cl)

                if (!missing(eval)) {
                    clusterEvalQ(cl, eval)
                }

                if (!missing(export)) {
                    clusterExport(cl, export, envir = envir)
                }

                ret$setSeed <- function(seed) {
                    clusterSetRNGStream(cl, iseed = seed)
                }

                ret$exportedObject <- function(obj) {
                    substitute(obj)
                }

                ret$lapply <- function(...) {
                    clusterApplyLB(cl, ...)
                }
            }

        }, error = function(...) {
            ret <- retlocal
        })
    }, error = function(e) {
        warning("Error during cluster setup: ", e)
    }, finally = {})

    return(ret)
}
