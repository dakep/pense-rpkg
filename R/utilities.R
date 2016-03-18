nameCoefVec <- function(coef, X) {
    dn <- dimnames(X)
    xnames <- paste("X", seq_len(ncol(X)), sep = "")

    if (!is.null(dn) && !is.null(dn[[2L]])) {
        xnames <- dn[[2L]]
    }

    names(coef) <- c("(Intercept)", xnames)
    return(coef)
}

orderOmitTies <- function(x, tol) {
    ord.x <- sort.list(x, na.last = NA, method = "quick")
    sorted <- x[ord.x]
    diffs <- diff(sorted)

    filtered.x <- c(sorted[1L], sorted[-1L][diffs > tol])
    filtered.ind <- c(ord.x[1L], ord.x[-1L][diffs > tol])

    return(list(
        clean = filtered.x,
        index = filtered.ind
    ))
}

facon <- function(delta) {
    23.9716 - 73.4391 * delta + 64.9480 * delta^2
}
