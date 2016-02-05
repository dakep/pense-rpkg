//
//  psc.c
//  penseinit
//
//  Created by David Kepplinger on 2016-01-30.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#include <Rinternals.h>
#include <float.h>
#include <stdlib.h>
#include "BLAS.h"

#include "psc.h"


static const int BLAS_M1L = -1;
static const int BLAS_1L = 1;
static const double BLAS_0F = 0.0;
static const double BLAS_1F = 1.0;
static const double BLAS_M1F = -1.0;

static const double LAPACK_EV_ABSTOL = 1e-18;
static const double LAPACK_EV_RANGE_LOWER = 1e-12;
static const double LAPACK_EV_RANGE_UPPER = DBL_MAX;

static const char * const BLAS_SIDE_LEFT = "L";
static const char * const BLAS_UPLO_UPPER = "U";
static const char * const BLAS_TRANS_NO = "N";
static const char * const BLAS_TRANS_TRANS = "T";
static const char * const BLAS_DIAG_NO = "N";

static int doEigenDecomposition(const char *const uplo, double *restrict matrix, const int n,
                                AuxMemory *const auxmem);

/**
 * Calculate the principal sensitivity components (PSCs)
 *
 * @param RXtr     numeric The (nvar by nobs) transposed X matrix
 * @param Ry	   numeric The (nobs) y vector
 * @param Rnobs    integer The number of observations in X and y
 * @param Rnvar    integer The number of variables in X
 * @param RnumPSCs integer The actual number of PSCs (OUTPUT!)
 *
 * @return Returns the (probably) too large (nobs by nvar) numeric matrix of PSCs
 */
SEXP C_pscs(SEXP RXtr, SEXP Ry, SEXP Rnobs, SEXP Rnvar, SEXP RnumPSCs)
{
    const int nobs = *INTEGER(Rnobs);
    const int nvar = *INTEGER(Rnvar);
    AuxMemory auxmem;
    initAuxMemory(&auxmem);

    SEXP ret = PROTECT(Rf_allocVector(REALSXP, nobs * nvar));

    *INTEGER(RnumPSCs) = calculatePSCs(REAL(ret), &auxmem, REAL(RXtr), REAL(Ry), nobs, nvar);

    UNPROTECT(1);
    return ret;
}



int calculatePSCs(double *restrict pscs, AuxMemory* auxmem,
                         const double *restrict Xtr, const double *restrict y,
                         const int nobs, const int nvar)
{
    int lapackInfo = 0;
    int i;
    int nevalues = 0;
    /*
     * G can point to the same address as auxmem->XtXinvX, because we don't need
     * auxmem->XtXinvX anymore when G is computed.
     * The same is true for Q and H, as well as Z and XtXinvX
     */
    double *restrict XtXinvX = pscs;
    double *restrict G = pscs;
    double *restrict Q;
    double *iterG, *iterH;
    double Wii;

    /**
     * Make sure we have enough space in the auxiliary memory
     */
    resizeAuxMemory(auxmem, nvar, nobs);

    Q = auxmem->H;

    memcpy(auxmem->XsqrtInvX, Xtr, nobs * nvar * sizeof(double));
    memcpy(auxmem->residuals, y, nobs * sizeof(double));

    /* Xsqrt = X %*% t(X) */
    BLAS_DGEMM(BLAS_TRANS_NO, BLAS_TRANS_TRANS, nvar, nvar, nobs,
        BLAS_1F, Xtr, nvar, Xtr, nvar,
        BLAS_0F, auxmem->Xsqrt, nvar);

    /* Xsqrt = chol(Xsqrt) */
    LAPACK_DPOTRF(BLAS_UPLO_UPPER, nvar, auxmem->Xsqrt, nvar, lapackInfo);

    if (lapackInfo != 0) {
        return -abs(lapackInfo);
    }

    /* t(XsqrtInvX) = t(inv(Xsqrt)) %*% t(X) */
    BLAS_DTRSM(BLAS_SIDE_LEFT, BLAS_UPLO_UPPER, BLAS_TRANS_TRANS, BLAS_DIAG_NO,
               nvar, nobs, BLAS_1F, auxmem->Xsqrt, nvar, auxmem->XsqrtInvX, nvar);

    memcpy(XtXinvX, auxmem->XsqrtInvX, nobs * nvar * sizeof(double));

    /* XtXinvX = inv(Xsqrt) %*% t(XsqrtInvX) */
    BLAS_DTRSM(BLAS_SIDE_LEFT, BLAS_UPLO_UPPER, BLAS_TRANS_NO, BLAS_DIAG_NO,
               nvar, nobs, BLAS_1F, auxmem->Xsqrt, nvar, XtXinvX, nvar);

    /* H = t(X) %*% XtXinvX */
    BLAS_DGEMM(BLAS_TRANS_TRANS, BLAS_TRANS_NO, nobs, nobs, nvar,
               BLAS_1F, Xtr, nvar, XtXinvX, nvar,
               BLAS_0F, auxmem->H, nobs);

    /* residuals = -H %*% y + y --> H is symmetric! */
    BLAS_DSYMV(BLAS_UPLO_UPPER, nobs,
               BLAS_M1F, auxmem->H, nobs, y, BLAS_1L,
               BLAS_1F, auxmem->residuals, BLAS_1L);

    memcpy(G, auxmem->XsqrtInvX, nobs * nvar * sizeof(double));

    /* t(G) = t(XsqrtInvX) %*% diag(W[i, i]) */
    for (i = 0, iterH = auxmem->H, iterG = G; i < nobs; ++i, iterH += nobs + 1, iterG += nvar) {
        /* W[i, i] = */
        Wii = auxmem->residuals[i] / (1 - *iterH);

        /* G[, i] = W[i, i] * XsqrtInvX[, i] */
        BLAS_DSCAL(nvar, Wii, iterG, BLAS_1L); // Do loop by hand!
    }

    /* Q = t(G) %*% G */
    BLAS_DGEMM(BLAS_TRANS_NO, BLAS_TRANS_TRANS, nvar, nvar, nobs,
               BLAS_1F, G, nvar, G, nvar,
               BLAS_0F, Q, nvar);

    nevalues = doEigenDecomposition(BLAS_UPLO_UPPER, Q, nvar, auxmem);

    if (nevalues > 0) {
        BLAS_DGEMM(BLAS_TRANS_TRANS, BLAS_TRANS_NO, nobs, nevalues, nvar,
                   BLAS_1F, auxmem->XsqrtInvX, nvar, auxmem->eigenvectors, nvar,
                   BLAS_0F, pscs, nobs);
    }

    return nevalues;
}


static int doEigenDecomposition(const char *const uplo, double *restrict matrix, const int n,
                                AuxMemory *const auxmem)
{
    int nevalues, lapackInfo;

    /* Compute needed size of working space */
    LAPACK_DSYEVR_RANGE(uplo, n, matrix, n,
                        LAPACK_EV_RANGE_LOWER, LAPACK_EV_RANGE_UPPER, LAPACK_EV_ABSTOL,
                        nevalues, auxmem->evalues,
                        auxmem->eigenvectors, n, auxmem->evectorsSupport,
                        auxmem->dblWorkMem, BLAS_M1L,
                        auxmem->intWorkMem, BLAS_M1L, lapackInfo);

    if (lapackInfo != 0) {
        return -abs(lapackInfo);
    }

    resizeDblWorkAuxMemory(auxmem, (int) auxmem->dblWorkMem[0]);
    resizeIntWorkAuxMemory(auxmem, auxmem->intWorkMem[0]);

    /* Actually perform eigenvalue decomposition */
    LAPACK_DSYEVR_RANGE(uplo, n, matrix, n,
                        LAPACK_EV_RANGE_LOWER, LAPACK_EV_RANGE_UPPER, LAPACK_EV_ABSTOL,
                        nevalues, auxmem->evalues,
                        auxmem->eigenvectors, n, auxmem->evectorsSupport,
                        auxmem->dblWorkMem, auxmem->dblWorkMemSize,
                        auxmem->intWorkMem, auxmem->intWorkMemSize, lapackInfo);


    if (lapackInfo != 0) {
        return -abs(lapackInfo);
    }

    return nevalues;
}
