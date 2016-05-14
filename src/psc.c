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
#include "olsreg.h"

#include "psc.h"


static BLAS_INT BLAS_M1L = -1;
static const double BLAS_0F = 0.0;
static const double BLAS_1F = 1.0;

static const double LAPACK_EV_ABSTOL = NUMERIC_EPS;
static const double LAPACK_EV_RANGE_LOWER = LAPACK_EV_MIN;
static const double LAPACK_EV_RANGE_UPPER = DBL_MAX;

static BLAS_CHAR BLAS_SIDE_LEFT = "L";
static BLAS_CHAR BLAS_UPLO_UPPER = "U";
static BLAS_CHAR BLAS_TRANS_NO = "N";
static BLAS_CHAR BLAS_TRANS_TRANS = "T";
static BLAS_CHAR BLAS_DIAG_NO = "N";

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
    SEXP ret;
    const int nobs = *INTEGER(Rnobs);
    const int nvar = *INTEGER(Rnvar);
    int * numPSCs = INTEGER(RnumPSCs);
    double *restrict Xtr = REAL(RXtr);
    double *restrict y = REAL(Ry);
    double *restrict coefs = (double *) malloc(nvar * sizeof(double));
    AuxMemory auxmem;

    initAuxMemory(&auxmem);
    resizeAuxMemory(&auxmem, nvar, nobs);

    *numPSCs = computeOLSCoefs(Xtr, y, nobs, nvar, coefs, auxmem.Xsqrt);

    if (*numPSCs == 0) {
        computeResiduals(Xtr, y, nobs, nvar, coefs, auxmem.residuals);

        ret = PROTECT(Rf_allocVector(REALSXP, nobs * nvar));

        *numPSCs = calculatePSCs(REAL(ret), &auxmem, Xtr, y, nobs, nvar);

        UNPROTECT(1);
        return ret;
    }

    return R_NilValue;
}

int calculatePSCs(double *restrict pscs, AuxMemory* auxmem,
                         const double *restrict Xtr, const double *restrict y,
                         const int nobs, const int nvar)
{
    int i, j;
    int nevalues = 0;
    /*
     * G can point to the same address as auxmem->XtXinvX, because we don't need
     * auxmem->XtXinvX anymore when G is computed.
     * The same is true for Q and H, as well as Z and XtXinvX
     */
    double *restrict XtXinvX = pscs;
    double *restrict G = pscs;
    const double *restrict iterX;
    const double *restrict iterXtXinvX;
    const double *restrict iterXsqrtInvX;
    double *iterG;
    double Wii;

    /**
     * Make sure we have enough space in the auxiliary memory
     */
    resizeAuxMemory(auxmem, nvar, nobs);

    memcpy(auxmem->XsqrtInvX, Xtr, nobs * nvar * sizeof(double));

    /* t(XsqrtInvX) = t(inv(Xsqrt)) %*% t(X) */
    BLAS_DTRSM(BLAS_SIDE_LEFT, BLAS_UPLO_UPPER, BLAS_TRANS_TRANS, BLAS_DIAG_NO,
               nvar, nobs, BLAS_1F, auxmem->Xsqrt, nvar, auxmem->XsqrtInvX, nvar);

    memcpy(XtXinvX, auxmem->XsqrtInvX, nobs * nvar * sizeof(double));

    /* XtXinvX = inv(Xsqrt) %*% t(XsqrtInvX) */
    BLAS_DTRSM(BLAS_SIDE_LEFT, BLAS_UPLO_UPPER, BLAS_TRANS_NO, BLAS_DIAG_NO,
               nvar, nobs, BLAS_1F, auxmem->Xsqrt, nvar, XtXinvX, nvar);

    /* t(G) = t(XsqrtInvX) %*% diag(W[i, i]) */
    iterG = G;
    iterX = Xtr;
    iterXtXinvX = XtXinvX;
    iterXsqrtInvX = auxmem->XsqrtInvX;

    for (i = 0; i < nobs; ++i) {
        /*
         * Set Wii to i-th diagonal element of H,
         * i.e., inner product of x_i and i-th column of XtXinvX
         */
        Wii = 0;
        for (j = 0; j < nvar; ++j, ++iterX, ++iterXtXinvX) {
            Wii += (*iterX) * (*iterXtXinvX);
        }

        /* W[i, i] = */
        Wii = auxmem->residuals[i] / (1 - Wii);

        /* G[ , i] = W[i, i] * XsqrtInvX[, i] */
        for (j = 0; j < nvar; ++j, ++iterG, ++iterXsqrtInvX) {
            (*iterG) = Wii * (*iterXsqrtInvX);
        }
    }


    /* Q = t(G) %*% G */
    BLAS_DGEMM(BLAS_TRANS_NO, BLAS_TRANS_TRANS, nvar, nvar, nobs,
               BLAS_1F, G, nvar, G, nvar,
               BLAS_0F, auxmem->Q, nvar);

    nevalues = doEigenDecomposition(BLAS_UPLO_UPPER, auxmem->Q, nvar, auxmem);

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
