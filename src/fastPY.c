//
//  fastPY.c
//  penseinit
//
//  Created by David Kepplinger on 2016-01-28.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>
#include <float.h>

#include "BLAS.h"
#include "Control.h"
#include "mscale.h"
#include "PartialSort.h"
#include "AuxMemory.h"
#include "psc.h"

#ifdef DEBUG
static void print_matf(int dr, int dc, const double *A, const char *header) {
    int r, c, i = 0;
    Rprintf("%s:\n", header);
    for (r = 0; r < dr; ++r) {
        for (c = 0, i = r; c < dc; ++c, i += dr) {
            Rprintf("%7.4f ", A[i]);
        }
        Rprintf("\n");
    }
}
#endif

static const int BLAS_1L = 1;
static const double BLAS_0F = 0.0;
static const double BLAS_1F = 1.0;
static const double BLAS_M1F = -1.0;

static const char * const BLAS_UPLO_UPPER = "U";
static const char * const BLAS_TRANS_NO = "N";
static const char * const BLAS_TRANS_TRANS = "T";
static const char * const BLAS_DIAG_NO = "N";

/**
 * NOTE: `estimates` must be able to hold nvar * (3 * nvar + 1) + nobs elements.
 * The first nvar * (3 * nvar + 1) elements are the initial estimates, the remaining
 * nobs elements can be discared!
 */
static int computeInitialEstimator(const double *restrict X, const double *restrict y,
                                   const int nobs, const int nvar, const Control* ctrl,
                                   double *restrict estimates);

static int estimateCoef(double *restrict X, const double *restrict y, const int nobs,
                        const int nvar, double *restrict coefEst, AuxMemory* auxmem);


static void calculateResiduals(const double *restrict X, const double *restrict y,
                               const double *restrict coefEst,
                               const int nobs, const int nvar, double *restrict residuals);

static int filterDataThreshold(const double *restrict X, const double *restrict y,
                                double *restrict newX, double *restrict newY,
                                const int nobs, const int nvar,
                                const double *restrict values, const double threshold,
                                CompareFunction compare);

/*
 * Compare functions
 */
static double lessThan(const double a, const double b);
static double greaterThan(const double a, const double b);
static double absoluteLessThan(const double a, const double b);


/**
 * Calculate the Pena-Yohai initial estimator
 *
 * @param RXtr   numeric The (nvar by nobs) transposed X matrix
 * @param Ry	 numeric The (nobs) y vector
 * @param Rnobs  integer The number of observations in X and y
 * @param Rnvar  integer The number of variables in X
 * @param RnumIt integer The maximum number of iterations
 * @param Reps   numeric The relative tolerance for convergence

 * @param RnumInits integer OUTPUT the number of actual initial estimators
 *
 * @return Returns the (probably) too large (3 * nvar + 2 by nvar) numeric
 *         matrix of initial estimators
 */
SEXP C_initpy(SEXP RXtr, SEXP Ry, SEXP Rnobs, SEXP Rnvar, SEXP RnumIt,
              SEXP Reps, SEXP RresidThreshold, SEXP RresidProportion,
              SEXP RpscProportion, SEXP RmscaleB, SEXP RmscaleCC,
              SEXP RmscaleMaxIt, SEXP RmscaleEPS, SEXP RmscaleRhoFun,
              SEXP RnumInits)
{
    Control ctrl = {
        .numIt = *INTEGER(RnumIt),
        .eps = *REAL(Reps),
        .residThreshold = *REAL(RresidThreshold),
        .residProportion = *REAL(RresidProportion),
        .pscProportion = *REAL(RpscProportion),
        .mscaleB = *REAL(RmscaleB),
        .mscaleCC = *REAL(RmscaleCC),
        .mscaleMaxIt = *INTEGER(RmscaleMaxIt),
        .mscaleEPS = *REAL(RmscaleEPS),
        .mscaleRhoFun = *INTEGER(RmscaleRhoFun),

        /* We don't need the elastic net parameters */
        .lambda1 = 0,
        .lambda2 = 0,
        .enMaxIt = 0,
        .enEPS = 0
    };

    const int nobs = *INTEGER(Rnobs);
    const int nvar = *INTEGER(Rnvar);

    /*
     * CAUTION: If the QR factorization is used to compute the estimates, we need extra
     * storage space!
     */
    /* SEXP ret = PROTECT(Rf_allocVector(REALSXP, (3 * nvar + 2) * nvar + nobs)); */
    SEXP ret = PROTECT(Rf_allocVector(REALSXP, (3 * nvar + 2) * nvar));

    *INTEGER(RnumInits) = computeInitialEstimator(REAL(RXtr), REAL(Ry), nobs, nvar, &ctrl, REAL(ret));

    UNPROTECT(1);
    return ret;
}


/***************************************************************************************************
 *
 * Main function to compute the initial estimator
 *
 **************************************************************************************************/
static int computeInitialEstimator(const double *restrict Xtr, const double *restrict y,
                                   const int nobs, const int nvar, const Control* ctrl,
                                   double *restrict estimates)
{
    AuxMemory auxmem;
    double *restrict bestCoefEst;
    double *restrict currentEst;
    double *restrict currentXtr = (double*) malloc(nobs * nvar * sizeof(double));
    double *restrict currentY = (double*) malloc(nobs * sizeof(double));
    double *restrict filteredXtr = (double*) malloc(nobs * nvar * sizeof(double));
    double *restrict filteredY = (double*) malloc(nobs * sizeof(double));
    double *restrict pscs = (double*) malloc(nobs * nvar * sizeof(double));
    double minObjective, tmpObjective;
    double scaledThreshold = 0;
    int currentNobs = nobs;
    int iter = 0;
    int j;
    int numPSCs = 0;
    int linalgError = 0;
    double *restrict currentPSC;
    int filteredNobs;
    double diff, normPrevBest = 0, normBest = 0;
    RhoFunction rhoFun = getRhoFunctionByName(ctrl->mscaleRhoFun);

    initAuxMemory(&auxmem);
    resizeAuxMemory(&auxmem, nvar, nobs);

//    this->setData(this->originalData);

    memcpy(currentXtr, Xtr, nobs * nvar * sizeof(double));
    memcpy(currentY, y, nobs * sizeof(double));

    /* 
     * At first, we will leave the front of the estimates-matrix empty -- this is the
     * place where the best estimate of the previous run will be stored
     */
    bestCoefEst = estimates;
    currentEst = estimates + nvar;
    minObjective = DBL_MAX;

    while(1) {
        /* 1. Calculate PSC for current work data */
        numPSCs = calculatePSCs(pscs, &auxmem, currentXtr, currentY, currentNobs, nvar);

        if (numPSCs < 0) {
            linalgError = -numPSCs;
            break;
        }

        memcpy(filteredXtr, currentXtr, currentNobs * nvar * sizeof(double));
        memcpy(filteredY, currentY, currentNobs * sizeof(double));

        /* 2. Estimate coefficients for work data (this could alter workData!) */
        linalgError = estimateCoef(filteredXtr, filteredY, currentNobs, nvar, currentEst, &auxmem);

        if (linalgError != 0) {
            break;
        }

        // Now evaluate this->coefEst on the
        calculateResiduals(currentXtr, currentY, currentEst, currentNobs, nvar, auxmem.residuals);
        tmpObjective = mscale(auxmem.residuals, currentNobs, ctrl->mscaleB, ctrl->mscaleEPS,
                              ctrl->mscaleMaxIt, rhoFun, ctrl->mscaleCC);

        if (tmpObjective < minObjective) {
            bestCoefEst = currentEst;
            minObjective = tmpObjective;
        }

        for(j = 0; j < numPSCs; ++j) {
            currentPSC = pscs + currentNobs * j;
            /* 4.1. Thin out X and y based on large values of PSCs */
            scaledThreshold = getQuantile(currentPSC, currentNobs, ctrl->pscProportion,
                                          lessThan);


            filteredNobs = filterDataThreshold(currentXtr, currentY, filteredXtr, filteredY,
                                               currentNobs, nvar, currentPSC,
                                               scaledThreshold, lessThan);


            /* 4.2. Estimate coefficients */
            currentEst += nvar;
            linalgError = estimateCoef(filteredXtr, filteredY, filteredNobs, nvar, currentEst,
                                       &auxmem);
            if (linalgError != 0) {
                break;
            }

            calculateResiduals(currentXtr, currentY, currentEst, currentNobs, nvar,
                               auxmem.residuals);
            tmpObjective = mscale(auxmem.residuals, currentNobs, ctrl->mscaleB, ctrl->mscaleEPS,
                                  ctrl->mscaleMaxIt, rhoFun, ctrl->mscaleCC);

            if (tmpObjective < minObjective) {
                minObjective = tmpObjective;
                bestCoefEst = currentEst;
            }

            /* 4.1. Thin out X and y based on large values of PSCs */
            scaledThreshold = getQuantile(currentPSC, currentNobs, ctrl->pscProportion,
                                          greaterThan);


            filteredNobs = filterDataThreshold(currentXtr, currentY, filteredXtr, filteredY,
                                               currentNobs, nvar, currentPSC,
                                               scaledThreshold, greaterThan);


            /* 4.2. Estimate coefficients */
            currentEst += nvar;
            linalgError = estimateCoef(filteredXtr, filteredY, filteredNobs, nvar, currentEst,
                                       &auxmem);
            if (linalgError != 0) {
                break;
            }

            calculateResiduals(currentXtr, currentY, currentEst, currentNobs, nvar,
                               auxmem.residuals);
            tmpObjective = mscale(auxmem.residuals, currentNobs, ctrl->mscaleB, ctrl->mscaleEPS,
                                  ctrl->mscaleMaxIt, rhoFun, ctrl->mscaleCC);

            if (tmpObjective < minObjective) {
                minObjective = tmpObjective;
                bestCoefEst = currentEst;
            }


            /* 4.1. Thin out X and y based on large values of PSCs */
            scaledThreshold = getQuantile(currentPSC, currentNobs, ctrl->pscProportion,
                                          absoluteLessThan);


            filteredNobs = filterDataThreshold(currentXtr, currentY, filteredXtr, filteredY,
                                               currentNobs, nvar, currentPSC,
                                               scaledThreshold, absoluteLessThan);


            /* 4.2. Estimate coefficients */
            currentEst += nvar;

            linalgError = estimateCoef(filteredXtr, filteredY, filteredNobs, nvar, currentEst,
                                       &auxmem);
            if (linalgError != 0) {
                break;
            }

            calculateResiduals(currentXtr, currentY, currentEst, currentNobs, nvar,
                               auxmem.residuals);
            tmpObjective = mscale(auxmem.residuals, currentNobs, ctrl->mscaleB, ctrl->mscaleEPS,
                                  ctrl->mscaleMaxIt, rhoFun, ctrl->mscaleCC);

            if (tmpObjective < minObjective) {
                minObjective = tmpObjective;
                bestCoefEst = currentEst;
            }
        }

        if (linalgError != 0) {
            break;
        }

        /*
         * Check if we converged to an best coef estimate
         */
        diff = 0;
        normBest = 0;
        for (j = 0; j < nvar; ++j) {
            diff += fabs(bestCoefEst[j] - estimates[j]);
            normBest += fabs(bestCoefEst[j]);
        }

        /* 5. Store best estimate for later at the beginning of estimates */
        memcpy(estimates, bestCoefEst, nvar * sizeof(double));
        bestCoefEst = estimates;
        currentEst = estimates + nvar;

        /*
         * Check if we reached the maximum number of iterations
         */
        if ((++iter >= ctrl->numIt) || (diff < ctrl->eps * normPrevBest)) {
            break;
        }

        normPrevBest = normBest;

        /* 6. Calculate residuals with best coefficient estimate for all observatins */
        calculateResiduals(Xtr, y, bestCoefEst, nobs, nvar, auxmem.residuals);

        if (ctrl->residThreshold < 0) {
            scaledThreshold = getQuantile(auxmem.residuals, nobs, ctrl->residProportion,
                                          absoluteLessThan);
        } else {
            scaledThreshold = ctrl->residThreshold * minObjective;
        }

        currentNobs = filterDataThreshold(Xtr, y, currentXtr, currentY, nobs, nvar, auxmem.residuals,
                                          scaledThreshold, absoluteLessThan);
        /* Update minObjective! */
        calculateResiduals(currentXtr, currentY, bestCoefEst, currentNobs, nvar, auxmem.residuals);
        minObjective = mscale(auxmem.residuals, currentNobs, ctrl->mscaleB, ctrl->mscaleEPS,
                              ctrl->mscaleMaxIt, rhoFun, ctrl->mscaleCC);
    }


    freeAuxMemory(&auxmem);

    free(pscs);
    free(currentXtr);
    free(currentY);
    free(filteredXtr);
    free(filteredY);

    if (linalgError != 0) {
        Rf_error("There was an error in one of the calls to LINPACK (%d)", linalgError);
    }

    return 3 * numPSCs + 2;
}

/***************************************************************************************************
 *
 * Little helper function to calculate the residuals
 *
 **************************************************************************************************/
static void calculateResiduals(const double *restrict Xtr, const double *restrict y,
                               const double *restrict coefEst,
                               const int nobs, const int nvar, double *restrict residuals)
{
    memcpy(residuals, y, nobs * sizeof(double));
    BLAS_DGEMV(BLAS_TRANS_TRANS, nvar, nobs,
               BLAS_M1F, Xtr, nvar,
               coefEst, BLAS_1L,
               BLAS_1F, residuals, BLAS_1L);
}

/***************************************************************************************************
 *
 * Function to compute the LS estimate
 *
 **************************************************************************************************/
static int estimateCoef(double *restrict Xtr, const double *restrict y, const int nobs,
                        const int nvar, double *restrict coefEst, AuxMemory* auxmem)
{
    int lapackInfo;
    /* We won't handle rank-deficient cases here!! */

    /* Xsqrt = t(X) %*% X */
    BLAS_DGEMM(BLAS_TRANS_NO, BLAS_TRANS_TRANS, nvar, nvar, nobs,
        BLAS_1F, Xtr, nvar, Xtr, nvar,
        BLAS_0F, auxmem->Xsqrt, nvar);

    /* Xsqrt = chol(Xsqrt) */
    LAPACK_DPOTRF(BLAS_UPLO_UPPER, nvar, auxmem->Xsqrt, nvar, lapackInfo);

    if (lapackInfo != 0) {
        return 0;
    }

    /* coefEst = t(X) %*% y */
    BLAS_DGEMV(BLAS_TRANS_NO, nvar, nobs, BLAS_1F, Xtr, nvar, y, BLAS_1L, BLAS_0F, coefEst,
               BLAS_1L);

    /* coefEst = inv(t(chol(Xsqrt))) %*% coefEst */
    BLAS_DTRSV(BLAS_UPLO_UPPER, BLAS_TRANS_TRANS, BLAS_DIAG_NO, nvar, auxmem->Xsqrt, nvar, coefEst,
               BLAS_1L);

    /* coefEst = inv(chol(Xsqrt)) %*% coefEst */
    BLAS_DTRSV(BLAS_UPLO_UPPER, BLAS_TRANS_NO, BLAS_DIAG_NO, nvar, auxmem->Xsqrt, nvar, coefEst,
               BLAS_1L);

    /*
     * Below code uses a QR factorization with LAPACK's DGELS
     * CAUTION: This needs extra storage space in the `coefEst` variable!!
     */
//    memcpy(coefEst, y, nobs * sizeof(double));
//
//    /* First query for optimal working memory size */
//    LAPACK_DGELS(BLAS_TRANS_TRANS, nvar, nobs, BLAS_1L, Xtr, nvar, coefEst,
//                 nobs, auxmem->dblWorkMem, BLAS_M1L, lapackInfo);
//
//    resizeDblWorkAuxMemory(auxmem, (int) auxmem->dblWorkMem[0]);
//
//    /* Actually calculate LS solution */
//    LAPACK_DGELS(BLAS_TRANS_TRANS, nvar, nobs, BLAS_1L, Xtr, nvar, coefEst,
//                 nobs, auxmem->dblWorkMem, auxmem->dblWorkMemSize, lapackInfo);

    return lapackInfo;
}

/***************************************************************************************************
 *
 * Helper functions to filter data based on a threshold
 *
 **************************************************************************************************/
static int filterDataThreshold(const double *restrict Xtr, const double *restrict y,
                                double *restrict newXtr, double *restrict newY,
                                const int nobs, const int nvar,
                                const double *restrict values, const double threshold,
                                CompareFunction compare)
{
    int i, toCount = 0;
    double *restrict toYIt = newY;
    const double *restrict fromYIt = y;

    int copyRows = 0;

    const double *restrict startFromX = Xtr;
    double *restrict startToX = newXtr;

    for (i = 0; i < nobs; ++i, ++fromYIt) {
        if (compare(values[i], threshold) < 0) {
            /* Copy value from y to newY */
            *toYIt = *fromYIt;
            ++toYIt;
            ++toCount;
            ++copyRows;
        } else {
            /*
             * This row (column in Xtr) should not be copied.
             * So copy everything so far
             */
            if (copyRows > 0) {
                memcpy(startToX, startFromX, copyRows * nvar * sizeof(double));
                startToX += copyRows * nvar;
            }

            startFromX += (copyRows + 1) * nvar;
            copyRows = 0;
        }
    }

    /* 
     * Copy last chunk of data
     */
    if (copyRows > 0) {
        memcpy(startToX, startFromX, copyRows * nvar * sizeof(double));
        startToX += copyRows * nvar;
    }

    return toCount;
}

/***************************************************************************************************
 *
 * Static compare functions
 *
 **************************************************************************************************/
static double lessThan(const double a, const double b)
{
    return a - b;
}

static double greaterThan(const double a, const double b)
{
    return b - a;
}

static double absoluteLessThan(const double a, const double b)
{
    return fabs(a) - fabs(b);
}
