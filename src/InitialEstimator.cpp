//
//  InitialEstimator.cpp
//  penseinit
//
//  Created by David Kepplinger on 2016-01-26.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#include <Rcpp.h>
#include <cstring>
#include <algorithm>
#include <stdexcept>
#include <Rmath.h>

#include "BLAS.h"
#include "LapackException.hpp"
#include "InitialEstimator.hpp"
#include "PartialSort.h"
#include "mscale.h"


#define MAX_NUM_PSCS(numVar) (3 * numVar + 2)

static const int BLAS_1L = 1;
static const double BLAS_0F = 0.0;
static const double BLAS_1F = 1.0;
static const double BLAS_M1F = -1.0;
static const char * const BLAS_DIAG_NO = "N";
static const char * const BLAS_TRANS_NO = "N";
static const char * const BLAS_TRANS_TRANS = "T";
static const char * const BLAS_UPLO_UPPER = "U";

#ifdef DEBUG
void print_matf(int dr, int dc, const double *A, const char *header) {
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

/**
 * Factor that corrects lambda to match RR-LS with -SE
 * It is estimated numerically based on (known?) ratios of scale^2/E(r^2)
 * for r~N(mu,sig^2) and different values of mu, sig, and b
 * where scale is the robust bisquare scale, and b is its BDP
 */
static inline double facon(const double delta);

/*
 * Compare functions
 */
static double lessThan(const double a, const double b);
static double greaterThan(const double a, const double b);
static double absoluteLessThan(const double a, const double b);


/**
 * Simple function to help filtering from one data to another data object
 *
 * CAUTION:
 * The memory of data `to` must be able to at least store everything from `from`!
 */
static inline void doFiltering(const Data& from, Data& to, const double *const values,
                               const double threshold, CompareFunction compare,
                               int nobs = 0);

/**
 * source.numObs() must give the actual number of observations, not the number of rows!!
 */
static inline void extendData(Data &dest, const Data& source, const double lambda, bool copy);

/**
 * Compute the residuals for the given data and coefficient estimate
 *
 */
static inline void computeResiduals(const Data& data, const double *RESTRICT coefEst,
                                    double *RESTRICT residuals);

/***************************************************************************************************
 *
 * Parent class for all initial estimators
 *
 **************************************************************************************************/

InitialEstimator::InitialEstimator(const Data& originalData, const Control& ctrl,
                                   PSC& psc,
                                   const int coefMemIncrease) :
        originalData(originalData),
		ctrl(ctrl),
		psc(psc)
{
    if (this->originalData.numVar() > 0) {
        /*
         * We have 3 * nvar + 1 parameter estimates
         * Due to the need for intial estimators that use OLS
         * for some additional storage space,
         * the memory is enlarged by `coefMemIncrease`
         */
        this->allCoefEstimates = new double[MAX_NUM_PSCS(this->originalData.numVar()) *
                                            this->originalData.numVar()
                                            + coefMemIncrease];

        this->coefObjFunScore = new double[MAX_NUM_PSCS(this->originalData.numVar())];
    }

    this->residuals = new double[this->originalData.numObs()];

    this->rhoFun = getRhoFunctionByName(ctrl.mscaleRhoFun);
}


InitialEstimator::~InitialEstimator()
{
    this->residualFilteredData.free();

    delete[] this->residuals;

    if (this->originalData.numVar() > 0) {
        delete[] this->allCoefEstimates;
        delete[] this->coefObjFunScore;
    }
}


int InitialEstimator::compute()
{
    const int origNvar = this->originalData.numVar();
    const int origNobs = this->originalData.numObs();
    int iter = 0;
    int j;
    int numPSCs = 0;
    double *const minObjective = this->coefObjFunScore;
    double * tmpObjective;
    double *RESTRICT bestCoefEst;
    double threshold;
    double diff, normPrevBest = 0, normBest = 0;
    const double *RESTRICT currentPSC;

    this->setData(this->originalData);

    bestCoefEst = this->allCoefEstimates;
    this->coefEst = this->allCoefEstimates + origNvar;
    *minObjective = DBL_MAX;

    while(1) {
        tmpObjective = this->coefObjFunScore + 1;

        /* 1. Estimate coefficients for residuals-filtered data */
        this->estimateCoefficients();

        // Now evaluate this->coefEst on the
        *tmpObjective = this->evaluateEstimate();
        if (*tmpObjective < *minObjective) {
            *minObjective = *tmpObjective;
            bestCoefEst = this->coefEst;
        }
        ++tmpObjective;

        /* 2. Calculate PSC for current work data */
        this->psc.setData(this->residualFilteredData);
        this->psc.setResiduals(this->residuals);
        numPSCs = this->psc.computePSC();

        for(j = 0; j < numPSCs; ++j) {
            currentPSC = this->psc.getPSC() + this->residualFilteredData.numObs() * j;
            /* 4.1. Thin out X and y based on large values of PSCs */
            threshold = getQuantile(currentPSC, this->residualFilteredNobs,
                                    this->ctrl.pscProportion, lessThan);

            this->filterDataPSC(currentPSC, threshold, lessThan);

            /* 4.2. Estimate coefficients */
            this->coefEst += origNvar;
            this->estimateCoefficients();
            *tmpObjective = this->evaluateEstimate();
            
            if (*tmpObjective < *minObjective) {
                *minObjective = *tmpObjective;
                bestCoefEst = this->coefEst;
            }
            ++tmpObjective;

            /* 4.1. Thin out X and y based on large values of PSCs */
            threshold = getQuantile(currentPSC, this->residualFilteredNobs,
                                    this->ctrl.pscProportion, greaterThan);

            this->filterDataPSC(currentPSC, threshold, greaterThan);

            /* 4.2. Estimate coefficients */
            this->coefEst += origNvar;
            this->estimateCoefficients();
            *tmpObjective = this->evaluateEstimate();

            if (*tmpObjective < *minObjective) {
                *minObjective = *tmpObjective;
                bestCoefEst = this->coefEst;
            }
            ++tmpObjective;

            /* 4.1. Thin out X and y based on large values of PSCs */
            threshold = getQuantile(currentPSC, this->residualFilteredNobs,
                                    this->ctrl.pscProportion, absoluteLessThan);

            this->filterDataPSC(currentPSC, threshold, absoluteLessThan);

            /* 4.2. Estimate coefficients */
            this->coefEst += origNvar;
            this->estimateCoefficients();
            *tmpObjective = this->evaluateEstimate();

            if (*tmpObjective < *minObjective) {
                *minObjective = *tmpObjective;
                bestCoefEst = this->coefEst;
            }
            ++tmpObjective;
        }

        /*
         * Check if we converged to the previous best coef estimate
         */
        diff = 0;
        normBest = 0;
        for (j = 0; j < origNvar; ++j) {
            diff += fabs(bestCoefEst[j] - this->allCoefEstimates[j]);
            normBest += fabs(bestCoefEst[j]);
        }

        /* 5. Store best estimate for later at the beginning of estimates */
        memcpy(this->allCoefEstimates, bestCoefEst, origNvar * sizeof(double));
        bestCoefEst = this->allCoefEstimates;

        /* Check for convergence or if we hit the maximum number of iterations */
        if ((++iter >= this->ctrl.numIt) || (diff < this->ctrl.eps * normPrevBest)) {
            break;
        }

        normPrevBest = normBest;

        /* 5. Calculate ALL residuals with best coefficient estimate */
        this->coefEst = bestCoefEst;
        computeResiduals(this->originalData, this->coefEst, this->residuals);

        if (this->ctrl.residThreshold < 0) {
            threshold = getQuantile(this->residuals, origNobs,
                                    this->ctrl.residProportion, absoluteLessThan);
        } else {
            threshold = this->ctrl.residThreshold * this->MscaleOfResiduals();
        }

        this->filterDataResiduals(threshold);

        /*
         * Now calculate the new minimum objective for the best estimate on the
         * filtered data
         */
        computeResiduals(this->residualFilteredData, this->coefEst, this->residuals);
        *minObjective = this->evaluateEstimate();

        this->coefEst = this->allCoefEstimates + origNvar;
    }

    return 3 * numPSCs + 2;
}

void InitialEstimator::setData(const Data &data)
{
    this->residualFilteredData.copy(data);
    this->residualFilteredNobs = data.numObs();
}

void InitialEstimator::filterDataResiduals(double threshold)
{
    doFiltering(this->originalData, this->residualFilteredData, this->residuals,
                threshold, absoluteLessThan);

    this->residualFilteredNobs = this->residualFilteredData.numObs();
}

double InitialEstimator::MscaleOfResiduals() const
{
    return mscale(this->residuals, this->residualFilteredNobs, this->ctrl.mscaleB,
                  this->ctrl.mscaleEPS, this->ctrl.mscaleMaxIt, this->rhoFun,
                  this->ctrl.mscaleCC);
}


/***************************************************************************************************
 *
 * OLS Initial estimator
 * This is the original one by Peña and Yohai
 *
 **************************************************************************************************/
OLS::OLS(const Data& originalData, const Control& ctrl) :
            InitialEstimator(originalData, ctrl, pscOls,
                             originalData.numObs() - originalData.numVar()),
            dataToUse(pscFilteredData)
{
    this->XtX = new double[originalData.numVar() * originalData.numVar()];
    this->pscOls.setXsqrtMemory(this->XtX);
}


OLS::~OLS()
{
    this->pscFilteredData.free();
    delete[] this->XtX;
}

void OLS::setData(const Data &data)
{
    this->pscFilteredData.free();
    this->pscFilteredData.setNumObs(data.numObs());
    this->pscFilteredData.setNumVar(data.numVar());
    this->pscFilteredData.resize();

    InitialEstimator::setData(data);
    dataToUse = this->residualFilteredData;
}

void OLS::filterDataResiduals(double threshold)
{
    InitialEstimator::filterDataResiduals(threshold);
    dataToUse = this->residualFilteredData;
}

void OLS::filterDataPSC(const double *RESTRICT values, const double threshold,
                        CompareFunction compare)
{
    doFiltering(this->residualFilteredData, this->pscFilteredData, values, threshold,
                compare);
    dataToUse = this->pscFilteredData;
}


void OLS::estimateCoefficients()
{
    const int nobs = this->dataToUse.numObs();
    const int nvar = this->dataToUse.numVar();
    int lapackInfo;

    /* XtX = t(X) %*% X */
    BLAS_DGEMM(BLAS_TRANS_NO, BLAS_TRANS_TRANS, nvar, nvar, nobs,
        BLAS_1F, this->dataToUse.getXtr(), nvar, this->dataToUse.getXtr(), nvar,
        BLAS_0F, this->XtX, nvar);

    /* XtX = chol(XtX) */
    LAPACK_DPOTRF(BLAS_UPLO_UPPER, nvar, this->XtX, nvar, lapackInfo);

    if (lapackInfo != 0) {
        throw LapackException("Could not solve least-squares equation", lapackInfo);
    }

    /* coefEst = t(X) %*% y */
    BLAS_DGEMV(BLAS_TRANS_NO, nvar, nobs, BLAS_1F, this->dataToUse.getXtr(), nvar,
               this->dataToUse.getY(), BLAS_1L, BLAS_0F, coefEst, BLAS_1L);

    /* coefEst = inv(t(chol(XtX))) %*% coefEst */
    BLAS_DTRSV(BLAS_UPLO_UPPER, BLAS_TRANS_TRANS, BLAS_DIAG_NO, nvar, this->XtX, nvar,
               this->coefEst, BLAS_1L);

    /* coefEst = inv(chol(XtX)) %*% coefEst */
    BLAS_DTRSV(BLAS_UPLO_UPPER, BLAS_TRANS_NO, BLAS_DIAG_NO, nvar, this->XtX, nvar,
               this->coefEst, BLAS_1L);

    /* Now update the residuals */
    computeResiduals(this->residualFilteredData, this->coefEst, this->residuals);
}

double OLS::evaluateEstimate() const
{
    return this->MscaleOfResiduals();
}

/***************************************************************************************************
 *
 * Elastic Net - Peña Yohai Initial estimator
 *
 **************************************************************************************************/
ENPY::ENPY(const Data& originalData, const Control& ctrl) :
           InitialEstimator(originalData, ctrl, pscOls),
           /* lambda1 for EN should not depend on N, but the given one does, so divide! */
           lambda1LS(ctrl.lambda1 / (facon(ctrl.mscaleB) * originalData.numObs())),
           lambda2LS(ctrl.lambda2 / facon(ctrl.mscaleB)),
           en(ctrl.enMaxIt, ctrl.enEPS, ctrl.enCentering),
           dataToUse(residualFilteredData)
{
    delete[] this->residuals;
    this->residuals = new double[this->originalData.numObs() + this->originalData.numVar()];
    if (this->lambda1LS == 0) {
        this->XtX = new double[this->originalData.numVar() * this->originalData.numVar()];
    }
}

ENPY::~ENPY()
{
    this->pscFilteredData.free();
    if (this->lambda1LS == 0) {
        delete[] this->XtX;
    }
}

void ENPY::estimateCoefficients()
{
    const int nobs = this->dataToUse.numObs();
    const int nvar = this->dataToUse.numVar();
    bool converged;
    int lapackInfo;
    int j;
    const double minusSqrtLambda2LS = -sqrt((this->lambda2LS * this->residualFilteredNobs) /
                                                this->originalData.numObs());

    if (this->lambda1LS == 0) {
        /* lambda1 is zero, i.e. we can do a simple OLS fit to the augmented X and Y */

        /* XtX = t(X) %*% X */
        BLAS_DGEMM(BLAS_TRANS_NO, BLAS_TRANS_TRANS, nvar, nvar, nobs,
            BLAS_1F, this->dataToUse.getXtr(), nvar, this->dataToUse.getXtr(), nvar,
            BLAS_0F, this->XtX, nvar);

        /* XtX = chol(XtX) */
        LAPACK_DPOTRF(BLAS_UPLO_UPPER, nvar, this->XtX, nvar, lapackInfo);

        if (lapackInfo != 0) {
            throw LapackException("Could not solve least-squares equation", lapackInfo);
        }

        /* coefEst = t(X) %*% y */
        BLAS_DGEMV(BLAS_TRANS_NO, nvar, nobs, BLAS_1F, this->dataToUse.getXtr(), nvar,
                   this->dataToUse.getY(), BLAS_1L, BLAS_0F, coefEst, BLAS_1L);

        /* coefEst = inv(t(chol(XtX))) %*% coefEst */
        BLAS_DTRSV(BLAS_UPLO_UPPER, BLAS_TRANS_TRANS, BLAS_DIAG_NO, nvar, this->XtX, nvar,
                   this->coefEst, BLAS_1L);

        /* coefEst = inv(chol(XtX)) %*% coefEst */
        BLAS_DTRSV(BLAS_UPLO_UPPER, BLAS_TRANS_NO, BLAS_DIAG_NO, nvar, this->XtX, nvar,
                   this->coefEst, BLAS_1L);
    } else {
        /* This already updates the residuals, but NOT for all observations */
        /* lambda1LS must not be adjusted for the number of observations! */
        this->en.setLambdas(this->lambda1LS, 0);
        converged = this->en.computeCoefs(this->dataToUse, this->coefEst, this->residuals);

        if (!converged) {
            throw std::runtime_error("LASSO did not converge. Either increase the number of "
                                     "iterations or lambda1.");
        }
    }

    /* Now update the residuals */
    computeResiduals(this->residualFilteredData, this->coefEst, this->residuals);

    if (this->lambda2LS > 0) {
        /* The last elements of the residuals are just - sqrt(lambda2) * beta[j] */
        for (j = 1; j < this->originalData.numVar(); ++j) {
            this->residuals[this->originalData.numObs() + j - 1] = minusSqrtLambda2LS * this->coefEst[j];
        }
    }
}

void ENPY::setData(const Data &data)
{
    this->residualFilteredNobs = data.numObs();

    this->residualFilteredData.free();
    this->residualFilteredData.setNumVar(data.numVar());

    if (this->lambda2LS > 0) {
        this->residualFilteredData.setNumObs(data.numObs() + data.numVar() - 1);
        this->residualFilteredData.resize();
        extendData(this->residualFilteredData, data, this->lambda2LS, true);
    } else {
        this->residualFilteredData.free();
        this->residualFilteredData.copy(data);
    }

    this->pscFilteredData.free();
    this->pscFilteredData.setNumObs(this->residualFilteredData.numObs());
    this->pscFilteredData.setNumVar(this->residualFilteredData.numVar());
    this->pscFilteredData.resize();
    this->dataToUse = this->residualFilteredData;
}


void ENPY::filterDataResiduals(double threshold)
{
    InitialEstimator::filterDataResiduals(threshold);

    if (this->lambda2LS > 0) {
        extendData(this->residualFilteredData, this->residualFilteredData,
                   (this->lambda2LS * this->residualFilteredNobs) / this->originalData.numObs(),
                   false);

        this->residualFilteredData.setNumObs(this->residualFilteredNobs +
                                             this->residualFilteredData.numVar() - 1);
    }

    dataToUse = this->residualFilteredData;
}

void ENPY::filterDataPSC(const double *RESTRICT values, const double threshold,
                         CompareFunction compare)
{
    doFiltering(this->residualFilteredData, this->pscFilteredData, values, threshold, compare,
                this->residualFilteredNobs);

    if (this->lambda2LS > 0) {
        extendData(this->pscFilteredData, this->pscFilteredData,
                   (this->lambda2LS * this->pscFilteredData.numObs()) / this->originalData.numObs(),
                   false);

        this->pscFilteredData.setNumObs(this->pscFilteredData.numObs() +
                                        this->originalData.numVar() - 1);
    }

    dataToUse = this->pscFilteredData;
}

double ENPY::evaluateEstimate() const
{
    int j;
    double penalty = 0;
    double scale = this->MscaleOfResiduals();

    for (j = 1; j < this->originalData.numVar(); ++j) {
        penalty += this->ctrl.lambda2 * this->coefEst[j] * this->coefEst[j] +
                   this->ctrl.lambda1 * fabs(this->coefEst[j]);
    }

    penalty *= ((double) this->residualFilteredNobs / (double) this->originalData.numObs());

    return this->residualFilteredNobs * scale * scale + penalty;
}

/***************************************************************************************************
 *
 * Ridge Regression Initial estimator
 * This is the one proposed by Maronna et al.
 *
 **************************************************************************************************/

//RRPY::RRPY(const Data& originalData, const Control& ctrl) :
//        InitialEstimator(originalData, ctrl, pscOls)
//{}
//
//RRPY::~RRPY()
//{
//    this->unextendedDataBackup.free();
//    this->unextendedData.free();
//}
//
//void RRPY::extendData(Data &dest, const Data &source)
//{
//    int j;
//    int srcNumObsBytes = source.numObs() * sizeof(double);
//    double *RESTRICT destPtr;
//    double *RESTRICT srcPtr;
//    double sqrtLambda = sqrt(this->ctrl.lambda1);
//
//    destPtr = dest.getX();
//    srcPtr = source.getX();
//
//    memset(destPtr, 0, dest.numObs() * dest.numVar() * sizeof(double));
//
//    /**
//     * First variable is 1's in source
//     * this will be followed by 0's in the extended column
//     */
//    memset(destPtr, 1, srcNumObsBytes);
//
//    destPtr += dest.numObs();
//    srcPtr += source.numObs();
//
//    for (j = 1; j < source.numVar(); ++j) {
//        memcpy(destPtr, srcPtr, srcNumObsBytes);
//        *(destPtr + j - 1) = sqrtLambda;
//        destPtr += dest.numObs();
//        srcPtr += source.numObs();
//    }
//
//    destPtr = dest.getY();
//    memcpy(destPtr, source.getY(), srcNumObsBytes);
//    memset(destPtr + source.numObs(), 0, source.numVar() * sizeof(double));
//}
//
//void RRPY::setData(const Data &data)
//{
//    /* Copy data into the unextended data and its backup */
//    this->unextendedData.free();
//    this->unextendedData.copy(data);
//
//    this->unextendedDataBackup.free();
//    this->unextendedDataBackup.copy(data);
//
//    /* Now extend the data for use in the PSC_OLS */
//    this->workData.free();
//    this->workData.setNumObs(data.numObs() + data.numVar());
//    this->workData.setNumVar(data.numVar());
//    this->workData.resize();
//    this->extendData(this->workData, data);
//}
//
//void RRPY::filterDataResiduals(double threshold)
//{
//    /* Filtered the unextended data by residuals */
//    doFiltering(this->originalData, this->unextendedData, this->getResiduals(),
//                IndexNumber(&threshold, 0, IndexNumber::ABSOLUTE_INCREASING));
//
//    this->currentUsefulNobs = this->unextendedData.numObs();
//
//    /* Backup the unextended data */
//    this->unextendedDataBackup.overrideMemory(this->unextendedData);
//
//    /* No that we filtered the unextended data, we have to extended this filtered data */
//    this->workData.setNumObs(this->unextendedData.numObs() + this->unextendedData.numVar());
//    this->workData.setNumVar(this->unextendedData.numVar());
//    this->extendData(this->workData, this->unextendedData);
//}
//
//void RRPY::filterDataPSC(std::vector<IndexNumber> &ordering)
//{
//    std::vector<IndexNumber>::iterator threshIt = ordering.begin() + this->workData.numObs();
//    std::partial_sort(ordering.begin(), threshIt + 1, ordering.end());
//
//    doFiltering(this->unextendedDataBackup, this->unextendedData, threshIt->getValues(), *threshIt);
//}
//
//
//void RRPY::estimateCoefficients()
//{
//    const int nobs = this->unextendedData.numObs();
//    const int nvar = this->unextendedData.numVar();
//    int i;
//    int lapackInfo;
//
//    double *RESTRICT K = new double[nobs * nobs];
//    double *RESTRICT alphavec = new double[nobs];
//    memset(K, 0, nvar * nvar * sizeof(double));
//
//    for (i = 0; i < nobs * nobs; i += nvar + 1) {
//        K[i] = this->ctrl.lambda1 * nobs;
//    }
//
//    BLAS_DGEMM(BLAS_TRANS_NO, BLAS_TRANS_TRANS, nobs, nobs, nvar,
//               BLAS_1F, this->unextendedData.getX(), nobs,
//               this->unextendedData.getX(), nvar,
//               BLAS_1F, K, nobs);
//
//    LAPACK_DPOTRF(BLAS_UPLO_UPPER, nvar, K, nvar, lapackInfo);
//
//    if (lapackInfo != 0) {
//        throw LapackException("Could not get Cholesky decomposition of X.t(X)", lapackInfo);
//    }
//
//    memcpy(alphavec, this->unextendedData.getY(), nobs * sizeof(double));
//
//    BLAS_DTRSV(BLAS_UPLO_UPPER, BLAS_TRANS_NO, BLAS_DIAG_NO, nobs, K, nobs, alphavec, BLAS_1L);
//
//    BLAS_DGEMV(BLAS_TRANS_TRANS, nobs, nvar,
//               BLAS_1F, this->unextendedData.getX(), nvar,
//               alphavec, BLAS_1L, BLAS_0F, this->coefEst, BLAS_1L);
//
//}
//
//double RRPY::evaluateEstimate()
//{
//    #warning Should we add the penalty here or not??
//    this->computeResiduals();
//    return this->MscaleOfResiduals();
//}


/**
 * Calculate lambda correction factor
 */
static inline double facon(const double delta) {
    return 23.9716 - 73.4391 * delta + 64.9480 * delta * delta;
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

/***************************************************************************************************
 *
 * Helper function to calculate the residuals
 *
 **************************************************************************************************/
static inline void computeResiduals(const Data& data, const double *RESTRICT coefEst,
                                    double *RESTRICT residuals)
{
    const int nvar = data.numVar();
    const int nobs = data.numObs();

    memcpy(residuals, data.getYConst(), nobs * sizeof(double));
    BLAS_DGEMV(BLAS_TRANS_TRANS, nvar, nobs,
               BLAS_M1F, data.getXtrConst(), nvar,
               coefEst, BLAS_1L,
               BLAS_1F, residuals, BLAS_1L);
}

/***************************************************************************************************
 *
 * Helper function to filter a data set
 *
 **************************************************************************************************/
static inline void doFiltering(const Data& from, Data& to, const double *const values,
                               const double threshold, CompareFunction compare, int nobs)
{
    int i, toCount = 0;
    double *RESTRICT toYIt = to.getY();
    const double *RESTRICT fromYIt = from.getYConst();

    int copyRows = 0;

    const double *RESTRICT startFromX = from.getXtrConst();
    double *RESTRICT startToX = to.getXtr();

    if (nobs <= 0) {
        nobs = from.numObs();
    }

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
                memcpy(startToX, startFromX, copyRows * from.numVar() * sizeof(double));
                startToX += copyRows * from.numVar();
            }

            startFromX += (copyRows + 1) * from.numVar();
            copyRows = 0;
        }
    }

    /* 
     * Copy last chunk of data
     */
    if (copyRows > 0) {
        memcpy(startToX, startFromX, copyRows * from.numVar() * sizeof(double));
        startToX += copyRows * from.numVar();
    }

    to.setNumObs(toCount);
}


/***************************************************************************************************
 *
 * Helper function to create the extended data set for ridge regression
 *
 **************************************************************************************************/
/**
 * source.numObs() should give the actual number of observations, not the number of rows!!
 */
static inline void extendData(Data &dest, const Data& source, const double lambda, bool copy)
{
    int j;
    double sqrtLambda = sqrt(lambda);
    double *RESTRICT destPtr = dest.getXtr();
    const double *RESTRICT srcPtr = source.getXtrConst();

    if (copy) {
        memcpy(destPtr, srcPtr, source.numObs() * source.numVar() * sizeof(double));
        memcpy(dest.getY(), source.getYConst(), source.numObs() * sizeof(double));
    }

    destPtr += source.numObs() * source.numVar();

    memset(dest.getY() + source.numObs(), 0, (dest.numVar() - 1) * sizeof(double));

    /* 
     * nvar * (nvar - 1) because nvar includes the intercept!
     */
    memset(destPtr, 0, source.numVar() * (source.numVar() - 1) * sizeof(double));

    for (j = 1, ++destPtr; j < source.numVar(); ++j, destPtr += source.numVar() + 1) {
        *destPtr = sqrtLambda;
    }
}
