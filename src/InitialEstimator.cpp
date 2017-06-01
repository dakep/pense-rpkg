//
//  InitialEstimator.cpp
//  pense
//
//  Created by David Kepplinger on 2016-01-26.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#include <RcppArmadillo.h>
#include <cstring>
#include <algorithm>
#include <stdexcept>
#include <Rmath.h>

#include "config.h"
#include "BLAS.h"
#include "LapackException.hpp"
#include "InitialEstimator.hpp"
#include "PartialSort.h"
#include "mscale.h"


#define MAX_NUM_PSCS(numVar) (3 * numVar + 2)

#define LAMBDA_1(lambda, alpha) lambda * alpha
#define LAMBDA_2(lambda, alpha, nobs) 0.5 * lambda * (1 - alpha) * nobs

/**
 * default options definitions
 */
static const RhoFunctionName DEFAULT_OPT_RHO_FUNCTION = BISQUARE;
static const int DEFAULT_OPT_MAXIT = 10;
static const double DEFAULT_OPT_RESID_THRESHOLD = 2;
static const double DEFAULT_OPT_RESID_PROPORTION = .8;
static const double DEFAULT_OPT_PSC_PROPORTION = .4;
static const double DEFAULT_OPT_LAMBDA = 1;
static const double DEFAULT_OPT_ALPHA = 1;
static const double DEFAULT_OPT_EPS = 1e-9;
static const double DEFAULT_OPT_MSCALE_B = .5;
static const double DEFAULT_OPT_MSCALE_CC = 1.5432;
static const int DEFAULT_OPT_MSCALE_MAXIT = 100;
static const double DEFAULT_OPT_MSCALE_EPS = 1e-9;

/**
 * BLAS constants definitions
 */
static BLAS_INT BLAS_1L = 1;

static const double BLAS_0F = 0.0;
static const double BLAS_1F = 1.0;
static const double BLAS_M1F = -1.0;

static BLAS_CHAR BLAS_DIAG_NO = "N";
static BLAS_CHAR BLAS_TRANS_NO = "N";
static BLAS_CHAR BLAS_TRANS_TRANS = "T";
static BLAS_CHAR BLAS_UPLO_UPPER = "U";

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

InitialEstimator::InitialEstimator(const Data& originalData, const Options& opts,
                                   PSC& psc, const int maxEstimators,
                                   const int coefMemIncrease) :
        originalData(originalData),
		ctrl(InitialEstimator::initIEControl(opts)),
		psc(psc)
{
    if (this->originalData.numVar() > 0) {
        /*
         * We have 3 * nvar + 1 parameter estimates
         * Due to the need for intial estimators that use OLS
         * for some additional storage space,
         * the memory is enlarged by `coefMemIncrease`
         */
        this->allCoefEstimates = new double[maxEstimators * this->originalData.numVar()
                                            + coefMemIncrease];

        this->coefObjFunScore = new double[maxEstimators];
    }

    this->residuals = new double[this->originalData.numObs()];

    this->rhoFun = getRhoFunctionByName(opts.get("rhoFunction", DEFAULT_OPT_RHO_FUNCTION));
}

InitialEstimator::IEControl InitialEstimator::initIEControl(const Options& opts)
{
    IEControl ctrl = {
        opts.get("maxIt", DEFAULT_OPT_MAXIT),
        opts.get("residThreshold", DEFAULT_OPT_RESID_THRESHOLD),
        opts.get("residProportion", DEFAULT_OPT_RESID_PROPORTION),
        opts.get("pscProportion", DEFAULT_OPT_PSC_PROPORTION),
        opts.get("lambda", DEFAULT_OPT_LAMBDA),
        opts.get("alpha", DEFAULT_OPT_ALPHA),
        opts.get("eps", DEFAULT_OPT_EPS),
        opts.get("mscaleB", DEFAULT_OPT_MSCALE_B),
        opts.get("mscaleCC", DEFAULT_OPT_MSCALE_CC),
        opts.get("mscaleMaxIt", DEFAULT_OPT_MSCALE_MAXIT),
        opts.get("mscaleEPS", DEFAULT_OPT_MSCALE_EPS)
    };
    return ctrl;
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
    int numPSCs = 0, newNumPSCs = 0;
    double *const minObjective = this->coefObjFunScore;
    double * tmpObjective;
    double *RESTRICT bestCoefEst;
    double threshold;
    double diff, normPrevBest = 0, normBest = 0;
    const double *RESTRICT currentPSC;

    this->resetData();

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
        computeResiduals(this->residualFilteredData, this->coefEst, this->residuals);
        this->psc.setData(this->residualFilteredData);
        this->psc.setResiduals(this->residuals);
        newNumPSCs = this->psc.computePSC();

        if (newNumPSCs < 2) {
            /*
             * Either no PSCs or only one have been found.
             * A single PSC is observed when all coefficients
             * are all zero.
             */
             break;
        }

        numPSCs = newNumPSCs;

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
        if ((++iter >= this->ctrl.maxIt) || (diff < this->ctrl.eps * normPrevBest)) {
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

        this->coefEst = this->allCoefEstimates + origNvar;
    }

    return 3 * numPSCs + 2;
}

void InitialEstimator::resetData()
{
    this->residualFilteredData.copy(this->originalData);
    this->residualFilteredNobs = this->originalData.numObs();
}

void InitialEstimator::filterDataResiduals(double threshold)
{
    doFiltering(this->originalData, this->residualFilteredData, this->residuals,
                threshold, absoluteLessThan);

    this->residualFilteredNobs = this->residualFilteredData.numObs();
}

double InitialEstimator::MscaleOfResiduals() const
{
    return mscale(this->residuals, this->originalData.numObs(), this->ctrl.mscaleB,
                  this->ctrl.mscaleEPS, this->ctrl.mscaleMaxIt, this->rhoFun,
                  this->ctrl.mscaleCC);
}


/***************************************************************************************************
 *
 * OLS Initial estimator
 * This is the original one by Peña and Yohai
 *
 **************************************************************************************************/
IEOls::IEOls(const Data& originalData, const Options& ctrl) :
             InitialEstimator(originalData, ctrl, pscOls, MAX_NUM_PSCS(originalData.numVar()),
                              originalData.numObs() - originalData.numVar()),
            dataToUse(pscFilteredData)
{
    this->XtX = new double[originalData.numVar() * originalData.numVar()];
    this->pscOls.setXsqrtMemory(this->XtX);
}


IEOls::~IEOls()
{
    this->pscFilteredData.free();
    delete[] this->XtX;
}

void IEOls::resetData()
{
    this->pscFilteredData.free();
    this->pscFilteredData.setNumObs(this->originalData.numObs());
    this->pscFilteredData.setNumVar(this->originalData.numVar());
    this->pscFilteredData.resize();

    InitialEstimator::resetData();
    dataToUse = this->residualFilteredData;
}

void IEOls::filterDataResiduals(double threshold)
{
    InitialEstimator::filterDataResiduals(threshold);
    dataToUse = this->residualFilteredData;
}

void IEOls::filterDataPSC(const double *RESTRICT values, const double threshold,
                        CompareFunction compare)
{
    doFiltering(this->residualFilteredData, this->pscFilteredData, values, threshold,
                compare);
    dataToUse = this->pscFilteredData;
}


void IEOls::estimateCoefficients()
{
    BLAS_INT nobs = this->dataToUse.numObs();
    BLAS_INT nvar = this->dataToUse.numVar();
    int lapackInfo;

    if (nobs == 0) {
        memset(this->coefEst, 0, nvar * sizeof(double));
        memcpy(this->residuals, this->originalData.getYConst(),
               this->originalData.numObs() * sizeof(double));
        return;
    }

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
    computeResiduals(this->originalData, this->coefEst, this->residuals);
}

double IEOls::evaluateEstimate() const
{
    return this->MscaleOfResiduals();
}



/***************************************************************************************************
 *
 * Elastic Net - Peña Yohai Initial estimator
 *
 **************************************************************************************************/
ENPY::ENPY(const Data& originalData, const Options& opts, const Options& enOpts) :
           InitialEstimator(originalData, opts, pscOls, MAX_NUM_PSCS(originalData.numVar())),

           lambdaLS(this->ctrl.lambda * 0.5),

           /* Select correct ElasticNet implementation */
           en(*getElasticNetImpl(enOpts, true)),

           dataToUse(residualFilteredData)
{
    /* Resize the residuals memory allocated from by base class */
    delete[] this->residuals;
    this->residuals = new double[this->originalData.numObs() + this->originalData.numVar()];

    if (this->ctrl.alpha == 0) {
        this->XtX = new double[this->originalData.numVar() * this->originalData.numVar()];
    }

    this->en.setAlphaLambda(1, LAMBDA_1(this->ctrl.alpha, this->lambdaLS));
}

ENPY::~ENPY()
{
    delete &this->en;
    this->pscFilteredData.free();
    if (this->ctrl.alpha == 0) {
        delete[] this->XtX;
    }
}

void ENPY::estimateCoefficients()
{
    BLAS_INT nobs = this->dataToUse.numObs();
    BLAS_INT nvar = this->dataToUse.numVar();
    int lapackInfo;
    int j;
    const double minusSqrtLambda2LS = -(1 + (this->ctrl.alpha > 0)) *
                                        sqrt(LAMBDA_2(this->lambdaLS, this->ctrl.alpha,
                                                      this->residualFilteredNobs));

    if (nobs == 0) {
        memset(this->coefEst, 0, nvar * sizeof(double));
        memcpy(this->residuals, this->originalData.getYConst(),
               this->originalData.numObs() * sizeof(double));
        return;
    }

    if (this->ctrl.alpha == 0) {
        /*
         * alpha is zero, i.e. we can do a simple OLS fit to the augmented X and Y
         */
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
        /*
         * We actually need to adjust lambda for the number of observations as
         * we operate on the *augmented* data, in order to get the same results
         * as from the elastic net on the true data
         */
        this->en.setAlphaLambda(1, LAMBDA_1(this->ctrl.alpha, this->lambdaLS) *
                                        this->residualFilteredNobs /
                                        (residualFilteredNobs + this->originalData.numVar() - 1));

        this->en.setData(this->dataToUse);
        this->en.computeCoefs(this->coefEst, this->residuals);

        if (this->en.getStatus() > 0) {
            Rcpp::warning("LASSO did not converge. Either increase the number of "
                          "iterations or lambda.");
        }
    }

    /* Now update the residuals */
    computeResiduals(this->originalData, this->coefEst, this->residuals);

    if (this->ctrl.alpha < 1) {
        /* The last elements of the residuals are just - sqrt(lambda2) * beta[j] */
        for (j = 1; j < this->originalData.numVar(); ++j) {
            this->residuals[this->originalData.numObs() + j - 1] = minusSqrtLambda2LS *
                                                                   this->coefEst[j];
        }
    }
}

void ENPY::resetData()
{
    this->residualFilteredNobs = this->originalData.numObs();

    this->residualFilteredData.free();
    this->residualFilteredData.setNumVar(this->originalData.numVar());

    if (this->ctrl.alpha < 1) {
        this->residualFilteredData.setNumObs(this->originalData.numObs() + this->originalData.numVar() - 1);
        this->residualFilteredData.resize();
        /* In case of EN, we need to multiply lambda_2 by 2 to account for the 1/2 in the EN objective */
        extendData(this->residualFilteredData, this->originalData,
                   (1 + (this->ctrl.alpha > 0)) *
                   LAMBDA_2(this->lambdaLS, this->ctrl.alpha, this->originalData.numObs()), true);
    } else {
        this->residualFilteredData.free();
        this->residualFilteredData.copy(this->originalData);
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

    if (this->ctrl.alpha < 1) {
        /* In case of EN, we need to multiply lambda_2 by 2 to account for the 1/2 in the EN objective */
        extendData(this->residualFilteredData, this->residualFilteredData,
                   (1 + (this->ctrl.alpha > 0)) *
                   LAMBDA_2(this->lambdaLS, this->ctrl.alpha, this->residualFilteredNobs),
                   false);

        this->residualFilteredData.setNumObs(this->residualFilteredNobs +
                                             this->residualFilteredData.numVar() - 1);
    }

    this->dataToUse = this->residualFilteredData;

}

void ENPY::filterDataPSC(const double *RESTRICT values, const double threshold,
                         CompareFunction compare)
{
    doFiltering(this->residualFilteredData, this->pscFilteredData, values, threshold, compare,
                this->residualFilteredNobs);

    if (this->ctrl.alpha < 1) {
        /* In case of EN, we need to multiply lambda_2 by 2 to account for the 1/2 in the EN objective */
        extendData(this->pscFilteredData, this->pscFilteredData,
                   (1 + (this->ctrl.alpha > 0)) *
                   LAMBDA_2(this->lambdaLS, this->ctrl.alpha, this->pscFilteredData.numObs()),
                   false);

        this->pscFilteredData.setNumObs(this->pscFilteredData.numObs() +
                                        this->originalData.numVar() - 1);
    }

    this->dataToUse = this->pscFilteredData;
}

double ENPY::evaluateEstimate() const
{
    int j;
    double penalty = 0;
    double tmp = 0.5 * (1 - this->ctrl.alpha);

    for (j = 1; j < this->originalData.numVar(); ++j) {
        penalty += this->ctrl.lambda * (
                        tmp * this->coefEst[j] * this->coefEst[j] +
                        this->ctrl.alpha * fabs(this->coefEst[j])
                   );
    }

    tmp = this->MscaleOfResiduals();
    return tmp * tmp + penalty;
}


/***************************************************************************************************
 *
 * Elastic Net - Peña Yohai Initial estimator with exact calculation
 * of the principal sensitivity components
 *
 **************************************************************************************************/
ENPY_Exact::ENPY_Exact(const Data& originalData, const Options& opts, const Options& enOpts) :
           InitialEstimator(originalData, opts, pscEn, MAX_NUM_PSCS(originalData.numObs())),

           /* Adjust lambda1 and lambda2 for LS objective */
           lambdaLS(this->ctrl.lambda * 0.5),

           /* Select correct ElasticNet implementation */
           en(*getElasticNetImpl(enOpts, true)),

           pscEn(en),
           dataToUse(residualFilteredData)
{
    this->en.setAlphaLambda(this->ctrl.alpha, this->lambdaLS);
}

ENPY_Exact::~ENPY_Exact()
{
    delete &this->en;
    this->pscFilteredData.free();
}

void ENPY_Exact::estimateCoefficients()
{
    if (this->dataToUse.numObs() == 0) {
        memset(this->coefEst, 0, this->dataToUse.numVar() * sizeof(double));
        memcpy(this->residuals, this->originalData.getYConst(),
               this->originalData.numObs() * sizeof(double));
        return;
    }

    this->en.setData(this->dataToUse);
    this->en.computeCoefs(this->coefEst, this->residuals);

    if (this->en.getStatus() > 0) {
        Rcpp::warning("EN did not converge. Either increase the number of "
                      "iterations or the penalty parameters.");
    }

    /* Now update the residuals */
    computeResiduals(this->originalData, this->coefEst, this->residuals);
}

void ENPY_Exact::resetData()
{
    this->residualFilteredNobs = this->originalData.numObs();

    this->residualFilteredData.free();
    this->residualFilteredData.setNumVar(this->originalData.numVar());


    this->residualFilteredData.free();
    this->residualFilteredData.copy(this->originalData);

    this->pscFilteredData.free();
    this->pscFilteredData.setNumObs(this->residualFilteredData.numObs());
    this->pscFilteredData.setNumVar(this->residualFilteredData.numVar());
    this->pscFilteredData.resize();
    this->dataToUse = this->residualFilteredData;
}


void ENPY_Exact::filterDataResiduals(double threshold)
{
    InitialEstimator::filterDataResiduals(threshold);
    this->dataToUse = this->residualFilteredData;
}

void ENPY_Exact::filterDataPSC(const double *RESTRICT values, const double threshold,
                         CompareFunction compare)
{
    doFiltering(this->residualFilteredData, this->pscFilteredData, values, threshold, compare);
    this->dataToUse = this->pscFilteredData;
}

double ENPY_Exact::evaluateEstimate() const
{
    int j;
    double penalty = 0;
    double tmp = 0.5 * (1 - this->ctrl.alpha);

    for (j = 1; j < this->originalData.numVar(); ++j) {
        penalty += this->ctrl.lambda * (
                        tmp * this->coefEst[j] * this->coefEst[j] +
                        this->ctrl.alpha * fabs(this->coefEst[j])
                   );
    }

    tmp = this->MscaleOfResiduals();
    return tmp * tmp + penalty;
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
