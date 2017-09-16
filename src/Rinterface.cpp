//
//  Rinterface.cpp
//  pense
//
//  Created by David Kepplinger on 2016-02-03.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#include <stdexcept>
#include <string>

#include "Rinterface.hpp"

#include "Options.hpp"
#include "Data.hpp"
#include "PSCxx.hpp"
#include "InitialEstimator.hpp"
#include "PENSEreg.hpp"
#include "MStep.hpp"
#include "olsreg.h"

using namespace Rcpp;
using namespace arma;

static inline Options listToOptions(SEXP list);
static inline void getMatDims(SEXP matrix, int* nrows, int* ncols);

RcppExport SEXP C_augtrans(SEXP RX)
{
    int nrow, ncol;
    getMatDims(RX, &nrow, &ncol);
    SEXP RXaug = PROTECT(Rf_allocMatrix(REALSXP, ncol + 1, nrow));
    const double *RESTRICT X = REAL(RX);
    const double *RESTRICT Xiter;
    double *RESTRICT Xaug = REAL(RXaug);
    int i, j;


    /* Augment X and transpose */
    for (i = 0; i < nrow; ++i) {
        (*Xaug) = 1;
        ++Xaug;

        Xiter = X + i;

        for (j = 0; j < ncol; ++j, ++Xaug, Xiter += nrow) {
            (*Xaug) = (*Xiter);
        }
    }

    UNPROTECT(1);

    return RXaug;
}

/***************************************************************************************************
 *
 * Elasitc Net with unweighted observations
 *
 **************************************************************************************************/
RcppExport SEXP C_elnet_sp(SEXP RXtr, SEXP Ry, SEXP Rcoefs, SEXP Ralpha,
                           SEXP Rlambda, SEXP Rintercept, SEXP Roptions,
                           SEXP RXtest)
{
    int nlambda = Rf_length(Rlambda);
    int nobs, nvar;
    getMatDims(RXtr, &nvar, &nobs);

    const Data data(REAL(RXtr), REAL(Ry), nobs, nvar);
    Options opts = listToOptions(Roptions);
    List retList;
    SEXP retResids = PROTECT(Rf_allocMatrix(REALSXP, data.numObs(), nlambda));
    SEXP retPreds = R_NilValue;
    sp_vec interceptSpVec(1);
    sp_vec currentBeta(nvar - 1);
    sp_mat coefEsts(nvar, nlambda);
    double adjFactor = 1;
    double intercept = 0;
    double* currentResidualsPtr = REAL(retResids);
    double* currentPredsPtr = NULL;
    const double *currentLambda = REAL(Rlambda);
    const double alpha = *REAL(Ralpha);
    const bool generatePredictions = Rf_isReal(RXtest);
    const bool estimate_intercept = (bool) *INTEGER(Rintercept);

    BEGIN_RCPP
    const mat Xtest = (generatePredictions ? as<mat>(RXtest) : mat());
    /* this does not need the column of 1's in the beginning */
    const mat XtrTrain = as<mat>(RXtr).tail_rows(nvar - 1);
    const vec yTrain = as< colvec >(Ry);


    if (generatePredictions) {
        retPreds = PROTECT(Rf_allocMatrix(REALSXP, Xtest.n_rows, nlambda));
        currentPredsPtr = REAL(retPreds);
    }

    if (opts.get("warmStart", true)) {
        sp_vec givenCoefs = as<sp_mat>(Rcoefs).col(0);
        interceptSpVec[0] = intercept = givenCoefs[0];
        currentBeta = givenCoefs.tail_rows(nvar - 1);
    }

    /*
     * We can always use a warm start, since the currentCoefs would be 0 anyways if it wasn't
     * specified otherwise.
     */
    opts.set("warmStart", true);

    ElasticNet *en = getElasticNetImpl(opts, estimate_intercept);
    sp_vec::const_iterator ccIt;
    en->setData(data);

    for (int i = 0; i < nlambda; ++i, currentResidualsPtr += nobs, ++currentLambda) {
        vec residuals(currentResidualsPtr, nobs, false, true);
        en->setAlphaLambda(alpha, *currentLambda);

        en->computeCoefs(intercept, currentBeta, residuals);
        interceptSpVec[0] = intercept;

        if (en->getStatus() != 0) {
            Rcpp::warning("EN algorithm had non-zero exit status for lambda=%g: %s",
                          *currentLambda, en->getStatusMessage());
        }

        if (!opts.get("naive", false) && alpha < 1) {
            adjFactor = sqrt(1.0 + (1 - alpha) * (*currentLambda));
            coefEsts.col(i) = join_cols(
                interceptSpVec,
                currentBeta * adjFactor
            );
            if (estimate_intercept) {
                coefEsts(0, i) = mean(yTrain - XtrTrain.t() * currentBeta * adjFactor);
            }
        } else {
            coefEsts.col(i) = join_cols(interceptSpVec, currentBeta);
        }

        if (generatePredictions) {
            vec predAlias(currentPredsPtr, Xtest.n_rows, false, true);
            predAlias = Xtest * currentBeta * adjFactor + coefEsts(0, i);
            currentPredsPtr += Xtest.n_rows;
        }
    }

    retList = List::create(
        Named("status") = en->getStatus(),
        Named("message") = en->getStatusMessage(),
        Named("coefficients") = coefEsts,
        Named("residuals") = retResids,
        Named("predictions") = retPreds
    );

    delete en;

    VOID_END_RCPP

    if (generatePredictions) {
        UNPROTECT(2);
    } else {
        UNPROTECT(1);
    }
    return wrap(retList);
}


/***************************************************************************************************
 *
 * Elasitc Net with weighted observations
 *
 **************************************************************************************************/
RcppExport SEXP C_elnet_weighted_sp(SEXP RXtr, SEXP Ry, SEXP Rweights, SEXP Rcoefs,
                                    SEXP Ralpha, SEXP Rlambda, SEXP Rintercept, SEXP Roptions,
                                    SEXP RXtest)
{
    int nlambda = Rf_length(Rlambda);
    int nobs, nvar;
    getMatDims(RXtr, &nvar, &nobs);

    const Data data(REAL(RXtr), REAL(Ry), nobs, nvar);
    const vec weights(REAL(Rweights), nobs, false, true);
    Options opts = listToOptions(Roptions);
    List retList;
    SEXP retResids = PROTECT(Rf_allocMatrix(REALSXP, data.numObs(), nlambda));
    SEXP retPreds = R_NilValue;
    sp_vec interceptSpVec(1);
    sp_vec currentBeta(nvar - 1);
    sp_mat coefEsts(nvar, nlambda);
    double intercept;
    double adjFactor = 1;
    double* currentResidualsPtr = REAL(retResids);
    double* currentPredsPtr = NULL;
    const double *currentLambda = REAL(Rlambda);
    const double alpha = *REAL(Ralpha);
    const bool generatePredictions = Rf_isReal(RXtest);
    const bool estimate_intercept = (bool) *INTEGER(Rintercept);

    BEGIN_RCPP
    const mat Xtest = (generatePredictions ? as<mat>(RXtest) : mat());
    /* this does not need the column of 1's in the beginning */
    const mat XtrTrain = as<mat>(RXtr).tail_rows(nvar - 1);
    const vec yTrain = as< colvec >(Ry);

    if (generatePredictions) {
        retPreds = PROTECT(Rf_allocMatrix(REALSXP, Xtest.n_rows, nlambda));
        currentPredsPtr = REAL(retPreds);
    }

    if (opts.get("warmStart", true)) {
        sp_vec givenCoefs = as<sp_mat>(Rcoefs).col(0);
        interceptSpVec[0] = intercept = givenCoefs[0];
        currentBeta = givenCoefs.tail_rows(nvar - 1);
    }

    /*
     * We can always use a warm start, since the currentCoefs would be 0 anyways if it wasn't
     * specified otherwise.
     */
    opts.set("warmStart", true);

    ElasticNet *en = getElasticNetImpl(opts, estimate_intercept);
    sp_vec::const_iterator ccIt;
    en->setData(data);

    for (int i = 0; i < nlambda; ++i, currentResidualsPtr += nobs, ++currentLambda) {
        vec residuals(currentResidualsPtr, nobs, false, true);
        en->setAlphaLambda(alpha, *currentLambda);

        en->computeCoefsWeighted(intercept, currentBeta, residuals, weights);

        interceptSpVec[0] = intercept;

        if (en->getStatus() != 0) {
            Rcpp::warning("EN algorithm had non-zero exit status for lambda=%g: %s",
                          *currentLambda, en->getStatusMessage());
        }

        if (!opts.get("naive", false) && alpha < 1) {
            adjFactor = sqrt(1.0 + (1 - alpha) * (*currentLambda));
            coefEsts.col(i) = join_cols(
                interceptSpVec,
                currentBeta * adjFactor
            );
            if (estimate_intercept) {
                coefEsts(0, i) = accu(weights % (yTrain - XtrTrain.t() * currentBeta * adjFactor)) /
                                    accu(weights);
            }
        } else {
            coefEsts.col(i) = join_cols(interceptSpVec, currentBeta);
        }

        if (generatePredictions) {
            vec predAlias(currentPredsPtr, Xtest.n_rows, false, true);
            predAlias = Xtest * currentBeta * adjFactor + coefEsts(0, i);
            currentPredsPtr += Xtest.n_rows;
        }
    }

    retList = List::create(
        Named("status") = en->getStatus(),
        Named("message") = en->getStatusMessage(),
        Named("coefficients") = coefEsts,
        Named("residuals") = retResids,
        Named("predictions") = retPreds
    );

    delete en;

    VOID_END_RCPP

    if (generatePredictions) {
        UNPROTECT(2);
    } else {
        UNPROTECT(1);
    }
    return wrap(retList);
}

/***************************************************************************************************
 *
 * Compute Principal Sensitivity Components for non-regularized regression
 *
 **************************************************************************************************/
RcppExport SEXP C_pscs_ols(SEXP RXtr, SEXP Ry)
{
    int nobs, nvar;
    getMatDims(RXtr, &nvar, &nobs);

    const Data data(REAL(RXtr), REAL(Ry), nobs, nvar);
    SEXP ret = R_NilValue;
    PSC_OLS psc;
    double *RESTRICT coefs = new double[data.numVar()];
    double *RESTRICT residuals = new double[data.numObs()];
    double *RESTRICT Xsqrt = new double[data.numVar() * data.numVar()];
    int npscs;

    BEGIN_RCPP

    npscs = computeOLSCoefs(data.getXtrConst(), data.getYConst(), data.numObs(), data.numVar(),
                            coefs, Xsqrt);

    if (npscs == 0) {
        computeResiduals(data.getXtrConst(), data.getYConst(), data.numObs(), data.numVar(), coefs,
                         residuals);

        psc.setData(data);
        psc.setXsqrtMemory(Xsqrt);
        psc.setResiduals(residuals);
        npscs = psc.computePSC();

        ret = PROTECT(Rf_allocVector(REALSXP, data.numObs() * npscs));
        memcpy(REAL(ret), psc.getPSC(), data.numObs() * npscs * sizeof(double));
        UNPROTECT(1);
    }

    VOID_END_RCPP

    delete[] coefs;
    delete[] residuals;
    delete[] Xsqrt;
    return ret;
}

/***************************************************************************************************
 *
 * Compute Principal Sensitivity Components for EN penalized regression
 *
 **************************************************************************************************/
RcppExport SEXP C_pscs_en(SEXP RXtr, SEXP Ry, SEXP Ralpha, SEXP Rlambda,
                          SEXP Rintercept, SEXP Roptions)
{
    int nobs, nvar;
    getMatDims(RXtr, &nvar, &nobs);

    const Data data(REAL(RXtr), REAL(Ry), nobs, nvar);
    const Options opts = listToOptions(Roptions);
    ElasticNet *en = getElasticNetImpl(opts, (bool) *INTEGER(Rintercept));

    SEXP ret = R_NilValue;

    double *RESTRICT coefs = new double[data.numVar()];
    double *RESTRICT residuals = new double[data.numObs()];
    int npscs;

    BEGIN_RCPP

    en->setAlphaLambda(*REAL(Ralpha), *REAL(Rlambda));
    en->setData(data);
    en->computeCoefs(coefs, residuals);

    if (en->getStatus() > 0) {
        throw std::runtime_error(en->getStatusMessage());
    }

    PSC_EN psc(*en);

    psc.setData(data);
    psc.setResiduals(residuals);
    npscs = psc.computePSC();

    ret = PROTECT(Rf_allocVector(REALSXP, data.numObs() * npscs));
    memcpy(REAL(ret), psc.getPSC(), data.numObs() * npscs * sizeof(double));

    UNPROTECT(1);
    delete en;

    VOID_END_RCPP

    delete[] coefs;
    delete[] residuals;
    return ret;
}


/***************************************************************************************************
 *
 * PY Initial Estimator for OLS (i.e., non-regularized) problems
 *
 **************************************************************************************************/
RcppExport SEXP C_py_ols(SEXP RXtr, SEXP Ry, SEXP RpyOptions)
{
    int nobs, nvar;
    getMatDims(RXtr, &nvar, &nobs);

    const Options opts = listToOptions(RpyOptions);
    const Data data(REAL(RXtr), REAL(Ry), nobs, nvar);

    IEOls ols(data, opts);
    SEXP coefs;
    SEXP objF;
    int niest;
    SEXP result = R_NilValue;

    BEGIN_RCPP

    niest = ols.compute();

    result = PROTECT(Rf_allocVector(VECSXP, 2));
    coefs = PROTECT(Rf_allocVector(REALSXP, niest * data.numVar()));
    objF = PROTECT(Rf_allocVector(REALSXP, niest));

    memcpy(REAL(coefs), ols.getInitialEstimators(), data.numVar() * niest * sizeof(double));
    memcpy(REAL(objF), ols.getObjectiveFunctionScores(), niest * sizeof(double));

    SET_VECTOR_ELT(result, 0, coefs);
    SET_VECTOR_ELT(result, 1, objF);

    UNPROTECT(3);

    VOID_END_RCPP

    return result;
}

/***************************************************************************************************
 *
 * PY Initial Estimator for EN regularized problems, using an Ridge approximation
 *
 **************************************************************************************************/
RcppExport SEXP C_enpy_rr(SEXP RXtr, SEXP Ry, SEXP Ralpha, SEXP Rlambda, SEXP RpyOptions,
                          SEXP RenOptions)
{
    int nobs, nvar;
    getMatDims(RXtr, &nvar, &nobs);

    const Options opts = listToOptions(RpyOptions);
    const Options enOpts = listToOptions(RenOptions);
    const Data data(REAL(RXtr), REAL(Ry), nobs, nvar);

    ENPY enpy(data, *REAL(Ralpha), *REAL(Rlambda), opts, enOpts);
    SEXP coefs;
    SEXP objF;
    int niest;
    SEXP result = R_NilValue;

    BEGIN_RCPP

    niest = enpy.compute();

    result = PROTECT(Rf_allocVector(VECSXP, 2));
    coefs = PROTECT(Rf_allocVector(REALSXP, niest * data.numVar()));
    objF = PROTECT(Rf_allocVector(REALSXP, niest));

    memcpy(REAL(coefs), enpy.getInitialEstimators(), data.numVar() * niest * sizeof(double));
    memcpy(REAL(objF), enpy.getObjectiveFunctionScores(), niest * sizeof(double));

    SET_VECTOR_ELT(result, 0, coefs);
    SET_VECTOR_ELT(result, 1, objF);

    UNPROTECT(3);

    VOID_END_RCPP

    return result;
}

/***************************************************************************************************
 *
 * PY Initial Estimator for EN regularized problems
 *
 **************************************************************************************************/
RcppExport SEXP C_enpy_exact(SEXP RXtr, SEXP Ry, SEXP Ralpha, SEXP Rlambda, SEXP RpyOptions,
                             SEXP RenOptions)
{
    int nobs, nvar;
    getMatDims(RXtr, &nvar, &nobs);

    const Options opts = listToOptions(RpyOptions);
    const Options enOpts = listToOptions(RenOptions);
    const Data data(REAL(RXtr), REAL(Ry), nobs, nvar);
    ENPY_Exact enpy(data, *REAL(Ralpha), *REAL(Rlambda), opts, enOpts);
    SEXP coefs, objF;
    int niest;
    SEXP result = R_NilValue;

    BEGIN_RCPP

    niest = enpy.compute();

    result = PROTECT(Rf_allocVector(VECSXP, 2));
    coefs = PROTECT(Rf_allocVector(REALSXP, niest * data.numVar()));
    objF = PROTECT(Rf_allocVector(REALSXP, niest));

    memcpy(REAL(coefs), enpy.getInitialEstimators(), data.numVar() * niest * sizeof(double));
    memcpy(REAL(objF), enpy.getObjectiveFunctionScores(), niest * sizeof(double));

    SET_VECTOR_ELT(result, 0, coefs);
    SET_VECTOR_ELT(result, 1, objF);

    UNPROTECT(3);

    VOID_END_RCPP

    return result;
}

/***************************************************************************************************
 *
 * Penalized Elastic Net S estimator for regression (PENSE)
 *
 **************************************************************************************************/
RcppExport SEXP C_pen_s_reg_sp(SEXP RXtr, SEXP Ry, SEXP Rintercept, SEXP Rcoefs,
                               SEXP Ralpha, SEXP Rlambda, SEXP RpenseOptions, SEXP RenOptions)
{
    int nobs, nvar;
    getMatDims(RXtr, &nvar, &nobs);

    const Options penseOpts = listToOptions(RpenseOptions);
    const Options enOpts = listToOptions(RenOptions);
    const Data data(REAL(RXtr), REAL(Ry), nobs, nvar);

    PENSEReg pr(data, *REAL(Ralpha), *REAL(Rlambda), penseOpts, enOpts);
    List retList;
    double intercept = *REAL(Rintercept);
    SEXP residuals = PROTECT(Rf_allocVector(REALSXP, data.numObs()));

    BEGIN_RCPP
    vec residVec(REAL(residuals), nobs, false, true);
    sp_vec beta = as<sp_mat>(Rcoefs).col(0);

    pr.compute(intercept, beta, residVec);

    retList = List::create(
        Named("intercept") = intercept,
        Named("beta") = sp_mat(beta),
        Named("residuals") = residuals,
        Named("scale") = pr.getScale(),
        Named("weights") = pr.getWeights(),
        Named("rel_change") = pr.relChange(),
        Named("iterations") = pr.iterations()
    );

    VOID_END_RCPP

    UNPROTECT(1);

    return wrap(retList);
}

/***************************************************************************************************
 *
 * Penalized Elastic Net M estimator for regression with initial scale (M-Step)
 *
 **************************************************************************************************/
RcppExport SEXP C_pen_mstep_sp(SEXP RXtr, SEXP Ry, SEXP Rintercept, SEXP Rcoefs, SEXP scale,
                               SEXP Ralpha, SEXP Rlambda, SEXP RmsOptions, SEXP RenOptions)
{
    int nobs, nvar;
    getMatDims(RXtr, &nvar, &nobs);

    const Options msOpts = listToOptions(RmsOptions);
    const Options enOpts = listToOptions(RenOptions);
    const Data data(REAL(RXtr), REAL(Ry), nobs, nvar);

    MStep ms(data, *REAL(Ralpha), *REAL(Rlambda), *REAL(scale), msOpts, enOpts);

    List retList;
    SEXP residuals = PROTECT(Rf_allocVector(REALSXP, data.numObs()));
    double intercept = *REAL(Rintercept);

    BEGIN_RCPP
    vec residVec(REAL(residuals), nobs, false, true);
    sp_vec beta = as<sp_mat>(Rcoefs).col(0);

    ms.compute(intercept, beta, residVec);

    retList = List::create(
        Named("intercept") = intercept,
        Named("beta") = sp_mat(beta),
        Named("residuals") = residuals,
        Named("weights") = ms.getWeights(),
        Named("rel_change") = ms.relChange(),
        Named("iterations") = ms.iterations()
    );

    VOID_END_RCPP

    UNPROTECT(1);

    return wrap(retList);
}


static inline Options listToOptions(SEXP Rlist)
{
    const List values = List(Rlist);
    const List names = values.attr("names");

    Options opts;
    List::const_iterator valuesIt = values.begin();
    List::const_iterator namesIt = names.begin();
    for (; valuesIt != values.end(); ++valuesIt, ++namesIt) {
        switch (TYPEOF(*valuesIt)) {
        case INTSXP:
            opts.set(as<std::string>(*namesIt), as<int>(*valuesIt));
            break;
        case REALSXP:
            opts.set(as<std::string>(*namesIt), as<double>(*valuesIt));
            break;
        case LGLSXP:
            opts.set(as<std::string>(*namesIt), as<bool>(*valuesIt));
            break;
        default:
            break;
        }
    }

    return opts;
}

static inline void getMatDims(SEXP matrix, int* nrows, int* ncols)
{
    SEXP Rdims;
    int* dims;
    PROTECT(Rdims = Rf_getAttrib(matrix, R_DimSymbol));
    dims = INTEGER(Rdims);
    *nrows = dims[0];
    *ncols = dims[1];
    UNPROTECT(1);
}
