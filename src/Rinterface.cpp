//
//  Rinterface.cpp
//  penseinit
//
//  Created by David Kepplinger on 2016-02-03.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#include <stdexcept>

#include "Rinterface.hpp"

#include "Control.h"
#include "Data.hpp"
#include "PSCxx.hpp"
#include "InitialEstimator.hpp"
#include "PENSEreg.hpp"
#include "MStep.hpp"
#include "olsreg.h"

using namespace Rcpp;

static inline Control parseControlList(SEXP control);

RcppExport SEXP C_augtrans(SEXP RX, SEXP Rnrow, SEXP Rncol)
{
    int nrow = *INTEGER(Rnrow);
    int ncol = *INTEGER(Rncol);
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

RcppExport SEXP C_pen_s_reg(SEXP RXtr, SEXP Ry, SEXP Rnobs, SEXP Rnvar, SEXP coefs,
                            SEXP Ralpha, SEXP Rlambda, SEXP Rcontrol)
{
    const Control ctrl = parseControlList(Rcontrol);
    const Data data(REAL(RXtr), REAL(Ry), *INTEGER(Rnobs), *INTEGER(Rnvar));

    PENSEReg pr(data, *REAL(Ralpha), *REAL(Rlambda), ctrl);
    SEXP newCoefs = PROTECT(Rf_allocVector(REALSXP, data.numVar()));
    SEXP residuals = PROTECT(Rf_allocVector(REALSXP, data.numObs()));
    SEXP relChange;
    SEXP scale;
    SEXP iterations;
    SEXP result = R_NilValue;
    double *RESTRICT newCoefsPtr = REAL(newCoefs);

    BEGIN_RCPP

    memcpy(newCoefsPtr, REAL(coefs), data.numVar() * sizeof(double));
    pr.compute(newCoefsPtr, REAL(residuals));

    result = PROTECT(Rf_allocVector(VECSXP, 5));
    scale = PROTECT(Rf_ScalarReal(pr.getScale()));
    relChange = PROTECT(Rf_ScalarReal(pr.getRelChange()));
    iterations = PROTECT(Rf_ScalarInteger(pr.getIterations()));

    SET_VECTOR_ELT(result, 0, newCoefs);
    SET_VECTOR_ELT(result, 1, residuals);
    SET_VECTOR_ELT(result, 2, scale);
    SET_VECTOR_ELT(result, 3, relChange);
    SET_VECTOR_ELT(result, 4, iterations);

    UNPROTECT(4);

    VOID_END_RCPP

    UNPROTECT(2);

    return result;
}

RcppExport SEXP C_pen_mstep(SEXP RXtr, SEXP Ry, SEXP Rnobs, SEXP Rnvar, SEXP coefs, SEXP scale,
                            SEXP Ralpha, SEXP Rlambda, SEXP Rcontrol)
{
    const Control ctrl = parseControlList(Rcontrol);
    const Data data(REAL(RXtr), REAL(Ry), *INTEGER(Rnobs), *INTEGER(Rnvar));

    MStep ms(data, *REAL(Ralpha), *REAL(Rlambda), ctrl);
    SEXP newCoefs = PROTECT(Rf_allocVector(REALSXP, data.numVar()));
    SEXP residuals = PROTECT(Rf_allocVector(REALSXP, data.numObs()));
    SEXP relChange;
    SEXP iterations;
    SEXP result = R_NilValue;
    double *RESTRICT newCoefsPtr = REAL(newCoefs);

    BEGIN_RCPP

    memcpy(newCoefsPtr, REAL(coefs), data.numVar() * sizeof(double));
    ms.compute(newCoefsPtr, *REAL(scale), REAL(residuals));

    result = PROTECT(Rf_allocVector(VECSXP, 4));
    relChange = PROTECT(Rf_ScalarReal(ms.getRelChange()));
    iterations = PROTECT(Rf_ScalarInteger(ms.getIterations()));

    SET_VECTOR_ELT(result, 0, newCoefs);
    SET_VECTOR_ELT(result, 1, residuals);
    SET_VECTOR_ELT(result, 2, relChange);
    SET_VECTOR_ELT(result, 3, iterations);

    UNPROTECT(3);

    VOID_END_RCPP

    UNPROTECT(2);

    return result;
}


RcppExport SEXP C_py_ols(SEXP RXtr, SEXP Ry, SEXP Rnobs, SEXP Rnvar, SEXP Rcontrol)
{
    const Control ctrl = parseControlList(Rcontrol);
    const Data data(REAL(RXtr), REAL(Ry), *INTEGER(Rnobs), *INTEGER(Rnvar));

    OLS ols(data, ctrl);
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


RcppExport SEXP C_enpy_rr(SEXP RXtr, SEXP Ry, SEXP Rnobs, SEXP Rnvar, SEXP Rcontrol)
{
    const Control ctrl = parseControlList(Rcontrol);
    const Data data(REAL(RXtr), REAL(Ry), *INTEGER(Rnobs), *INTEGER(Rnvar));

    ENPY enpy(data, ctrl);
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


RcppExport SEXP C_elnet(SEXP RXtr, SEXP Ry, SEXP Rcoefs, SEXP Rnobs, SEXP Rnvar, SEXP Ralpha,
                        SEXP Rlambda, SEXP RmaxIt, SEXP Reps, SEXP Rcentering, SEXP Rwarm,
                        SEXP RenAlgorithm)
{
    const Data data(REAL(RXtr), REAL(Ry), *INTEGER(Rnobs), *INTEGER(Rnvar));
    ElasticNet *en = getElasticNetImpl((ENAlgorithm) *INTEGER(RenAlgorithm),
                                       *REAL(Reps), (bool) *INTEGER(Rcentering),
                                       *INTEGER(RmaxIt));
    bool warm = (*INTEGER(Rwarm) == 1);
    SEXP result = R_NilValue;
    SEXP retCoef = PROTECT(Rf_allocVector(REALSXP, data.numVar()));
    SEXP retResid = PROTECT(Rf_allocVector(REALSXP, data.numObs()));
    SEXP converged = PROTECT(Rf_allocVector(LGLSXP, 1));
    double *retCoefPtr = REAL(retCoef);

    BEGIN_RCPP

    if (warm) {
        memcpy(retCoefPtr, REAL(Rcoefs), data.numVar() * sizeof(double));
    }

    en->setAlphaLambda(*REAL(Ralpha), *REAL(Rlambda));
    *LOGICAL(converged) = en->computeCoefs(data, retCoefPtr, REAL(retResid), warm);

    result = PROTECT(Rf_allocVector(VECSXP, 3));

    SET_VECTOR_ELT(result, 0, converged);
    SET_VECTOR_ELT(result, 1, retCoef);
    SET_VECTOR_ELT(result, 2, retResid);

    delete en;
    UNPROTECT(1);

    VOID_END_RCPP

    UNPROTECT(3);
    return result;
}


RcppExport SEXP C_elnet_ll(SEXP RXtr, SEXP Ry, SEXP Rcoefs, SEXP Rnobs, SEXP Rnvar, SEXP Rlambda1,
                           SEXP Rlambda2, SEXP RmaxIt, SEXP Reps, SEXP Rcentering, SEXP Rwarm,
                           SEXP RenAlgorithm)
{
    const Data data(REAL(RXtr), REAL(Ry), *INTEGER(Rnobs), *INTEGER(Rnvar));
    ElasticNet *en = getElasticNetImpl((ENAlgorithm) *INTEGER(RenAlgorithm),
                                       *REAL(Reps), (bool) *INTEGER(Rcentering),
                                       *INTEGER(RmaxIt));
    bool warm = (*INTEGER(Rwarm) == 1);
    SEXP result = R_NilValue;
    SEXP retCoef = PROTECT(Rf_allocVector(REALSXP, data.numVar()));
    SEXP retResid = PROTECT(Rf_allocVector(REALSXP, data.numObs()));
    SEXP converged = PROTECT(Rf_allocVector(LGLSXP, 1));
    double *retCoefPtr = REAL(retCoef);

    BEGIN_RCPP

    if (warm) {
        memcpy(retCoefPtr, REAL(Rcoefs), data.numVar() * sizeof(double));
    }

    en->setLambdas(*REAL(Rlambda1), *REAL(Rlambda2));
    *LOGICAL(converged) = en->computeCoefs(data, retCoefPtr, REAL(retResid), warm);

    result = PROTECT(Rf_allocVector(VECSXP, 3));

    SET_VECTOR_ELT(result, 0, converged);
    SET_VECTOR_ELT(result, 1, retCoef);
    SET_VECTOR_ELT(result, 2, retResid);

    delete en;
    UNPROTECT(1);

    VOID_END_RCPP

    UNPROTECT(3);
    return result;
}


RcppExport SEXP C_enpy_exact(SEXP RXtr, SEXP Ry, SEXP Rnobs, SEXP Rnvar, SEXP Rcontrol)
{
    const Control ctrl = parseControlList(Rcontrol);
    const Data data(REAL(RXtr), REAL(Ry), *INTEGER(Rnobs), *INTEGER(Rnvar));
    ENPY_Exact enpy(data, ctrl);
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

RcppExport SEXP C_pscs_ols(SEXP RXtr, SEXP Ry, SEXP Rnobs, SEXP Rnvar)
{
    const Data data(REAL(RXtr), REAL(Ry), *INTEGER(Rnobs), *INTEGER(Rnvar));
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

RcppExport SEXP C_pscs_en(SEXP RXtr, SEXP Ry, SEXP Rnobs, SEXP Rnvar, SEXP Ralpha, SEXP Rlambda,
                          SEXP RmaxIt, SEXP Reps, SEXP Rcentering, SEXP RenAlgorithm)
{
    const Data data(REAL(RXtr), REAL(Ry), *INTEGER(Rnobs), *INTEGER(Rnvar));
    ElasticNet *en = getElasticNetImpl((ENAlgorithm) *INTEGER(RenAlgorithm),
                                       *REAL(Reps), (bool) *INTEGER(Rcentering),
                                       *INTEGER(RmaxIt));

    SEXP ret = R_NilValue;

    double *RESTRICT coefs = new double[data.numVar()];
    double *RESTRICT residuals = new double[data.numObs()];
    int npscs;
    bool converged;

    BEGIN_RCPP

    en->setAlphaLambda(*REAL(Ralpha), *REAL(Rlambda));
    converged = en->computeCoefs(data, coefs, residuals);

    if (!converged) {
        throw std::runtime_error("Elastic Net did not converge");
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


static inline Control parseControlList(SEXP Rcontrol)
{
    List control = List(Rcontrol);
    Control tmp = {
        as<double>(control["lambda1"]),
        as<double>(control["lambda2"]),
        as<int>(control["numIt"]),
        as<double>(control["eps"]),
        as<double>(control["resid.threshold"]),
        as<double>(control["resid.proportion"]),
        as<double>(control["psc.proportion"]),

        as<int>(control["en.maxit"]),
        as<double>(control["en.tol"]),
        as<int>(control["en.centering"]),
        (ENAlgorithm) as<int>(control["en.algorithm"]),

        as<double>(control["mscale.delta"]),
        as<double>(control["mscale.cc"]),
        as<int>(control["mscale.maxit"]),
        as<double>(control["mscale.tol"]),
        (RhoFunctionName) as<int>(control["mscale.rho.fun"])
    };

    return tmp;
}

