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
#include "olsreg.h"

using namespace Rcpp;

static inline Control parseControlList(SEXP control);


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

RcppExport SEXP C_elnet(SEXP RXtr, SEXP Ry, SEXP Rnobs, SEXP Rnvar, SEXP Ralpha, SEXP Rlambda,
                        SEXP RmaxIt, SEXP Reps, SEXP Rcentering)
{
    const Data data(REAL(RXtr), REAL(Ry), *INTEGER(Rnobs), *INTEGER(Rnvar));
    ElasticNet en(*INTEGER(RmaxIt), *REAL(Reps), (bool) *INTEGER(Rcentering));
    SEXP result = R_NilValue;
    SEXP retCoef = PROTECT(Rf_allocVector(REALSXP, data.numVar()));
    SEXP retResid = PROTECT(Rf_allocVector(REALSXP, data.numObs()));
    SEXP converged = PROTECT(Rf_allocVector(LGLSXP, 1));

    BEGIN_RCPP

    en.setAlphaLambda(*REAL(Ralpha), *REAL(Rlambda));
    *LOGICAL(converged) = en.computeCoefs(data, REAL(retCoef), REAL(retResid));

    result = PROTECT(Rf_allocVector(VECSXP, 3));

    SET_VECTOR_ELT(result, 0, converged);
    SET_VECTOR_ELT(result, 1, retCoef);
    SET_VECTOR_ELT(result, 2, retResid);

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
                          SEXP RmaxIt, SEXP Reps, SEXP Rcentering)
{
    const Data data(REAL(RXtr), REAL(Ry), *INTEGER(Rnobs), *INTEGER(Rnvar));
    ElasticNet en(*INTEGER(RmaxIt), *REAL(Reps), (bool) *INTEGER(Rcentering));
    PSC_EN psc(en);

    SEXP ret = R_NilValue;

    double *RESTRICT coefs = new double[data.numVar()];
    double *RESTRICT residuals = new double[data.numObs()];
    int npscs;
    bool converged;

    BEGIN_RCPP

    en.setAlphaLambda(*REAL(Ralpha), *REAL(Rlambda));
    converged = en.computeCoefs(data, coefs, residuals);

    if (!converged) {
        throw std::runtime_error("Elastic Net did not converge");
    }

    psc.setData(data);
    psc.setResiduals(residuals);
    npscs = psc.computePSC();

    ret = PROTECT(Rf_allocVector(REALSXP, data.numObs() * npscs));
    memcpy(REAL(ret), psc.getPSC(), data.numObs() * npscs * sizeof(double));
    UNPROTECT(1);

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
        as<double>(control["residThreshold"]),
        as<double>(control["residProportion"]),
        as<double>(control["pscProportion"]),

        as<int>(control["enMaxIt"]),
        as<double>(control["enEPS"]),
        as<int>(control["enCentering"]),

        as<double>(control["mscaleB"]),
        as<double>(control["mscaleCC"]),
        as<int>(control["mscaleMaxIt"]),
        as<double>(control["mscaleEPS"]),
        (RhoFunctionName) as<int>(control["mscaleRhoFun"])
    };

    return tmp;
}

