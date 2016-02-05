//
//  Rinterface.cpp
//  penseinit
//
//  Created by David Kepplinger on 2016-02-03.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#include "Rinterface.hpp"

#include "Control.h"
#include "Data.hpp"
#include "InitialEstimator.hpp"

using namespace Rcpp;

static inline Control parseControlList(SEXP control);

RcppExport SEXP C_enpy_rr(SEXP RXtr, SEXP Ry, SEXP Rnobs, SEXP Rnvar, SEXP Rcontrol)
{
    Control ctrl = parseControlList(Rcontrol);
    Data data(REAL(RXtr), REAL(Ry), *INTEGER(Rnobs), *INTEGER(Rnvar));

    ENPY enpy(data, ctrl);
    SEXP coefs, objF;
    int niest;
    SEXP result;

    BEGIN_RCPP

    try {
        niest = enpy.compute();
    } catch (std::exception &e) {
        throw;
    }

    result = PROTECT(Rf_allocVector(VECSXP, 2));
    coefs = PROTECT(Rf_allocVector(REALSXP, niest * data.numVar()));
    objF = PROTECT(Rf_allocVector(REALSXP, niest));

    memcpy(REAL(coefs), enpy.getInitialEstimators(), data.numVar() * niest * sizeof(double));
    memcpy(REAL(objF), enpy.getObjectiveFunctionScores(), niest * sizeof(double));

    SET_VECTOR_ELT(result, 0, coefs);
    SET_VECTOR_ELT(result, 1, objF);

    UNPROTECT(3);
    return result;

    END_RCPP
}

RcppExport SEXP C_enpy_Mn(SEXP RXtr, SEXP Ry, SEXP Rnobs, SEXP Rnvar, SEXP Rcontrol)
{
    Control ctrl = parseControlList(Rcontrol);
    Data data(REAL(RXtr), REAL(Ry), *INTEGER(Rnobs), *INTEGER(Rnvar));

    ENPY enpy(data, ctrl);
    SEXP coefs, objF;
    int niest;
    SEXP result;

    BEGIN_RCPP

    try {
        niest = enpy.compute();
    } catch (std::exception &e) {
        throw;
    }

    result = PROTECT(Rf_allocVector(VECSXP, 2));
    coefs = PROTECT(Rf_allocVector(REALSXP, niest * data.numVar()));
    objF = PROTECT(Rf_allocVector(REALSXP, niest));

    memcpy(REAL(coefs), enpy.getInitialEstimators(), data.numVar() * niest * sizeof(double));
    memcpy(REAL(objF), enpy.getObjectiveFunctionScores(), niest * sizeof(double));

    SET_VECTOR_ELT(result, 0, coefs);
    SET_VECTOR_ELT(result, 1, objF);

    UNPROTECT(3);
    return result;

    END_RCPP
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

