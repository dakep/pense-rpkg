//
//  r_interface.cc
//  pense
//
//  Created by David Kepplinger on 2019-04-03.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifdef HAVE_RCPP
#include <R_ext/Rdynload.h>

#include "rcpp_integration.hpp"
#include "r_en_regression.hpp"
#include "r_pense_regression.hpp"
#include "r_mesten_regression.hpp"
#include "r_robust_utils.hpp"
#include "r_enpy.hpp"
#include "r_utilities.hpp"

extern "C" SEXP run_testthat_tests() noexcept;

//! R initialzing function (must be in the global namespace).
extern "C" void R_init_pense(DllInfo *dll) noexcept;

using namespace pense::r_interface;

namespace {
//! Exported methods
const R_CallMethodDef kExportedCallMethods[] = {
  // {"C_run_testthat_tests", (DL_FUNC) &run_testthat_tests, 0},
  {"C_tau_size", (DL_FUNC) &TauSize, 1},
  {"C_approx_match", (DL_FUNC) &ApproximateMatch, 3},
  {"C_mscale", (DL_FUNC) &MScale, 2},
  {"C_mloc", (DL_FUNC) &MLocation, 3},
  {"C_mlocscale", (DL_FUNC) &MLocationScale, 3},
  {"C_lsen_regression", (DL_FUNC) &LsEnRegression, 5},
  {"C_pense_regression", (DL_FUNC) &PenseEnRegression, 7},
  {"C_pense_max_lambda", (DL_FUNC) &PenseMaxLambda, 4},
  {"C_mesten_regression", (DL_FUNC) &MestEnRegression, 6},
  {"C_mesten_max_lambda", (DL_FUNC) &MestEnMaxLambda, 5},
  {"C_penpy", (DL_FUNC) &PenPyInitialEstimator, 6},
  {"C_pscs", (DL_FUNC) &PrincipalSensitivityComponents, 5},
  {NULL, NULL, 0}
};
}  // namespace

extern "C" void R_init_pense(DllInfo *dll) noexcept {
    R_registerRoutines(dll, NULL, kExportedCallMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

#endif  // HAVE_RCPP
