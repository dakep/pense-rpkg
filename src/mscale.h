//
//  mscale.h
//  pense
//
//  Created by David Kepplinger on 2016-01-30.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#ifndef mscale_h
#define mscale_h

#include "config.h"
#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef double (*RhoFunction)(double, const double);

/**
 * Get the rho function by name
 */
RhoFunction getRhoFunctionByName(RhoFunctionName name);

/**
 * Calculate the M-Scale of a vector of numbers
 *
 * @param Rvalues numeric The vector of REAL's
 * @param Rlength integer The length of Rvalues
 * @param Rb      numeric The target average (right side) in the M-scale equation
 * @param Rcc     numeric The constant for the rho function
 * @param RmaxIt  integer The maximum number of iterations
 * @param Reps    numeric The relative tolerance for convergence
 * @param Rrhofun integer The rho function to use (see Control.h for possible values).
 *                        If the selected rho function is unknown, the bisquare function
 *                        will be used by default.
 *
 * @return numeric The M-scale
 */
SEXP C_mscale(SEXP Rvalues, SEXP Rlength, SEXP Rb, SEXP Rcc, SEXP RmaxIt, SEXP Reps, SEXP Rrhofun);

/**
 * Calculate the robust M-scale for the given values `x` of length `n`
 */
double mscale(const double *RESTRICT x, const int n, const double b, const double eps,
              const int maxIt, RhoFunction rho, const double cc);

#ifdef __cplusplus
}
#endif

#endif /* mscale_h */
