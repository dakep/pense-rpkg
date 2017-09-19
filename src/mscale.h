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
 * Calculate the robust M-scale for the given values `x` of length `n`
 */
double mscale(const double *RESTRICT x, const int n, const double b, const double eps,
              const int maxIt, RhoFunction rho, const double cc);

#ifdef __cplusplus
}
#endif

#endif /* mscale_h */
