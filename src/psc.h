//
//  psc.h
//  penseinit
//
//  Created by David Kepplinger on 2016-01-30.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#ifndef psc_h
#define psc_h

#include "AuxMemory.h"

/**
 * Calculate the principal sensitivity components (PSCs)
 *
 * @param pscs	 The (nobs by nvar) memory where the PSCs will be computed
 * @param Xtr	 The (nvar by nobs) transposed X matrix
 * @param y		 The (nobs) y vector
 * @param nobs	 The number of observations in X and y
 * @param nvar	 The number of variables in X
 * @param auxmem Auxilliary memory used during computation
 * @return Returns the number of PSCs or values less than 0 if an error
 *		   occured.
 */
int calculatePSCs(double *restrict pscs, AuxMemory* auxmem,
                  const double *restrict Xtr, const double *restrict y,
                  const int nobs, const int nvar);

#endif /* psc_h */
