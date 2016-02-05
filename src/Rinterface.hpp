//
//  Rinterface.hpp
//  penseinit
//
//  Created by David Kepplinger on 2016-02-03.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#ifndef Rinterface_hpp
#define Rinterface_hpp

#include <Rcpp.h>

/**
 * The control list must have following entries:
 *	lambda1			numeric
 *	lambda2			numeric
 *	numIt			integer
 *	eps				numeric
 *	residThreshold  numeric
 *	residProportion numeric
 *	pscProportion	numeric
 *
 *	enMaxIt			integer
 *	enEPS			numeric
 *
 *	mscaleB			numeric
 *	mscaleCC		numeric
 *	mscaleMaxIt		integer
 *	mscaleEPS		numeric
 *	mscaleRhoFun	integer
 */


/**
 * @param Xtr     numeric The transpose of the numeric X matrix (size `nvar` x `nobs`)
 * @param y       numeric The numeric y vector (size `nobs`)
 * @param nobs    integer The number of observations
 * @param nvar    integer The number of variables (including the intercept)
 * @param control List    The control list as described above
 *
 * @return numeric Returns the numeric matrix of size `nvar` x (3 * `nvar` + 2)
 */
RcppExport SEXP C_enpy_rr(SEXP Xtr, SEXP y, SEXP nobs, SEXP nvar, SEXP control);

/**
 * @param Xtr     numeric The transpose of the numeric X matrix (size `nvar` x `nobs`)
 * @param y       numeric The numeric y vector (size `nobs`)
 * @param nobs    integer The number of observations
 * @param nvar    integer The number of variables (including the intercept)
 * @param control List    The control list as described above
 *
 * @return numeric Returns the numeric matrix of size `nvar` x (3 * `nvar` + 2)
 */
RcppExport SEXP C_enpy_Mn(SEXP Xtr, SEXP y, SEXP nobs, SEXP nvar, SEXP control);

#endif /* Rinterface_hpp */
