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
 *  enCentering     integer (0/1)
 *
 *	mscaleB			numeric
 *	mscaleCC		numeric
 *	mscaleMaxIt		integer
 *	mscaleEPS		numeric
 *	mscaleRhoFun	integer
 */


/**
 * Solve following minimzation problem:
 * (1 / (2*N)) * RSS + lambda * ( ((1 - alpha)/2)*L1(beta) + alpha*L2(beta) )
 *
 *
 * @param Xtr       numeric The transpose of the numeric X matrix (size `nvar` x `nobs`)
 * @param y         numeric The numeric y vector (size `nobs`)
 * @param nobs      integer The number of observations
 * @param nvar      integer The number of variables (including the intercept)
 * @param alpha     numeric The alpha parameter for the penalization
 * @param lambda    numeric The lambda parameter for the penalization
 * @param maxIt     integer The maximum number of iterations allowed
 * @param eps       numeric The relative tolerance for convergence
 * @param centering integer Should centering be done (1=yes, 0=no) for rows with a leading
                            1.
 *
 * @return List Returns a list with two elements:
 *		item 1: Logical telling if the algorithm converged
 *      item 1: Numeric vector with the cofficient estimates
 *      item 2: Numeric vector with the residuals
 */
RcppExport SEXP C_elnet(SEXP Xtr, SEXP y, SEXP nobs, SEXP nvar, SEXP alpha, SEXP lambda,
                        SEXP maxIt, SEXP eps, SEXP centering);

/**
 * @param Xtr     numeric The transpose of the numeric X matrix (size `nvar` x `nobs`)
 * @param y       numeric The numeric y vector (size `nobs`)
 * @param nobs    integer The number of observations
 * @param nvar    integer The number of variables (including the intercept)
 *
 * @return Returns a numeric matrix of size `nobs` x `npscs`, where `npscs` is at
 *		   most `nvar`.
 */
RcppExport SEXP C_pscs2(SEXP Xtr, SEXP y, SEXP nobs, SEXP nvar);

/**
 * @param Xtr     numeric The transpose of the numeric X matrix (size `nvar` x `nobs`)
 * @param y       numeric The numeric y vector (size `nobs`)
 * @param nobs    integer The number of observations
 * @param nvar    integer The number of variables (including the intercept)
 * @param control List    The control list as described above
 *
 * @return List Returns a list with two elements:
 *      item 1: The numeric matrix of size `nvar` x (3 * `nvar` + 2)
 *      item 2: The value of the objective function for each initial estimate
 */
RcppExport SEXP C_enpy_rr(SEXP Xtr, SEXP y, SEXP nobs, SEXP nvar, SEXP control);

/**
 * @param Xtr     numeric The transpose of the numeric X matrix (size `nvar` x `nobs`)
 * @param y       numeric The numeric y vector (size `nobs`)
 * @param nobs    integer The number of observations
 * @param nvar    integer The number of variables (including the intercept)
 * @param control List    The control list as described above
 *
 * @return List Returns a list with two elements:
 *      item 1: The numeric matrix of size `nvar` x (3 * `nvar` + 2)
 *      item 2: The value of the objective function for each initial estimate
 */
RcppExport SEXP C_enpy_Mn(SEXP Xtr, SEXP y, SEXP nobs, SEXP nvar, SEXP control);

#endif /* Rinterface_hpp */
