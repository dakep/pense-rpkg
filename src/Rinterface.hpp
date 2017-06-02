//
//  Rinterface.hpp
//  pense
//
//  Created by David Kepplinger on 2016-02-03.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#ifndef Rinterface_hpp
#define Rinterface_hpp

#include <RcppArmadillo.h>

/**
 * Add a column of 1's to the matrix X and transpose the augmented
 * matrix
 *
 *
 * @param X    numeric The numeric matrix
 *
 * @return numeric matrix The matrix X with a column of 1's prepended and transposed.
						  The matrix is of size (ncol + 1) x nrow.
 */
RcppExport SEXP C_augtrans(SEXP X);

/**
 * Solve following minimzation problem:
 * (1 / (2*N)) * L2(y - beta0 - X . beta)^2 + lambda * ( ((1 - alpha)/2)*L2(beta)^2 + alpha*L1(beta) )
 *
 *
 * @param Xtr         numeric The transpose of the numeric X matrix (size `nvar` x `nobs`)
 * @param y           numeric The numeric y vector (size `nobs`)
 * @param coefs		  numeric The inital coefficients, will be copied if `warm` is 1, otherwise
 *							  not referenced.
 * @param alpha       numeric The alpha parameter for the penalization
 * @param lambda      numeric The lambda parameter for the penalization
 * @param intercept   integer Should an intercept be estimated be done (1=yes, 0=no)
 *                            considering only rows with a leading 1.
 * @param options     list    A list with options for the specific EN algorithm
 *
 * @return List Returns a list with two elements:
 *		item 1: Integer telling if the status of the algorithm
 *		item 2: Error message explaining the status of the algorithm
 *      item 3: Numeric vector with the cofficient estimates
 *      item 4: Numeric vector with the residuals
 */
RcppExport SEXP C_elnet(SEXP Xtr, SEXP y, SEXP coefs, SEXP alpha,
						SEXP lambda, SEXP intercept, SEXP options);

/**
 * Solve following minimzation problem:
 * (1 / (2*N)) * L2(weights (y - beta0 - X . beta))^2 + lambda * ( ((1 - alpha)/2)*L2(beta)^2 + alpha*L1(beta) )
 *
 *
 * @param Xtr         numeric The transpose of the numeric X matrix (size `nvar` x `nobs`)
 * @param y           numeric The numeric y vector (size `nobs`)
 * @param weights	  numeric The numeric weights vector (size `nobs`)
 * @param coefs		  numeric The inital coefficients, will be copied if `warm` is 1, otherwise
 *							  not referenced.
 * @param alpha       numeric The alpha parameter for the penalization
 * @param lambda      numeric The lambda parameter for the penalization
 * @param intercept   integer Should an intercept be estimated be done (1=yes, 0=no)
 *                            considering only rows with a leading 1.
 * @param enOptions     list    A list with options for the specific EN algorithm
 *
 * @return List Returns a list with two elements:
 *		item 1: Integer telling if the status of the algorithm
 *		item 2: Error message explaining the status of the algorithm
 *      item 3: Numeric vector with the cofficient estimates
 *      item 4: Numeric vector with the residuals
 */
RcppExport SEXP C_elnet_weighted(SEXP Xtr, SEXP y, SEXP weights, SEXP coefs,
								 SEXP alpha, SEXP lambda, SEXP intercept, SEXP enOptions);

/**
 * @param Xtr     numeric The transpose of the numeric X matrix (size `nvar` x `nobs`)
 * @param y       numeric The numeric y vector (size `nobs`)
 *
 * @return Returns a numeric matrix of size `nobs` x `npscs`, where `npscs` is at
 *		   most `nvar`.
 */
RcppExport SEXP C_pscs_ols(SEXP Xtr, SEXP y);

/**
 * @param Xtr           numeric The transpose of the numeric X matrix (size `nvar` x `nobs`)
 * @param y             numeric The numeric y vector (size `nobs`)
 * @param alpha         numeric The alpha parameter for the penalization
 * @param lambda        numeric The lambda parameter for the penalization
 * @param intercept     integer Should an intercept be estimated be done (1=yes, 0=no)
 *                              considering only rows with a leading 1.
 * @param enOptions     list    A list with options for the specific EN algorithm
 *
 * @return Returns a numeric matrix of size `nobs` x `npscs`, where `npscs` is at
 *		   most `nobs`.
 */
RcppExport SEXP C_pscs_en(SEXP Xtr, SEXP y, SEXP alpha, SEXP lambda,
                          SEXP intercept, SEXP enOptions);

/**
 * @param Xtr       numeric The transpose of the numeric X matrix (size `nvar` x `nobs`)
 * @param y         numeric The numeric y vector (size `nobs`)
 * @param pyOptions List    The control list for the PY algorithm
 *
 * @return List Returns a list with two elements:
 *      item 1: The numeric matrix of size `nvar` x (3 * `nvar` + 2)
 *      item 2: The value of the objective function for each initial estimate
 */
RcppExport SEXP C_py_ols(SEXP Xtr, SEXP y, SEXP pyOptions);

/**
 * @param Xtr       numeric The transpose of the numeric X matrix (size `nvar` x `nobs`)
 * @param y         numeric The numeric y vector (size `nobs`)
 * @param alpha         numeric The alpha parameter for the penalization
 * @param lambda        numeric The lambda parameter for the penalization
 * @param pyOptions List    The control list for the PY algorithm
 * @param enOptions List    The options for the EN algorithm
 *
 * @return List Returns a list with two elements:
 *      item 1: The numeric matrix of size `nvar` x (3 * `nvar` + 2)
 *      item 2: The value of the objective function for each initial estimate
 */
RcppExport SEXP C_enpy_rr(SEXP Xtr, SEXP y, SEXP alpha, SEXP lambda,
                          SEXP pyOptions, SEXP enOptions);

/**
 * @param Xtr       numeric The transpose of the numeric X matrix (size `nvar` x `nobs`)
 * @param y         numeric The numeric y vector (size `nobs`)
 * @param alpha         numeric The alpha parameter for the penalization
 * @param lambda        numeric The lambda parameter for the penalization
 * @param pyOptions List    The control list for the PY algorithm
 * @param enOptions List    The options for the specific EN algorithm
 *
 * @return List Returns a list with two elements:
 *      item 1: The numeric matrix of size `nvar` x (3 * `nvar` + 2)
 *      item 2: The value of the objective function for each initial estimate
 */
RcppExport SEXP C_enpy_exact(SEXP Xtr, SEXP y, SEXP alpha, SEXP lambda,
                             SEXP pyOptions, SEXP enOptions);


/**
 * @param Xtr           numeric The transpose of the numeric X matrix (size `nvar` x `nobs`)
 * @param y             numeric The numeric y vector (size `nobs`)
 * @param alpha         numeric The alpha parameter for the EN penalty
 * @param lambda        numeric The lambda parameter for the EN penalty
 * @param coefs         numeric The vector of inital coeffiecients (including the intercept)
 * @param penseOptions  List    The options for PENSE
 * @param enOptions     List    The options for the specific EN algorithm
 *
 * @return List Returns a list with 5 elements:
 *      item 1: coefficient vector (first element is the intercept)
 *      item 2: residuals vector
 *      item 3: estimate scale of the residuals
 *      item 4: relative change in the last iteration
 *      item 5: number of iterations
 */
RcppExport SEXP C_pen_s_reg(SEXP Xtr, SEXP y, SEXP coefs, SEXP alpha, SEXP lambda,
                            SEXP penseOptions, SEXP enOptions);


/**
 * @param Xtr           numeric The transpose of the numeric X matrix (size `nvar` x `nobs`)
 * @param y             numeric The numeric y vector (size `nobs`)
 * @param coefs         numeric The vector of inital coeffiecients (including the intercept)
 * @param scale         numeric The (fixed) scale estimate
 * @param alpha         numeric The alpha parameter for penalization
 * @param lambda        numeric The lambda parameter for penalization
 * @param msOptions     List    The options for the M-Step
 * @param enOptions     List    The options for the specific EN algorithm
 *
 * @return List Returns a list with four elements:
 *      item 1: coefficient vector (first element is the intercept)
 *      item 2: residuals vector
 *      item 3: relative change in the last iteration
 *      item 4: number of iterations
 */
RcppExport SEXP C_pen_mstep(SEXP Xtr, SEXP y, SEXP coefs, SEXP scale,
                            SEXP alpha, SEXP lambda, SEXP msOptions, SEXP enOptions);

#endif /* Rinterface_hpp */
