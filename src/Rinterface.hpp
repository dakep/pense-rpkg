//
//  Rinterface.hpp
//  pense
//
//  Created by David Kepplinger on 2016-02-03.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#ifndef Rinterface_hpp
#define Rinterface_hpp

#include "config.h"
#include <RcppArmadillo.h>

extern "C" void R_init_pense(DllInfo *dll);

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
 * Compute the robust *size* of the vector x
 *
 * @param x numeric vector
 * @return numeric scalar robust size of vector x
 */
RcppExport SEXP C_tau_size(SEXP x);

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
RcppExport SEXP C_mscale(SEXP Rvalues, SEXP Rlength, SEXP Rb, SEXP Rcc, SEXP RmaxIt, SEXP Reps,
                         SEXP Rrhofun);

/**
 * Solve following minimzation problem:
 * (1 / (2*N)) * L2(y - beta0 - X . beta)^2 + lambda * ( ((1 - alpha)/2)*L2(beta)^2 + alpha*L1(beta) )
 *
 * It is the caller's responsibility to ensure that at least one predictor is present
 * in Xtr. The function will NOT check the arguments and will crash if Xtr does not
 * contain enough values.
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
 * @param Xtest       numeric numeric matrix (size `nvar` x `nobs-test`) WITHOUT intercept column
 *                            used to generate predictions. Can be NULL.
 *
 * @return List Returns a list with two elements:
 *		item 1: Integer telling if the status of the algorithm
 *		item 2: Error message explaining the status of the algorithm
 *      item 3: sparse matrix (dgCMatrix) with coefficient estimates
 *      item 4: Numeric matrix with the residuals
 *      item 5: numeric matrix with predictions for Xtest
 */
RcppExport SEXP C_elnet_sp(SEXP Xtr, SEXP y, SEXP coefs, SEXP alpha,
						   SEXP lambda, SEXP intercept, SEXP options,
                           SEXP Xtest);

/**
 * Solve following minimzation problem:
 * (1 / (2*N)) * L2(weights (y - beta0 - X . beta))^2 + lambda * ( ((1 - alpha)/2)*L2(beta)^2 + alpha*L1(beta) )
 *
 * It is the caller's responsibility to ensure that at least one predictor is present
 * in Xtr. The function will NOT check the arguments and will crash if Xtr does not
 * contain enough values.
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
 * @param Xtest       numeric numeric matrix (size `nvar` x `nobs-test`) WITHOUT intercept column
 *                            used to generate predictions. Can be NULL.
 *
 * @return List Returns a list with two elements:
 *		item 1: Integer telling if the status of the algorithm
 *		item 2: Error message explaining the status of the algorithm
 *      item 3: sparse matrix (dgCMatrix) with the cofficient estimates
 *      item 4: numeric matrix with the residuals
 *      item 5: numeric matrix with predictions for Xtest
 */
RcppExport SEXP C_elnet_weighted_sp(SEXP Xtr, SEXP y, SEXP weights, SEXP coefs,
								    SEXP alpha, SEXP lambda, SEXP intercept, SEXP enOptions,
                                    SEXP Xtest);


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
 * @param Xtr           numeric   transpose of the numeric X matrix (size `nvar` x `nobs`)
 * @param y             numeric   numeric y vector (size `nobs`)
 * @param alpha         numeric   alpha parameter for the EN penalty
 * @param lambda        numeric   lambda parameter for the EN penalty
 * @param intercept     numeric   initial intercept
 * @param coefs         dgCMatrix sparse vector of inital coefficients
 * @param penseOptions  List      options for PENSE
 * @param enOptions     List      options for the specific EN algorithm
 *
 * @return List Returns a list with 6 elements:
 *      item 1: intercept
 *      item 2: sparse coefficient vector
 *      item 3: residuals vector
 *      item 4: estimated scale of the residuals
 *      item 5: relative change in the last iteration
 *      item 6: number of iterations
 */
RcppExport SEXP C_pen_s_reg_sp(SEXP Xtr, SEXP y, SEXP intercept, SEXP coefs, SEXP alpha,
                               SEXP lambda, SEXP penseOptions, SEXP enOptions);

/**
 * @param Xtr           numeric   transpose of the numeric X matrix (size `nvar` x `nobs`)
 * @param y             numeric   numeric y vector (size `nobs`)
 * @param intercept     numeric   initial intercept
 * @param coefs         dgCMatrix sparse vector of inital coefficients
 * @param scale         numeric   (fixed) scale estimate
 * @param alpha         numeric   alpha parameter for penalization
 * @param lambda        numeric   lambda parameter(s) for penalization
 * @param msOptions     List      options for the M-Step
 * @param enOptions     List      options for the specific EN algorithm
 *
 * @return List Returns a list with four elements:
 *      item 1: intercept
 *      item 2: sparse coefficient vector (first element is the intercept)
 *      item 3: residuals vector
 *      item 4: relative change in the last iteration
 *      item 5: number of iterations
 */
RcppExport SEXP C_pen_mstep_sp(SEXP Xtr, SEXP y, SEXP intercept, SEXP coefs, SEXP scale,
                               SEXP alpha, SEXP lambda, SEXP msOptions, SEXP enOptions);

#endif /* Rinterface_hpp */
