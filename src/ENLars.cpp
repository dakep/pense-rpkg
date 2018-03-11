//
//  ENLars.cpp
//  pense
//
//  Created by David Kepplinger on 2016-01-31.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//  Parts of the augmented LARS algorithm are Copyright (C) Andreas Alfons <alfons@ese.eur.nl>
#include "config.h"

#include <cfloat>
#include <RcppArmadillo.h>
#include <Rmath.h>

#include "ENLars.hpp"
#include "olsreg.h"

/**
 *
 */
static const ENLars::UseGram DEFAULT_OPT_GRAM_MODE = ENLars::AUTO;
static const double DEFAULT_OPT_EPS = 1e-8;


using namespace arma;

/**
 * Static inline functions only used in this compile unit
 */
static inline sword sign(const double x);
static inline double findStepGram(const double corActiveY, const vec& corY, const uvec& inactive,
                                  const double corActiveEqui, const vec& corEquiV,
                                  const double eps);
static inline double findStep(const double corActiveY, const vec& corY, const uvec& inactive,
                              const double corActiveEqui, const vec& corInactiveEquiV,
                              const double eps);
static inline uvec findDrops(const vec& beta, const uvec& active, const vec& w,
                             const double eps, double& step);

static void downdateCholesky(const uvec& drops, mat *const L, uword& nrActive, uword& rank);


/****************************************************************************
 *
 * Elastic Net using the Least Angular Regression method
 *
 * This augments the data matrix with a diagonal matrix and then computes
 * the LASSO solution.
 *
 ***************************************************************************/
ENLars::ENLars(const bool intercept) : ElasticNet(intercept),
            lambda1(0), sqrtLambda2(0),
            eps(DEFAULT_OPT_EPS),
            gramMode(DEFAULT_OPT_GRAM_MODE)
{
}

ENLars::ENLars(const bool intercept, const Options& options) : ElasticNet(intercept),
            lambda1(0), sqrtLambda2(0),
            eps(DEFAULT_OPT_EPS),
            gramMode(DEFAULT_OPT_GRAM_MODE)
{
    this->setOptions(options);
}


ENLars::~ENLars()
{
}

void ENLars::setOptions(const Options& options)
{
    this->gramMode = options.get("gramMode", this->gramMode);
    this->eps = options.get("eps", this->eps);
}

void ENLars::setLambdas(const double lambda1, const double lambda2)
{
    this->lambda1 = lambda1;
    this->sqrtLambda2 = sqrt(lambda2);
}

void ENLars::setAlphaLambda(const double alpha, const double lambda)
{
    this->setLambdas(lambda * alpha, lambda * (1 - alpha));
}

void ENLars::computeCoefsWeighted(double& intercept, sp_vec &coefs, vec &residuals, const vec &weights)
{
    const uword trueNobs = this->yOrig.n_elem;
    const uword nvar = this->XtrAug.n_rows;
    const vec sqrtWeights = sqrt(weights);
    const double recipSumWeights = 1 / accu(weights);

    /*
     * First augment data
     */
    this->augmentData();

    /*
     * Second weight the observations
     *
     * If we assume the data is centered already -- we only need to weight the observations
     * Otherwise, we also have make the design orthogonal to (y - beta0)
     */
    const mat XtrAugOrig = this->XtrAug;

    this->yAug.head(trueNobs) %= sqrtWeights;

    if (this->intercept) {
        mat tmp = this->XtrAug.head_cols(trueNobs).each_row() % sqrtWeights.t();
        this->XtrAug.head_cols(trueNobs) = tmp -
                    tmp * sqrtWeights * sqrtWeights.t() * recipSumWeights;
    } else {
        this->XtrAug.head_cols(trueNobs).each_row() %= sqrtWeights.t();
    }

    this->XtrAug.row(0).zeros();

    /*
     * Then perform LASSO on the augmented data (the algorithm is aware of the augmentation)
     * We don't have an intercept (due to the orthogonalization)
     */
    vec coefsDense(coefs.n_elem + 1, fill::zeros);

    /*
     * Then perform LASSO on the augmented data (the algorithm is aware of the augmentation)
     */
    if (this->lambda1 > 0) {
        this->augmentedLASSO(coefsDense, residuals, trueNobs, false);
    } else {
        this->augmentedOLS(coefsDense, residuals, trueNobs, false);
    }

    intercept = coefsDense[0];
    coefs = sp_vec(coefsDense.tail_rows(nvar - 1));

    /* Restore original data */
    this->XtrAug = XtrAugOrig;

    residuals = this->yOrig - conv_to< colvec >::from(
        coefsDense.t() * this->XtrAug.head_cols(trueNobs)
    );

    /*
     * Recover intercept if requested
     */
    if (this->intercept) {
        residuals -= dot(this->meanX, coefsDense);
        intercept = accu(weights % residuals) * recipSumWeights;
        residuals -= intercept;
    }
}

void ENLars::computeCoefsWeighted(double *RESTRICT coefs, double *RESTRICT resids,
                                  const double *RESTRICT weights)
{
    const uword trueNobs = this->yOrig.n_elem;
    const uword nvar = this->XtrAug.n_rows;
    uword i;
    vec sqrtWeights(trueNobs);
    double recipSumWeights = 0; /* This is also the squared L2 norm of the sqrtWeights */

    /*
     * First augment data
     */
    this->augmentData();

    /*
     * Second weight the observations
     *
     * If we assume the data is centered already -- we only need to weight the observations
     * Otherwise, we also have make the design orthogonal to (y - beta0)
     */
    mat origDataCopy = this->XtrAug;
    mat tmpMat;
    mat& XtrWeighted = (this->intercept ? tmpMat : this->XtrAug);

    if (this->intercept) {
        XtrWeighted.set_size(nvar, trueNobs);
    }

    recipSumWeights = 0;
    for (i = 0; i < trueNobs; ++i) {
        sqrtWeights[i] = sqrt(weights[i]);
        recipSumWeights += weights[i];

        this->yAug[i] *= sqrtWeights[i];
        XtrWeighted.unsafe_col(i) = sqrtWeights[i] * this->XtrAug.unsafe_col(i);
    }

    if (this->intercept) {
        recipSumWeights = 1 / recipSumWeights;
        this->XtrAug.head_cols(trueNobs) = XtrWeighted -
                    XtrWeighted * sqrtWeights * sqrtWeights.t() * recipSumWeights;
    }

    this->XtrAug.row(0).zeros();

    /*
     * Then perform LASSO on the augmented data (the algorithm is aware of the augmentation)
     * We don't have an intercept (due to the orthogonalization)
     */
    vec beta(coefs, nvar, false, true);
    vec residuals(resids, trueNobs, false, true);

    beta.zeros();

    if (this->lambda1 > 0) {
        this->augmentedLASSO(beta, residuals, trueNobs, false);
    } else {
        this->augmentedOLS(beta, residuals, trueNobs, false);
    }

    /* Restore original data */
    this->XtrAug = origDataCopy;

    residuals = this->yOrig - conv_to< colvec >::from(
        beta.t() * this->XtrAug.head_cols(trueNobs)
    );

    /*
     * Recover intercept if requested
     */
    if (this->intercept) {
        residuals -= dot(this->meanX, beta);

        coefs[0] = 0;

        for (i = 0; i < trueNobs; ++i) {
            coefs[0] += weights[i] * resids[i];
        }

        coefs[0] *= recipSumWeights;

        residuals -= coefs[0];
    }
}

void ENLars::computeCoefs(double& intercept, sp_vec& coefs, vec& residuals)
{
    const uword trueNobs = this->yOrig.n_elem;
    const uword nvar = this->XtrAug.n_rows;

    if (trueNobs == 0 || nvar == 0) {
        coefs.zeros();
        residuals = this->yOrig;
    } else if (nvar == 1) {
        /* Data has only the first column of 1's – compute mean and return */
        coefs.zeros();
        coefs[0] = mean(residuals);
        residuals -= coefs[0];
    } else {
        /*
         * First augment data
         */
        this->augmentData();

        vec coefsDense(coefs.n_elem + 1, fill::zeros);

        /*
         * Then perform LASSO on the augmented data (the algorithm is aware of the augmentation)
         */
        if (this->lambda1 > 0) {
            this->augmentedLASSO(coefsDense, residuals, trueNobs, this->intercept);
        } else {
            this->augmentedOLS(coefsDense, residuals, trueNobs, this->intercept);
        }

        intercept = coefsDense[0];
        coefs = sp_vec(coefsDense.tail_rows(nvar - 1));
    }
}


void ENLars::computeCoefs(double *RESTRICT coefs, double *RESTRICT resids)
{
    const uword trueNobs = this->yOrig.n_elem;
    const uword nvar = this->XtrAug.n_rows;
    vec beta(coefs, nvar, false, true);
    vec residuals(resids, trueNobs, false, true);

    if (trueNobs == 0 || nvar == 0) {
        memset(coefs, 0, nvar * sizeof(double));
        memcpy(resids, this->yOrig.memptr(), trueNobs * sizeof(double));
    } else if (nvar == 1) {
        /* Data has only the first column of 1's – compute mean and return */
        memcpy(residuals.memptr(), this->yOrig.memptr(), trueNobs * sizeof(double));
        beta.zeros();
        beta[0] = mean(residuals);
        residuals -= beta[0];
    } else {
        /*
         * First augment data
         */
        this->augmentData();

        /*
         * Then perform LASSO on the augmented data (the algorithm is aware of the augmentation)
         */
        if (this->lambda1 > 0) {
            this->augmentedLASSO(beta, residuals, trueNobs, this->intercept);
        } else {
            this->augmentedOLS(beta, residuals, trueNobs, this->intercept);
        }
    }
}


void ENLars::augmentedOLS(vec& coefs, vec& residuals, const uword nobs,
                          const bool intercept)
{
    if (!intercept) {
        coefs[0] = 0;
        coefs.tail_rows(coefs.n_elem - 1) = arma::solve(
            this->XtrAug.tail_rows(this->XtrAug.n_rows - 1).t(), this->yAug,
            arma::solve_opts::fast + arma::solve_opts::no_approx
        );
    } else {
        arma::solve(coefs, this->XtrAug.t(), this->yAug,
                    arma::solve_opts::fast + arma::solve_opts::no_approx);
    }

    residuals = this->yAug.head_rows(nobs) - this->XtrAug.head_cols(nobs).t() * coefs;

    if (intercept) {
        coefs[0] -= dot(this->meanX, coefs);
    }
}


/**
 * This algorithm is an extended version of the "fastLasso" algorithm from 
 * the R package robustHD-0.5.1: Copyright (C) Andreas Alfons <alfons@ese.eur.nl>
 */
void ENLars::augmentedLASSO(vec& beta, vec& residuals, const uword nobs, const bool intercept)
{
    /*
     * START of fastLasso algorithm from Andreas Alfons from
     * R package robustHD
     */
    const double rescaledLambda = nobs * this->lambda1;

    uword i, j;
    sword ri;

	double meanY = 0;

	uvec inactive(this->XtrAug.n_rows - 1),
         ignores;

    /*
     * Check if we should use the Gram matrix
     */
    bool useGram = (this->gramMode == YES ||
                    (this->gramMode == AUTO && this->XtrAug.n_rows <= MAX_PREDICTORS_GRAM));

    uword nrInactive = this->XtrAug.n_rows;
    uword nrActive = 0;	// number of active predictors

    /*
     * Center response if requested
     * (the predictors are already centered)
     */
	if(intercept) {
		meanY = mean(this->yAug.head_rows(nobs));
		this->yAug.head_rows(nobs) -= meanY;
	}

    /*
     * Initialize inactive set. At the beginng, all predictors are inactive,
     * but column 0 (intercept) is always active!
     */
    for (j = 1; j < this->XtrAug.n_rows; ++j) {
        inactive[j - 1] = j;
    }

	if(useGram) {
		this->gramMat = this->XtrAug * this->XtrAug.t();
	}

    /*
     * Compute current correlations
     */
    this->corY = this->XtrAug * this->yAug;

    /*
     * Reset coefficients
     */
    beta.zeros();


    if(this->XtrAug.n_rows == 2) {
        // special case of only one variable (with sufficiently large norm)
        // set maximum step size in the direction of that variable
        double maxStep = this->corY(1);
        if(maxStep < 0) {
            maxStep = -maxStep; // absolute value
        }
        // compute coefficients for least squares solution
        vec betaLS = solve(this->XtrAug.row(1).t(), this->yAug);

        // compute lasso coefficients
        if(rescaledLambda < maxStep) {
            // interpolate towards least squares solution
            beta(1) = (maxStep - rescaledLambda) * betaLS[0] / maxStep;
        }
    } else if (this->XtrAug.n_rows > 2) {
        /*
         * further initializations for iterative steps
         */

        /*
         * previous and current regression coefficients
         */
        vec previousBeta = zeros(this->XtrAug.n_rows);

        /* previous and current penalty parameter */
        double previousLambda = R_PosInf,
               currentLambda = R_PosInf;

        /* indicates variables to be dropped */
        uvec drops;

        /*
         * Cholesky L of Gram matrix of active variables
         * and rank of L
         */
        mat L;
        uword rank = 0;

        /* maximum number of variables to be sequenced */
        uword maxActive = std::min(this->XtrAug.n_cols - intercept, this->XtrAug.n_rows - 1);

        uword usableVariables = this->XtrAug.n_rows;
        uword newPred;
        vec xNewPred;

        uvec active;        // active predictors

        /*
         * keep track of sign of correlations for the active variables
         * (double precision is necessary for solve())
         */
        vec signs;

        double maxCor, tmp;

        /*
         * Modified LARS algorithm for lasso solution
         */
        while(nrActive < maxActive) {
            /*
             * Find maximum absolute correlation with Y for the inactive
             */
            maxCor = 0;
            for (i = 0; i < inactive.n_elem; ++i) {
                tmp = fabs(this->corY[inactive[i]]);
                if (maxCor < tmp) {
                    maxCor = tmp;
                }
            }

            if(nrActive == 0) {	// no active variables
                previousLambda = maxCor;
            } else {
                previousLambda = currentLambda;
            }
            currentLambda = maxCor;

            /*
             * Check if we are already at the requested lambda
             */
            if(currentLambda <= rescaledLambda) {
                break;
            }

            if(drops.n_elem == 0) {
                for(j = 0; j < inactive.n_elem; ++j) {
                    newPred = inactive[j];
                    if (maxCor - eps > fabs(this->corY[newPred])) {
                        /*
                         * This one won't be activated -- skip
                         */
                        continue;
                    }

                    /*
                     * update Cholesky L of Gram matrix of active variables
                     * this cannot be put into its own void function since
                     * insert_rows() doesn't work with referenced matrices
                     */

                    double newX, normNewX;
                    if(useGram) {
                        xNewPred = this->gramMat.unsafe_col(newPred);	// reuses memory
                        newX = xNewPred[newPred];
                    } else {
                        xNewPred = conv_to<colvec>::from(this->XtrAug.row(newPred));
                        newX = dot(xNewPred, xNewPred);
                    }

                    normNewX = sqrt(newX);

                    if(nrActive == 0) {	// no active variables, L is empty
                        L.set_size(1, 1);
                        L[0] = normNewX;
                        rank = 1;
                    } else {
                        vec oldX;
                        if(useGram) {
                            oldX = xNewPred.elem(active);
                        } else {
                            oldX = this->XtrAug.rows(active) * xNewPred;
                        }
                        vec l = solve(trimatl(L), oldX);
                        tmp = newX - dot(l, l);

                        /*
                         * check if L is machine singular
                         */
                        if(tmp > eps) {
                            /*
                             * no singularity: update Cholesky L
                             */
                            tmp = sqrt(tmp);
                            ++rank;

                            /*
                             * add new row and column to Cholesky L
                             * this is a little trick: sometimes we need
                             * lower triangular matrix, sometimes upper
                             * hence we define a quadratic matrix and use
                             * triangularView() to interpret the matrix in
                             * correct way
                             * insert row and column without initializing memory
                             * (set_size() and reshape() have strange behavior)
                             */
                            L.insert_rows(nrActive, 1, false);
                            L.insert_cols(nrActive, 1, false);
                            // fill new parts of the matrix
                            for(i = 0; i < nrActive; ++i) {
                                L(nrActive, i) = l(i);
                                L(i, nrActive) = l(i);
                            }
                            L(nrActive, nrActive) = tmp;
                        }
                    }

                    /*
                     * add new variable to active set or drop it for good
                     * in case of singularity
                     */
                    if(rank == nrActive) {
                        /*
                         * singularity: drop variable for good
                         * --> - increase number of ignored variables
                         *     - decrease number of usable variables
                         *     - adjust maximum number of active variables
                         */
                        ignores.insert_rows(ignores.n_elem, 1, false);
                        ignores(ignores.n_elem - 1) = newPred;
                        --usableVariables;
                        if(usableVariables < maxActive) {
                            maxActive = usableVariables;
                        }
                    } else {
                        /*
                         * no singularity: add variable to active set
                         */
                        active.insert_rows(nrActive, 1, false);
                        active(nrActive) = newPred;

                        /*
                         * keep track of sign of correlation for new active variable
                         */
                        signs.insert_rows(nrActive, 1, false);
                        signs(nrActive) = sign(this->corY(newPred));

                        ++nrActive;
                    }
                }

                for(ri = inactive.n_elem - 1; ri >= 0; --ri) {
                    i = inactive(ri);
                    if (maxCor - eps <= fabs(this->corY(i))) {
                        /*
                         * This one was activated -- remove from inactivated list
                         */
                        inactive.shed_row(ri);
                    }
                }

                nrInactive = inactive.n_elem;	// update number of inactive variables
            }


            /*
             * Prepare for computation of step size
             */
            vec b = solve(trimatl(L), signs); // (here double prec. of signs is required)
            vec w = solve(trimatu(L), b);
            vec equiangularVec;
            vec corEquiV;
            // correlations of active variables with equiangular vector
            double corActiveEqui = 1 / sqrt(dot(w, signs));

            /*
             * Coefficients of active variables in linear combination forming the
             * equiangular vector.
             *
             * Note that this has the right signs
             */
            w *= corActiveEqui;

            if(useGram) {
                /*
                 * If we use the Gram matrix, we need the correlations between the
                 * equiangular vector and ALL variables
                 */
                corEquiV = this->gramMat.cols(active) * w;
            } else {
                /*
                 * we only need equiangular vector if we don't use the precomputed
                 * Gram matrix, otherwise we can compute the correlations directly
                 * from the Gram matrix
                 */
                equiangularVec = this->XtrAug.rows(active).t() * w;
            }

            // compute step size in equiangular direction
            double step;
            if(nrActive < maxActive) {
                if(useGram) {
                    /*
                     * compute step size in the direction of the equiangular vector
                     * Here we only need the correlations for the INACTIVE variables
                     */
                    step = findStepGram(maxCor, this->corY, inactive, corActiveEqui, corEquiV,
                                        this->eps);
                } else {
                    corEquiV = this->XtrAug.rows(inactive) * equiangularVec;

                    // compute step size in the direction of the equiangular vector
                    step = findStep(maxCor, this->corY, inactive, corActiveEqui, corEquiV, this->eps);
                }
            } else {
                // last step: take maximum possible step
                step = maxCor / corActiveEqui;
            }

            /*
             * Adjust step size if any sign changes and drop corresponding variables
             */
            drops = findDrops(beta, active, w, this->eps, step);

            /*
             * Update regression coefficients
             */
            previousBeta = beta;
            beta.elem(active) += step * w;

            /*
             * Update correlations
             */
            if(useGram) {
                /*
                 * Here we also need the correlation for active variables, since they may be
                 * dropped at a later stage
                 */
                this->corY -= step * corEquiV;
            } else {
                this->yAug -= step * equiangularVec;	// take step in equiangular direction
                this->corY = this->XtrAug * this->yAug;
            }

            /*
             * Drop variables if necessary
             */
            if(drops.n_elem > 0) {
                uword newInactive, drop;

                /*
                 * Downdate Cholesky decomposition `L`
                 */
                downdateCholesky(drops, &L, nrActive, rank);

                // mirror lower triangular part
                L = symmatu(L);

                /*
                 * add dropped variables to inactive set and make sure
                 * coefficients are 0
                 */
                inactive.insert_rows(nrInactive, drops.n_elem, false);

                for(j = 0; j < drops.n_elem; ++j) {
                    newInactive = active(drops(j));
                    inactive(nrInactive + j) = newInactive;
                    beta(newInactive) = 0;	// make sure coefficient is 0
                }

                nrInactive = inactive.n_elem;	// update number of inactive variables

                /*
                 * drop variables from active set and sign vector
                 * number of active variables is already updated above
                 */
                for(ri = drops.n_elem - 1; ri >= 0; --ri) {
                    drop = drops(ri);	// index with respect to active set
                    /*
                     * drop variables from active set and sign vector
                     * number of active variables is already updated above
                     */
                    active.shed_row(drop);
                    signs.shed_row(drop);
                }
            }
        }

        /*
         * Interpolate coefficients for given lambda
         */
        if(nrActive == 0) {
            /*
             * lambda larger than largest lambda from steps of LARS algorithm
             */
            beta.zeros();
        } else {
            /*
             * penalty parameter within two steps
             */
            if(nrActive == maxActive) {
                /*
                 * Current coefficients are the least squares solution (in the
                 * high-dimensional case, as far along the solution path as possible)
                 * current and previous values of the penalty parameter need to be
                 * reset for interpolation
                 */
                previousLambda = currentLambda;
                currentLambda = 0;
            }

            /*
             * interpolate coefficients
             */
            beta = ((rescaledLambda - currentLambda) * previousBeta +
                    (previousLambda - rescaledLambda) * beta) /
                        (previousLambda - currentLambda);
        }
    }

    /*
     * Compute residuals and intercept
     */
    residuals = this->yAug.head_rows(nobs) - this->XtrAug.head_cols(nobs).t() * beta;

    if(intercept) {
        beta(0) = meanY - dot(beta, this->meanX);
    }
}

void ENLars::setData(const Data& data)
{
    uword reqNobs = this->yOrig.n_elem,
          reqNvar = this->XtrAug.n_rows;

    if (((uword) data.numVar() != reqNvar) || ((uword) data.numObs() != reqNobs)) {
        /* The size of the data matrix changed --> we need to change the buffer as well */
        reqNobs = (uword) data.numObs();
        reqNvar = (uword) data.numVar();
    }

    /* Always copy the original y */
    this->yOrig = vec(data.getYConst(), data.numObs());

    if (this->sqrtLambda2 > 0) {
        /* A ridge penalty term is requested --> make room for the augmented data */
        reqNobs = data.numObs() + data.numVar() - 1;
    }

    if (reqNobs != this->XtrAug.n_cols) {
        this->XtrAug.set_size(reqNvar, reqNobs);
        this->yAug.set_size(reqNobs);

        if (this->sqrtLambda2 > 0 && reqNvar > 0) {
            this->XtrAug.tail_cols(reqNvar - 1).zeros();
       }
    }

    /*
     * Copy the X matrix
     */
    memcpy(this->XtrAug.memptr(), data.getXtrConst(),
            data.numObs() * data.numVar() * sizeof(double));

    /*
     * Center predictors if requested
     */
	if(this->intercept) {
		this->meanX = mean(this->XtrAug.head_cols(data.numObs()), 1);
        this->meanX[0] = 0; // Don't change the intercept-column
        this->XtrAug.head_cols(data.numObs()).each_col() -= this->meanX;
	}
}

void ENLars::augmentData()
{
    uword reqNobs = this->yOrig.n_elem,
          reqNvar = this->XtrAug.n_rows;
    double *RESTRICT XtrAugIter;
    double tmp;

    if (this->sqrtLambda2 > 0) {
        /* A ridge penalty term is requested --> make room for the augmented data */
        reqNobs += reqNvar - 1;
    }

    if (reqNobs != this->XtrAug.n_cols) {
        /*
         * This is only necessary if the L2 penalty changed between calls
         */
        this->XtrAug.resize(reqNvar, reqNobs);
        this->yAug.set_size(reqNobs);

        if (this->sqrtLambda2 > 0 && reqNvar > 0) {
            this->XtrAug.tail_cols(reqNvar - 1).zeros();
       }
    }

    /*
     * Copy first n elements of y and set the remaining to 0
     * The augmented part of Xtr won't be changed, so we don't have to re-initialize
     * it.
     */
    this->yAug.zeros();
    memcpy(this->yAug.memptr(), this->yOrig.memptr(), this->yOrig.n_elem * sizeof(double));

    /*
     * Set the diagonal to the penalty
     */
    if (this->sqrtLambda2 > 0) {
        tmp = sqrt((double) this->yOrig.n_elem) * this->sqrtLambda2;

        XtrAugIter = this->XtrAug.memptr() + this->yOrig.n_elem * this->XtrAug.n_rows + 1;
        for (uword i = 1; i < this->XtrAug.n_rows; ++i) {
            (*XtrAugIter) = tmp;
            XtrAugIter += this->XtrAug.n_rows + 1;
        }
    }
}

/**
 * get the sign of the double x (either -1 if x < 0, +1 if x > 0, 0 if x == 0)
 *
 * x ... the number to compute the sign from
 */
static inline sword sign(const double x) {
	return (x > 0) - (x < 0);
}

/**
 * compute step size in the direction of the equiangular vector
 *
 * corActiveY.. ... correlations of active variables with current response
 * corInactiveY ... correlations of inactive variables with current response
 * corActiveU ..... correlations of active variables with equiangular vector
 * corInactiveU ... correlations of inactive variables with equiangular vector
 * eps ............ small numerical value (effective zero)
 */
static inline double findStepGram(const double corActiveY, const vec& corY, const uvec& inactive,
                       const double corActiveEqui, const vec& corEquiV,
                       const double eps) {
    double step = corActiveY / corActiveEqui;      // maximum possible step;
    uword i;
    double corWY, corWEquiv;
    double tmp;

    for (i = 0; i < inactive.n_elem; ++i) {
        corWY = corY(inactive(i));
        corWEquiv = corEquiV(inactive(i));

        tmp = (corActiveY - corWY) / (corActiveEqui - corWEquiv);

        if ((tmp > eps) && (tmp < step)) {
            step = tmp;
        }


        tmp = (corActiveY + corWY) / (corActiveEqui + corWEquiv);

        if ((tmp > eps) && (tmp < step)) {
            step = tmp;
        }
    }

	return step;
}

/**
 * compute step size in the direction of the equiangular vector
 *
 * corActiveY.. ... correlations of active variables with current response
 * corInactiveY ... correlations of inactive variables with current response
 * corActiveU ..... correlations of active variables with equiangular vector
 * corInactiveU ... correlations of inactive variables with equiangular vector
 * eps ............ small numerical value (effective zero)
 */
static inline double findStep(const double corActiveY, const vec& corY, const uvec& inactive,
                       const double corActiveEqui, const vec& corInactiveEquiV,
                       const double eps) {
    double step = corActiveY / corActiveEqui;      // maximum possible step;
    uword i;
    double corWY;
    double tmp;

    for (i = 0; i < inactive.n_elem; ++i) {
        corWY = corY(inactive(i));

        tmp = (corActiveY - corWY) / (corActiveEqui - corInactiveEquiV(i));

        if ((tmp > eps) && (tmp < step)) {
            step = tmp;
        }


        tmp = (corActiveY + corWY) / (corActiveEqui + corInactiveEquiV(i));

        if ((tmp > eps) && (tmp < step)) {
            step = tmp;
        }
    }

	return step;
}


/**
 * adjust step size if any sign changes before the designated step size,
 * and return the corresponding variables to be dropped
 *
 * beta   ... current regression coefficients
 * active ... indices of inactive variables
 * w ........ coefficients of active variables in linear combination forming
 *            the equiangular vector
 * eps ...... small numerical value (effective zero)
 * step ..... step size in direction of equiangular vector
 */
static inline uvec findDrops(const vec& beta, const uvec& active, const vec& w,
                      const double eps, double& step) {
	/*
     * for each variable, compute step size where sign change would take place,
	 * and keep track of indices of variables that are potentially dropped
     */
	vec steps = -beta.elem(active) / w;
	uvec drops = find(steps > eps);
	if(drops.n_elem > 0) {
		// check if sign change occurs before the designated step size
		// if so, adjust step size and find variables to be dropped
		steps = steps.elem(drops);
		double smallestPositive = steps.min();
		if(smallestPositive < step) {
			step = smallestPositive;
			drops = drops.elem(find(steps == smallestPositive));
		} else drops.reset();
	}

	/*
     * if there are no sign changes or sign change would occur after the
	 * designated step size, an empty vector is returned
     */
	return drops;
}

/*
 * downdate Cholesky L
 * this cannot be put into its own void function since
 * shed_col() and shed_row() don't work with referenced matrices
 */
static inline void downdateCholesky(const uvec& drops, mat *const L, uword& nrActive, uword& rank)
{
    double a, b, tau, s, c;
    sword ri;
    uword i, j, drop;

    for(ri = drops.n_elem - 1; ri >= 0; --ri) {	// reverse order
        // variables need to be dropped in descending order
        drop = drops(ri);	// index with respect to active set

        /*
         * modify upper triangular part as in R package 'lars'
         * simply deleting columns is not enough, other computations
         * necessary but complicated due to Fortran code
         */
        L->shed_col(drop);
        --nrActive; // decrease number of active variables
        for(i = drop; i < nrActive; ++i) {
            a = L->at(i, i);
            b = L->at(i + 1, i);

            if(b != 0.0) {
                /* compute the rotation */
                if(std::abs(b) > std::abs(a)) {
                    tau = -a / b;
                    s = 1.0 / sqrt(1.0 + tau * tau);
                    c = s * tau;
                } else {
                    tau = -b / a;
                    c = 1.0 / sqrt(1.0 + tau * tau);
                    s = c * tau;
                }

                /*
                 * update 'L'
                 */
                L->at(i,i) = c * a - s * b;
                for(j = i + 1; j < nrActive; ++j) {
                    a = L->at(i, j);
                    b = L->at(i + 1, j);
                    L->at(i, j) = c * a - s * b;
                    L->at(i + 1, j) = s * a + c * b;
                }
            }
        }

        L->shed_row(nrActive);
        --rank;
    }
}

