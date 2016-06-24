//
//  ElasticNet.cpp
//  pense
//
//  Created by David Kepplinger on 2016-01-31.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//
#include "config.h"

#include <cfloat>
#include <Rmath.h>

#include <RcppArmadillo.h>

#include "Control.h"
#include "ElasticNet.hpp"
#include "olsreg.h"

using namespace arma;

/**
 * Static inline functions only used in this compile unit
 */
static inline double softThreshold(const double z, const double gamma);
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
 * Convenience functions to choose one of the two possible
 * Elastic Net implementations
 *
 ***************************************************************************/
ElasticNet* getElasticNetImpl(const ENAlgorithm enAlgorithm, const double eps, const bool center,
							  const int maxIt)
{
    ElasticNet* en = NULL;
    ElasticNetLARS::UseGram useGram = ElasticNetLARS::AUTO;

    switch (enAlgorithm) {
        case GRADIENT_DESCENT:
            en = new ElasticNetGDESC(maxIt, eps, center);
            break;
        case AUGMENTED_LARS_GRAM:
            useGram = ElasticNetLARS::YES;
            break;
        case AUGMENTED_LARS_NOGRAM:
            useGram = ElasticNetLARS::NO;
            break;
        case AUGMENTED_LARS_AUTO:
        default:
            useGram = ElasticNetLARS::AUTO;
            break;
    }

    if (enAlgorithm != GRADIENT_DESCENT) {
        en = new ElasticNetLARS(eps, center, useGram);
    }

    return en;
}

ElasticNet* getElasticNetImpl(const Control& ctrl)
{
    return getElasticNetImpl(ctrl.enAlgorithm, ctrl.enEPS, ctrl.enCentering, ctrl.enMaxIt);
}


/****************************************************************************
 *
 * Elastic Net using the Coordinatewise Gradient Descent method
 *
 *
 ***************************************************************************/

ElasticNetGDESC::ElasticNetGDESC(const int maxIt, const double eps, const bool center) :
			ElasticNet(eps, center), maxIt(maxIt), XtrSize(0), XmeansSize(0)
{}

ElasticNetGDESC::~ElasticNetGDESC()
{
    if(this->XtrSize > 0) {
        delete[] this->Xtr;
    }

    if(this->XmeansSize > 0) {
        delete[] this->Xmeans;
        delete[] this->Xvars;
    }
}

void ElasticNetGDESC::setLambdas(const double lambda1, const double lambda2)
{
    this->lambda = lambda2 + lambda1;
    this->alpha = lambda1 / (lambda2 + lambda1);
    if (lambda2 == 0 && lambda1 == 0) {
        this->alpha = 0;
    }
}

void ElasticNetGDESC::setAlphaLambda(const double alpha, const double lambda)
{
    this->alpha = alpha;
    this->lambda = lambda;
}


bool ElasticNetGDESC::computeCoefsWeighted(const Data& data, double *RESTRICT coefs,
                                           double *RESTRICT residuals,
                                           const double *RESTRICT weights,
                                           const bool warm)
{

    throw Rcpp::exception("Weighted EN is not yet implemented with coordinate-descent");

    return true;
}

bool ElasticNetGDESC::computeCoefs(const Data& data, double *RESTRICT coefs,
                                   double *RESTRICT residuals, const bool warm)
{
    const double la = (this->lambda * this->alpha);
    const double updateDenom = this->lambda * (1 - this->alpha);

    double yMean;
    double tmp;
    double coefChange;
    double totalChange;
    const double *RESTRICT yPtr = data.getYConst();
    const double *RESTRICT XtrConstIter = data.getXtrConst();
    const double *RESTRICT weight;
    const double *RESTRICT actualXtr;
    int i, j, iter;
    double centerDenom = 0;
    double *RESTRICT XtrIter;
    double norm = 0;
    double newNorm;

    this->resizeBuffer(data);

    /**
     * We can take the given coefficients as warm start.
     */
    if (!warm) {
        memset(coefs, 0, data.numVar() * sizeof(double));
    }

    /*
     * First we calculate the mean of y and the X variables
     */
    if (this->center) {
        yMean = 0;
        memset(this->Xmeans, 0, (data.numVar() - 1) * sizeof(double));

        weight = data.getXtrConst();
        for (i = 0; i < data.numObs(); ++i, weight += data.numVar()) {
            if (*weight > 0) {
                yMean += yPtr[i] * (*weight);

                /*
                 * If we deal with augmented data, the first column is 1's and 0's.
                 * Count only those rows which have a leading 1.
                 */
                centerDenom += (*weight);

                ++XtrConstIter; /* Skip intercept */
                for (j = 0; j < data.numVar() - 1; ++j, ++XtrConstIter) {
                    this->Xmeans[j] += (*XtrConstIter) * (*weight);
                }
            } else {
                XtrConstIter += data.numVar();
            }
        }

        if (centerDenom > 0) {
            yMean /= centerDenom;

            for (j = 0; j < data.numVar() - 1; ++j) {
                this->Xmeans[j] /= centerDenom;
            }
        } else {
            yMean = 0;
            memset(this->Xmeans, 0, (data.numVar() - 1) * sizeof(double));
        }

        /*
         * We start by "setting" all beta_j's to 0, i.e. the residuals
         * are the deviation from the average.
         *
         * And by centering the variables (at least for those observations with weight > 0)
         */
        XtrConstIter = data.getXtrConst();
        XtrIter = this->Xtr;
        weight = data.getXtrConst();
        for (i = 0; i < data.numObs(); ++i, weight += data.numVar()) {
            residuals[i] = yPtr[i] - yMean * (*weight);

            /* Skip the first row (intercept!) */
            ++XtrConstIter;
            ++XtrIter;
            for (j = 0; j < data.numVar() - 1; ++j, ++XtrConstIter, ++XtrIter) {
                *XtrIter = (*XtrConstIter) - this->Xmeans[j] * (*weight);
            }
        }

        actualXtr = this->Xtr;
    } else {
        yMean = 0;
        memset(this->Xmeans, 0, (data.numVar() - 1) * sizeof(double));
        memcpy(residuals, yPtr, data.numObs() * sizeof(double));
        actualXtr = data.getXtrConst();
    }

    /*
     * Residuals are now either y or y - mean(y).
     * If we use a warm start, we have to update the residuals
     * and calculate the norm of the coefficient vector
     */
    if (warm) {
        XtrConstIter = actualXtr;
        for (i = 0; i < data.numObs(); ++i) {
            /* Skip the first row (intercept!) */
            ++XtrConstIter;
            for (j = 1; j < data.numVar(); ++j, ++XtrConstIter) {
                residuals[i] -= (*XtrConstIter) * coefs[j];
            }
        }

        /* Compute norm */
        for (j = 1; j < data.numVar(); ++j) {
            norm += fabs(coefs[j]);
        }
    }

    /*
     * Compute length of the vectors of variables (= N * Var(X_j))
     */
    memset(this->Xvars, 0, (data.numVar() - 1) * sizeof(double));
    XtrConstIter = actualXtr;
    for (i = 0; i < data.numObs(); ++i) {
        ++XtrConstIter;
        for (j = 0; j < data.numVar() - 1; ++j, ++XtrConstIter) {
            this->Xvars[j] += (*XtrConstIter) * (*XtrConstIter);
        }
    }

    for (j = 0; j < data.numVar() - 1; ++j, ++XtrConstIter) {
        this->Xvars[j] /= data.numObs();
    }

    iter = 0;

    while(1) {
        totalChange = 0;
        coefChange = 0;
        newNorm = 0;

        /* Start at j = 1 because j = 0 is the intercept which will be calculated at the end */
        for (j = 1; j < data.numVar(); ++j) {

            /* Update coefficient beta_j */
            XtrConstIter = actualXtr + j;
            tmp = 0;

            for (i = 0; i < data.numObs(); ++i, XtrConstIter += data.numVar()) {
                tmp += *XtrConstIter * residuals[i];
            }


            tmp = softThreshold(tmp / data.numObs() + coefs[j] * this->Xvars[j - 1], la) /
                        (this->Xvars[j - 1] + updateDenom);


            coefChange = coefs[j] - tmp;
            coefs[j] = tmp;

            newNorm += fabs(tmp);

            /* Update residuals if the coefficient has changed */
            if (coefChange != 0) {
                XtrConstIter = actualXtr + j;
                for (i = 0; i < data.numObs(); ++i, XtrConstIter += data.numVar()) {
                    residuals[i] += (*XtrConstIter) * coefChange;
                }

                totalChange += fabs(coefChange);
            }
        }

        /*
         * Check for max. iterations
         */
        if ((++iter > this->maxIt) || (totalChange <= this->eps * norm)) {
            break;
        }

        norm = newNorm;
    }

    /* Update intercept */
    coefs[0] = yMean;
    for (j = 1; j < data.numVar(); ++j) {
        coefs[0] -= coefs[j] * this->Xmeans[j - 1];
    }

    /* Residuals are already adjusted! */
    return (iter <= this->maxIt);
}


void ElasticNetGDESC::resizeBuffer(const Data& data)
{
    if (this->center && ((data.numObs() * data.numVar()) > (this->XtrSize))) {
        if(this->XtrSize > 0) {
            delete[] this->Xtr;
        }
        this->XtrSize = data.numObs() * data.numVar();
        this->Xtr = new double[this->XtrSize];
    }

    if (data.numVar() > this->XmeansSize) {
        if(this->XmeansSize > 0) {
            delete[] this->Xmeans;
            delete[] this->Xvars;
        }
        this->XmeansSize = data.numVar();
        this->Xmeans = new double[this->XmeansSize - 1];
        this->Xvars = new double[this->XmeansSize - 1];
    }
}


static inline double softThreshold(const double z, const double gamma)
{
    if (fabs(z) <= gamma) {
        return 0.;
    } else if (z < 0) {
        return z + gamma;
    }
    return z - gamma;
}



/****************************************************************************
 *
 * Elastic Net using the Least Angular Regression method
 * 
 * This augments the data matrix with a diagonal matrix and then computes
 * the LASSO solution.
 *
 ***************************************************************************/
ElasticNetLARS::ElasticNetLARS(const double eps, const bool center, const UseGram useGram) :
        ElasticNet(eps, center), lambda1(0), sqrtLambda2(0), gramMode(useGram), augNobs(0)
{
}


ElasticNetLARS::~ElasticNetLARS()
{
}

void ElasticNetLARS::setLambdas(const double lambda1, const double lambda2)
{
    this->lambda1 = lambda1;
    this->sqrtLambda2 = sqrt(lambda2);
}

void ElasticNetLARS::setAlphaLambda(const double alpha, const double lambda)
{
    this->setLambdas(lambda * alpha, lambda * (1 - alpha));
}

void ElasticNetLARS::setThreshold(const double eps)
{
    /* Don't allow adjusting the threshold */
}


bool ElasticNetLARS::computeCoefsWeighted(const Data& data, double *RESTRICT coefs,
                                          double *RESTRICT resids,
                                          const double *RESTRICT weights,
                                          const bool warm)
{
    uword i;
    vec sqrtWeights(data.numObs());
    double recipSumWeights = 0; /* This is also the squared L2 norm of the sqrtWeights */

    /*
     * First augment data
     */
    this->augmentData(data);

    /*
     * Second weight the observations
     *
     * If we assume the data is centered already -- we only need to weight the observations
     * Otherwise, we also have make the design orthogonal to (y - beta0)
     */
    mat tmpMat;
    mat& XtrWeighted = (this->center ? tmpMat : this->XtrAug);

    if (this->center) {
        XtrWeighted.set_size(data.numVar(), data.numObs());
    }

    recipSumWeights = 0;
    for (i = 0; i < (uword) data.numObs(); ++i) {
        sqrtWeights[i] = sqrt(weights[i]);
        recipSumWeights += weights[i];

        this->yAug[i] *= sqrtWeights[i];
        XtrWeighted.unsafe_col(i) = sqrtWeights[i] * this->XtrAug.unsafe_col(i);
    }

    if (this->center) {
        recipSumWeights = 1 / recipSumWeights;
        this->XtrAug.cols(0, data.numObs() - 1) = XtrWeighted -
                            XtrWeighted * sqrtWeights * sqrtWeights.t() * recipSumWeights;

        XtrWeighted.reset();
    }

    /*
     * Then perform LASSO on the augmented data (the algorithm is aware of the augmentation)
     * We don't have an intercept (due to the orthogonalization)
     */
    vec beta(coefs, data.numVar(), false, true);
    vec residuals(resids, data.numObs(), false, true);

    this->augmentedLASSO(beta, residuals, data.numObs(), false);


    /*
     * Recover intercept if requested
     */
    if (this->center) {
        coefs[0] = 0;

        computeResiduals(data.getXtrConst(), data.getYConst(), data.numObs(), data.numVar(), coefs,
                         resids);

        for (i = 0; i < (uword) data.numObs(); ++i) {
            coefs[0] += weights[i] * resids[i];
        }

        coefs[0] *= recipSumWeights;

        residuals -= coefs[0];
    }

    return true;
}


bool ElasticNetLARS::computeCoefs(const Data& data, double *RESTRICT coefs,
								  double *RESTRICT resids, const bool warm)
{
    vec beta(coefs, data.numVar(), false, true);
    vec residuals(resids, data.numObs(), false, true);

    if (data.numObs() == 0 || data.numVar() == 0) {
        memset(coefs, 0, data.numVar() * sizeof(double));
        memcpy(resids, data.getYConst(), data.numObs() * sizeof(double));
    } else if (data.numVar() == 1) {
        /* Data has only the first column of 1's – compute mean and return */
        memcpy(residuals.memptr(), data.getYConst(), data.numObs() * sizeof(double));
        beta.zeros();
        beta[0] = mean(residuals);
        residuals -= beta[0];
        return true;
    }

    /*
     * First augment data
     */
    this->augmentData(data);

    /*
     * Then perform LASSO on the augmented data (the algorithm is aware of the augmentation)
     */

    this->augmentedLASSO(beta, residuals, data.numObs(), this->center);

	return true;
}


/**
 * This algorithm is an extended version of the "fastLasso" algorithm written by
 * Andreas Alfons for the R package robustHD-0.5.1
 */
void ElasticNetLARS::augmentedLASSO(vec& beta, vec& residuals, const uword nobs,
                                    const bool intercept)
{
    /*
     * START of fastLasso algorithm from Andreas Alfons from
     * R package robustHD
     */
    const double rescaledLambda = this->augNobs * this->lambda1;

//    vec beta(coefs, this->XtrAug.n_rows, false, true);
//    vec residuals(resids, nobs, false, true);

    uword i, j;
    sword ri;

	double meanY;

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
     * Center data if requested
     */
	if(intercept) {
		this->meanX = mean(this->XtrAug.cols(0, nobs - 1), 1);
        this->meanX[0] = 0; // Don't change the intercept-column
        this->XtrAug.cols(0, nobs - 1).each_col() -= this->meanX;
		meanY = mean(this->yAug.rows(0, nobs - 1));
		this->yAug.rows(0, nobs - 1) -= meanY;
	} else {
		meanY = 0;		  // just to avoid warning, this is never used
//		intercept = 0;	  // zero intercept
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
        j = inactive(0);
        // set maximum step size in the direction of that variable
        double maxStep = this->corY(j);
        if(maxStep < 0) {
            maxStep = -maxStep; // absolute value
        }
        // compute coefficients for least squares solution
        vec betaLS = solve(this->XtrAug.unsafe_col(j), this->yAug);

        // compute lasso coefficients
        if(rescaledLambda < maxStep) {
            // interpolate towards least squares solution
            beta(j) = (maxStep - rescaledLambda) * betaLS[0] / maxStep;
        }
    } else if (this->XtrAug.n_rows > 2) {
        /*
         * further initializations for iterative steps
         */

        /*
         * previous and current regression coefficients
         */
        vec previousBeta = zeros(this->XtrAug.n_rows);
//            currentBeta = zeros(this->XtrAug.n_rows);



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
    residuals = this->yAug.rows(0, nobs - 1) - this->XtrAug.cols(0, nobs - 1).t() * beta;

    if(intercept) {
        beta(0) = meanY - dot(beta, this->meanX);
    }
}


void ElasticNetLARS::augmentData(const Data& data)
{
    uword i;
    bool resize = false;
    uword newNobs = this->XtrAug.n_cols,
          newNvar = this->XtrAug.n_rows;
    double *RESTRICT XtrAugIter;
    double tmp;

    if (((uword) data.numVar() != this->XtrAug.n_rows) || ((uword) data.numObs() != this->XtrAug.n_cols)) {
        newNobs = (uword) data.numObs();
        newNvar = (uword) data.numVar();
        this->augNobs = newNobs;
        resize = true;
    }

    if ((this->augNobs == newNobs) && (this->sqrtLambda2 > 0)) {
        /* Resize matrix to make room for augmented data */
        newNobs = newNobs + newNvar - 1;
        resize = true;
    } else if (this->sqrtLambda2 == 0) {
        /* No need to resize allocated memory */
        newNobs = this->augNobs;
        resize = true;
    }


    if (resize) {
        this->XtrAug.set_size(newNvar, newNobs);
        this->yAug.set_size(newNobs);

        if (newNobs > this->augNobs) {
            memset(this->XtrAug.memptr() + this->augNobs * this->XtrAug.n_rows, 0,
                   (this->XtrAug.n_rows - 1) * this->XtrAug.n_rows * sizeof(double));
       }
    }

    /*
     * Copy first n elements of y and set the remaining to 0
     * (this has to be done everytime because it will be overwritten by
     * LARS algorithm)
     * The augmented part of Xtr won't be changed, so we don't have to re-initialize
     * it.
     */
    memcpy(this->yAug.memptr(), data.getYConst(), this->augNobs * sizeof(double));

    if (newNobs > this->augNobs) {
        memset(this->yAug.memptr() + this->augNobs, 0, (this->XtrAug.n_rows - 1) * sizeof(double));
    }

    memcpy(this->XtrAug.memptr(), data.getXtrConst(),
           this->augNobs * this->XtrAug.n_rows * sizeof(double));

    if (this->sqrtLambda2 > 0) {
        tmp = sqrt(this->augNobs) * this->sqrtLambda2;

        XtrAugIter = this->XtrAug.memptr() + this->augNobs * this->XtrAug.n_rows + 1;
        for (i = 1; i < this->XtrAug.n_rows; ++i) {
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

