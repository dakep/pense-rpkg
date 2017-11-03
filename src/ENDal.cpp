//
//  ENDal.cpp
//  pense
//
//  Created by David Kepplinger on 2017-05-21.
//  Copyright © 2017 David Kepplinger. All rights reserved.
//
#include "config.h"

#include <cfloat>
#include <RcppArmadillo.h>
#include <Rmath.h>

#include "ENDal.hpp"

using namespace arma;

/**
 * Default argument values
 */
static const int    DEFAULT_OPT_VERBOSITY = 0;
static const int    DEFAULT_OPT_MAXIT = 100;
static const double DEFAULT_OPT_EPS = 1e-5;
static const bool   DEFAULT_OPT_WARM_START = true;
static const double DEFAULT_OPT_ETA_START = -1;
static const double DEFAULT_OPT_ETA_START_NUMERATOR = 0.01;
static const double DEFAULT_OPT_ETA_MULTIPLIER = 2;
static const ENDal::Preconditioner DEFAULT_OPT_PRECONDITIONER = ENDal::NONE;


/**
 * constants defining the decrease in the step size in each line search iteration
 * and the fraction of the decrease in the objective function predicted by
 * linear extrapolation that we will accept in the line search
 */
static const double LINESEARCH_STEPSIZE_MULTIPLIER = 0.8;  // 0 < x < 1
static const double LINESEARCH_STEP_CONDITION = 0.3;      // 0 < x < 0.5
static const int LINESEARCH_MAX_STEP = 200;
static const int SOLVE_PCG_MAXIT = 100;
static const int SOLVE_PCG_MAXIT_UPDATE = 35;

static const uword VEC_SOFTTHRESH_SMALL_SPARSE = 50;
static const double VEC_SOFTTHRESH_SPARSE = 0.1;


/**
 * Static inline functions only used in this compile unit
 */

/**
 * (Vectorized) soft threshold functions, computing
 * sign(z) * max(0, |z| - gamma)
 */
static inline double softThreshold(const double z, const double gamma);
static inline void vecSoftThreshold(vec& z, const double gamma);

/**
 * Optimized soft threshold function to compute
 * sign(z1 + c * z2) * max(0, |z1 + c * z2| - gamma)
 * where z1 is a sparse vector
 */
static inline void vecSoftThresholdInplace(sp_vec& zsoft, const sp_vec& z1, const double c, const vec&z2, const double gamma);
static inline void vecSoftThresholdInplaceSmallSparse(sp_vec& zsoft, const sp_vec& sp_z1, const double c, const vec& z2, const double gamma);
static inline void vecSoftThresholdInplaceLargeSparse(sp_vec& zsoft, const sp_vec& sp_z1, const double c, const vec& z2, const double gamma);
static inline void vecSoftThresholdInplaceNonSparse(sp_vec& zsoft, const sp_vec& sp_z1, const double c, const vec& z2, const double gamma);

static int solve_pcg_diag(const mat &A, vec &x, const vec &b, const double eps);
static int solve_pcg(const mat &A, vec &x, const vec &b, const mat &precond, const double eps);

static inline double squaredL2Norm(const vec& x);
static inline double spSquaredL2Norm(const sp_vec& x);

/**
 * Dual of the squared loss function 1/(2) * ||a - y||_2^2
 * @param a dual vector
 * @param y response vector
 * @param aNeg should `-a` be used instead of `a`
 */
static inline double lossDual(const vec& a, const vec& y, const bool aNeg);


/****************************************************************************
 *
 * Elastic Net using the Coordinatewise Gradient Descent method
 *
 *
 ***************************************************************************/

ENDal::ENDal(const bool intercept) :
			ElasticNet(intercept),
            verbosity(DEFAULT_OPT_VERBOSITY),
            maxIt(DEFAULT_OPT_MAXIT),
            eps(DEFAULT_OPT_EPS),
            warmStart(DEFAULT_OPT_WARM_START),
            etaStart(DEFAULT_OPT_ETA_START),
            etaStartNumerator(DEFAULT_OPT_ETA_START_NUMERATOR),
            etaMultiplier(DEFAULT_OPT_ETA_MULTIPLIER),
            precondType(DEFAULT_OPT_PRECONDITIONER),
            bufferSizeNobs(0), bufferSizeNvar(0),
            y(NULL), Xtr(NULL),
            useWeights(false)
{
    this->sqrtWeights = ones(1);
}

ENDal::ENDal(const bool intercept, const Options& options) :
			ElasticNet(intercept),
            verbosity(DEFAULT_OPT_VERBOSITY),
            maxIt(DEFAULT_OPT_MAXIT),
            eps(DEFAULT_OPT_EPS),
            warmStart(DEFAULT_OPT_WARM_START),
            etaStart(DEFAULT_OPT_ETA_START),
            etaStartNumerator(DEFAULT_OPT_ETA_START_NUMERATOR),
            etaMultiplier(DEFAULT_OPT_ETA_MULTIPLIER),
            precondType(DEFAULT_OPT_PRECONDITIONER),
            bufferSizeNobs(0), bufferSizeNvar(0),
            y(NULL), Xtr(NULL),
            useWeights(false)
{
    this->sqrtWeights = ones(1);
    this->setOptions(options);
}

ENDal::~ENDal()
{
    if (this->y) {
        delete this->y;
    }
    if (this->Xtr) {
        delete this->Xtr;
    }
}

void ENDal::setOptions(const Options& options)
{
    this->verbosity = options.get("verbosity", this->verbosity);
    this->maxIt = options.get("maxit", this->maxIt);
    this->eps = options.get("eps", this->eps);
    this->warmStart = options.get("warmStart", this->warmStart);
    this->etaStart = options.get("etaStart", this->etaStart);
    this->etaStartNumerator = options.get("etaStartNumerator", this->etaStartNumerator);
    this->etaMultiplier = options.get("etaMultiplier", this->etaMultiplier);
    this->precondType = (Preconditioner) options.get("preconditioner", (int) this->precondType);
}

void ENDal::setLambdas(const double lambda1, const double lambda2)
{
    this->lambda = lambda2 + lambda1;
    this->alpha = lambda1 / (lambda2 + lambda1);
    if (lambda2 == 0 && lambda1 == 0) {
        this->alpha = 0;
    }
}

void ENDal::setAlphaLambda(const double alpha, const double lambda)
{
    this->alpha = alpha;
    this->lambda = lambda;
}

inline double ENDal::fullObjectiveFun(const double intercept, const arma::sp_vec& beta)
{
    double objf = this->nLambda * (
        0.5 * (1 - this->alpha) * spSquaredL2Norm(beta) +
        this->alpha * norm(beta, 1)
    );

    if (this->intercept) {
        if (this->useWeights) {
            objf += 0.5 * squaredL2Norm(this->Xtr->t() * beta + intercept * this->sqrtWeights - (*this->y));
        } else {
            objf += 0.5 * squaredL2Norm(this->Xtr->t() * beta + intercept - (*this->y));
        }
    } else {
        objf += 0.5 * squaredL2Norm(this->Xtr->t() * beta - (*this->y));
    }

    return objf;
}

void ENDal::setData(const Data& data)
{
    /*
     * @TODO Don't copy the data here. Only if necessary!
     * sth. along the lines
     * this->Xtr = new mat(const_cast<double *>(data.getXtrConst()), data.numVar(), data.numObs());
     * this->y = new mat(const_cast<double *>(data.getYConst()), data.numObs());
     */
    /* Remove previous data */
    if (this->y) {
        delete this->y;
        this->y = NULL;
    }
    if (this->Xtr) {
        delete this->Xtr;
        this->Xtr = NULL;
    }

    this->bufferSizeNvar = data.numVar();

    /* Initialize new data */
    if (data.numObs() > 0 || data.numVar() > 0) {
        this->y = new vec(const_cast<double *>(data.getYConst()), data.numObs(), false, false);
        this->Xtr = new mat(const_cast<double *>(data.getXtrConst()), data.numVar(), data.numObs(),
                            false, false);
        if (data.numVar() > 0) {
            /* This actually will bind the data to a new memory anyways */
            this->Xtr->shed_row(0); // The intercept is handled differently!
        }
    } else {
        this->y = new vec();
        this->Xtr = new mat();
    }

    if (data.numObs() != this->bufferSizeNobs) {
        this->bufferSizeNobs = data.numObs();
        this->precond.reset();
    }

    this->hessBuffKeep.reset();
    this->hessBuff.zeros(this->bufferSizeNobs, this->bufferSizeNobs);
}

void ENDal::computeCoefsWeighted(double *RESTRICT coefs, double *RESTRICT resids,
                                 const double *RESTRICT weights)
{
    /* First check the data if something has to be done at all */
    if (this->bufferSizeNvar == 0) {
        if (this->bufferSizeNobs > 0) {
            memcpy(resids, this->y->memptr(), this->y->n_elem * sizeof(double));
        }
        return;
    }

    const vec weightsVec(const_cast<double *>(weights), this->bufferSizeNobs, false, true);
    vec residuals(resids, this->bufferSizeNobs, false, true);
    double* intercept = coefs;
    vec coefsRetDense(coefs + 1, this->bufferSizeNvar - 1, false, true);
    sp_vec coefsSparse(coefsRetDense.n_elem);

    if (this->warmStart) {
        coefsSparse = sp_vec(coefsRetDense);
    }

    this->computeCoefsWeighted(*intercept, coefsSparse, residuals, weightsVec);

    /* update returned beta */
    coefsRetDense = vec(coefsSparse);
}

void ENDal::computeCoefsWeighted(double& intercept, sp_vec& beta, vec& residuals,
                                 const vec& weights)
{
    /* First check the data if something has to be done at all */
    if (this->bufferSizeNvar == 0) {
        if (this->bufferSizeNobs > 0) {
            residuals = *(this->y);
        }
        return;
    }

    if (this->bufferSizeNobs > 0) {
        this->sqrtWeights = sqrt(weights);
    }

    if (this->bufferSizeNvar == 1 || this->bufferSizeNobs == 0) {
        intercept = 0;
        if (this->bufferSizeNobs > 0) {
            intercept = mean(this->sqrtWeights % *this->y);
        }
        residuals = (*this->y) - intercept;
        return;
    }

    /* We have at least an intercept plus one real predictor... */
    vec *const origY = this->y;
    mat *const origXtr = this->Xtr;

    if (!this->warmStart) {
        beta.zeros();
    }

    this->status = 0;
    this->statusMessage = "";

    this->sqrtWeightsOuter = this->sqrtWeights * this->sqrtWeights.t();
    this->useWeights = true;

    vec weightedY = (*this->y) % this->sqrtWeights;
    mat weightedXtr = (*this->Xtr).each_row() % this->sqrtWeights.t();

    this->y = &weightedY;
    this->Xtr = &weightedXtr;

    this->hessBuffKeep.reset();
    this->hessBuff.zeros();

    /*
     * Theoretically, we would also need to reset the preconditioner.
     * However, if the data does not change than the old one might still
     * be useful
     */

    this->dal(intercept, beta);

    if (this->intercept) {
        residuals = (*origY) - (intercept) - origXtr->t() * beta;
    } else {
        residuals = (*origY) - origXtr->t() * beta;
    }

    /* restore original data */
    this->y = origY;
    this->Xtr = origXtr;
    this->useWeights = false;
}

void ENDal::computeCoefs(double *RESTRICT coefs, double *RESTRICT resids)
{
    /* First check the data if something has to be done at all */
    if (this->bufferSizeNvar == 0) {
        if (this->bufferSizeNobs > 0) {
            memcpy(resids, this->y->memptr(), this->y->n_elem * sizeof(double));
        }
        return;
    }

    vec residuals(resids, this->bufferSizeNobs, false, true);
    double* intercept = coefs;
    vec coefsRetDense(coefs + 1, this->bufferSizeNvar - 1, false, true);
    sp_vec coefsSparse(coefsRetDense.n_elem);

    if (this->warmStart) {
        coefsSparse = sp_vec(coefsRetDense);
    }

    this->computeCoefs(*intercept, coefsSparse, residuals);

    /* update returned beta */
    coefsRetDense = vec(coefsSparse);
}

void ENDal::computeCoefs(double& intercept, arma::sp_vec& beta, arma::vec& residuals)
{
    /* First check the data if something has to be done at all */
    if (this->bufferSizeNvar == 0) {
        if (this->bufferSizeNobs > 0) {
            residuals = *(this->y);
        }
        return;
    } else if (beta.n_elem == 0 || this->bufferSizeNobs == 0) {
        if (this->bufferSizeNobs > 0) {
            intercept = mean(*this->y);
        }
        residuals = (*this->y) - intercept;
        return;
    }

    /* We have at least an intercept plus one real predictor... */
    if (!this->warmStart) {
        beta.zeros();
    }

    this->status = 0;
    this->statusMessage = "";

    this->useWeights = false;

    this->dal(intercept, beta);

    if (this->intercept) {
        residuals = (*this->y) - (intercept) - this->Xtr->t() * beta;
    } else {
        residuals = (*this->y) - this->Xtr->t() * beta;
    }
    this->useWeights = false;

    /* Update residuals */
    if (this->intercept) {
        residuals = (*this->y) - (intercept) - this->Xtr->t() * beta;
    } else {
        residuals = (*this->y) - this->Xtr->t() * beta;
    }
}

inline void ENDal::dal(double& intercept, arma::sp_vec& beta)
{
    const int nobs = this->y->n_elem;
    this->nLambda = nobs * this->lambda;
    const double la = (this->nLambda * this->alpha);
    const double updateDenomMult = 1 / (this->nLambda * (1 - this->alpha));

    vec dualVec(nobs);
    vec tmpInnerProd(nobs);

    int iter;
    double dualFunVal, dualFunValPrev;
    double primalFunVal;
    double relativeDualityGap;
    double aL1, aL1Prev;

    vec a;

    /* Variables for inner minimization problem */
    vec stepDir(nobs, fill::zeros);
    vec XtrStepDir;
    vec phiGradient;
    vec Xtra;

    int innerIter = 0, lineSearchIter = 0;
    int pcgSteps = 0;
    double stepSize, stepSizePrev;
    double phiVal, phiValStep;
    double decr;
    double threshold;
    double normGradient;
    double normDiffInt;

    /*
     * Initialize variables
     */
    if (this->etaStart > 0) {
        this->eta[0] = this->etaStart;
    } else {
        this->eta[0] = this->etaStartNumerator / this->nLambda;
        if (this->eta[0] > 1e6) {
            this->eta[0] = 1e6;
        }
    }

    this->eta[1] = this->eta[0];

    if (!this->intercept) {
        intercept = 0;
    }

    if (this->warmStart) {
        a = (*this->y) - this->Xtr->t() * beta;

        switch (this->intercept + 2 * this->useWeights)
        {
        case 1:
            /* We have an intercept but no weights */
            intercept = mean(a);
            a -= intercept;
            break;
        case 3:
            /* We have an intercept and weights */
            intercept = mean(this->sqrtWeights % a);
            a -= this->sqrtWeights * intercept;
            break;
        default:
            /* we have no intercept */
            intercept = 0;
            break;
        }
    } else {
        beta.zeros();

        switch (this->intercept + 2 * this->useWeights)
        {
        case 1:
            /* We have an intercept but no weights */
            intercept = mean(*this->y);
            a = (*this->y) - intercept;
            break;
        case 3:
            /* We have an intercept and weights */
            intercept = mean(this->sqrtWeights % *this->y);
            a = (*this->y) - this->sqrtWeights * intercept;
            break;
        default:
            /* we have no intercept */
            a = *this->y;
            break;
        }
    }

    /*
     * Outer minimization
     */
    iter = 0;
    dualFunValPrev = dualFunVal = 0;
    primalFunVal = 0;
    aL1 = aL1Prev = 0;

    Xtra = (*this->Xtr) * a;

    while (1) {
        /*
         * Check if relative duality gap (rdg) is below the threshold
         */
        dualFunValPrev = dualFunVal;
        if (this->intercept) {
            if (this->useWeights) {
                dualVec = a - this->sqrtWeights * mean(this->sqrtWeights % a);
            } else {
                dualVec = a - mean(a);
            }
        } else {
            dualVec = a;
        }

        tmpInnerProd = (*this->Xtr) * dualVec;

        if (this->alpha < 1) {
            vecSoftThreshold(tmpInnerProd, la);
            dualFunVal = lossDual(dualVec, *this->y, true) + 0.5 * squaredL2Norm(tmpInnerProd) * updateDenomMult;
        } else {
            dualVec *= fmin(this->nLambda / max(abs(tmpInnerProd)), 1);
            dualFunVal = lossDual(dualVec, *this->y, true);
        }

        if (iter > 0 && dualFunVal > dualFunValPrev) {
            dualFunVal = dualFunValPrev;
        }

        primalFunVal = this->fullObjectiveFun(intercept, beta);
        relativeDualityGap = (primalFunVal + dualFunVal) / primalFunVal;

#ifdef DEBUG
        if (this->verbosity > 0) {
            Rcpp::Rcout << "[[" << iter << "]]" <<
                " fval=" << primalFunVal <<
                "; dval=" << dualFunVal <<
                "; rdg=" << relativeDualityGap <<
                "; eta=" << this->eta[0] <<
                "; eta(int)=" << this->eta[1] <<
            std::endl;
        }
#endif

        /*
         * Check for max. iterations
         */
        if ((relativeDualityGap < this->eps) || (++iter > this->maxIt) || this->status != 0) {
            break;
        }

        /*========================================================================================*/
        /* Minimize phi --> updates `a` and `beta` & `intercept`                                  */
        /*========================================================================================*/
        const double cutoff = this->nLambda * this->eta[0] * this->alpha;
        const double multFact = 1 / (1 + this->nLambda * this->eta[0] * (1 - this->alpha));
        const sp_vec betaOrig(beta);
        const double interceptOrig = intercept;

        innerIter = 0;
        stepSize = -1;

        while (this->status == 0) {
            vecSoftThresholdInplace(beta, betaOrig, this->eta[0], Xtra, cutoff);

            phiVal = lossDual(a, *this->y, true) +
                (0.5 / this->eta[0]) * multFact * spSquaredL2Norm(beta);

            if (this->intercept) {
                intercept = interceptOrig + this->eta[1] * (this->useWeights ?
                        accu(this->sqrtWeights % a) :
                        accu(a)
                    );

                phiVal += (0.5 / this->eta[1]) * intercept * intercept;
            }

            this->evalPhiGrad(a, beta, intercept, multFact, phiGradient);

            /*
             * Check for max. iterations
             */
            if (++innerIter > this->maxIt) {
                break;
            }

            normGradient = squaredL2Norm(phiGradient);

            threshold = (1 / this->eta[0]) * spSquaredL2Norm(betaOrig - multFact * beta);
            if (this->intercept) {
                normDiffInt = (interceptOrig - intercept);
                threshold += (1 / this->eta[1]) * normDiffInt * normDiffInt;
            }

            /* It is not necessary to go far below the actual precision */
            if (threshold < this->eps) {
                threshold = 0.5 * this->eps;
            }

#ifdef DEBUG
            if (this->verbosity > 1) {
                Rcpp::Rcout << "[" << innerIter << "]" <<
                    " fval=" << phiVal <<
                    "; norm(gg)=" << sqrt(normGradient) <<
                    "; thresh=" << sqrt(threshold) <<
                    "; step=" << stepSize <<
                    "; #pcg=" << pcgSteps <<
                std::endl;
            }
#endif

            /*
             * Check for convergence to desired tolerance
             */
            if ((innerIter > 1) && (normGradient < threshold)) {
                break;
            }

            /*
             * Gradient is not small enough -- continue Newton steps
             */
            pcgSteps = this->getPhiStepDir(stepDir, phiGradient, a, beta, intercept, multFact);
            decr = dot(stepDir, phiGradient);
            XtrStepDir = (*this->Xtr) * stepDir;

            /*
             * Backtracking line search for step size to update `a`
             */
            stepSizePrev = stepSize = 1;
            a -= stepDir;
            Xtra -= XtrStepDir;
            lineSearchIter = 0;

            while (++lineSearchIter <= LINESEARCH_MAX_STEP) {
                vecSoftThresholdInplace(beta, betaOrig, this->eta[0], Xtra, cutoff);

                phiValStep = lossDual(a, *this->y, true) +
                    (0.5 / this->eta[0]) * multFact * spSquaredL2Norm(beta);

                if (this->intercept) {
                    intercept = interceptOrig + this->eta[1] * (this->useWeights ?
                            accu(this->sqrtWeights % a) :
                            accu(a)
                        );

                    phiValStep += (0.5 / this->eta[1]) * intercept * intercept;
                }

                if (phiValStep < phiVal - LINESEARCH_STEP_CONDITION * stepSize * decr) {
                    break;
                }

                stepSizePrev = stepSize;
                stepSize *= LINESEARCH_STEPSIZE_MULTIPLIER;
                a -= (stepSize - stepSizePrev) * stepDir;
                Xtra -= (stepSize - stepSizePrev) * XtrStepDir;
            }

            if (lineSearchIter > LINESEARCH_MAX_STEP) {
                this->status = 2;
                this->statusMessage = "Newton-step did not converge. Gradient is too steep.";
            }
        }

        beta *= multFact;

        /*========================================================================================*/
        /* END minimize phi                                                                       */
        /*========================================================================================*/

        /* Update `eta` */
        this->eta[0] *= this->etaMultiplier;

        /* Update `eta (intercept)` */
        if (this->intercept) {
            aL1 = (this->useWeights ? accu(this->sqrtWeights % a) : accu(a));
            if ((iter > 1) && (aL1 > this->eps) && (aL1 > 0.5 * aL1Prev)) {
                this->eta[1] *= 10 * this->etaMultiplier;
            } else {
                this->eta[1] *= this->etaMultiplier;
            }
            aL1Prev = aL1;
        }
    }

    if (iter > this->maxIt) {
        this->status = 1;
        this->statusMessage = "Algorithm did not converge in maxit steps";
    }
}


inline void ENDal::evalPhiGrad(const arma::vec &a, const arma::sp_vec& beta, const double intercept, const double multFact, arma::vec &grad)
{
    grad = a - (*this->y) + multFact * this->Xtr->t() * beta;

    switch (this->intercept + 2 * this->useWeights)
    {
    case 1:
        /* We have an intercept but no weights */
        grad += intercept;
        break;
    case 3:
        /* We have an intercept and weights */
        grad += intercept * this->sqrtWeights;
        break;
    default:
        /* we have no intercept */
        break;
    }
}


inline int ENDal::getPhiStepDir(arma::vec &stepDir, const arma::vec &grad, const arma::vec &a, const arma::sp_vec& beta, const double intercept, const double multFact)
{
    /* sp_vec does not have an overload for the find method, so it has to be done by hand */
    // uvec active = find(beta);
    uvec active(beta.n_nonzero);
    mat hess;

    sp_vec::const_iterator indexIt = beta.begin();
    for (uword i = 0; indexIt != beta.end(); ++indexIt, ++i) {
        active[i] = indexIt.row();
    }

    if (beta.n_nonzero == 0) {
        /* all beta are 0, so the solution can be easily computed */
        switch (this->intercept + 2 * this->useWeights)
        {
        case 1:
            /* We have an intercept but no weights */
            stepDir = grad / (this->eta[1] + 1);
            return 0;
            break;
        case 3:
            /* We have an intercept and weights */
            hess = this->eta[1] * this->sqrtWeightsOuter;
            hess.diag() += 1;
            break;
        default:
            /* we have no intercept */
            stepDir = grad;
            return 0;
            break;
        }
    } else {
        hess = this->eta[0] * multFact * this->getHessBuff(active);
        hess.diag() += 1;

        switch (this->intercept + 2 * this->useWeights)
        {
        case 1:
            /* We have an intercept but no weights */
            hess += this->eta[1];
            break;
        case 3:
            /* We have an intercept and weights */
            hess += this->eta[1] * this->sqrtWeightsOuter;
            break;
        default:
            /* we have no intercept */
            break;
        }

    }

#ifdef DEBUG
    if (this->verbosity > 2) {
        active.t().raw_print(Rcpp::Rcout, "active coefficients:");
        Rcpp::Rcout << "frobenius norm(hess)=" << norm(hess, "fro") <<
            "; L2 norm(beta)=" << norm(beta) <<
        std::endl;
    }
#endif

    int iters = 0;
    switch (this->precondType)
    {
    case APPROX:
        if (this->precond.n_elem == 0) {
#ifdef DEBUG
        if (this->verbosity > 1) {
            Rcpp::Rcout << "=================== Re-compute preconditioner.... ================" << std::endl;
        }
#endif
            this->precond = inv_sympd(hess);
            stepDir = this->precond * grad;
            return 0;
        }

        iters = solve_pcg(hess, stepDir, grad, this->precond, this->eps);
        if (iters > SOLVE_PCG_MAXIT_UPDATE) {
            this->precond.reset();
        }
        return iters;
    case DIAG:
        return solve_pcg_diag(hess, stepDir, grad, this->eps);
    default:
    case NONE:
        return 1 - solve(stepDir, hess, grad, solve_opts::fast + solve_opts::no_approx);
    }
}

const arma::mat& ENDal::getHessBuff(const uvec& keep)
{
    if (size(keep) != size(this->hessBuffKeep) || any(keep - this->hessBuffKeep)) {
        this->hessBuff = this->Xtr->rows(keep).t() * this->Xtr->rows(keep);
        this->hessBuffKeep = keep;
    }
    return this->hessBuff;
}


static inline double softThreshold(const double z, const double gamma)
{
    if (fabs(z) <= gamma) {
        return 0.;
    } else if (z < 0) {
        return (z + gamma);
    }
    return (z - gamma);
}

static inline void vecSoftThreshold(vec& z, const double gamma)
{
    for (vec::iterator elIterator = z.begin(); elIterator != z.end(); ++elIterator) {
        (*elIterator) = softThreshold(*elIterator, gamma);
    }
}

static inline double lossDual(const vec& a, const vec& y, const bool aNeg)
{
    if (aNeg) {
        return as_scalar(0.5 * dot(a, a) - dot(a, y));
    } else {
        return as_scalar(0.5 * dot(a, a) + dot(a, y));
    }
}

static inline double squaredL2Norm(const vec& x)
{
    double tmp = norm(x, 2);
    return tmp * tmp;
}
static inline double spSquaredL2Norm(const sp_vec& x)
{
    double tmp = norm(x, 2);
    return tmp * tmp;
}


/*****************************************************************
 * Iterative routine -- CG
 *
 * CG solves the symmetric positive definite linear
 * system Ax=b using the Preconditioned Conjugate Gradient method.
 * The (reciprocal) diagonal of A is used as the preconditioner.
 *
 * CG follows the algorithm described on p. 15 in the
 * SIAM Templates book.
 *
 * Upon successful return, output arguments have the following values:
 *
 *        x  --  approximate solution to Ax = b
 *
 * original source code from http://math.nist.gov/iml++/cg.h.txt
 * as well as from the Eigen C++ template library for linear algebra.
 * Copyright (C) 2011-2014 Gael Guennebaud <gael.guennebaud@inria.fr>
 *
 * This Source Code Form is subject to the terms of the Mozilla
 * Public License v. 2.0. If a copy of the MPL was not distributed
 * with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
 ****************************************************************/
static int solve_pcg_diag(const mat &A, vec &x, const vec &b, const double eps)
{
    int iter = 0;
    double alpha, beta;
    double normb = norm(b);
    vec r = b - A * x;

    if (normb == 0.0) {
        normb = 1;
    }

    const double thresh = eps * normb;

    if (norm(r) <= thresh) {
        return iter;
    }

    vec z(A.n_cols);
    vec tmp(A.n_cols);
    vec p = r / A.diag();
    double absOld, absNew = dot(r, p);

    do {
        tmp = A * p;
        alpha = absNew / dot(p, tmp);

        x += alpha * p;
        r -= alpha * tmp;

        if (norm(r) <= thresh) {
            return iter;
        }

        z = r / A.diag();

        absOld = absNew;
        absNew = dot(r, z);

        beta = absNew / absOld;
        p = z + beta * p;
    } while (++iter <= SOLVE_PCG_MAXIT);

    return iter;
}


/*****************************************************************
 * Iterative routine -- CG
 *
 * CG solves the symmetric positive definite linear
 * system Ax=b using the Preconditioned Conjugate Gradient method.
 * The matrix T is used as the preconditioner.
 *
 * CG follows the algorithm described on p. 15 in the
 * SIAM Templates book.
 *
 * Upon successful return, output arguments have the following values:
 *
 *        x  --  approximate solution to Ax = b
 *
 * original source code from http://math.nist.gov/iml++/cg.h.txt
 * as well as from the Eigen C++ template library for linear algebra.
 * Copyright (C) 2011-2014 Gael Guennebaud <gael.guennebaud@inria.fr>
 *
 * This Source Code Form is subject to the terms of the Mozilla
 * Public License v. 2.0. If a copy of the MPL was not distributed
 * with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
 ****************************************************************/
static int solve_pcg(const mat &A, vec &x, const vec &b, const mat &precond, const double eps)
{
    int iter = 0;
    double alpha, beta;
    double normb = norm(b);
    vec r = b - A * x;

    if (normb == 0.0) {
        normb = 1;
    }

    const double thresh = eps * normb;

    if (norm(r) <= thresh) {
        return iter;
    }

    vec z(A.n_cols);
    vec tmp(A.n_cols);
    vec p = precond * r;
    double absOld, absNew = dot(r, p);

    do {
        tmp = A * p;
        alpha = absNew / dot(p, tmp);

        x += alpha * p;
        r -= alpha * tmp;

        if (norm(r) <= thresh) {
            return iter;
        }

        z = precond * r;

        absOld = absNew;
        absNew = dot(r, z);

        beta = absNew / absOld;
        p = z + beta * p;
    } while (++iter <= SOLVE_PCG_MAXIT);

    return iter;
}


/**
 * Soft-thresholding functions for sparse vectors
 */
static inline void vecSoftThresholdInplace(sp_vec& zsoft, const sp_vec& sp_z1, const double c, const vec& z2, const double gamma)
{
    if (sp_z1.n_nonzero < VEC_SOFTTHRESH_SMALL_SPARSE) {
        /* This branch is better for small sparse vectors */
        vecSoftThresholdInplaceSmallSparse(zsoft, sp_z1, c, z2, gamma);
    } else if (sp_z1.n_nonzero < VEC_SOFTTHRESH_SPARSE * sp_z1.n_elem) {
        /* This branch is better for large sparse vectors */
        vecSoftThresholdInplaceLargeSparse(zsoft, sp_z1, c, z2, gamma);
    } else {
        /* This branch is better for not-so-sparse vectors */
        vecSoftThresholdInplaceNonSparse(zsoft, sp_z1, c, z2, gamma);
    }
}


static inline void vecSoftThresholdInplaceNonSparse(sp_vec& zsoft, const sp_vec& sp_z1, const double c, const vec& z2, const double gamma)
{
    const vec z1 = vec(sp_z1);
    vec zsoft_tmp(z1.n_elem);
    vec::const_iterator readIter1 = z1.begin();
    vec::const_iterator readIter2 = z2.begin();
    vec::iterator writeIter = zsoft_tmp.begin();

    for (; writeIter != zsoft_tmp.end(); ++writeIter, ++readIter1, ++readIter2) {
        (*writeIter) = *readIter1 + c * (*readIter2);
        if ((*writeIter) > gamma) {
            (*writeIter) -= gamma;
        } else if ((*writeIter) < -gamma) {
            (*writeIter) += gamma;
        } else {
            (*writeIter) = 0;
        }
    }
    zsoft = sp_vec(zsoft_tmp);
}

static inline void vecSoftThresholdInplaceSmallSparse(sp_vec& zsoft, const sp_vec& z1, const double c, const vec&z2, const double gamma)
{
    double tmp;
    uword i = 0, upper_bound = 0;
    zsoft.zeros(size(z1));

    sp_vec::const_iterator z1Iter = z1.begin();
    vec::const_iterator z2Iter = z2.begin();

    do {
        upper_bound = (z1Iter == z1.end() ? z1.n_elem : z1Iter.row());
        while (i < upper_bound) {
            tmp = c * (*z2Iter);
            if (tmp > gamma) {
                zsoft[i] = tmp - gamma;
            } else if (tmp < -gamma) {
                zsoft[i] = tmp + gamma;
            }
            ++z2Iter;
            ++i;
        }

        if (i < z1.n_elem) {
            tmp = (*z1Iter) + c * (*z2Iter);
            if (tmp > gamma) {
                zsoft[i] = tmp - gamma;
            } else if (tmp < -gamma) {
                zsoft[i] = tmp + gamma;
            }

            ++z1Iter;
            ++z2Iter;
        }
    } while (++i < z1.n_elem);
}

static inline void vecSoftThresholdInplaceLargeSparse(sp_vec& zsoft, const sp_vec& z1, const double c, const vec&z2, const double gamma)
{
    double tmp;
    std::vector<uword> ind;
    std::vector<double> val;
    uword i = 0;
    uword upper_bound = 0;

    ind.reserve(std::min(z1.n_elem, (uword) 1.4 * z1.n_nonzero));
    val.reserve(ind.capacity());

    sp_vec::const_iterator z1Iter = z1.begin();
    vec::const_iterator z2Iter = z2.begin();

    do {
        upper_bound = (z1Iter == z1.end() ? z1.n_elem : z1Iter.row());
        while (i < upper_bound) {
            tmp = c * (*z2Iter);
            if (tmp > gamma) {
                ind.push_back(i);
                val.push_back(tmp - gamma);
            } else if (tmp < -gamma) {
                ind.push_back(i);
                val.push_back(tmp + gamma);
            }
            ++z2Iter;
            ++i;
        }

        if (i < z1.n_elem) {
            tmp = (*z1Iter) + c * (*z2Iter);
            if (tmp > gamma) {
                ind.push_back(i);
                val.push_back(tmp - gamma);
            } else if (tmp < -gamma) {
                ind.push_back(i);
                val.push_back(tmp + gamma);
            }
        }
        ++z1Iter;
        ++z2Iter;
    } while (++i < z1.n_elem);

    umat indmat(2, (uword) ind.size(), fill::zeros);
    indmat.row(0) = Row<uword>(ind);
    zsoft = sp_mat(indmat, vec(val), z1.n_elem, 1, false, false).col(0);
}
