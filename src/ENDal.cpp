//
//  ENDal.cpp
//  pense
//
//  Created by David Kepplinger on 2017-05-21.
//  Copyright © 2017 David Kepplinger. All rights reserved.
//
#include "config.h"

#include <cfloat>
#include <Rmath.h>

#include <RcppArmadillo.h>

#include "ENDal.hpp"

/**
 * Default argument values
 */
static const int    DEFAULT_OPT_MAXIT = 100;
static const double DEFAULT_OPT_EPS = 1e-5;
static const bool   DEFAULT_OPT_WARM_START = true;
static const double DEFAULT_OPT_ETA_START = -1;
static const double DEFAULT_OPT_ETA_START_NUMERATOR = 0.01;
static const double DEFAULT_OPT_ETA_MULTIPLIER = 2;
static const bool   DEFAULT_OPT_USE_BUFFER = true;


/**
 * constants defining the decrease in the step size in each line search iteration
 * and the fraction of the decrease in the objective function predicted by
 * linear extrapolation that we will accept in the line search
 */
static const double LINESEARCH_STEPSIZE_MULTIPLIER = 0.8;  // 0 < x < 1
static const double LINESEARCH_STEP_CONDITION = 0.3;      // 0 < x < 0.5
static const int LINESEARCH_MAX_STEP = 20;

using namespace arma;

/**
 * Static inline functions only used in this compile unit
 */

/**
 * (Vectorized) soft threshold functions, computing
 * mult * sign(z) * max(0, |z| - gamma)
 */
static inline double softThreshold(const double z, const double gamma, const double mult);
static inline void vecSoftThreshold(vec& z, const double gamma, const double mult);

static inline double squaredL2Norm(const vec& x);

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
            maxIt(DEFAULT_OPT_MAXIT),
            eps(DEFAULT_OPT_EPS),
            warmStart(DEFAULT_OPT_WARM_START),
            etaStart(DEFAULT_OPT_ETA_START),
            etaStartNumerator(DEFAULT_OPT_ETA_START_NUMERATOR),
            etaMultiplier(DEFAULT_OPT_ETA_MULTIPLIER),
            useHessBuffer(DEFAULT_OPT_USE_BUFFER),
            bufferSizeNobs(0), bufferSizeNvar(0),
            y(NULL), Xtr(NULL),
            useWeights(false)
{
    this->sqrtWeights = ones(1);
}

ENDal::ENDal(const bool intercept, const Options& options) :
			ElasticNet(intercept),
            maxIt(DEFAULT_OPT_MAXIT),
            eps(DEFAULT_OPT_EPS),
            warmStart(DEFAULT_OPT_WARM_START),
            etaStart(DEFAULT_OPT_ETA_START),
            etaStartNumerator(DEFAULT_OPT_ETA_START_NUMERATOR),
            etaMultiplier(DEFAULT_OPT_ETA_MULTIPLIER),
            useHessBuffer(DEFAULT_OPT_USE_BUFFER),
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
    this->maxIt = options.get("maxit", this->maxIt);
    this->eps = options.get("eps", this->eps);
    this->warmStart = options.get("warmStart", this->warmStart);
    this->etaStart = options.get("etaStart", this->etaStart);
    this->etaStartNumerator = options.get("etaStartNumerator", this->etaStartNumerator);
    this->etaMultiplier = options.get("etaMultiplier", this->etaMultiplier);
    this->useHessBuffer = options.get("useBuffer", this->useHessBuffer);
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

inline double ENDal::fullObjectiveFun(const double intercept, const arma::vec& beta)
{
    double objf = this->nLambda * (
        0.5 * (1 - this->alpha) * squaredL2Norm(beta) +
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
        this->y = new vec(data.getYConst(), data.numObs());
        this->Xtr = new mat(data.getXtrConst(), data.numVar(), data.numObs());
        if (data.numVar() > 0) {
            this->Xtr->shed_row(0); // The intercept is handled differently!
        }
    } else {
        this->y = new vec();
        this->Xtr = new mat();
    }

    if (data.numObs() != this->bufferSizeNobs) {
        this->a.resize(data.numObs());
        this->a = -(*this->y);
        this->bufferSizeNobs = data.numObs();
    }


    if (this->useHessBuffer) {
        this->hessBuffKeep.reset();
        this->hessBuff.zeros(this->bufferSizeNobs, this->bufferSizeNobs);
    }
}

void ENDal::computeCoefsWeighted(double *RESTRICT coefs, double *RESTRICT resids,
                                 const double *RESTRICT weights)
{
    vec *const origY = this->y;
    mat *const origXtr = this->Xtr;
    double *intercept = coefs;
    vec residuals(resids, this->bufferSizeNobs, false, true);
    vec beta(coefs + 1, this->bufferSizeNvar - 1, false, true);

    /* First check the data if something has to be done at all */
    if (this->bufferSizeNvar == 0) {
        if (this->bufferSizeNobs > 0) {
            residuals = *(this->y);
        }
        return;
    }

    this->sqrtWeights = sqrt(vec(weights, this->bufferSizeNobs));
    this->useWeights = true;

    if (this->bufferSizeNvar == 1 || this->bufferSizeNobs == 0) {
        *intercept = ((this->bufferSizeNobs > 0) ? mean(this->sqrtWeights % *this->y) : 0);
        beta.zeros();
        residuals = (*this->y) - (*intercept);
        return;
    }

    vec weightedY = (*this->y) % this->sqrtWeights;
    mat weightedXtr = (*this->Xtr).each_row() % this->sqrtWeights.t();

    this->y = &weightedY;
    this->Xtr = &weightedXtr;

    this->dal(*intercept, beta);

    if (this->intercept) {
        residuals = (*origY) - (*intercept) - origXtr->t() * beta;
    } else {
        residuals = (*origY) - origXtr->t() * beta;
    }

    this->y = origY;
    this->Xtr = origXtr;
    this->useWeights = false;
}

void ENDal::computeCoefs(double *RESTRICT coefs, double *RESTRICT resids)
{
    double *intercept = coefs;
    vec residuals(resids, this->bufferSizeNobs, false, true);
    vec beta(coefs + 1, this->bufferSizeNvar - 1, false, true);

    /* First check the data if something has to be done at all */
    if (this->bufferSizeNvar == 0) {
        if (this->bufferSizeNobs > 0) {
            residuals = *(this->y);
        }
        return;
    } else if (this->bufferSizeNvar == 1 || this->bufferSizeNobs == 0) {
        *intercept = ((this->bufferSizeNobs > 0) ? mean(*this->y) : 0);
        beta.zeros();
        residuals = (*this->y) - (*intercept);
        return;
    }

    this->useWeights = false;

    this->dal(*intercept, beta);

    if (this->intercept) {
        residuals = (*this->y) - (*intercept) - this->Xtr->t() * beta;
    } else {
        residuals = (*this->y) - this->Xtr->t() * beta;
    }
    this->useWeights = false;

    /* Update residuals */
    if (this->intercept) {
        residuals = (*this->y) - (*intercept) - this->Xtr->t() * beta;
    } else {
        residuals = (*this->y) - this->Xtr->t() * beta;
    }
}

inline void ENDal::dal(double& intercept, arma::vec& beta)
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

    if (this->etaStart > 0) {
        this->eta[0] = this->etaStart;
    } else {
        this->eta[0] = this->etaStartNumerator / this->nLambda;
        if (this->eta[0] > 1e6) {
            this->eta[0] = 1e6;
        }
    }

    this->eta[1] = this->eta[0];

    if (this->intercept) {

    } else {
        intercept = 0;
    }

    if (this->warmStart) {
        this->a = (*this->y) - this->Xtr->t() * beta;
        if (this->intercept) {
            if (this->useWeights) {
                intercept = mean(this->sqrtWeights % this->a);
                this->a -= this->sqrtWeights * intercept;
            } else {
                intercept = mean(this->a);
                this->a -= intercept;
            }
        } else {
            intercept = 0;
        }
    } else {
        if (this->intercept) {
            if (this->useWeights) {
                intercept = mean(this->sqrtWeights % *this->y);
                this->a = (*this->y) - this->sqrtWeights * intercept;
            } else {
                intercept = mean(*this->y);
                this->a = (*this->y) - intercept;
            }
        } else {
            this->a = *this->y;
        }
        beta.zeros();

    }

    iter = 0;
    dualFunValPrev = dualFunVal = 0;
    primalFunVal = 0;
    aL1 = aL1Prev = 0;

    while (1) {
        /*
         * Check if relative duality gap (rdg) is below the threshold
         */
        dualFunValPrev = dualFunVal;
        if (this->intercept) {
            if (this->useWeights) {
                dualVec = this->a - this->sqrtWeights * mean(this->sqrtWeights % this->a);
            } else {
                dualVec = this->a - mean(this->a);
            }
        } else {
            dualVec = this->a;
        }

        tmpInnerProd = (*this->Xtr) * dualVec;

        if (this->alpha < 1) {
            vecSoftThreshold(tmpInnerProd, la, 1);
            dualFunVal = lossDual(dualVec, *this->y, true) + 0.5 * squaredL2Norm(tmpInnerProd) * updateDenomMult;
        } else {
            dualVec *= fmin(nLambda / max(abs(tmpInnerProd)), 1);
            dualFunVal = lossDual(dualVec, *this->y, true);
        }

        if (iter > 0 && dualFunVal > dualFunValPrev) {
            dualFunVal = dualFunValPrev;
        }

        primalFunVal = this->fullObjectiveFun(intercept, beta);
        relativeDualityGap = (primalFunVal + dualFunVal) / primalFunVal;

#ifdef DEBUG
        Rcpp::Rcout << "[[" << iter << "]]" <<
            " fval=" << primalFunVal <<
            "; dval=" << dualFunVal <<
            "; rdg=" << relativeDualityGap <<
            "; eta=" << this->eta[0] <<
            "; eta(int)=" << this->eta[1] <<
            std::endl;
#endif

        if (relativeDualityGap < this->eps) {
            break;
        }

        /*
         * Check for max. iterations
         */
        if (++iter > this->maxIt) {
            break;
        }

        /* Minimize phi --> updates `a` and `beta` & `intercept` */
        this->minimizePhi(beta, intercept);

        /* Update `eta` */
        this->eta[0] *= this->etaMultiplier;

        /* Update `eta (intercept)` */
        if (this->intercept) {
            aL1 = (this->useWeights ? accu(this->sqrtWeights % this->a) : accu(this->a));
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
        this->statusMessage = "algorithm did not converge";
    }
}


inline bool ENDal::minimizePhi(vec& beta, double& intercept)
{
    const vec betaOrig(beta);
    const double interceptOrig = intercept;
    int iter = 0, lineSearchIter = 0;
    double stepSize = 1;
    double phiVal, phiValStep;
    double decr;
    double threshold;
    double normGradient;
    double normDiffInt;
    vec candidateSolution(this->a);
    vec stepDir;
    vec phiGradient;
    mat phiHessian;

    while (1) {
        beta = betaOrig;
        intercept = interceptOrig;

        phiVal = this->evalPhi(this->a, beta, intercept, phiGradient, phiHessian, true);

        /*
         * Check for max. iterations
         */
        if (++iter > this->maxIt) {
            break;
        }

        stepDir = solve(phiHessian, phiGradient, solve_opts::fast + solve_opts::no_approx);
        decr = dot(stepDir, phiGradient);

        normGradient = squaredL2Norm(phiGradient);

        threshold = (1 / this->eta[0]) * squaredL2Norm(betaOrig - beta);
        if (this->intercept) {
            normDiffInt = (interceptOrig - intercept);
            threshold += (1 / this->eta[1]) * normDiffInt * normDiffInt;
        }

        /* It is not necessary to go far below the actual precision */
        if (threshold < this->eps) {
            threshold = 0.5 * this->eps;
        }

#ifdef DEBUG
        Rcpp::Rcout << "[" << iter << "]" <<
            " fval=" << phiVal <<
            "; norm(gg)=" << sqrt(normGradient) <<
            "; thresh=" << sqrt(threshold) <<
            "; decr=" << decr <<
            "; step=" << stepSize <<
            std::endl;
#endif

        /*
         * Check for convergence to desired tolerance
         */
        if ((iter > 1) && (normGradient < threshold)) {
            break;
        }

        /*
         * Backtracking line search for step size to update `a`
         */
        lineSearchIter = 0;
        stepSize = 1;
        while (++lineSearchIter <= LINESEARCH_MAX_STEP) {
            candidateSolution = this->a - stepSize * stepDir;
            beta = betaOrig;
            intercept = interceptOrig;
            phiValStep = this->evalPhi(candidateSolution, beta, intercept, phiGradient,
                                       phiHessian, false);

            if (phiValStep < phiVal - LINESEARCH_STEP_CONDITION * stepSize * decr) {
                break;
            }

            stepSize *= LINESEARCH_STEPSIZE_MULTIPLIER;
        }

        this->a = candidateSolution;
    }

    return (iter <= this->maxIt);
}

inline double ENDal::evalPhi(const arma::vec& a, arma::vec& beta, double& intercept,
                             arma::vec &grad, arma::mat& hess, bool evalGrad)
{
    const double cutoff = this->nLambda * this->eta[0] * this->alpha;
    const double multFact = 1 / (1 + this->nLambda * this->eta[0] * (1 - this->alpha));
    beta += this->eta[0] * (*this->Xtr) * a;

    if (this->intercept) {
        intercept += this->eta[1] * (this->useWeights ? accu(this->sqrtWeights % a) : accu(a));
    }

    vecSoftThreshold(beta, cutoff, 1);

    const double phiMoreauEnv = 0.5 * multFact * squaredL2Norm(beta);
    const double interceptMorauEnv = (this->intercept ? 0.5 * intercept * intercept : 0 );

    const double phiVal = lossDual(a, *this->y, true) + (1 / this->eta[0]) * phiMoreauEnv +
        (1 / this->eta[1]) * interceptMorauEnv;

    if (evalGrad) {
        const uvec keep = find(beta);
        grad = a - (*this->y) + multFact * this->Xtr->t() * beta;
        if (this->useHessBuffer) {
            hess = this->eta[0] * multFact * this->getHessBuff(keep);
        } else {
            hess = this->eta[0] * multFact * this->Xtr->rows(keep).t() * this->Xtr->rows(keep);
        }

        hess.diag() += 1;

        switch (this->intercept + 2 * this->useWeights)
        {
        case 1:
            /* We have an intercept but no weights */
            grad += intercept;
            hess += this->eta[1];
            break;
        case 3:
            /* We have an intercept and weights */
            grad += intercept * this->sqrtWeights;
            hess += this->eta[1] * this->sqrtWeights * this->sqrtWeights.t();
            break;
        default:
            /* we have no intercept */
            break;
        }
    }

    beta *= multFact;

    return phiVal;
}

const arma::mat& ENDal::getHessBuff(const uvec& keep)
{
    if (size(keep) != size(this->hessBuffKeep) || any(keep - this->hessBuffKeep)) {
        this->hessBuff = this->Xtr->rows(keep).t() * this->Xtr->rows(keep);
        this->hessBuffKeep = keep;
    }
    return this->hessBuff;
}


static inline double softThreshold(const double z, const double gamma, const double mult)
{
    if (fabs(z) <= gamma) {
        return 0.;
    } else if (z < 0) {
        return mult * (z + gamma);
    }
    return mult * (z - gamma);
}

static inline void vecSoftThreshold(vec& z, const double gamma, const double mult)
{
    for (vec::iterator elIterator = z.begin(); elIterator != z.end(); ++elIterator) {
        (*elIterator) = softThreshold(*elIterator, gamma, mult);
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
