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


/**
 * constants defining the decrease in the step size in each line search iteration
 * and the fraction of the decrease in the objective function predicted by
 * linear extrapolation that we will accept in the line search
 */
static const double LINESEARCH_STEPSIZE_MULTIPLIER = 0.9;  // 0 < x < 1
static const double LINESEARCH_STEP_CONDITION = 0.45;      // 0 < x < 0.5

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
			ElasticNet(intercept), bufferSizeNobs(0), y(), Xtr()
{
    this->setOptions(Options());
}

ENDal::ENDal(const bool intercept, const Options& options) :
			ElasticNet(intercept), bufferSizeNobs(0), y(), Xtr()
{
    this->setOptions(options);
}

ENDal::~ENDal()
{
}

void ENDal::setOptions(const Options& options)
{
    this->maxIt = options.get("maxit", DEFAULT_OPT_MAXIT);
    this->eps = options.get("eps", DEFAULT_OPT_EPS);
    this->warmStart = options.get("warmStart", DEFAULT_OPT_WARM_START);
    this->etaStart = options.get("etaStart", DEFAULT_OPT_ETA_START);
    this->etaStartNumerator = options.get("etaStartNumerator", DEFAULT_OPT_ETA_START_NUMERATOR);
    this->etaMultiplier = options.get("etaMultiplier", DEFAULT_OPT_ETA_MULTIPLIER);
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
    return 0.5 * squaredL2Norm(this->Xtr.t() * beta + intercept - this->y) + this->nLambda * (
        0.5 * (1 - this->alpha) * squaredL2Norm(beta) +
        this->alpha * norm(beta, 1)
    );

}

void ENDal::setData(const Data& data)
{
    if (data.numObs() > 0 && data.numVar() > 0) {
        this->y = vec(data.getYConst(), data.numObs());
        this->Xtr = mat(data.getXtrConst(), data.numVar(), data.numObs());
        this->Xtr.shed_row(0); // We need to take care of the intercept in a different way!
    } else {
        this->y.reset();
        this->Xtr.reset();
    }

    if (data.numObs() != this->bufferSizeNobs) {
        this->a.resize(data.numObs());
        this->a = -this->y;
        this->bufferSizeNobs = data.numObs();
    }

    this->nLambda = this->y.n_elem * this->lambda;
}

void ENDal::computeCoefs(double *RESTRICT coefs, double *RESTRICT resids)
{
    const double la = (this->nLambda * this->alpha);
    const double updateDenomMult = 1 / (this->nLambda * (1 - this->alpha));
    const int nobs = this->y.n_elem;
    const int nvar = this->Xtr.n_rows;

    double *intercept = coefs;
    vec dualVec(nobs);
    vec tmpInnerProd(nobs);
    vec beta(coefs + 1, nvar, false, true);
    vec residuals(resids, nobs, false, true);

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

    if (!this->warmStart) {
        *intercept = 0;
        beta.zeros();
        this->a = -this->y;
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
        dualVec = this->a - mean(this->a); // Now we have the duality vector
        tmpInnerProd = this->Xtr * dualVec;

        if (this->alpha < 1) {
            vecSoftThreshold(tmpInnerProd, la, 1);
            dualFunVal = lossDual(dualVec, this->y, true) + 0.5 * squaredL2Norm(tmpInnerProd) * updateDenomMult;
        } else {
            dualVec *= fmin(nLambda / max(abs(tmpInnerProd)), 1);
            dualFunVal = lossDual(dualVec, this->y, true);
        }

        if (iter > 0 && dualFunVal > dualFunValPrev) {
            dualFunVal = dualFunValPrev;
        }

        primalFunVal = this->fullObjectiveFun(*intercept, beta);
        relativeDualityGap = (primalFunVal + dualFunVal) / primalFunVal;

        Rcpp::Rcout << "[[" << iter << "]]" <<
            " fval=" << primalFunVal <<
            "; dval=" << dualFunVal <<
            "; rdg=" << relativeDualityGap <<
            "; eta=" << this->eta[0] <<
            "; eta(int)=" << this->eta[1] <<
            std::endl;

        if (relativeDualityGap < this->eps) {
            break;
        }

        /*
         * Check for max. iterations
         */
        if (++iter > this->maxIt) {
            break;
        }

        /* Minimize phi --> updates `a` and `beta` + `intercept` */
        this->minimizePhi(beta, *intercept);

        /* Update `eta` */
        this->eta[0] *= this->etaMultiplier;

        /* Update `eta (intercept)` */
        aL1 = accu(this->a);
        if ((iter > 1) && (aL1 > this->eps) && (aL1 > 0.5 * aL1Prev)) {
            this->eta[1] *= 10 * this->etaMultiplier;
        } else {
            this->eta[1] *= this->etaMultiplier;
        }
        aL1Prev = aL1;
    }

    /* Update residuals */
    residuals = this->y - (*intercept) - this->Xtr.t() * beta;

    if (iter > this->maxIt) {
        this->status = 1;
        this->statusMessage = "Algorithm did not converge";
//        throw Error("Algorithm did not converge");
    }
}

inline bool ENDal::minimizePhi(vec& beta, double& intercept)
{
    const vec betaOrig(beta);
    const double interceptOrig = intercept;
    int iter = 0, lineSearchIter = 0;
    double stepSize = 1;
    double phiVal;
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

        stepDir = solve(phiHessian, phiGradient);
        decr = dot(stepDir, phiGradient);

        normDiffInt = (interceptOrig - intercept);
        normGradient = squaredL2Norm(phiGradient);

        Rcpp::Rcout << "[" << iter << "]" <<
            " fval=" << phiVal <<
            "; norm(gg)=" << normGradient <<
            "; decr=" << decr <<
            "; step=" << stepSize <<
            std::endl;

        /*
         * Check for convergence to desired tolerance
         */
        threshold = (1 / this->eta[0]) * squaredL2Norm(betaOrig - beta) +
            (1 / this->eta[1]) * normDiffInt * normDiffInt;
        if ((iter > 1) && (threshold < normGradient)) {
            break;
        }

        /*
         * Backtracking line search for step size to update `a`
         */
        lineSearchIter = 0;
        stepSize = 1;
        while (1) {
            candidateSolution = this->a - stepSize * stepDir;
            beta = betaOrig;
            intercept = interceptOrig;
            phiVal = this->evalPhi(candidateSolution, beta, intercept, phiGradient,
                                   phiHessian, false);

            if (phiVal < phiVal - LINESEARCH_STEP_CONDITION * stepSize * decr) {
                break;
            }

            if (++lineSearchIter > this->maxIt) {
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
    beta += this->eta[0] * this->Xtr * a;
    intercept += this->eta[1] * accu(a);

    vecSoftThreshold(beta, cutoff, 1);

    const double phiMoreauEnv = 0.5 * multFact * squaredL2Norm(beta);
    const double interceptMorauEnv = 0.5 * intercept * intercept;

    const double phiVal = lossDual(a, this->y, true) + (1 / this->eta[0]) * phiMoreauEnv +
        (1 / this->eta[1]) * interceptMorauEnv;

    if (evalGrad) {
        mat tmp = this->Xtr.rows(find(beta));
        grad = a - this->y + multFact * (this->Xtr.t() * beta) + intercept;
        hess = this->eta[0] * multFact * (tmp.t() * tmp) + this->eta[1];
        hess.diag() += 1;
    }

    beta *= multFact;

    return phiVal;
}

void ENDal::computeCoefsWeighted(double *RESTRICT coefs, double *RESTRICT resids,
                                 const double *RESTRICT weights)
{
    throw Rcpp::exception("Weighted EN is not yet implemented with DAL");
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
