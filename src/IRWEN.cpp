//
//  IRWEN.cpp
//  pense
//
//  Created by David Kepplinger on 2017-05-30.
//  Copyright © 2017 David Kepplinger. All rights reserved.
//
#include "config.h"
#include <RcppArmadillo.h>

#include "IRWEN.hpp"

#include <Rmath.h>

#include "ElasticNet.hpp"
#include "olsreg.h"

using namespace arma;

static const int DEFAULT_OPT_VERBOSITY = 0;
static const int DEFAULT_OPT_MAXIT = 1000;
static const double DEFAULT_OPT_EPS = 1e-6;

static const double NUMERICAL_TOLERANCE = NUMERIC_EPS;

static const Options WARM_START_OPTION = Options::createSimple("warmStart", true);

IRWEN::IRWEN(const Data& data, const double alpha, const double lambda, const Options& opts, const Options &enOpts) :
    Xtr(const_cast<double *>(data.getXtrConst()), data.numVar(), data.numObs(), false, true),
    y(const_cast<double *>(data.getYConst()), data.numObs(), false, true),
    weights(data.numObs()),
    verbosity(opts.get("verbosity", DEFAULT_OPT_VERBOSITY)),
    maxIt(opts.get("maxit", DEFAULT_OPT_MAXIT)),
    eps(opts.get("eps", DEFAULT_OPT_EPS)),
    lambda(lambda)
{
    this->en = getElasticNetImpl(enOpts, true);
    this->en->setOptions(WARM_START_OPTION);
    this->en->setAlphaLambda(alpha, lambda);
    this->en->setData(data);
}

IRWEN::~IRWEN()
{
    delete this->en;
}

void IRWEN::compute(double& intercept, arma::sp_vec& beta, arma::vec& residuals)
{
    double oldIntercept;
    sp_vec oldBeta(beta);
    double normOldBeta;

    this->iteration = 0;
    residuals = this->y - this->Xtr.tail_rows(beta.n_elem).t() * beta - intercept;

    do {
        ++this->iteration;

        this->updateWeights(residuals);

        oldIntercept = intercept;
        oldBeta = beta;

        /*
         * Perform EN using current coefficients as warm start (only applicable for the coordinate
         * descent algorithm)
         */
        en->computeCoefsWeighted(intercept, beta, residuals, this->weights);

        if (en->getStatus() > 0) {
            std::ostringstream stringStream;
            stringStream << "Weighted elastic net had non-successful status for lambda=" <<
                this->lambda << ": " << en->getStatusMessage();
            Rcpp::warning(stringStream.str());
        }

        /*
         * Compute relative change
         */
        normOldBeta = norm(oldBeta) + sqrt(oldIntercept * oldIntercept);
        this->relChangeVal = norm(beta - oldBeta) + sqrt((intercept - oldIntercept) * (intercept - oldIntercept));

        if (normOldBeta < NUMERICAL_TOLERANCE) {
            if (this->relChangeVal < NUMERICAL_TOLERANCE) {
                /* We barely moved away from the zero-vector --> "converged" */
                this->relChangeVal = 0;
            } else {
                /* We moved away from the zero-vector --> continue */
                this->relChangeVal = 2 * this->eps;
            }
        }

#ifdef DEBUG
        if (this->verbosity > 0) {
            Rcpp::Rcout << "IRWEN [[" << this->iteration << "]] "
                "lambda=" << this->lambda <<
                "; norm(old coefs)=" << normOldBeta <<
                "; norm(weights)=" << norm(this->weights) <<
                "; rel. change=" << this->relChangeVal / normOldBeta <<
            std::endl;
        }
#endif

    } while((this->iteration < this->maxIt) && (this->relChangeVal > normOldBeta * this->eps));

    this->relChangeVal /= normOldBeta;
}
