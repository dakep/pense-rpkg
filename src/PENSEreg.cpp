//
//  PENSEreg.cpp
//  pense
//
//  Created by David Kepplinger on 2016-04-22.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//
#include "config.h"
#include <RcppArmadillo.h>

#include "PENSEreg.hpp"

#include <Rmath.h>

#include "Data.hpp"
#include "ElasticNet.hpp"
#include "olsreg.h"
#include "mscale.h"

using namespace arma;

static const RhoFunction rhoBisquare2 = getRhoFunctionByName(BISQUARE);
static const double DEFAULT_OPT_BDP = 0.5;
static const double DEFAULT_OPT_CC = 1.5476445356;
static const double DEFAULT_OPT_MSCALE_EPS = 1e-8;
static const int DEFAULT_OPT_MSCALE_MAX_IT = 200;

static inline double wgtBisquare2(double x, double c);

PENSEReg::PENSEReg(const Data& data, const double alpha, const double lambda, const Options& opts, const Options &enOpts) :
        IRWEN(data, alpha, lambda, opts, enOpts),
        bdp(opts.get("bdp", DEFAULT_OPT_BDP)),
        cc(opts.get("cc", DEFAULT_OPT_CC)),
        mscaleEps(opts.get("mscaleEps", DEFAULT_OPT_MSCALE_EPS)),
        mscaleMaxIt(opts.get("mscaleMaxit", DEFAULT_OPT_MSCALE_MAX_IT)),
        scale(1)
{
}

PENSEReg::~PENSEReg()
{
}

void PENSEReg::updateObjective(const vec& residuals, const double betaENPenalty)
{
    this->objectiveVal = this->scale * this->scale + betaENPenalty;
}

void PENSEReg::updateWeights(const vec& residuals)
{
    double tmp;
    this->scale = mscale(residuals.memptr(), residuals.n_elem, this->bdp, this->mscaleEps,
                         this->mscaleMaxIt, rhoBisquare2, this->cc);

    tmp = 0;
    for (uword i = 0; i < residuals.n_elem; ++i) {
        this->weights[i] = wgtBisquare2(residuals[i], this->scale * this->cc);
        tmp += this->weights[i];
    }

    /*
     * Normalize weights to sum to n --> just as in the unweighted case
     */
    this->weights *= residuals.n_elem / tmp;
}

static inline double wgtBisquare2(double x, double c)
{
    if (fabs(x) > (c)) {
        return(0.);
    }

    x /= c;
    x = (1 - x) * (1 + x);
    return x * x; // * 6 / (c * c);
}
