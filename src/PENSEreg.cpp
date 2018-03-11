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

static inline double wgtBisquare2(double x, const double c);

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
    /*
     * The weights computed result in a quadratic majorizer for the true objective
     * at the current estimates.
     */
    double accu = 0;
    this->scale = mscale(
        residuals.memptr(),
        residuals.n_elem,
        this->bdp,
        this->mscaleEps,
        this->mscaleMaxIt,
        rhoBisquare2,
        this->cc
    );

    const double cc_scaled = this->scale * this->cc;

    for (uword i = 0; i < residuals.n_elem; ++i) {
        this->weights[i] = wgtBisquare2(residuals[i], cc_scaled);
        accu += this->weights[i] * residuals[i] * residuals[i];
    }

    this->weights *= 2 * residuals.n_elem * this->scale * this->scale / accu;
}

/**
 * Weight function rho'(x, c) / x for the UNSTANDARDIZED bisquare rho function.
 * This means, that rho(inf, c) != 1
 */
static inline double wgtBisquare2(double x, const double c)
{
    if (fabs(x) > (c)) {
        return(0.);
    }

    x /= c;
    x = (1 - x) * (1 + x);
    return x * x; /* this is missing the "* 6 / (c * c)" part which would be for the standardized rho function */
}
