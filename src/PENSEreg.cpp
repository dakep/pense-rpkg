//
//  PENSEreg.cpp
//  pense
//
//  Created by David Kepplinger on 2016-04-22.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//
#include "PENSEreg.hpp"

#include <RcppArmadillo.h>
#include <Rmath.h>

#include "Data.hpp"
#include "ElasticNet.hpp"
#include "olsreg.h"
#include "mscale.h"

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

void PENSEReg::updateWeights(const double *RESTRICT residuals)
{
    double tmp;
    int i;
    this->scale = mscale(residuals, this->data.numObs(), this->bdp, this->mscaleEps,
                         this->mscaleMaxIt, rhoBisquare2, this->cc);

    tmp = 0;
    for (i = 0; i < this->data.numObs(); ++i) {
        this->weights[i] = wgtBisquare2(residuals[i] / this->scale, this->cc);
        tmp += this->weights[i];
    }

    /*
     * Normalize weights to sum to n --> just as in the unweighted case
     */
    tmp = this->data.numObs() / tmp;
    for (i = 0; i < this->data.numObs(); ++i) {
        this->weights[i] *= tmp;
    }
}

static inline double wgtBisquare2(double x, double c)
{
    if (fabs(x) > (c)) {
        return(0.);
    }

    x /= c;
    x = (1 - x) * (1 + x);
    return x * x * 6 / (c * c);
}
