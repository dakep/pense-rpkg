//
//  MStep.cpp
//  pense
//
//  Created by David Kepplinger on 2016-05-13.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//
#include "config.h"
#include <RcppArmadillo.h>

#include "MStep.hpp"

#include <Rmath.h>

#include "Data.hpp"
#include "ElasticNet.hpp"
#include "olsreg.h"
#include "mscale.h"

using namespace arma;

static const double DEFAULT_OPT_CC = 3.44;
static const RhoFunction rhoBisquare2 = getRhoFunctionByName(BISQUARE);

static inline double wgtBisquare2(double x, double c);

MStep::MStep(const Data& data, const double alpha, const double lambda, const double scale, const Options& opts, const Options &enOpts) :
    IRWEN(data, alpha, lambda, opts, enOpts),
    ccScaled(scale * opts.get("cc", DEFAULT_OPT_CC))
{
}

MStep::~MStep()
{
}

void MStep::updateObjective(const vec& residuals, const double betaENPenalty)
{
    this->objectiveVal = betaENPenalty;
    for (uword i = 0; i < residuals.n_elem; ++i) {
        this->objectiveVal += rhoBisquare2(residuals[i], this->ccScaled);
    }
}

void MStep::updateWeights(const vec& residuals)
{
    double tmp = 0;
    for (uword i = 0; i < residuals.n_elem; ++i) {
        this->weights[i] = wgtBisquare2(residuals[i], this->ccScaled);
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

