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

static inline double wgtBisquareSTD2(double x, const double c);

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
    this->objectiveVal /= residuals.n_elem;
}

void MStep::updateWeights(const vec& residuals)
{
    for (uword i = 0; i < residuals.n_elem; ++i) {
        this->weights[i] = wgtBisquareSTD2(residuals[i], this->ccScaled);
    }
}

/**
 * Weight function rho'(x, c) / x for the STANDARDIZED bisquare rho function.
 * This means, that rho(inf, c) == 1 and thus wgt(x, s*c) == wgt(x/s, c) / s^2!
 */
static inline double wgtBisquareSTD2(double x, const double c)
{
    if (fabs(x) > (c)) {
        return(0.);
    }

    x /= c;
    x = (1 - x) * (1 + x);
    return x * x * 6 / (c * c);
}
