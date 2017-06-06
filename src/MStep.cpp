//
//  MStep.cpp
//  pense
//
//  Created by David Kepplinger on 2016-05-13.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#include "MStep.hpp"

#include <RcppArmadillo.h>
#include <Rmath.h>


#include "Data.hpp"
#include "ElasticNet.hpp"
#include "olsreg.h"
#include "mscale.h"

static const double DEFAULT_OPT_CC = 3.44;

static inline double wgtBisquare2(double x, double c);

MStep::MStep(const Data& data, const double alpha, const double lambda, const double scale, const Options& opts, Options &enOpts) :
    IRWEN(data, alpha, lambda, opts, enOpts),
    cc(opts.get("cc", DEFAULT_OPT_CC)),
    scale(scale)
{
}


MStep::~MStep()
{
}

void MStep::updateWeights(const double *RESTRICT residuals)
{
    int i;
    double tmp = 0;
    for (i = 0; i < this->data.numObs(); ++i) {
        this->weights[i] = wgtBisquare2(residuals[i], this->scale * this->cc);
        tmp += weights[i];
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
    return x * x; // * 6 / (c * c);
}

