//
//  MStep.cpp
//  penseinit
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

static const double NUMERICAL_TOLERANCE = NUMERIC_EPS;

static inline double wgtBisquare2(double x, double c);


MStep::MStep(const Data& data, const double alpha, const double lambda,
			 const Control& ctrl) : data(data), alpha(alpha), lambda(lambda), ctrl(ctrl),
                                    rhoBisquare(getRhoFunctionByName(BISQUARE))

{
}


MStep::~MStep()
{
}


void MStep::compute(double *RESTRICT currentCoef, const double scale, double *RESTRICT residuals)
{
    ElasticNet *en = getElasticNetImpl(this->ctrl);

    bool enConverged;
    double *RESTRICT oldCoef = new double[data.numVar()];
    double *RESTRICT weights = new double[data.numObs()];
    double tmp;
    double norm2Old;
    int i, j;

    this->iteration = 0;

    en->setAlphaLambda(this->alpha, this->lambda);

    computeResiduals(data.getXtrConst(), data.getYConst(), data.numObs(), data.numVar(),
                     currentCoef, residuals);


    /*
     * For the coordinate descent EN algorithm we have to perform
     * some additional steps:
     *  - adjust convergency threshold
     */
    if (this->ctrl.enAlgorithm == GRADIENT_DESCENT) {
        tmp = mscale(this->data.getYConst(), this->data.numObs(), 0.5, 1e-8, 200, rhoBisquare,
                     1.54764);
        en->setThreshold(this->ctrl.enEPS * this->data.numObs() * tmp * tmp);
    }

    do {
        ++this->iteration;

        /*
         * Compute weights
         */
        for (i = 0; i < this->data.numObs(); ++i) {
            weights[i] = wgtBisquare2(residuals[i] / scale, this->ctrl.mscaleCC);
        }

        /*
         * Copy current coefficients to check for convergence later
         */
        memcpy(oldCoef, currentCoef, this->data.numVar() * sizeof(double));

        /*
         * Perform EN using current coefficients as warm start (only applicable for the coordinate
         * descent algorithm)
         */
        enConverged = en->computeCoefsWeighted(this->data, currentCoef, residuals, weights, true);

        if (!enConverged) {
            Rcpp::warning("Weighted elastic net did not converge");
        }

        /*
         * Compute relative change
         */
        this->relChange = 0;
        norm2Old = 0;
        for (j = 0; j < data.numVar(); ++j) {
            tmp = oldCoef[j] - currentCoef[j];
            this->relChange += tmp * tmp;
            norm2Old += oldCoef[j] * oldCoef[j];
        }

        if (norm2Old < NUMERICAL_TOLERANCE) {
            if (this->relChange < NUMERICAL_TOLERANCE) {
                /* We barely moved away from the zero-vector --> "converged" */
                this->relChange = 0;
            } else {
                /* We moved away from the zero-vector --> continue */
                this->relChange = 2 * this->ctrl.eps;
            }
        } else {
            this->relChange /= norm2Old;
        }

    } while((this->iteration < this->ctrl.numIt) && (this->relChange > this->ctrl.eps));

    delete en;
    delete[] oldCoef;
    delete[] weights;
}


static inline double wgtBisquare2(double x, double c)
{
    if (fabs(x) > (c)) {
        return(0.);
    }

    x /= c;
    x = (1 - x) * (1 + x);
    return x * x; // * 6 / (c * c)?
}

