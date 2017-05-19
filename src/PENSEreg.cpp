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

static const double NUMERICAL_TOLERANCE = NUMERIC_EPS;

static inline double wgtBisquare2(double x, double c);


PENSEReg::PENSEReg(const Data& data, const double alpha, const double lambda,
			 const Control& ctrl) : data(data), alpha(alpha), lambda(lambda), ctrl(ctrl),
                                    rhoBisquare(getRhoFunctionByName(BISQUARE))

{
}

PENSEReg::~PENSEReg()
{
}


void PENSEReg::compute(double *RESTRICT currentCoef, double *RESTRICT residuals)
{
    ElasticNet *en = getElasticNetImpl(this->ctrl);
    RhoFunction rhoFun = getRhoFunctionByName(this->ctrl.mscaleRhoFun);

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

    this->scale = mscale(residuals, data.numObs(), this->ctrl.mscaleB, this->ctrl.mscaleEPS,
                         this->ctrl.mscaleMaxIt, rhoFun, this->ctrl.mscaleCC);

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
        tmp = 0;
        for (i = 0; i < this->data.numObs(); ++i) {
            weights[i] = wgtBisquare2(residuals[i] / this->scale, this->ctrl.mscaleCC);
            tmp += weights[i];
        }

        /*
         * Normalize weights to sum to n --> just as in the unweighted case
         */
        tmp = this->data.numObs() / tmp;
        for (i = 0; i < this->data.numObs(); ++i) {
            weights[i] *= tmp;
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
         * Compute scale
         */
        this->scale = mscale(residuals, this->data.numObs(), this->ctrl.mscaleB,
                             this->ctrl.mscaleEPS, this->ctrl.mscaleMaxIt, rhoFun,
                             this->ctrl.mscaleCC);

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
    return x * x * 6 / (c * c);
}

