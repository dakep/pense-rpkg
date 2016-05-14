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
    const double *RESTRICT yIter;
    const double *RESTRICT XtrIter;

    ElasticNet *en = getElasticNetImpl(this->ctrl);
    Data wgtData;

    bool enConverged;
    double *RESTRICT oldCoef = new double[data.numVar()];
    double *RESTRICT XtrWgtIter;
    double *RESTRICT yWgtIter;
    double *RESTRICT weightBeta = new double[data.numObs()];
    double weightBetaSum;
    double tmp;
    double norm2Old;
    int i, j;

    this->iteration = 0;

    en->setAlphaLambda(this->alpha, this->lambda);

    wgtData.setNumObs(data.numObs());
    wgtData.setNumVar(data.numVar());
    wgtData.resize();

    computeResiduals(data.getXtrConst(), data.getYConst(), data.numObs(), data.numVar(),
                     currentCoef, residuals);

    do {
        ++this->iteration;

        /*
         * Compute weights and weight observations
         */

        XtrWgtIter = wgtData.getXtr();
        yWgtIter = wgtData.getY();
        yIter = this->data.getYConst();
        XtrIter = this->data.getXtrConst();

        weightBetaSum = 0;
        for (i = 0; i < this->data.numObs(); ++i, ++yIter, ++yWgtIter) {
            weightBeta[i] = wgtBisquare2(residuals[i] / scale, this->ctrl.mscaleCC);
            tmp = sqrt(weightBeta[i]);
            weightBetaSum += weightBeta[i];

            *yWgtIter = ((*yIter) - currentCoef[0]) * tmp;

            /* Skip first row of 1's! */
            ++XtrWgtIter;
            ++XtrIter;

            for (j = 1; j < data.numVar(); ++j, ++XtrWgtIter, ++XtrIter) {
                *XtrWgtIter = (*XtrIter) * tmp;
            }
        }

        /*
         * For the coordinate descent EN algorithm we have to perform
         * some additional steps:
         *  - adjust convergency threshold
         */
        if (this->ctrl.enAlgorithm == GRADIENT_DESCENT) {
            tmp = mscale(wgtData.getY(), this->data.numObs(), 0.5, 1e-8, 200, rhoBisquare, 1.54764);
            en->setThreshold(this->ctrl.enEPS * this->data.numObs() * tmp * tmp);
        }

        /*
         * Copy current coefficients to check for convergence later
         */
        memcpy(oldCoef, currentCoef, this->data.numVar() * sizeof(double));

        /*
         * Perform EN using current coefficients as warm start (only applicable for the coordinate
         * descent algorithm)
         */
        enConverged = en->computeCoefs(wgtData, currentCoef, residuals, true);

        if (!enConverged) {
            Rcpp::warning("Weighted elastic net did not converge");
        }

        /*
         * Compute residuals for ORIGINAL data (not weighted)
         */
        computeResiduals(data.getXtrConst(), data.getYConst(), data.numObs(), data.numVar(),
                         currentCoef, residuals);

        /*
         * Compute intercept and adjust residuals
         */
        currentCoef[0] = 0;
        for (i = 0; i < data.numObs(); ++i) {
            currentCoef[0] += residuals[i] * weightBeta[i];
        }
        currentCoef[0] /= weightBetaSum;

        for (i = 0; i < data.numObs(); ++i) {
            residuals[i] -= currentCoef[0];
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

    wgtData.free();
    delete en;
    delete[] oldCoef;
    delete[] weightBeta;
}


static inline double wgtBisquare2(double x, double c)
{
    if (fabs(x) > (c)) {
        return(0.);
    }

    x /= c;
    x = (1 - x * x);
    return 6. * x * x / (c * c);
}

