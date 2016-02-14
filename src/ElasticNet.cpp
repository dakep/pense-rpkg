//
//  ElasticNet.cpp
//  penseinit
//
//  Created by David Kepplinger on 2016-01-31.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#include <cfloat>
#include <Rmath.h>

#include "ElasticNet.hpp"

static double softThreshold(const double z, const double gamma);

ElasticNet::ElasticNet(const int maxIt, const double eps) : maxIt(maxIt), eps(eps),
            XtrSize(0), XmeansSize(0)
{}

ElasticNet::~ElasticNet()
{
    if(this->XtrSize > 0) {
        delete[] this->Xtr;
    }

    if(this->XmeansSize > 0) {
        delete[] this->Xmeans;
        delete[] this->Xvars;
    }
}

void ElasticNet::setLambdas(const double lambda1, const double lambda2)
{
    this->lambda = lambda2 + lambda1 / 2;
    this->alpha = lambda1 / (2 * lambda2 + lambda1);
}

void ElasticNet::setAlphaLambda(const double alpha, const double lambda)
{
    this->alpha = alpha;
    this->lambda = lambda;
}


bool ElasticNet::computeCoefs(const Data& data, double *RESTRICT coefs, double *RESTRICT residuals,
                             const bool center)
{
    /*
     * NOTE:
     * We have to divide the lambdas by N in order to get the same results as
     * from lars
     */
    const double la = (this->lambda * this->alpha);
    const double updateDenom = this->lambda * (1 - this->alpha);

    double yMean;
    double tmp;
    double coefChange;
    double totalChange;
    const double *RESTRICT yPtr = data.getYConst();
    const double *RESTRICT XtrConstIter = data.getXtrConst();
    const double *RESTRICT weight;
    const double *RESTRICT actualXtr;
    int i, j, iter;
    double centerDenom = 0;
    double *RESTRICT XtrIter;
    double norm = 0;
    double newNorm;

    this->resizeBuffer(data, center);

    memset(coefs, 0, data.numVar() * sizeof(double));

    /*
     * First we calculate the mean of y and the X variables
     */
    if (center) {
        yMean = 0;
        memset(this->Xmeans, 0, (data.numVar() - 1) * sizeof(double));

        weight = data.getXtrConst();
        for (i = 0; i < data.numObs(); ++i, weight += data.numVar()) {
            if (*weight > 0) {
                yMean += yPtr[i] * (*weight);

                /*
                 * If we deal with augmented data, the first column is 1's and 0's.
                 * Count only those rows which have a leading 1.
                 */
                centerDenom += (*weight);

                ++XtrConstIter; /* Skip intercept */
                for (j = 0; j < data.numVar() - 1; ++j, ++XtrConstIter) {
                    this->Xmeans[j] += (*XtrConstIter) * (*weight);
                }
            } else {
                XtrConstIter += data.numVar();
            }
        }

        if (centerDenom > 0) {
            yMean /= centerDenom;

            for (j = 0; j < data.numVar() - 1; ++j) {
                this->Xmeans[j] /= centerDenom;
            }
        } else {
            yMean = 0;
            memset(this->Xmeans, 0, (data.numVar() - 1) * sizeof(double));
        }

        /*
         * We start by "setting" all beta_j's to 0, i.e. the residuals
         * are the deviation from the average.
         *
         * And by centering the variables (at least for those observations with weight > 0)
         */
        XtrConstIter = data.getXtrConst();
        XtrIter = this->Xtr;
        weight = data.getXtrConst();
        for (i = 0; i < data.numObs(); ++i, weight += data.numVar()) {
            residuals[i] = yPtr[i] - yMean * (*weight);

            /* Skip the first row (intercept!) */
            ++XtrConstIter;
            ++XtrIter;
            for (j = 0; j < data.numVar() - 1; ++j, ++XtrConstIter, ++XtrIter) {
                *XtrIter = (*XtrConstIter) - this->Xmeans[j] * (*weight);
            }
        }

        actualXtr = this->Xtr;
    } else {
        yMean = 0;
        memset(this->Xmeans, 0, (data.numVar() - 1) * sizeof(double));
        memcpy(residuals, yPtr, data.numObs() * sizeof(double));
        actualXtr = data.getXtrConst();
    }


    /*
     * Compute length of the vectors of variables (= N * Var(X_j))
     */
    memset(this->Xvars, 0, (data.numVar() - 1) * sizeof(double));
    XtrConstIter = actualXtr;
    for (i = 0; i < data.numObs(); ++i) {
        ++XtrConstIter;
        for (j = 0; j < data.numVar() - 1; ++j, ++XtrConstIter) {
            this->Xvars[j] += (*XtrConstIter) * (*XtrConstIter);
        }
    }

    for (j = 0; j < data.numVar() - 1; ++j, ++XtrConstIter) {
        this->Xvars[j] /= data.numObs();
    }

    iter = 0;

    while(1) {
        totalChange = 0;
        coefChange = 0;
        newNorm = 0;

        /* Start at j = 1 because j = 0 is the intercept which will be calculated at the end */
        for (j = 1; j < data.numVar(); ++j) {

            /* Update coefficient beta_j */
            XtrConstIter = actualXtr + j;
            tmp = 0;

            for (i = 0; i < data.numObs(); ++i, XtrConstIter += data.numVar()) {
                tmp += *XtrConstIter * residuals[i];
            }


            tmp = softThreshold(tmp / data.numObs() + coefs[j] * this->Xvars[j - 1], la) /
                        (this->Xvars[j - 1] + updateDenom);


            coefChange = coefs[j] - tmp;
            coefs[j] = tmp;

            newNorm += fabs(tmp);

            /* Update residuals if the coefficient has changed */
            if (coefChange != 0) {
                XtrConstIter = actualXtr + j;
                for (i = 0; i < data.numObs(); ++i, XtrConstIter += data.numVar()) {
                    residuals[i] += (*XtrConstIter) * coefChange;
                }

                totalChange += fabs(coefChange);
            }
        }

        /*
         * Check for max. iterations
         */
        if ((++iter > this->maxIt) || (totalChange <= this->eps * norm)) {
            break;
        }

        norm = newNorm;
    }

    /* Update intercept */
    coefs[0] = yMean;
    for (j = 1; j < data.numVar(); ++j) {
        coefs[0] -= coefs[j] * this->Xmeans[j - 1];
    }

    /* Residuals are already adjusted! */
    return (iter <= this->maxIt);
}


void ElasticNet::resizeBuffer(const Data& data, const bool center)
{
    if (center && ((data.numObs() * data.numVar()) > (this->XtrSize))) {
        if(this->XtrSize > 0) {
            delete[] this->Xtr;
        }
        this->XtrSize = data.numObs() * data.numVar();
        this->Xtr = new double[this->XtrSize];
    }

    if (data.numVar() > this->XmeansSize) {
        if(this->XmeansSize > 0) {
            delete[] this->Xmeans;
            delete[] this->Xvars;
        }
        this->XmeansSize = data.numVar();
        this->Xmeans = new double[this->XmeansSize - 1];
        this->Xvars = new double[this->XmeansSize - 1];
    }
}


static inline double softThreshold(const double z, const double gamma)
{
    if (fabs(z) <= gamma) {
        return 0.;
    } else if (z < 0) {
        return z + gamma;
    }
    return z - gamma;
}
