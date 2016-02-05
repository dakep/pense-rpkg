//
//  mscale.c
//  penseinit
//
//  Created by David Kepplinger on 2016-01-30.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "mscale.h"
#include "PartialSort.h"

static const double MAD_SCALE_CONSTANT = 1.4826;

static double rhoBisquare(double x, const double c);
static double rhoHuber(double x, const double c);
static double rhoGaussWeight(double x, const double c);
static double absoluteLessThan(const double a, const double b);

/**
 * Calculate the M-Scale of a vector of numbers
 *
 * @param Rvalues numeric The vector of REAL's
 * @param Rlength integer The length of Rvalues
 * @param Rb      numeric The target average (right side) in the M-scale equation
 * @param Rcc     numeric The constant for the rho function
 * @param RmaxIt  integer The maximum number of iterations
 * @param Reps    numeric The relative tolerance for convergence
 * @param Rrhofun integer The rho function to use (see Control.h for possible values).
 *                        If the selected rho function is unknown, the bisquare function
 *                        will be used by default.
 *
 * @return numeric The M-scale
 */
SEXP C_mscale(SEXP Rvalues, SEXP Rlength, SEXP Rb, SEXP Rcc, SEXP RmaxIt, SEXP Reps, SEXP Rrhofun)
{
    SEXP Rscale = PROTECT(Rf_allocVector(REALSXP, 1));
    RhoFunction rhoFun = getRhoFunctionByName(*INTEGER(Rrhofun));

    *REAL(Rscale) = mscale(REAL(Rvalues), *INTEGER(Rlength), *REAL(Rb), *REAL(Reps),
                           *INTEGER(RmaxIt), rhoFun, *REAL(Rcc));

    UNPROTECT(1);
    return Rscale;
}




RhoFunction getRhoFunctionByName(RhoFunctionName name)
{
    switch (name) {
        case HUBER:
            return rhoHuber;
            break;
        case GAUSS_WEIGHT:
            return rhoGaussWeight;
            break;
        case BISQUARE:
        default:
            break;
    }

    return rhoBisquare;
}

double mscale(const double * x, const int n, const double b, const double eps, const int maxIt,
              RhoFunction rho, const double cc)
{
    int medianPos = (n / 2);
    int takeAvgForMedian = (n % 2);
    double scale;
    double tmpScale;
    double rhoSum;
    double rhoDenomInv = 1 / (b * n);
    int it, i;
    double err = 1 + eps;

    /*
     * Initialize the M-Scale with the MAD
     */
    double *restrict xcpy = (double*) malloc((n + 1) * sizeof(double));

    memcpy(xcpy, x, n * sizeof(double));
    xcpy[n] = DBL_MAX;

    partialQsort(xcpy, 0, medianPos + 2, n - 1, absoluteLessThan);

    scale = fabs(xcpy[medianPos]);

    if (takeAvgForMedian == 0) {
        scale = 0.5 * (scale + fabs(xcpy[medianPos - 1]));
    }

    scale *= MAD_SCALE_CONSTANT; // == MAD

    /*
     * Start iterations
     */

    it = 0;
    do {
        rhoSum = 0;
        for (i = 0; i < n; ++i) {
            rhoSum += rho(x[i] / scale, cc);
        }

        tmpScale = scale * sqrt(rhoSum * rhoDenomInv);
        err = fabs(tmpScale / scale - 1);
        scale = tmpScale;
    } while ((++it < maxIt) && (err > eps));

    free(xcpy);

    return scale;
}



/**
 * The rho functions are taken from R package robustbase_0.92-5
 * Copyright
 * Peter Rousseeuw and Christophe Croux, see file 'Copyrights';
 * Valentin Todorov <valentin.todorov@chello.at>,
 * Andreas Ruckstuhl <andreas.ruckstuhl@zhaw.ch>,
 * Matias Salibian-Barrera <matias@stat.ubc.ca>,
 * Tobias Verbeke <tobias.verbeke@openanalytics.eu>,
 * Manuel Koller <mkoller@ispm.unibe.ch>,
 * Martin Maechler
 */
static double rhoBisquare(double x, const double c)
{
    if (fabs(x) > (c)) {
        return(1.);
    }
	x /= c;
	x *= x;
	return (x * (3. + x * (-3. + x)));
}

static double rhoHuber(double x, const double c)
{
    return (fabs(x) <= c) ? (x * x * 0.5) : (c * (fabs(x) - c * 0.5));
}

static double rhoGaussWeight(double x, const double c)
{
    x /= c;
    return -expm1(-(x * x) * 0.5);
}

static double absoluteLessThan(const double a, const double b)
{
    return fabs(a) - fabs(b);
}
