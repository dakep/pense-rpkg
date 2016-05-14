//
//  InitialEstimator.hpp
//  penseinit
//
//  Created by David Kepplinger on 2016-01-26.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#ifndef InitialEstimator_h
#define InitialEstimator_h

#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>


#include "Control.h"
#include "Data.hpp"
#include "PSCxx.hpp"
#include "PartialSort.h"
#include "mscale.h"
#include "ElasticNet.hpp"

class InitialEstimator
{
public:
	virtual ~InitialEstimator();

	int compute();

	const double* getInitialEstimators() const
	{
		return this->allCoefEstimates;
	}

	const double* getObjectiveFunctionScores() const
	{
		return this->coefObjFunScore;
	}

protected:
	InitialEstimator(const Data& originalData, const Control& ctrl, PSC& psc,
					 const int maxEstimators, const int coefMemIncrease = 0);

	const Data& originalData;
	const Control& ctrl;

    /**
     * Reset to the original data
     */
    virtual void resetData();


	/**
	 * Evaluate the current estimate on the
	 * objective function
	 */
	virtual double evaluateEstimate() const = 0;

	/**
	 * Temporarily filter the `complete` data based on the supplied ordering
	 * and the configured threshold or proportion for PSCs
     */
	virtual void filterDataPSC(const double *RESTRICT values, const double threshold,
							   CompareFunction compare) = 0;


	/**
	 * Permanently filter the data based on the supplied ordering
     * Thus, the filtered data becomes the `complete` data
	 */
	virtual void filterDataResiduals(double threshold);

	/**
	 * Estimate the regression coefficients using the current working data
	 * (i.e. the data filtered by residuals and PSCs)
	 * AND update the residuals vector with the residuals for all
	 * observations in the residual-filtered data
	 */
	virtual void estimateCoefficients() = 0;

	PSC& psc;

	Data residualFilteredData;
	/**
	 * The number of *useful* observations in workData after
	 * the original data was filtered by residuals
	 */
	int residualFilteredNobs;

	double *RESTRICT coefEst;
	double *RESTRICT residuals;

	double MscaleOfResiduals() const;
	RhoFunction rhoFun;

private:
	double *RESTRICT allCoefEstimates;
	double *RESTRICT coefObjFunScore;
};


class OLS : public InitialEstimator
{
public:
    OLS(const Data& originalData, const Control& ctrl);
    virtual ~OLS();

protected:

    virtual void resetData();
    virtual void filterDataResiduals(double threshold);
    virtual void filterDataPSC(const double *RESTRICT values, const double threshold,
							   CompareFunction compare);
    virtual void estimateCoefficients();
	virtual double evaluateEstimate() const;

private:
    PSC_OLS pscOls;
    Data pscFilteredData;
	Data dataToUse;

    double *RESTRICT XtX;
};

class ENPY : public InitialEstimator
{
public:
    ENPY(const Data& originalData, const Control& ctrl);
    virtual ~ENPY();

protected:
    virtual void resetData();
    virtual void filterDataResiduals(double threshold);
    virtual void filterDataPSC(const double *RESTRICT values, const double threshold,
							   CompareFunction compare);
    virtual void estimateCoefficients();
	virtual double evaluateEstimate() const;

private:
    const double lambdaLS;

	ElasticNet &en;

    PSC_OLS pscOls;
    Data pscFilteredData;
	Data dataToUse;

	double currentNullDeviance;

	double *RESTRICT XtX;
};


class ENPY_Exact : public InitialEstimator
{
public:
    ENPY_Exact(const Data& originalData, const Control& ctrl);
    virtual ~ENPY_Exact();

protected:
    virtual void resetData();
    virtual void filterDataResiduals(double threshold);
    virtual void filterDataPSC(const double *RESTRICT values, const double threshold,
							   CompareFunction compare);
    virtual void estimateCoefficients();
	virtual double evaluateEstimate() const;

private:
    const double lambdaLS;

	ElasticNet &en;

    PSC_EN pscEn;
    Data pscFilteredData;
	Data dataToUse;

	double currentNullDeviance;
};


#endif /* InitialEstimator_h */
