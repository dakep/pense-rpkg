//
//  PSCxx.hpp
//  pense
//
//  Created by David Kepplinger on 2016-01-24.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#ifndef PSCxx_hpp
#define PSCxx_hpp

#include "config.h"
#include "Data.hpp"
#include "ElasticNet.hpp"


class PSC
{
public:
    virtual ~PSC();

	virtual void setData(const Data &data);

    /**
     * Set the residuals for the next computation of the PSCs
     */
    virtual void setResiduals(const double *RESTRICT residuals) = 0;

    /**
     * Compute the PSC given the previously set residuals
     */
	virtual int computePSC() = 0;

    /**
     * Get the most recently computed PSCs
     */
	virtual const double* getPSC() = 0;

protected:
	PSC(bool nobsEigenvalues);

	int doEigenDecomposition(const char* const uplo, double *RESTRICT matrix,
							 const int n);

	const double* getEigenvalues() const
	{
		return this->evalues;
	}

	const double* getEigenvectors() const
	{
		return this->evectors;
	}


	Data data;
private:
	const bool nobsEigenvalues;

    double *RESTRICT EVDworkMem;
    int *RESTRICT EVIworkMem;
    int EVDworkMemSize;
    int EVIworkMemSize;

    double *RESTRICT evalues;
    double *RESTRICT evectors;
    int *RESTRICT evectorsSupport;
};

class PSC_OLS : public PSC {
public:
    PSC_OLS();
    ~PSC_OLS();

    virtual void setData(const Data &data);

    virtual void setResiduals(const double *RESTRICT residuals);

	/**
	 * Make sure to set the residuals before calling
	 * this method!
	 */
    virtual int computePSC();

	virtual const double* getPSC()
    {
		return this->Z;
	}

	/**
	 * Use an external memory for Xsqrt that will be
	 * updated by the caller before calling computePSC
	 */
	void setXsqrtMemory(double *RESTRICT Xsqrt);

private:
	bool XsqrtProvided;
    bool initialized;

    double *RESTRICT Z;

    /* We can keep a few memory regions as their sizes will not change */
    double *RESTRICT Xsqrt;
    double *RESTRICT XsqrtInvX;
    double *RESTRICT Q;
    double *RESTRICT residuals;
};

class PSC_EN : public PSC {
public:
    PSC_EN(ElasticNet &en);
    ~PSC_EN();

    virtual void setData(const Data &data);

    virtual void setResiduals(const double *RESTRICT residuals);

    virtual int computePSC();

	virtual const double* getPSC()
    {
		return this->Z;
	}

private:
    bool initialized;

	ElasticNet &en;

    double *RESTRICT Z;
    double *RESTRICT residMat;
	double *RESTRICT buffer;
    double *RESTRICT residuals;
};

#endif /* PSCxx_hpp */
