//
//  PSC.hpp
//  penseinit
//
//  Created by David Kepplinger on 2016-01-24.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#ifndef PenaYohai_hpp
#define PenaYohai_hpp

#include "config.h"
#include "Data.hpp"


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
	PSC();

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

    virtual int computePSC();

	virtual const double* getPSC()
    {
		return this->Z;
	}

private:
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
    PSC_EN();
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

    double *RESTRICT Z;
    double *RESTRICT residMat;
    double *RESTRICT residuals;
};

#endif /* PenaYohai_hpp */
