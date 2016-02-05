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
    virtual void setResiduals(const double *RESTRICT residuals) = 0;
	virtual int computePSC() = 0;
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

	virtual void free();


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
	bool gotResiduals;
    bool initialized;
    void free();

    double *RESTRICT Z;

    /* We can keep a few memory regions as their sizes will not change */
    double *RESTRICT Xsqrt;
    double *RESTRICT XsqrtInvX;
    double *RESTRICT H;
    double *RESTRICT residuals;
};

//class PSC_Lasso : public PSC
//{
//public:
//    void computePSC();
//    const double *const getPSC();
//};

#endif /* PenaYohai_hpp */
