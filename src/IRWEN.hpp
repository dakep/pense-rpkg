//
//  IRWEN.hpp
//  pense
//
//  Created by David Kepplinger on 2017-05-31.
//  Copyright © 2017 David Kepplinger. All rights reserved.
//

#ifndef IRWEN_hpp
#define IRWEN_hpp

#include "Options.hpp"
#include "ElasticNet.hpp"
#include "Data.hpp"

class IRWEN {
public:
	IRWEN(const Data& data, const double alpha, const double lambda, const Options& opts, const Options &enOpts);
	~IRWEN();

	void compute(double *RESTRICT coefficients, double *RESTRICT residuals);

	int iterations() const
	{
		return this->iteration;
	}

	double relChange() const
	{
		return this->relChangeVal;
	}

protected:
    virtual void updateWeights(const double *RESTRICT residuals) = 0;

	const Data& data;
    double* weights;

private:
    const int verbosity;
    const int maxIt;
    const double eps;
    const double lambda;

    ElasticNet* en;

	int iteration;
	double relChangeVal;
};

#endif /* IRWEN_hpp */
