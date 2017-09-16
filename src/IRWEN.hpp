//
//  IRWEN.hpp
//  pense
//
//  Created by David Kepplinger on 2017-05-31.
//  Copyright © 2017 David Kepplinger. All rights reserved.
//

#ifndef IRWEN_hpp
#define IRWEN_hpp

#include "config.h"
#include <RcppArmadillo.h>


#include "Options.hpp"
#include "ElasticNet.hpp"
#include "Data.hpp"


class IRWEN {
public:
	IRWEN(const Data& data, const double alpha, const double lambda, const Options& opts, const Options &enOpts);
	~IRWEN();

	void compute(double& intercept, arma::sp_vec& beta, arma::vec& residuals);

    void setLambda(const double lambda);

    arma::vec getWeights() const
    {
        return this->weights;
    }

	int iterations() const
	{
		return this->iteration;
	}

	double relChange() const
	{
		return this->relChangeVal;
	}

    double getObjective() const
    {
        return this->objectiveVal;
    }

protected:
    virtual void updateWeights(const arma::vec& residuals) = 0;
    virtual void updateObjective(const arma::vec& residuals, const double betaENPenalty) = 0;

    const arma::mat Xtr;
    const arma::vec y;
    arma::vec weights;

    double objectiveVal;

private:
    const int verbosity;
    const int maxIt;
    const double eps;
    const double alpha;
    double lambda;

    ElasticNet* en;

	int iteration;
	double relChangeVal;
};

#endif /* IRWEN_hpp */
