//
//  ENLars.hpp
//  pense
//
//  Created by David Kepplinger on 2016-01-31.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#ifndef ENLars_hpp
#define ENLars_hpp

#include "config.h"

#include <RcppArmadillo.h>

#include "Data.hpp"
#include "ElasticNet.hpp"

class ENLars : public ElasticNet {
public:
	enum UseGram {
		AUTO = 0,
		YES,
		NO
	};

    ENLars(const bool intercept);
    ENLars(const bool intercept, const Options& options);
    ~ENLars();

    void setOptions(const Options& options);
    void setData(const Data& data);
	void setLambdas(const double lambda1, const double lambda2);
	void setAlphaLambda(const double alpha, const double lambda);

    void computeCoefs(double *RESTRICT coefs, double *RESTRICT residuals);
    void computeCoefs(double& intercept, arma::sp_vec& coefs, arma::vec& residuals);
	void computeCoefsWeighted(double *RESTRICT coefs,
							  double *RESTRICT residuals,
							  const double *RESTRICT weights);
	void computeCoefsWeighted(double& intercept,
                              arma::sp_vec &coefs,
							  arma::vec &residuals,
							  const arma::vec &weights);

private:
    /**
     * Automatically switch to non-Gram when more than the following
     * number of predictors are present.
     * 1400 ~ 15 MiByte of memory -- should be okay for most systems
     */
	static const int MAX_PREDICTORS_GRAM = 1400;

	double lambda1;
	double sqrtLambda2;

    double eps;
	UseGram gramMode;

    arma::vec yOrig;
	arma::mat XtrAug;
	arma::vec yAug;

	/**
	 * Certain matrices/vectors which are useful to store
	 * between calls if the size of the data doesn't change too often.
	 */
	arma::mat gramMat;
	arma::vec corY;
	arma::vec meanX;

	void augmentData();

	void augmentedLASSO(arma::vec& coefs, arma::vec& residuals, const arma::uword nobs,
						const bool intercept);

	void augmentedOLS(arma::vec& coefs, arma::vec& residuals, const arma::uword nobs,
					  const bool intercept);
};

#endif /* ENLars_hpp */
