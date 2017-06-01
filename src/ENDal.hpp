//
//  ENDal_hpp
//  pense
//
//  Created by David Kepplinger on 2017-05-21.
//  Copyright © 2017 David Kepplinger. All rights reserved.
//

#ifndef ENDal_hpp
#define ENDal_hpp

#include "config.h"

#include <RcppArmadillo.h>

#include "Data.hpp"
#include "ElasticNet.hpp"


/**
 * Solve the EN problem
 *
 * (1 / 2N) * RSS + lambda * (((1 - alpha) / 2) * L2(beta)^2 + alpha * L1(beta))
 *
 */
class ENDal : public ElasticNet
{
public:
	ENDal(const bool intercept);
	ENDal(const bool intercept, const Options& options);
	~ENDal();

    void setData(const Data& data);
    void setOptions(const Options& options);
	void setLambdas(const double lambda1, const double lambda2);
	void setAlphaLambda(const double alpha, const double lambda);
    void computeCoefs(double *RESTRICT coefs, double *RESTRICT residuals);
    void computeCoefsWeighted(double *RESTRICT coefs, double *RESTRICT residuals, const double *RESTRICT weights);

private:
    int maxIt;
    double etaMultiplier;
    double eps;
    double etaStart;
    double etaStartNumerator;
    bool warmStart;

    double lambda;
    double alpha;

    double nLambda;
    double eta[2];

    int bufferSizeNobs;

    arma::vec a;
    arma::vec y;
    arma::mat Xtr;

    double fullObjectiveFun(const double intercept, const arma::vec& beta);

    bool minimizePhi(arma::vec& beta, double& intercept);

    double evalPhi(const arma::vec& a, arma::vec& beta, double& intercept, arma::vec &grad, arma::mat& hess, bool evalGrad);
};


#endif /* ENDal_hpp */
