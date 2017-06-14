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
	enum Preconditioner {
		NONE = 0,
		DIAG,
		APPROX
	};

	ENDal(const bool intercept);
	ENDal(const bool intercept, const Options& options);
	~ENDal();

    void setData(const Data& data);
    void setOptions(const Options& options);
	void setLambdas(const double lambda1, const double lambda2);
	void setAlphaLambda(const double alpha, const double lambda);
    void computeCoefs(double *RESTRICT coefs, double *RESTRICT residuals);
    void computeCoefs(double& intercept, arma::sp_vec& beta, arma::vec& residuals);
    void computeCoefsWeighted(double *RESTRICT coefs, double *RESTRICT residuals, const double *RESTRICT weights);
    void computeCoefsWeighted(double& intercept, arma::sp_vec& beta, arma::vec& residuals,
                              const arma::vec& weights);

private:
    int verbosity;
    int maxIt;
    double eps;
    bool warmStart;
    double etaStart;
    double etaStartNumerator;
    double etaMultiplier;
    Preconditioner precondType;

    double lambda;
    double alpha;

    double nLambda;
    double eta[2];

    int bufferSizeNobs;
    int bufferSizeNvar;

    arma::vec* y;
    arma::mat* Xtr;
    arma::vec sqrtWeights;
    arma::mat sqrtWeightsOuter;
    bool useWeights;

    void dal(double& intercept, arma::sp_vec& beta);

    double fullObjectiveFun(const double intercept, const arma::sp_vec& beta);

    void evalPhiGrad(const arma::vec &a, const arma::sp_vec& beta, const double intercept, const double multFact, arma::vec &grad);

    int getPhiStepDir(arma::vec &stepDir, const arma::vec &grad, const arma::vec &a, const arma::sp_vec& beta, const double intercept, const double multFact);

    const arma::mat& getHessBuff(const arma::uvec& keep);
    arma::uvec hessBuffKeep;
    arma::mat hessBuff;
    arma::mat precond;
};


#endif /* ENDal_hpp */
