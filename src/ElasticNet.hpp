//
//  ElasticNet.hpp
//  pense
//
//  Created by David Kepplinger on 2016-01-31.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#ifndef ElasticNet_hpp
#define ElasticNet_hpp

#include "config.h"
#include <string>
#include <stdexcept>

#include <RcppArmadillo.h>

#include "Data.hpp"
#include "Options.hpp"

/**
 * Solve the EN problem
 * argmin_{beta0, beta} (1 / 2N) * L2(y - beta0 - X . beta)^2 +
 *						 lambda * (((1 - alpha) / 2) * L2(beta)^2 + alpha * L1(beta))
 *
 * or the weighted EN problem
 * argmin_{beta0, beta} (1 / 2N) * L2(weights * (y - beta0 - X . beta))^2 +
 *						 lambda * (((1 - alpha) / 2) * L2(beta)^2 + alpha * L1(beta))
 *
 *
 */
class ElasticNet
{
public:
    class Error : public std::runtime_error
    {
    public:
        Error(const std::string& what) : std::runtime_error(what)
        {}
    };

	ElasticNet(const bool intercept) : intercept(intercept), status(0), statusMessage("")
	{
	}

	virtual ~ElasticNet()
    {
    }

    /*
     * Set algorithm-specific options
     */
    virtual void setOptions(const Options& options) = 0;

    virtual int getStatus()
    {
        return this->status;
    }

    virtual const std::string& getStatusMessage()
    {
        return this->statusMessage;
    }

    /*
     * Set the data for the upcoming operations
     *
     * @param data Is assumed to have the leading column of 1's for the intercept (even if
     *             it won't be estimated, or assumed to be 0).
     */
    virtual void setData(const Data& data) = 0;

	/**
	 * Set the regularization based on one parameter for the L1 penalization
     * and one for the L2 penalization, i.e.,
     * lambda = lambda2 + lambda1
     * alpha = lambda1 / (lambda2 + lambda1)
	 *
	 * NOTE: The values for lambda1 and lambda2 are INDEPENDENT
	 *		 of the number of observations!
	 */
	virtual void setLambdas(const double lambda1, const double lambda2) = 0;

	/**
	 * Set the regularization based on alpha and lambda
	 *
	 * NOTE: The values for lambda and alpha are INDEPENDENT
	 *		 of the number of observations!
	 */
	virtual void setAlphaLambda(const double alpha, const double lambda) = 0;

    /**
     * Solve the EN problem
     *
     * argmin_{beta0, beta} (1 / 2N) * L2(y - beta0 - X . beta)^2 +
	 *						 lambda * (((1 - alpha) / 2) * L2(beta)^2 + alpha * L1(beta))
     *
     *
     * @param coefs Is assumed to be at least data.numVar() long! This can be taken as
     *          the starting point for the coordinate-descend algorithm. (see argument warm)
     * @param residuals A vector of residuals with as many observations as in data. The
     *          residuals will be recalculated for the given (warm) coefficients
     *
     * NOTE: The leading column of X is used as weight for the row in
     *		 in centering the data.
     *
     * @returns TRUE if the algorithm converged, FALSE otherwise
     */
    virtual void computeCoefs(double *RESTRICT coefs, double *RESTRICT residuals) = 0;
    virtual void computeCoefs(double& intercept, arma::sp_vec& beta, arma::vec& residuals) = 0;

    /**
     * Solve the weighted EN problem
     *
     * argmin_{beta0, beta} (1 / 2N) * L2(weights * (y - beta0 - X . beta))^2 +
	 *						 lambda * (((1 - alpha) / 2) * L2(beta)^2 + alpha * L1(beta))
     *
     *
     * @param coefs Is assumed to be at least data.numVar() long! This can be taken as
     *          the starting point for the coordinate-descend algorithm. (see argument warm)
     * @param residuals A vector of residuals with as many observations as in data. The
     *          residuals will be recalculated for the given (warm) coefficients
     * @param weights a vector of weights for each observation. NOTE: the weights will not
     *          be normalized!
     *
     * NOTE: The leading column of X is used as weight for the row in centering the data.
     *
     * @returns TRUE if the algorithm converged, FALSE otherwise
     */
    virtual void computeCoefsWeighted(double *RESTRICT coefs, double *RESTRICT residuals,
									  const double *RESTRICT weights) = 0;
    virtual void computeCoefsWeighted(double& intercept, arma::sp_vec& coefs, arma::vec& residuals,
                                      const arma::vec& weights) = 0;

protected:
    const bool intercept;

    int status;
    std::string statusMessage;
};

/**
 * Convenience functions to choose one of the above Elastic Net implementation.
 */
ElasticNet* getElasticNetImpl(const Options& options, const bool intercept);


#endif /* ElasticNet_hpp */
