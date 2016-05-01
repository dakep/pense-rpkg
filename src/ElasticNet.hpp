//
//  ElasticNet.hpp
//  penseinit
//
//  Created by David Kepplinger on 2016-01-31.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#ifndef ElasticNet_hpp
#define ElasticNet_hpp

#include "config.h"

#include <RcppArmadillo.h>

#include "Data.hpp"


class ElasticNet
{
public:
	ElasticNet(const double eps, const bool center) : center(center), eps(eps)
	{
	}

	virtual ~ElasticNet()
    {}

	/*
	 * Set the regularization based on two independent
	 * lambda values.
	 * The conversion formula is:
	 * lambda = lambda2 + lambda1 / 2
	 * alpha = lambda1 / (2 * lambda)
	 *
	 * NOTE: The values for lambda1 and lambda2 are INDEPENDENT
	 *		 of the number of observations! This is different
	 *		 from the R package lars!
	 *		 lambda1 = lambda_lars / numObs
	 *
	 */
	virtual void setLambdas(const double lambda1, const double lambda2) = 0;

	virtual void setAlphaLambda(const double alpha, const double lambda) = 0;

	void setThreshold(const double eps)
	{
		this->eps = eps;
	}

	/*
	 * Solve the EN problem
	 *
	 * (1 / 2N) * RSS + lambda * (((1 - alpha) / 2) * L2(beta) + alpha * L1(beta))
	 *
	 * which is equivalent to solving
	 *
	 * (1 / N) * RSS + lambda2 * L2(beta) + lambda1 * L1(beta))
	 *
	 * @param data Is assumed to have the leading column of 1's for the intercept
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
    virtual bool computeCoefs(const Data& data, double *RESTRICT coefs,
							  double *RESTRICT residuals, const bool warm = FALSE) = 0;

protected:
    const bool center;
    double eps;
};






class ElasticNetGDESC : public ElasticNet
{
public:
	ElasticNetGDESC(const int maxIt, const double eps, const bool center);
	~ElasticNetGDESC();

	/*
	 * Set the regularization based on two independent
	 * lambda values.
	 * The conversion formula is:
	 * lambda = lambda2 + lambda1 / 2
	 * alpha = lambda1 / (2 * lambda)
	 *
	 * NOTE: The values for lambda1 and lambda2 are INDEPENDENT
	 *		 of the number of observations! This is different
	 *		 from the R package lars!
	 *		 lambda1 = lambda_lars / numObs
	 *
	 */
	void setLambdas(const double lambda1, const double lambda2);

	void setAlphaLambda(const double alpha, const double lambda);

	/*
	 * Solve the EN problem
	 *
	 * (1 / 2N) * RSS + lambda * (((1 - alpha) / 2) * L2(beta) + alpha * L1(beta))
	 *
	 * which is equivalent to solving
	 *
	 * (1 / N) * RSS + lambda2 * L2(beta) + lambda1 * L1(beta))
	 *
	 * @param data Is assumed to have the leading column of 1's for the intercept
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
    bool computeCoefs(const Data& data, double *RESTRICT coefs,
							  double *RESTRICT residuals, const bool warm = FALSE);

private:
	const int maxIt;
	double alpha;
	double lambda;

	void resizeBuffer(const Data& data);

	double *RESTRICT Xtr;
	double *RESTRICT Xmeans;
	double *RESTRICT Xvars;
	int XtrSize;
	int XmeansSize;
};





class ElasticNetLARS : public ElasticNet {
public:
    ElasticNetLARS(const double eps, const bool center);
    ~ElasticNetLARS();

	void setLambdas(const double lambda1, const double lambda2);

	void setAlphaLambda(const double alpha, const double lambda);

    bool computeCoefs(const Data& data, double *RESTRICT coefs, double *RESTRICT residuals,
                      const bool warm = FALSE);


private:
	double lambda1;
	double sqrtLambda2;

	arma::mat XtrAug;
	arma::vec yAug;
	arma::uword augNobs;

	void augmentData(const Data& data);

};

#endif /* ElasticNet_hpp */
