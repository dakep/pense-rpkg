//
//  PENSEreg.hpp
//  penseinit
//
//  Created by David Kepplinger on 2016-04-22.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#ifndef PENSEreg_hpp
#define PENSEreg_hpp

#include "mscale.h"
#include "Control.h"
#include "Data.hpp"

class PENSEReg {
public:
	PENSEReg(const Data& data, const double alpha, const double lambda,
			 const Control& ctrl);
	~PENSEReg();

	void compute(double *RESTRICT coefficients, double *RESTRICT residuals);


	int getIterations() const
	{
		return this->iteration;
	}

	double getScale() const
	{
		return this->scale;
	}

	double getRelChange() const
	{
		return this->relChange;
	}

private:
	const Data& data;
	const double alpha;
	const double lambda;
	const Control& ctrl;
	const RhoFunction rhoBisquare;

	int iteration;
	double scale;
	double relChange;
};

#endif /* PENSEreg_hpp */
