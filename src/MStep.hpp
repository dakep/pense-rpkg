//
//  MStep.hpp
//  pense
//
//  Created by David Kepplinger on 2016-05-13.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#ifndef MStep_hpp
#define MStep_hpp

#include "Control.h"
#include "Data.hpp"
#include "mscale.h"

class MStep {
public:
	MStep(const Data& data, const double alpha, const double lambda,
          const Control& ctrl);
	~MStep();

	void compute(double *RESTRICT coefficients, const double scale, double *RESTRICT residuals);

	int getIterations() const
	{
		return this->iteration;
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
	double relChange;
};

#endif /* MStep_hpp */
