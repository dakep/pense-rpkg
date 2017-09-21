//
//  MStep.hpp
//  pense
//
//  Created by David Kepplinger on 2016-05-13.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#ifndef MStep_hpp
#define MStep_hpp

#include "config.h"
#include <RcppArmadillo.h>

#include "IRWEN.hpp"
#include "Options.hpp"
#include "Data.hpp"
#include "mscale.h"

class MStep : public IRWEN {
public:
	MStep(const Data& data, const double alpha, const double lambda, const double scale, const Options& opts, const Options &enOpts);
	~MStep();

protected:
    void updateWeights(const arma::vec& residuals);
    void updateObjective(const arma::vec& residuals, const double betaENPenalty);

private:
    const double ccScaled;    // tuning constant times the initial scale
};

#endif /* MStep_hpp */
