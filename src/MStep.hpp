//
//  MStep.hpp
//  pense
//
//  Created by David Kepplinger on 2016-05-13.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#ifndef MStep_hpp
#define MStep_hpp

#include "IRWEN.hpp"
#include "Options.hpp"
#include "Data.hpp"
#include "mscale.h"

class MStep : public IRWEN {
public:
	MStep(const Data& data, const double alpha, const double lambda, const double scale, const Options& opts, Options &enOpts);
	~MStep();

protected:
    void updateWeights(const double *RESTRICT residuals);

private:
    const double cc;    // tuning constant
    const double scale; // initial scale
};

#endif /* MStep_hpp */
