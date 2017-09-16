//
//  PENSEreg.hpp
//  pense
//
//  Created by David Kepplinger on 2016-04-22.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#ifndef PENSEreg_hpp
#define PENSEreg_hpp

#include "config.h"
#include <RcppArmadillo.h>

#include "IRWEN.hpp"
#include "mscale.h"
#include "Data.hpp"

class PENSEReg : public IRWEN {
public:
	PENSEReg(const Data& data, const double alpha, const double lambda, const Options& opts, const Options &enOpts);
	~PENSEReg();

    double getScale() const
    {
        return this->scale;
    }

protected:
    void updateWeights(const arma::vec& residuals);
    void updateObjective(const arma::vec& residuals, const double betaENPenalty);

private:
    const double bdp;       // breakdown point
    const double cc;        // tuning constant to achieve the breakdown point
    const double mscaleEps; // tolerance for mscale
    const int mscaleMaxIt;  // max iterations for mscale

    double scale;
};

#endif /* PENSEreg_hpp */
