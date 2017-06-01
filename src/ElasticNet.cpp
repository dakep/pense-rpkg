//
//  ElasticNet.cpp
//  pense
//
//  Created by David Kepplinger on 2016-01-31.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//
#include "config.h"
#include "ElasticNet.hpp"
#include "ENLars.hpp"
#include "ENDal.hpp"
#include "Options.hpp"

static const ENAlgorithm DEFAULT_OPT_EN_ALGORITHM = AUGMENTED_LARS;

/****************************************************************************
 *
 * Convenience functions to choose one of the two possible
 * Elastic Net implementations
 *
 ***************************************************************************/
ElasticNet* getElasticNetImpl(const Options& options, const bool intercept)
{
    ElasticNet* en = NULL;
    ENAlgorithm enAlgorithm = options.get("algorithm", DEFAULT_OPT_EN_ALGORITHM);

    switch (enAlgorithm) {
        case DAL:
            en = new ENDal(intercept, options);
            break;
        case AUGMENTED_LARS:
        default:
            en = new ENLars(intercept, options);
            break;
    }

    return en;
}
