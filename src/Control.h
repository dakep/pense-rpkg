//
//  Control.h
//  pense
//
//  Created by David Kepplinger on 2016-01-30.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#ifndef Control_h
#define Control_h

typedef enum RhoFunctionNameTag {
    BISQUARE = 0,
    HUBER = 1,
    GAUSS_WEIGHT = 2
} RhoFunctionName;


typedef enum ENAlgorithmTag {
    GRADIENT_DESCENT = 0,
    AUGMENTED_LARS_GRAM = 1,
    AUGMENTED_LARS_NOGRAM = 2,
    AUGMENTED_LARS_AUTO = 3
} ENAlgorithm;

typedef struct ControlTag {
    const double lambda;
    const double alpha;
    const int numIt;
    const double eps;
    const double residThreshold;
    const double residProportion;
    const double pscProportion;

    const int enMaxIt;
    const double enEPS;
    const int enCentering;
    const ENAlgorithm enAlgorithm;

    const double mscaleB;
    const double mscaleCC;
    const int mscaleMaxIt;
    const double mscaleEPS;
    const RhoFunctionName mscaleRhoFun;
} Control;

#endif /* Control_h */
