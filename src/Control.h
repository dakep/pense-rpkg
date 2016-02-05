//
//  Control.h
//  penseinit
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

typedef struct ControlTag {
    /*
     * NOTE: Adjust lambdas for cc.scale BUT NOT FOR LS!
     * Also the EN solver used solves
     *  (1/N) * RSS + lambda1 || beta ||_1 + lambda2 || beta ||_2^2
     * so make sure the lambda values are ADJUSTED for this (i.e. divided by N).
     *
     * e.g.: lambda1_C = lambda1_R * cc.scale^2 / N
     */
    const double lambda1;
    const double lambda2;
    const int numIt;
    const double eps;
    const double residThreshold;
    const double residProportion;
    const double pscProportion;

    const int enMaxIt;
    const double enEPS;
    const int enCentering;

    const double mscaleB;
    const double mscaleCC;
    const int mscaleMaxIt;
    const double mscaleEPS;
    const RhoFunctionName mscaleRhoFun;
} Control;

#endif /* Control_h */
