//
//  PartialSort.h
//  pense
//
//  Created by David Kepplinger on 2016-01-30.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#ifndef PartialSort_h
#define PartialSort_h

#ifdef __cplusplus
extern "C" {
#endif

typedef double (*CompareFunction)(const double, const double);

double getQuantile(const double * values, const int length, const double quantile,
                   CompareFunction compare);

void partialQsort(double *values, const int lower, const int middle, const int upper,
                  CompareFunction compare);

#ifdef __cplusplus
}
#endif

#endif /* PartialSort_h */
