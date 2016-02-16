//
//  PartialSort.c
//  penseinit
//
//  Created by David Kepplinger on 2016-01-30.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#include <Rmath.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "PartialSort.h"

double getQuantile(const double *values, const int length, const double quantile,
                          CompareFunction compare)
{
    double retVal;
    int quantIndex = (int) ceil(length * quantile);
    double *restrict valuescpy = (double*) malloc((length + 1) * sizeof(double));

    memcpy(valuescpy, values, length * sizeof(double));

    valuescpy[length] = compare(DBL_MAX, 0.0);

    partialQsort(valuescpy, 0, quantIndex + 2, length - 1, compare);

    retVal = valuescpy[quantIndex];
    free(valuescpy);

    return retVal;
}

/**
 * NOTE: `values` must be one longer than upper and must hold the maximum possible value!
 */
/*
 *  (C) Copyright by Ariel Faigon, 1987
 *  Released under the GNU GPL (General Public License) version 2
 *  or any later version (http://www.gnu.org/licenses/licenses.html)
 */
void partialQsort(double *values, const int lower, const int middle, const int upper,
                  CompareFunction compare)
{
    int i, j;
    double tmp;
	double pivot;

    if (lower < upper) {
        tmp = values[lower];
        values[lower] = values[(upper + lower) / 2];
        values[(upper + lower) / 2] = tmp;

        i = lower;
        j = upper + 1;
        pivot = values[lower];

        while (1) {
            do ++i; while (compare(values[i], pivot) < 0);
            do --j; while (compare(values[j], pivot) > 0);


            if (j < i) break;

            tmp = values[i];
            values[i] = values[j];
            values[j] = tmp;
        }

        tmp = values[lower];
        values[lower] = values[j];
        values[j] = tmp;

        partialQsort(values, lower, middle, j - 1, compare);

        if (i < middle) {
            partialQsort(values, i, middle, upper, compare);
        }
    }
}
