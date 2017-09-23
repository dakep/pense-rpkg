//
//  PartialSort.c
//  pense
//
//  Created by David Kepplinger on 2016-01-30.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//  the partial quicksort algorithm is (C) Copyright by Ariel Faigon, 1987
//

#include "config.h"

#include <Rmath.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "PartialSort.h"

double getQuantile(const double *values, const int length, const double quantile,
                          CompareFunction compare)
{
    double retVal = 0;
    int quantIndex = (int) ceil(length * quantile);
	int start = 0;
	int end = quantIndex;
	int less = 0;
	int i;
    double *restrict valuescpy = (double*) malloc((length + 1) * sizeof(double));
	double cmp;

    memcpy(valuescpy, values, length * sizeof(double));

    valuescpy[length] = compare(DBL_MAX, 0.0);

	while (less < quantIndex && end <= length) {
		partialQsort(valuescpy, start, end + 2, length - 1, compare);

		retVal = valuescpy[end];

		less = 0;
		for (i = 0; i < end && less < quantIndex; ++i) {
			less += (compare(valuescpy[i], retVal) < 0);
		}

		end = ((end + quantIndex < length) ? (end + quantIndex) : length);
		start += quantIndex;
	}

	if (start > 0) {
		less = 0;
		retVal = valuescpy[quantIndex];
		for (i = 0; i < end && less < quantIndex; ++i) {
			cmp = compare(valuescpy[i], retVal);
			if (cmp < 0) {
				++less;
			} else if (cmp > 0) {
				retVal = valuescpy[i];
				less = i;
			}
		}
	}

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
