/*
 * config.h
 * pense
 *
 * Created by David Kepplinger on 2016-01-22.
 * Copyright (c) 2015 David Kepplinger. All rights reserved.
 *
 * This file is part of the R package pense.
 *
 * pense is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * pense is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with R. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef pense_config_h
#define pense_config_h

#include "autoconfig.h"

#ifdef HAVE_INTTYPES_H
#   include <inttypes.h>
#endif

#ifdef HAVE_STDINT_H
#   include <stdint.h>
#endif

#ifndef __cplusplus
#   undef RESTRICT
#   define RESTRICT restrict
#elif !defined(RESTRICT)
#   define RESTRICT
#endif

#ifndef ARMA_DONT_USE_FORTRAN_HIDDEN_ARGS
#   define ARMA_DONT_USE_FORTRAN_HIDDEN_ARGS
#endif

#if !defined DEBUG && !defined ARMA_NO_DEBUG
#   define ARMA_NO_DEBUG
#endif

#if !defined HAVE_OPENMP_CXX && !defined ARMA_DONT_USE_OPENMP
#   define ARMA_DONT_USE_OPENMP
#endif

#define STRICT_R_HEADERS

#define NUMERIC_EPS 1e-32
#define LAPACK_EV_MIN 1e-12

typedef enum RhoFunctionNameTag {
    HUBER = 0,
    BISQUARE = 1,
    GAUSS_WEIGHT = 5
} RhoFunctionName;


typedef enum ENAlgorithmTag {
    AUGMENTED_LARS = 0,
    DAL
} ENAlgorithm;


#endif


