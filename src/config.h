/*
 * config.h
 * penseinit
 *
 * Created by David Kepplinger on 2016-01-22.
 * Copyright (c) 2015 David Kepplinger. All rights reserved.
 *
 * This file is part of the R package penseinit.
 *
 * penseinit is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * penseinit is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with R. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef penseinit_config_h
#define penseinit_config_h

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


#ifndef DEBUG
#   define ARMA_NO_DEBUG
#endif


#define NUMERIC_EPS 1e-16
#define LAPACK_EV_MIN 1e-12

#endif


