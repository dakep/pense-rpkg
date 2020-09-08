//
//  r_robust_utils.hpp
//  pense
//
//  Created by David Kepplinger on 2019-04-03.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef R_ROBUST_UTILS_HPP_
#define R_ROBUST_UTILS_HPP_

#include "nsoptim_forward.hpp"

namespace pense {
namespace r_interface {
//! Compute the tau-Scale of Centered Values
//!
//! @param x numeric values.
//! @return the tau-scale of the centered values.
SEXP TauSize(SEXP x) noexcept;

//! Compute the M-scale of Centered Values
//!
//! @param x numeric values.
//! @param mscale_opts a list of options for the M-scale equation.
//! @return the M-scale of the centered values.
SEXP MScale(SEXP x, SEXP mscale_opts) noexcept;

//! Compute the M-location
//!
//! @param x numeric values.
//! @param scale the scale of the values.
//! @param opts a list of options for the M-estimating equation.
//! @return the M-estimate of location.
SEXP MLocation(SEXP x, SEXP scale, SEXP opts) noexcept;

//! Compute the M-estimate of the Location and Scale
//!
//! @param x numeric values.
//! @param mscale_opts a list of options for the M-estimating equation.
//! @param location_opts a list of options for the location rho-function
//! @return a vector with 2 elements: the location and the scale estimate.
SEXP MLocationScale(SEXP x, SEXP mscale_opts, SEXP location_opts) noexcept;
}  // namespace r_interface
}  // namespace pense

#endif  // R_ROBUST_UTILS_HPP_
