//
//  r_utilities.hpp
//  pense
//
//  Created by David Kepplinger on 2019-05-12.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef R_UTILITIES_HPP_
#define R_UTILITIES_HPP_

#include "nsoptim_forward.hpp"

namespace pense {
namespace r_interface {
//! Approximate value matching.
//!
//! Returns a vector of 1-based positions of the (first) matches of `x` in `table`.
//!
//! @param x numeric predictor matrix with `n` rows and `p` columns.
//! @param y numeric response vector with `n` elements.
//! @return a vector the same lenght of `x` with integers giving the position in `table` of the first match
//!         if there is a match, or `NA_integer_` otherwise.
SEXP ApproximateMatch(SEXP x, SEXP table, SEXP eps) noexcept;

//! @brief Extract the robustness weights from minima of the robust penalized objective function
//!
//! @param r_solutions a nested list of solutions (lambda > solution)
//! @param r_nobs the total number of observations for estimating the global solutions
//! @param r_cc the cutoff constant used for Tukey's bisquare rho function
//! @return a list of the same length as r_solutions, with each element being a matrix
//!   of weights.
SEXP RobustnessWeight (SEXP r_solutions, SEXP r_nobs, SEXP r_cc) noexcept;

}  // namespace r_interface
}  // namespace pense

#endif  // R_UTILITIES_HPP_
