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

//! @brief Find the best matches between the global solutions and the CV solutions.
//!
//! @param r_solutions_cv a nested list of solutions (cv fold > lambda > solution)
//! @param r_solutions_global a list of solutions (lambda > solution)
//! @param r_cv_k the number of CV folds in each replication. `r_solutions_cv` is expected to
//!   be of length `cv_k * cv_repl`.
//! @param r_nobs the total number of observations for estimating the global solutions
//! @param r_cc the cutoff constant used for Tukey's bisquare rho function
//! @param r_ncores number of cores to use in parallel.
//! @return a nested list of the form lambda > global solution index > {distances, wmspe},
//!   where `distances` is a vector of distances to the CV solutions used for prediction and
//!   `wmspe` the weighted means-squared prediction error from these predictions.
SEXP MatchSolutionsByWeight (SEXP r_solutions_cv, SEXP r_solutions_global, SEXP r_cv_k, SEXP r_nobs,
                             SEXP r_cc, SEXP r_ncores) noexcept;

}  // namespace r_interface
}  // namespace pense

#endif  // R_UTILITIES_HPP_
