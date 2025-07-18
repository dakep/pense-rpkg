//
//  r_utilities.cc
//  pense
//
//  Created by David Kepplinger on 2019-05-12.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#include "r_utilities.hpp"

#include <cmath>
#include <limits>
#include <iterator>

#include "alias.hpp"
#include "rho.hpp"
#include "r_interface_utils.hpp"
#include "omp_utils.hpp"
#include "robust_scale_location.hpp"

namespace {
template<class T>
using FwdList = pense::alias::FwdList<T>;
using pense::r_interface::MakeVectorView;
using pense::r_interface::MakeUIntVector;
using arma::uword;
using arma::uvec;

//! @brief Transform a single solution into a weight vector.
//! @param solution a solution list (containing residuals and scale)
//! @param rho the rho function defining the weight
//! @param nobs the number of observations to expect
//! @return a weight vector
arma::vec ExtractWeightsVector (const Rcpp::List& solution, const pense::RhoBisquare& rho) {
  const double scale = Rcpp::as<double>(solution["scale"]);
  auto residuals = MakeVectorView(solution["residuals"]);
  const auto wgts = rho.Weight(*residuals, scale);
  return wgts / arma::mean(wgts % arma::square(*residuals / scale));
}

//! @brief Transform a list of solutions into a weight matrix.
//! @param solutions a list of solutions (level1 > solution {residuals, scale})
//! @param rho the rho function defining the weights
//! @param nobs the number of observations to expect
//! @return a weight matrix, with `nobs` rows and `solutions.size()` columns
arma::mat ExtractWeights (const Rcpp::List& solutions, const pense::RhoBisquare& rho, const int nobs) {
  arma::mat wgts(nobs, solutions.size(), arma::fill::none);

  for (int i = 0; i < solutions.size(); ++i) {
    wgts.col(i) = ExtractWeightsVector(solutions[i], rho);
  }

  return wgts;
}

} // namespace

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
SEXP ApproximateMatch(SEXP r_x, SEXP r_table, SEXP r_eps) noexcept {
  const R_xlen_t len_x = Rf_xlength(r_x);
  const int len_table = Rf_length(r_table);
  SEXP r_matches = PROTECT(Rf_allocVector(INTSXP, len_x));
  int* matches = INTEGER(r_matches);
  double const * x = REAL(r_x);
  double const * table = REAL(r_table);
  const double eps = *REAL(r_eps);

  for (R_xlen_t i = 0; i < len_x; ++i) {
    matches[i] = NA_INTEGER;
    for (int j = 0; j < len_table; ++j) {
      if (std::abs(x[i] - table[j]) < eps) {
        matches[i] = j + 1;
        break;
      }
    }
  }

  UNPROTECT(1);
  return r_matches;
}

//! @brief Extract the robustness weights from minima of the robust penalized objective function
//!
//! @param r_solutions a nested list of solutions (lambda > solution)
//! @param r_nobs the total number of observations for estimating the global solutions
//! @param r_cc the cutoff constant used for Tukey's bisquare rho function
//! @return a list of the same length as r_solutions, with each element being a matrix
//!   of weights.
SEXP RobustnessWeight (SEXP r_solutions, SEXP r_nobs, SEXP r_cc) noexcept {
  using namespace pense;
  using Rcpp::List;
  using Rcpp::as;
  using Rcpp::Named;

  BEGIN_RCPP

  const List solutions = as<List>(r_solutions);
  const int nobs = as<int>(r_nobs);
  const double cc = as<double>(r_cc);
  List all_weights(solutions.size());
  RhoBisquare rho(cc);

  // Loop over each lambda.
  for (int lambda_ind = 0; lambda_ind < solutions.size(); ++lambda_ind) {
    all_weights[lambda_ind] = ExtractWeights(solutions[lambda_ind], rho, nobs);
  }

  return Rcpp::wrap(all_weights);

  END_RCPP
}

}  // namespace r_interface
}  // namespace pense
