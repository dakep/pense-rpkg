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

class BestMatch {
 public:
  BestMatch (const int n_global_sol, const int n_cv_folds) :
    similarity(n_global_sol, n_cv_folds, arma::fill::value(std::numeric_limits<double>::min())),
    sol_ind(n_global_sol, n_cv_folds, arma::fill::zeros) {}

  void Update(const uword row, const uword col, const double similarity_, const uword sol_ind_) {
    similarity.at(row, col) = similarity_;
    sol_ind.at(row, col) = sol_ind_;
  }

  arma::mat similarity;
  arma::umat sol_ind;
};

BestMatch FindBestMatch (const arma::mat& global_wgts, const Rcpp::List& solutions_cv, const int lambda_ind,
                         const pense::RhoBisquare& rho) {
  using Rcpp::List;
  using namespace pense;

  BestMatch best_match(global_wgts.n_cols, solutions_cv.size());

  // In each CV fold, compare the global weights with the weights of the CV solutions for the current lambda
  for (int cv_fold_ind = 0; cv_fold_ind < solutions_cv.size(); ++cv_fold_ind) {
    const List& cv_fold = solutions_cv[cv_fold_ind];
    const List& cv_fold_estimates = cv_fold["estimates"];
    const List& cv_lambda_solutions = cv_fold_estimates[lambda_ind];
    const uvec train_ind = MakeUIntVector(cv_fold["train_ind"]) - 1;

    for (int sol_ind = 0; sol_ind < cv_lambda_solutions.size(); ++sol_ind) {
      const List& sol = cv_lambda_solutions[sol_ind];
      const auto cv_wgts = ExtractWeightsVector(sol, rho);
      for (int global_sol_ind = 0; global_sol_ind < global_wgts.n_cols; ++global_sol_ind) {
        // Compute correlation between the global weights and the CV solution weights
        const double cor = arma::as_scalar(
          arma::cor(global_wgts.unsafe_col(global_sol_ind).elem(train_ind), cv_wgts));

        if (cor > best_match.similarity(global_sol_ind, cv_fold_ind)) {
          best_match.Update(global_sol_ind, cv_fold_ind, cor, sol_ind);
        }
      }
    }
  }
  return best_match;
};

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
                             SEXP r_cc) noexcept {
  using namespace pense;
  using Rcpp::List;
  using Rcpp::as;
  using Rcpp::Named;

  BEGIN_RCPP

  const List solutions_cv = as<List>(r_solutions_cv);
  const List solutions_global = as<List>(r_solutions_global);
  const int cv_k = as<int>(r_cv_k);
  const int nobs = as<int>(r_nobs);
  const double cc = as<double>(r_cc);
  const int cv_repl = solutions_cv.size() / cv_k;
  List results(solutions_global.size());
  RhoBisquare rho(cc);

  // Loop over each lambda.
  for (int lambda_ind = 0; lambda_ind < solutions_global.size(); ++lambda_ind) {
    // For this lambda, extract the weights of all global solutions,
    // and determine the solution with smallest distance for each global solution in each fold.
    const auto global_wgts = ExtractWeights(solutions_global[lambda_ind], rho, nobs);
    const auto best_match = FindBestMatch(global_wgts, solutions_cv, lambda_ind, rho);

    // Collect the correlations and weighted mean squared prediction errors
    // for all global solutions at the current lambda index.
    List lambda_result(global_wgts.n_cols);

    for (int global_sol_ind = 0; global_sol_ind < global_wgts.n_cols; ++global_sol_ind) {
      arma::vec pred_wmse(cv_repl, arma::fill::zeros);
      arma::vec pred_tau_size(cv_repl, arma::fill::zeros);
      arma::mat similarities(cv_k, cv_repl, arma::fill::none);

      int cv_repl_ind = -1;
      arma::vec all_test_resids(nobs, arma::fill::none);
      int insert_index = 0;

      for (int cv_fold_ind = 0; cv_fold_ind < solutions_cv.size(); ++cv_fold_ind) {
        const List& cv_fold = solutions_cv[cv_fold_ind];
        const List& cv_fold_solutions = cv_fold["estimates"];
        const List& cv_lambda_solutions = cv_fold_solutions[lambda_ind];
        const List& selected_sol = cv_lambda_solutions[best_match.sol_ind(global_sol_ind, cv_fold_ind)];

        auto test_residuals = MakeVectorView(selected_sol["test_residuals"]);
        const arma::uvec test_ind = MakeUIntVector(cv_fold["test_ind"]) - 1;

        similarities(cv_fold_ind) = best_match.similarity(global_sol_ind, cv_fold_ind);

        if (cv_fold_ind % cv_k == 0) {
          if (cv_repl_ind >= 0) {
            pred_tau_size(cv_repl_ind) = pense::TauSize(all_test_resids);
            insert_index = 0;
          }
          ++cv_repl_ind;
        }
        // Divide each chunk by the sum of the weights to get the overall weighted mean in the end
        pred_wmse(cv_repl_ind) += arma::dot(global_wgts.unsafe_col(global_sol_ind).elem(test_ind),
                                            arma::square(*test_residuals)) / nobs;

        // Copy the test residuals from each chunk to compute the tau-size afterwards.
        const int upper_index = insert_index + test_residuals->n_elem - 1;
        all_test_resids.subvec(insert_index, upper_index) = *test_residuals;
        insert_index += test_residuals->n_elem;
      }

      // Process the tau-size in the last CV replication
      pred_tau_size(cv_repl_ind) = pense::TauSize(all_test_resids);

      lambda_result[global_sol_ind] = List::create(Named("rankcorr") = similarities,
                                                   Named("wmspe") = pred_wmse,
                                                   Named("tau_size") = pred_tau_size);
    }

    results[lambda_ind] = lambda_result;
  }

  return Rcpp::wrap(results);

  END_RCPP
}

}  // namespace r_interface
}  // namespace pense
