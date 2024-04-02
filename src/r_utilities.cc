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

//! @brief Sort a vector of **1-based** indices according to one or more vectors of values of equal or greater length.
class IndexSort {
 public:
  //! @brief Initialize the index sorter for the given indices.
  //! @param indices a vector of (0-based!) indices.
  explicit IndexSort(const uvec indices) : ind_(indices) {}

  //! @brief Initialize the index sorter for a regular index from 0 to the given length.
  //! @param length length of the indices.
  explicit IndexSort(const uword length)
    : ind_(arma::regspace<uvec>(0, length - 1)) {}

  //! @brief Initialize the index sorter as a copy of the given object.
  //! @param other an instance of IndexSort of which the indices should be copied.
  explicit IndexSort(const IndexSort& other) : ind_(other.ind_) {}

  //! @brief Sort the 1-based indices according to the given values.
  //! @param values vector of values. Must be at least as long as the largest index.
  //! @param max_swaps maximum number of swaps allowed. If < 0, no maximum is enforced.
  //! @return number of swaps performed.
  int operator()(const arma::vec& values, const int max_swaps = -1) {
    values_ = &values;
    max_swaps_ = max_swaps;
    auto buf = ind_;
    int swaps = SplitMerge(0, ind_.n_elem, &buf, &ind_);
    return swaps;
  }

  //! @brief Access the indices, sorted according to the last sorting operation.
  //! @return sorted indices.
  const uvec& SortedIndices() const {
    return ind_;
  }

 private:
  int max_swaps_;
  uvec ind_;
  arma::vec const * values_;

  int Merge (const uvec& sorted_halves, const int start, const int middle, const int end, uvec* output) {
    int i = start;
    int j = middle;
    int swaps = 0;

    for (int k = start; k < end; ++k) {
      if (i < middle && (j >= end || (*values_)(sorted_halves(i)) <= (*values_)(sorted_halves(j)))) {
        (*output)(k) = sorted_halves(i);
        ++i;
      } else {
        (*output)(k) = sorted_halves(j);
        swaps += middle - i;
        ++j;
      }
    }
    return swaps;
  }

  int SplitMerge (const int start, const int end, uvec* from, uvec* to) {
    if (end - start <= 1) {
      return 0;
    }

    int middle = (start + end) / 2;
    int swaps = SplitMerge(start, middle, to, from);

    // Check if we should continue with the other half or if we already reached the maximum number of swaps.
    if (max_swaps_ >= 0 && swaps > max_swaps_) {
      return swaps;
    }

    swaps += SplitMerge(middle, end, to, from);

    // Check if we should continue merging the two halves or if we already reached the maximum number of swaps.
    if (max_swaps_ >= 0 && swaps > max_swaps_) {
      return swaps;
    }

    swaps += Merge(*from, start, middle, end, to);

    return swaps;
  }
};


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

struct DuplicateCount {
  uword in_x = 0;
  uword in_y = 0;
  uword in_both = 0;
};

///! Count the number of zeros in x and y, with x being restricted to a subset of indices.
DuplicateCount DuplicateZeros (const arma::vec& x, const arma::vec& y, const uvec& x_subset) {
  DuplicateCount cnt;
  for (uword i = 0; i < x_subset.n_elem; ++i) {
    const bool x_zero = x[x_subset[i]] == 0;
    const bool y_zero = y[i] == 0;
    if (x_zero && y_zero) {
      ++cnt.in_both;
    }
    if (x_zero) {
      ++cnt.in_x;
    }
    if (y_zero) {
      ++cnt.in_y;
    }
  }

  return cnt;
}

class BestMatch {
 public:
  BestMatch (const int n_global_sol, const int n_cv_folds) :
    swaps(n_global_sol, n_cv_folds, arma::fill::value(std::numeric_limits<uword>::max() / 2)),
    kendall_tau(n_global_sol, n_cv_folds, arma::fill::zeros),
    sol_ind(n_global_sol, n_cv_folds, arma::fill::zeros) {}

  void Update(const uword row, const uword col, const uword swaps_,
              const double kendall_tau_, const uword sol_ind_) {
    swaps.at(row, col) = swaps_;
    kendall_tau.at(row, col) = kendall_tau_;
    sol_ind.at(row, col) = sol_ind_;
  }

  arma::umat swaps;
  arma::mat kendall_tau;
  arma::umat sol_ind;
};

void FindBestMatchesForFold (const arma::uword cv_fold_ind,
                             const arma::uword cv_sol_ind,
                             const arma::mat& global_wgts,
                             const arma::uvec& sorted_train_ind,
                             const arma::uvec& unsorted_train_ind,
                             const arma::vec& cv_sol_weights,
                             BestMatch* best_match) {
  for (int global_sol_ind = 0; global_sol_ind < global_wgts.n_cols; ++global_sol_ind) {
    // Now sort the training index according to the global weights.
    // The number of swaps is the number of discordant pairs between the CV weights and the global weights.
    IndexSort sol_sorter(sorted_train_ind);

    const auto zeros_dupl = DuplicateZeros(global_wgts.col(global_sol_ind), cv_sol_weights, unsorted_train_ind);
    const uword numerator_adjustment = zeros_dupl.in_x * (zeros_dupl.in_x - 1) / 2 -
      zeros_dupl.in_both * (zeros_dupl.in_both - 1) / 2 +
      zeros_dupl.in_y * (zeros_dupl.in_y - 1) / 2;

    const uword n_pairs = unsorted_train_ind.n_elem * (unsorted_train_ind.n_elem - 1) / 2;
    // Allow up to 2x as many swaps as for the maximum kendall's tau to
    // ensure the numerator adjustment for the ties doesn't make the sorting
    // stop too early.
    const uword max_swaps = 2 * best_match->swaps(global_sol_ind, cv_fold_ind);
    const uword swaps = sol_sorter(global_wgts.col(global_sol_ind), max_swaps);
    if (swaps < max_swaps) {
      const uword kendall_num = numerator_adjustment + 2 * swaps;
      const double kendall_denom = std::sqrt(n_pairs - zeros_dupl.in_x * (zeros_dupl.in_x - 1) / 2) *
        std::sqrt(n_pairs - zeros_dupl.in_y * (zeros_dupl.in_y - 1) / 2);
      const double kendall_tau = (n_pairs > kendall_num) ? (n_pairs - kendall_num) / kendall_denom : 0;

      if (kendall_tau > best_match->kendall_tau(global_sol_ind, cv_fold_ind)) {
        best_match->Update(global_sol_ind, cv_fold_ind, swaps, kendall_tau, cv_sol_ind);
      }
    }
  }
}

BestMatch FindBestMatchMT (const arma::mat& global_wgts, const Rcpp::List& solutions_cv, const int lambda_ind,
                           const pense::RhoBisquare& rho, const int num_threads) {
  using Rcpp::List;
  using namespace pense;
  using pense::alias::FwdList;

  const int n_cv_folds = solutions_cv.size();
  BestMatch best_match(global_wgts.n_cols, n_cv_folds);

  // Rcpp::List is not thread-safe. Hence we extract the necessary data before spawning multiple threads.
  FwdList<uvec> cv_train_indices;
  FwdList<arma::mat> cv_sol_weights;

  auto cv_train_indices_it = cv_train_indices.before_begin();
  auto cv_sol_weights_it = cv_sol_weights.before_begin();
  for (int cv_fold_ind = 0; cv_fold_ind < n_cv_folds; ++cv_fold_ind) {
    const List& cv_fold = solutions_cv[cv_fold_ind];
    const List& cv_fold_estimates = cv_fold["estimates"];
    const List& cv_lambda_solutions = cv_fold_estimates[lambda_ind];

    cv_train_indices_it = cv_train_indices.insert_after(cv_train_indices_it, MakeUIntVector(cv_fold["train_ind"]) - 1);
    cv_sol_weights_it = cv_sol_weights.insert_after(cv_sol_weights_it, ExtractWeights(cv_lambda_solutions, rho,
                                                                                      cv_train_indices_it->n_elem));
  }

  // In each CV fold, compare the global weights with the weights of the CV solutions for the current lambda
  cv_train_indices_it = cv_train_indices.begin();
  cv_sol_weights_it = cv_sol_weights.begin();

  #pragma omp parallel \
              num_threads(num_threads) \
              default(none) \
              shared(best_match, cv_sol_weights_it, cv_train_indices_it) \
              const_local_shared(global_wgts, n_cv_folds, lambda_ind)
  #pragma omp single
  for (int cv_fold_ind = 0; cv_fold_ind < n_cv_folds; ++cv_fold_ind) {
    const auto cv_sol_weights_fold_it = cv_sol_weights_it++;
    const auto cv_train_indices_fold_it = cv_train_indices_it++;

    #pragma omp task default(none) \
      firstprivate(cv_fold_ind, cv_sol_weights_fold_it, cv_train_indices_fold_it) \
      shared(best_match) \
      const_local_shared(global_wgts, lambda_ind)
    for (int sol_ind = 0; sol_ind < cv_sol_weights_fold_it->n_cols; ++sol_ind) {
      // First sort the indices according to the CV solution weights ...
      IndexSort sol_sorter(cv_sol_weights_fold_it->n_rows);
      sol_sorter(cv_sol_weights_fold_it->col(sol_ind));

      // Then sort the actual training indices based on the sorted index of the weights (-1 to account
      // for R's 1-based index)
      const uvec cv_sol_sorted_train_ind = cv_train_indices_fold_it->elem(sol_sorter.SortedIndices());
      FindBestMatchesForFold(cv_fold_ind, sol_ind, global_wgts, cv_sol_sorted_train_ind, *cv_train_indices_fold_it,
                             cv_sol_weights_fold_it->col(sol_ind), &best_match);
    }
  }

  return best_match;
}

BestMatch FindBestMatchST (const arma::mat& global_wgts, const Rcpp::List& solutions_cv, const int lambda_ind,
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
      IndexSort sol_sorter(train_ind.n_elem);
      // First sort the indices according to the CV solution weights ...
      const auto cv_wgts = ExtractWeightsVector(sol, rho);
      sol_sorter(cv_wgts);
      // Then sort the actual training indices based on the sorted index of the weights (-1 to account
      // for R's 1-based index)
      const uvec cv_sol_sorted_train_ind = train_ind.elem(sol_sorter.SortedIndices());
      FindBestMatchesForFold(cv_fold_ind, sol_ind, global_wgts, cv_sol_sorted_train_ind, train_ind, cv_wgts,
                             &best_match);
    }
  }
  return best_match;
}

BestMatch FindBestMatch (const arma::mat& global_wgts, const Rcpp::List& solutions_cv, const int lambda_ind,
                         const pense::RhoBisquare& rho, const int num_threads) {
  if (pense::omp::Enabled(num_threads)) {
    return FindBestMatchMT(global_wgts, solutions_cv, lambda_ind, rho, num_threads);
  } else {
    return FindBestMatchST(global_wgts, solutions_cv, lambda_ind, rho);
  }
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
                             SEXP r_cc, SEXP r_ncores) noexcept {
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
  const double ncores = as<int>(r_ncores);
  const int cv_repl = solutions_cv.size() / cv_k;
  List results(solutions_global.size());
  RhoBisquare rho(cc);

  // Loop over each lambda.
  for (int lambda_ind = 0; lambda_ind < solutions_global.size(); ++lambda_ind) {
    // For this lambda, extract the weights of all global solutions,
    // and determine the solution with smallest distance for each global solution in each fold.
    const auto global_wgts = ExtractWeights(solutions_global[lambda_ind], rho, nobs);
    const auto best_match = FindBestMatch(global_wgts, solutions_cv, lambda_ind, rho, ncores);

    // Collect the correlations and weighted mean squared prediction errors
    // for all global solutions at the current lambda index.
    List lambda_result(global_wgts.n_cols);

    for (int global_sol_ind = 0; global_sol_ind < global_wgts.n_cols; ++global_sol_ind) {
      arma::vec pred_wmse(cv_repl, arma::fill::zeros);
      arma::vec pred_tau_size(cv_repl, arma::fill::zeros);
      arma::mat kendall_taus(cv_k, cv_repl, arma::fill::none);
      const double wgt_sum = arma::accu(global_wgts.col(global_sol_ind));

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

        kendall_taus(cv_fold_ind) = best_match.kendall_tau(global_sol_ind, cv_fold_ind);

        if (cv_fold_ind % cv_k == 0) {
          if (cv_repl_ind >= 0) {
            pred_tau_size(cv_repl_ind) = pense::TauSize(all_test_resids);
            insert_index = 0;
          }
          ++cv_repl_ind;
        }
        // Divide each chunk by the sum of the weights to get the overall weighted mean in the end
        pred_wmse(cv_repl_ind) += arma::dot(global_wgts.unsafe_col(global_sol_ind).elem(test_ind),
                                            arma::square(*test_residuals)) / wgt_sum;

        // Copy the test residuals from each chunk to compute the tau-size afterwards.
        const int upper_index = insert_index + test_residuals->n_elem - 1;
        all_test_resids.subvec(insert_index, upper_index) = *test_residuals;
        insert_index += test_residuals->n_elem;
      }

      // Process the tau-size in the last CV replication
      pred_tau_size(cv_repl_ind) = pense::TauSize(all_test_resids);

      lambda_result[global_sol_ind] = List::create(Named("rankcorr") = kendall_taus,
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
