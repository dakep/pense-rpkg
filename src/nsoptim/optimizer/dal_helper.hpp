//
//  dal_helper.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-11-30.
//  Copyright © 2018 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_OPTIMIZER_DAL_HELPER_HPP_
#define NSOPTIM_OPTIMIZER_DAL_HELPER_HPP_

#include <type_traits>
#include <memory>
#include <algorithm>
#include <limits>

#include "../armadillo.hpp"
#include "../traits/traits.hpp"
#include "../objective/ls_regression_loss.hpp"
#include "../container/metrics.hpp"
#include "../container/data.hpp"

namespace nsoptim {
namespace _optim_dal_internal {
constexpr double kPhiStepDirMinEps = 1e-8;  //!< Minimum convergence tolerance for PCG

//! Constants defining the decrease in the step size in each line search iteration
//! and the fraction of the decrease in the objective function predicted by
//! linear extrapolation that we will accept in the line search
constexpr double kLinesearchStepsizeMultiplier = 0.8;  //!< 0 < x < 1
constexpr double kLinesearchStepCondition = 0.3;       //!< 0 < x < 0.5
constexpr int kLinesearchMaxSteps = 50;

//! Maximum relative change in weights tolerated without recomputing the preconditioner.
//! A smaller number results in more "fresh starts", i.e., matrix inversions,
//! while a larger number results in more PCG iterations.
constexpr double kMaximumRelativeWeightChange = 1e-3;

//! If the norm of the gradient increases by this much, re-compute the exact inverse of the Hessian.
constexpr double kMaximumGradientNormIncrease = 1.1;

//! Possible negative return values for `Hessian::SolveFor`
constexpr int kPhiStepDirFullInversion = -1;  //< Computed the inverse of the Hessian. No iterations required.
constexpr int kPhiStepDirInversionFailed = -2;  //< Could not compute the inverse of the Hessian.

//! Possible negative return values for `DalEnOptimizer::MinimizePhi`
constexpr int kMinimizePhiGradientTooSteep = -1;  //< Could not minimize phi because the gradient is too steep!
constexpr int kMinimizePhiHessianSingular = -2;  //< Could not compute the inverse of the Hessian.

constexpr double kMaximumEtaStart = 1e6;  //< Maximum start value for eta.

//! POD structure for the proximity and regularization parameters.
struct ProximityParameters {
  double nxlambda = -1;
  double slope;
  double intercept;
};

//! Enumeration of how data changed.
struct DataChanges {
  bool data_changed;
  int weights_changed;  //!< 0 ... no change, 1 ... small change, 2 ... large change
};

//! Compute the dual of the LS loss.
//!
//! @param a (negative of the) dual vector.
//! @param y response vector.
//! @return the value of the dual loss.
inline double DualLoss(const arma::vec& a, const arma::vec& y);

//! Compute the solution `x` to the linear system of equations `A . x = b` using the preconditioned conjuage gradient
//! method and the given matrix as preconditioner.
//!
//! @param A the matrix `A` in the equation.
//! @param b the vector `b` in the equation.
//! @param precond the preconditioner matrix.
//! @param eps the convergence threshold.
//! @param x a pointer to the solution vector (`x` in the equation). Is used as the starting point for the PCG method.
//! @return the number of iterations or `kPhiStepDirInversionFailed` if not converged.
inline int SolvePcg(const arma::mat &A, const arma::vec &b, const arma::mat& precond,
                    const double eps, arma::vec* x) noexcept;

//! Proxy class to manage the actual data.
//! The default implementation simply proxies to the data used by the loss function.
template<typename LossFunction, typename = typename traits::is_weighted<LossFunction>::type>
class DataProxy {
 public:
  DataProxy() = default;
  // Create a new DataProxy for the given loss function.
  explicit DataProxy(LossFunction const * const loss) noexcept : data_(loss ? &(loss->data()) : nullptr) {}

  //! Copying is allowed.
  DataProxy(const DataProxy& other) = default;
  DataProxy& operator=(const DataProxy& other) = default;
  //! Moving is allowed.
  DataProxy(DataProxy&& other) = default;
  DataProxy& operator=(DataProxy&& other) = default;

  //! Update the data proxy with a new loss function.
  //!
  //! @param the new loss function.
  //! @return information on what data changed.
  DataChanges Update(const LossFunction& loss) noexcept {
    // if (data_ != &loss.data() ||
    //     (data_ && (data_->n_obs() != loss.data().n_obs() || data_->n_pred() != loss.data().n_pred()))) {
    //   data_ = &loss.data();
    //   return DataChanges {true, 0};
    // }
    data_ = &loss.data();
    return DataChanges {true, 0};
  }

  //! Access the data set.
  //!
  //! @return a reference to the data set.
  const PredictorResponseData& operator*() const {
    return *data_;
  }

  //! Get the pointer to the data set.
  //!
  //! @return a reference to the data set.
  PredictorResponseData const * get() const {
    return data_;
  }

  //! Get the average weight.
  //!
  //! @return average weight.
  constexpr double mean_weight() const noexcept {
    return 1.;
  }

  //! Access the data set.
  //!
  //! @return a pointer to the data set.
  PredictorResponseData const * operator->() const {
    return data_;
  }

 private:
  PredictorResponseData const * data_ = nullptr;
};

//! Proxy class to manage the actual data.
//! This implementation for the weighted LS loss proxies to the weighted data.
template<typename LossFunction>
class DataProxy<LossFunction, std::true_type> {
 public:
  DataProxy() = default;

  //! Important: the DataProxy only retains a reference to the loss' data and weights. Therefore, it is the client's
  //! responsibility to ensure that the passed-in loss is available until *after* the `Update` method is called with
  //! a new loss!
  explicit DataProxy(LossFunction const * const loss)
    : data_(loss ? &(loss->data()) : nullptr), sqrt_weights_(loss ? &(loss->sqrt_weights()) : nullptr),
      mean_weight_(loss ? loss->mean_weight() : 1.),
      weighted_data_(loss ?
                     PredictorResponseData(data_->cx().each_col() % *sqrt_weights_, data_->cy() % *sqrt_weights_) :
                     PredictorResponseData()) {}

  //! Copying is allowed.
  DataProxy(const DataProxy& other) = default;
  DataProxy& operator=(const DataProxy& other) = default;
  //! Moving is allowed.
  DataProxy(DataProxy&& other) = default;
  DataProxy& operator=(DataProxy&& other) = default;

  //! Update the data proxy with a new loss function.
  //! Important: the DataProxy only retains a reference to the loss' data and weights. Therefore, it is the client's
  //! responsibility to ensure that the passed-in loss is available until *after* the `Update` method is called with
  //! a new loss!
  //!
  //! @param the new loss function.
  //! @return information on what data changed.
  DataChanges Update(const LossFunction& loss) {
    // Check if it's actually new data and/or new weights
    // DataChanges changes {data_ != &loss.data() || (data_ && (data_->n_obs() != loss.data().n_obs() ||
    //                                                          data_->n_pred() != loss.data().n_pred())), 0};
    // if (sqrt_weights_ != &loss.sqrt_weights()) {
    //   if (sqrt_weights_ && (sqrt_weights_->n_elem == loss.sqrt_weights().n_elem)) {
    //     // Check by how much the weights changed. If the change is small enough, the caller does not need to know.
    //     const double diff_norm = arma::norm(loss.sqrt_weights() - *sqrt_weights_);

    //     if (diff_norm * diff_norm < sqrt_weights_->n_elem * _optim_dal_internal::kMaximumRelativeWeightChange) {
    //       changes.weights_changed = 1;
    //     } else {
    //       changes.weights_changed = 2;
    //     }
    //     sqrt_weights_ = &loss.sqrt_weights();
    //   } else {
    //     // The dimensions of the weights vector changed.
    //     changes.weights_changed = 2;
    //   }
    // }
    DataChanges changes {true, 2};

    if (changes.data_changed || changes.weights_changed) {
      sqrt_weights_ = &loss.sqrt_weights();
      data_ = &loss.data();
      mean_weight_ = loss.mean_weight();
      weighted_data_ = PredictorResponseData(data_->cx().each_col() % *sqrt_weights_, data_->cy() % *sqrt_weights_);
    }
    return changes;
  }

  //! Get the square-root of the observation weights.
  //!
  //! @return vector with the sqrt of the observation weights.
  const arma::vec& sqrt_weights() const noexcept {
    return *sqrt_weights_;
  }

  //! Get the average weight.
  //!
  //! @return average weight.
  double mean_weight() const noexcept {
    return mean_weight_;
  }

  //! Access the weighted data set.
  //!
  //! @return a reference to the weighted data set.
  const PredictorResponseData& operator*() const noexcept {
    return weighted_data_;
  }

  //! Get the pointer to the data set.
  //!
  //! @return a reference to the data set.
  PredictorResponseData const * get() const {
    return &weighted_data_;
  }

  //! Access the weighted data set.
  //!
  //! @return a pointer to the weighted data set.
  PredictorResponseData const * operator->() const noexcept {
    return &weighted_data_;
  }

 private:
  PredictorResponseData const * data_ = nullptr;
  arma::vec const * sqrt_weights_ = nullptr;
  double mean_weight_;
  PredictorResponseData weighted_data_;
};


//! Proxy class to manage the Hessian of the phi function.
//!
//! Storage requirements are $2 n_obs^2 + n_pred + n_obs^2 if weighted$
template<class LossFunction>
class Hessian {
  using WeightsTag = typename traits::is_weighted<LossFunction>::type;
  using WeightsOuter = typename std::conditional<WeightsTag::value, arma::mat, int>::type;

 public:
  Hessian() = default;
  explicit Hessian(const DataProxy<LossFunction>& data) noexcept : Hessian(data, WeightsTag{}) {}

  // //! Copy the other Hessian, but replace the data pointer.
  Hessian(const Hessian& other, const DataProxy<LossFunction>& data) noexcept :
    data_(data.get()), pcg_convergence_tol_(other.pcg_convergence_tol_),
    sqrt_weights_outer_(other.sqrt_weights_outer_), outer_product_(other.outer_product_),
    preconditioner_(other.preconditioner_), active_(other.active_) {}

  //! Copying is allowed.
  Hessian(const Hessian& other) = default;
  Hessian& operator=(const Hessian& other) = default;
  //! Moving is allowed.
  Hessian(Hessian&& other) = default;
  Hessian& operator=(Hessian&& other) = default;

  //! Set the convergence tolerance for the PCG algorithm.
  //! This should be smaller than the desired convergence tolerance for the main DAL algorithm!
  void convergence_tolerance(double convergence_tolerance) noexcept {
    if (convergence_tolerance < kPhiStepDirMinEps) {
      pcg_convergence_tol_ = convergence_tolerance;
    }
  }

  //! Update the data pointer, without changing the state of the Hessian.
  //! This is useful if the data changed slightly, but the preconditioner and the active set
  //! remain the same.
  void UpdateData(const DataProxy<LossFunction>& data) {
    data_ = &(*data);
    UpdateWeights(data, WeightsTag{});
    if (active_.n_elem != 0) {
      auto active_view = data_->cx().cols(active_);
      outer_product_ = active_view * active_view.t();
    }
  }

  //! Solve the equation
  //! H . s = x
  //! for the step direction `s` given the matrix
  //! `H = eta.slope * moreau_factor * (X_a . X_a') + eta.intercept * w . w'` and the gradient `x`.
  //! using preconditioned conjugate gradient if possible.
  //!
  //! @param beta the current slope coefficients.
  //! @param gradient the gradient at the current coefficients value.
  //! @param moreau_factor the Moreau factor in the current iteration.
  //! @param state the current state of the DAL algorithm.
  //! @param include_intercept if the intercept should be included or not.
  //! @param convergence_tolerance the convergence tolerance to solve the equation.
  //! @param metric pointer to the Metrics to store metrics about the solution.
  //! @param step_dir pointer to the vector where the step direction will be stored at.
  //! @return if >= 0, the number of iterations required to compute the step direction. If < 0, the following status:
  //!         `kPhiStepDirFullInversion` ... fully inverted the Hessian.
  //!         `kPhiStepDirInversionFailed` ... failed to solve the system of equations.
  int SolveFast(const arma::sp_vec& beta, const arma::vec& gradient, const double moreau_factor,
                const ProximityParameters& eta, const bool include_intercept, Metrics* metric, arma::vec* step_dir) {
    if (beta.n_nonzero == 0) {
      // No predictors are active. We can compute the inverse of the Hessian and the step direction directly.
      if (include_intercept) {
        *step_dir = StepDirNoPredictors(gradient, eta, WeightsTag{});
        return 0;
      } else {
        // If there's no intercept, compute the step direction directly. Don't touch the internal state at all.
        *step_dir = gradient;
      }
    }

    // sp_vec does not have an overload for the find method, so it has to be done in a hacky way.
    // Instead of the easy way `const uvec active = find(beta);` we create a view of the row-indices.
    beta.sync();
    const arma::uvec active(const_cast<arma::uword*>(beta.row_indices), beta.n_nonzero, false, true);
    bool solve_directly = false;

    // Compute/update the outer product.
    if (active_.n_elem == 0) {
      // The current outer product is empty, but there are new active predictors. Compute from scratch.
      active_ = active;
      auto active_view = data_->cx().cols(active);
      outer_product_ = active_view * active_view.t();
      if (active.n_elem < data_->n_pred()) {
        preconditioner_.eye(arma::size(outer_product_));
        // If not all predictors are activated at once, perform rank-1 updates sequentially.
        for (auto new_active_it = active.cbegin(), new_active_end = active.cend(); new_active_it != new_active_end;
             ++new_active_it) {
          auto var_view = data_->cx().col(*new_active_it);
          arma::vec tmp = preconditioner_ * var_view;
          preconditioner_ -= tmp * tmp.t() / (1 + arma::as_scalar(var_view.t() * tmp));
        }
        metric->AddDetail("step_dir_rank1_updates", static_cast<int>(active.n_elem));
      } else {
        // All predictors are activated at once. We will compute the preconditioner by inverting the hessian below.
        solve_directly = true;
      }
    } else {
      // There is a change in active variables. Perform rank-1 updates to the outer product and the preconditioner.
      int rank1_updates = 0;
      const auto old_active_end = active_.end();
      auto old_active_it = active_.begin();
      const auto new_active_end = active.end();
      auto new_active_it = active.begin();
      while (new_active_it != new_active_end || old_active_it != old_active_end) {
        if ((old_active_it != old_active_end) && (new_active_it == new_active_end || *new_active_it > *old_active_it)) {
          // All remaining old active predictors must be removed or a predictor is removed from the active set.
          auto var_view = data_->cx().col(*old_active_it);
          outer_product_ -= var_view * var_view.t();

          // Update the preconditioner according to the Sherman–Morrison-Woodbury formula
          arma::vec tmp = preconditioner_ * var_view;
          preconditioner_ += tmp * tmp.t() / arma::as_scalar(1 - var_view.t() * tmp);
          ++old_active_it;
          ++rank1_updates;
        } else if (old_active_it == old_active_end || *new_active_it < *old_active_it) {
          // All remaining new active predictors must be added or a predictor is added to the new active set.
          auto var_view = data_->cx().col(*new_active_it);
          outer_product_ += var_view * var_view.t();

          // Update the preconditioner according to the Sherman–Morrison-Woodbury formula
          arma::vec tmp = preconditioner_ * var_view;
          preconditioner_ -= tmp * tmp.t() / arma::as_scalar(1 + var_view.t() * tmp);
          ++new_active_it;
          ++rank1_updates;
        } else {
          // The predictor is already in the old active set.
          if (old_active_it != old_active_end) {
            ++old_active_it;
          }
          if (new_active_it != new_active_end) {
            ++new_active_it;
          }
        }
      }
      active_ = active;
      metric->AddDetail("step_dir_rank1_updates", rank1_updates);
    }

    // Compute the Hessian.
    arma::mat hessian;
    if (include_intercept) {
      hessian = eta.slope * moreau_factor * outer_product_ + eta.intercept * sqrt_weights_outer_;
    } else {
      hessian = eta.slope * moreau_factor * outer_product_;
    }
    hessian.diag() += 1;  //< this is usuall faster than adding a diagonal matrix.

    if (solve_directly) {
      // If we solve the equation directly, invert the Hessian and store the result as the future preconditioner.
      const bool success = arma::inv_sympd(preconditioner_, hessian);
      metric->AddDetail("step_dir_invert_hessian", 1);
      *step_dir = preconditioner_ * gradient;
      return success ? kPhiStepDirFullInversion : kPhiStepDirInversionFailed;
    }

    // Use PCG with the current preconditioner to solve for the step direction.
    if (step_dir->n_elem == 0) {
      step_dir->zeros(gradient.n_elem);
    }
    const auto pcg_iters = SolvePcg(hessian, gradient, preconditioner_, pcg_convergence_tol_, step_dir);
    const int max_tolerated_pcg_iters = static_cast<int>(std::sqrt(static_cast<double>(gradient.n_elem))) + 1;
    if (pcg_iters < 0 || pcg_iters > max_tolerated_pcg_iters) {
      // PCG did not converge. Recompute the preconditiioner.
      metric->AddDetail("step_dir_invert_hessian", 1);
      metric->AddDetail("step_dir_pcg_iter", (pcg_iters < 0) ? static_cast<int>(gradient.n_elem) : pcg_iters);
      if (!arma::inv_sympd(preconditioner_, hessian)) {
        return kPhiStepDirInversionFailed;
      }
      *step_dir = preconditioner_ * gradient;
      return kPhiStepDirFullInversion;
    } else {
      metric->AddDetail("step_dir_pcg_iter", static_cast<int>(pcg_iters));
    }
    return pcg_iters;
  }

 private:
  Hessian(const DataProxy<LossFunction>& data, std::true_type) noexcept
      : data_(&(*data)), sqrt_weights_outer_(data.sqrt_weights() * data.sqrt_weights().t()) {}

  Hessian(const DataProxy<LossFunction>& data, std::false_type) noexcept
      : data_(&(*data)), sqrt_weights_outer_(1) {}

  //! Update the outer-product of the sqrt-weights.
  void UpdateWeights(const DataProxy<LossFunction>& data, std::true_type) noexcept {
    sqrt_weights_outer_ = data.sqrt_weights() * data.sqrt_weights().t();
  }

  //! Update the outer-product of the sqrt-weights.
  void UpdateWeights(const DataProxy<LossFunction>&, std::false_type) const noexcept {}

  //! Solve the equation as in `SolveFor` if no predictors are active and a weighted LS loss.
  arma::vec StepDirNoPredictors(const arma::vec& gradient, const ProximityParameters& eta, std::true_type) noexcept {
    // If there's an intercept and weights, the preconditioner is set to the inverse of the Hessian.
    preconditioner_ = sqrt_weights_outer_ * (-eta.intercept / (1 + eta.intercept * arma::trace(sqrt_weights_outer_)));
    preconditioner_.diag() += 1.;
    // Clean the active coefficients and the outer product to reflect that.
    active_.reset();
    // The preconditioner is now the actual inverse itself and can be used directly to compute the step direction.
    return preconditioner_ * gradient;
  }

  //! Solve the equation as in `SolveFor` if no predictors are active and an un-weighted LS loss.
  arma::vec StepDirNoPredictors(const arma::vec& gradient, const ProximityParameters& eta,
                                std::false_type) const noexcept {
    return gradient / (1 + eta.intercept);
  }

  PredictorResponseData const * data_;
  double pcg_convergence_tol_ = kPhiStepDirMinEps;
  WeightsOuter sqrt_weights_outer_;
  arma::mat outer_product_;
  arma::mat preconditioner_;
  arma::uvec active_;
};

//! Compute the dual of the LS loss.
//!
//! @param a (negative of the) dual vector.
//! @param y response vector.
//! @return the value of the dual loss.
inline double DualLoss(const arma::vec& a, const arma::vec& y) {
  return 0.5 * arma::dot(a, a) - arma::dot(a, y);
}

//! Compute the solution `x` to the linear system of equations `A . x = b` using the preconditioned conjuage gradient
//! method and the given matrix as preconditioner.
//!
//! @param A the matrix `A` in the equation.
//! @param b the vector `b` in the equation.
//! @param precond the preconditioner matrix.
//! @param eps the convergence threshold.
//! @param x a pointer to the solution vector (`x` in the equation). Is used as the starting point for the PCG method.
//! @return the number of iterations or `kPhiStepDirInversionFailed` if not converged.
inline int SolvePcg(const arma::mat &A, const arma::vec &b, const arma::mat& precond, const double eps,
                    arma::vec* x) noexcept {
  using arma::norm;
  using arma::dot;

  arma::uword iter = 0;
  auto normb = norm(b);
  arma::vec r = b - A * (*x);

  if (normb < std::numeric_limits<double>::epsilon()) {
    normb = 1;
  }

  const double thresh = eps * normb;

  if (norm(r) <= thresh) {
    return 0;
  }

  arma::vec p = precond * r;
  double abs = arma::dot(r, p);

  do {
    const arma::vec tmp = A * p;
    const auto alpha = abs / arma::dot(p, tmp);

    *x += alpha * p;
    const arma::vec new_r = r - alpha * tmp;

    if (norm(new_r) <= thresh) {
      return static_cast<int>(iter);
    }

    const arma::vec z = precond * new_r;
    const auto beta = dot(z, new_r - r) / abs;
    abs = dot(new_r, z);
    p = z + beta * p;
    r = new_r;
  } while (++iter <= b.n_elem);

  return kPhiStepDirInversionFailed;
}
}  // namespace _optim_dal_internal
}  // namespace nsoptim

#endif  // NSOPTIM_OPTIMIZER_DAL_HELPER_HPP_
