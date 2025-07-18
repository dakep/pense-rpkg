//
//  robust_scale_location.hpp
//  pense
//
//  Created by David Kepplinger on 2019-01-30.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef ROBUST_SCALE_LOCATION_HPP_
#define ROBUST_SCALE_LOCATION_HPP_

#include <memory>
#include <cmath>

#include "nsoptim.hpp"
#include "rho.hpp"
#include "constants.hpp"
#include "rcpp_utils.hpp"

namespace pense {
namespace robust_scale_location {
//! Default breakdown point for the M-scale equation.
constexpr double kDefaultMscaleDelta = 0.5;
//! Default number of iterations for the M-scale algorithm.
constexpr int kDefaultMscaleMaxIt = 100;

//! Compute an initial estimate of the scale of `values`.
//!
//! @param values a vector of values.
//! @param delta the right-hand-side in the M-scale equation.
//! @param eps numerical tolerance.
//! @return an initial, inaccurate, estimate of the scale of `values`.
double InitialScaleEstimate(const arma::vec& values, const double delta, const double eps);

}  // namespace robust_scale_location

//! The result of simultaneous estimation of the M-location and M-scale.
struct LocationScaleEstimate {
  double location;
  double scale;
};

//! Exception class for when all the weights become zero.
class ZeroWeightsException : public std::runtime_error {
 public:
  ZeroWeightsException() : std::runtime_error("all weights are zero") {}
};

//! Compute the tau-scale of the uncentered values.
//! This function fixes the value of C_2 to 3!
//!
//! @param values values
double TauSize(const arma::vec& values) noexcept;

//! Functor to compute the Mscale using the specified `rho` function.
//! The M-scale `s` is defined by the M-estimation equation:
//!  (1/n) sum_{i = 1}^n rho(value[i] / s) == delta
class Mscale {
  std::shared_ptr<const RhoFunction> rho_;
  double delta_;
  int max_it_;
  int it_ = -1;
  double eps_;
  double scale_;

 public:
  //! Construct the M-scale function.
  //!
  //! @param user_options an R list of user options.
  explicit Mscale(const Rcpp::List& user_options) noexcept
    : rho_(RhoFactory(user_options)),
      delta_(GetFallback(user_options, "delta",
        robust_scale_location::kDefaultMscaleDelta)),
      max_it_(GetFallback(user_options, "max_it",
        robust_scale_location::kDefaultMscaleMaxIt)),
      eps_(GetFallback(user_options, "eps", kDefaultConvergenceTolerance)),
      scale_(-1) {}

  //! Construct the M-scale function.
  //!
  //! @param rho rho function to use.
  //! @param delta right-hand side of the M-estimation equation.
  //! @param max_it maximum number of iterations.
  //! @param eps numerical tolerance for convergence.
  Mscale(std::shared_ptr<const RhoFunction> rho, const double delta, const int max_it,
         const double eps) noexcept
      : rho_(rho), delta_(delta), max_it_(max_it), eps_(eps), scale_(-1) {}

  Mscale(const Mscale&) = default;
  Mscale& operator=(const Mscale&) = default;
  Mscale(Mscale&&) = default;
  Mscale& operator=(Mscale&&) = default;

  //! Reset the M-scale estimate.
  void Reset() noexcept {
    scale_ = -1;
  }

  //! Compute the M-scale of the given values. The initial guess is the median of the absolute values (MAD).
  //!
  //! @param values a vector of values.
  //! @return the M-scale of the given values.
  double operator()(const arma::vec& values) const {
    return ComputeMscale(values, InitialEstimate(values));
  }

  //! Compute the M-scale of the given values. The initial guess is either the given scale (if positive),
  //! the previous scale estimate (if availalbe), or the median of the absolute values (MAD).
  //!
  //! @param values a vector of values.
  //! @return the M-scale of the given values.
  double operator()(const arma::vec& values, double scale) const {
    return ComputeMscale(values, scale);
  }

  //! Compute the M-scale of the given values. The previous scale estimate (if present) is used as an initial guess.
  //! If no previous estimate is present (or the M-scale object was reset), the median of the absolute values (MAD)
  //! is used as an initial guess.
  //!
  //! @param values a vector of values.
  //! @return the M-scale of the given values.
  double operator()(const arma::vec& values) {
    scale_ = ComputeMscale(values, InitialEstimate(values));
    return scale_;
  }

  //! Set the initial estimate of the scale.
  //!
  //! @param scale initial scale estimate (> 0).
  void SetInitial(const double scale) noexcept {
    scale_ = scale;
  }

  //! Get the number of iterations required for the last scale estimate.
  //!
  //! @return number of iterations, or -1 if none have been performed yet.
  int LastIterations() const noexcept {
    return it_;
  }

  //! Compute the 1st derivative of the M-scale function with respect to each element.
  //!
  //! @param values vector of values
  //! @return a vector of derivatives, one for each element in `values`. If the scale is 0 or the M-scale equation
  //!   is violated, an empty vector is returned.
  arma::vec Derivative(const arma::vec& values) const;

  //! Compute the gradient and Hessian of the M-scale
  //! function evaluated at the given vector.
  //!
  //! @param values vector of values
  //! @return a matrix of dimension n x n + 1, where the 1st column is the
  //!    gradient and the other columns are the Hessian matrix.
  //!    If the scale is 0 or the M-scale equation
  //!    is violated, an empty matrix is returned.
  arma::mat GradientHessian(const arma::vec& values) const;

  //! Compute the maximum of the 1st and 2nd derivatives of the M-scale
  //! function evaluated at all elements in the given vector.
  //!
  //! @param values vector of values
  //! @return a vector with 3 elements: the M-scale,
  //!    the maximum element of the gradient
  //!    and the maximum element in the Hessian.
  //!    If the scale is 0 or the M-scale equation
  //!    is violated, an empty vector is returned.
  arma::vec::fixed<3> MaxGradientHessian(const arma::vec& values) const;

  //! Get the rho function object.
  std::shared_ptr<const RhoFunction> rho() const noexcept {
    return rho_;
  }

  //! Get the maximum number of iterations.
  int max_it() const noexcept {
    return max_it_;
  }

  //! Get the numerical tolerance.
  double eps() const noexcept {
    return eps_;
  }

  //! Get delta, i.e., the right-hand-side of the M-scale equation
  double delta() const noexcept {
    return delta_;
  }

 private:
  double ComputeMscale(const arma::vec& values, const double init_scale) const;
  double ComputeMscale(const arma::vec& values, const double init_scale);

  //! The Newton iterations are unstable if outliers have a strong effect on the scale estimate.
  //! In these cases, the sublinear method works more reliably, but it is slow.
  double ComputeMscaleFallback(const arma::vec& values, const int max_it, double scale) const;

  double InitialEstimate(const arma::vec& values) const {
    // If the internal scale is already set, use it as initial estimate.
    if (scale_ > eps_) {
      return scale_;
    }
    return robust_scale_location::InitialScaleEstimate(values, delta_, eps_);
  }

  double HessianElementUnscaled(const int i, const int k,
                                const arma::vec& rho_1st,
                                const arma::vec& rho_2nd,
                                const arma::vec& values,
                                const double sum_2nd,
                                const double diag_offset) const {
    return diag_offset +
      rho_1st[i] * rho_1st[k] * sum_2nd -
      rho_1st[i] * rho_2nd[k] * values[k] -
      rho_1st[k] * rho_2nd[i] * values[i];
  }
};

//! Computation of the M-location of the given vector.
//!
//! @param values values to compute the location and scale from.
//! @param rho rho-function for the M-location.
//! @param scale the scale of the values.
//! @param convergence_tol numeric convergence tolerance.
//! @param max_it maximum number of iterations.
//! @return location of the given values.
double MLocation(const arma::vec& values, const RhoFunction& rho, const double scale,
                 const double convergence_tol, const int max_it);

//! Simultaneous computation of the M-location and M-scale of the given vector.
//!
//! @param values values to compute the location and scale from.
//! @param mscale M-scale definition.
//! @param location_rho rho-function for the M-location.
//! @return location (first) and scale (second) of the given values.
LocationScaleEstimate MLocationScale(const arma::vec& values, const Mscale& mscale,
                                     const RhoFunction& location_rho);
}  // namespace pense

#endif  // ROBUST_SCALE_LOCATION_HPP_
