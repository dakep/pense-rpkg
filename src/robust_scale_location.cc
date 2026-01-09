//
//  robust_scale_location.cc
//  pense
//
//  Created by David Kepplinger on 2019-01-30.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#include "nsoptim.hpp"
#include "robust_scale_location.hpp"
#include "constants.hpp"

using arma::vec;
using arma::median;
using arma::abs;
using arma::uword;

namespace {
constexpr double kTauSizeC2Squared = 9.;
constexpr double kTauSizeConsistencyConstant = 1. / 0.961;

constexpr double kMadScaleConsistencyConstant = 1.4826;
}  // namespace

namespace pense {
double TauSize(const vec& values) noexcept {
  const vec abs_values(abs(values));
  const double sigma_0 = median(abs_values);

  if (sigma_0 < kNumericZero) {
    return 0.;
  }

  const double tau_size = arma::mean(arma::clamp(arma::square(abs_values / sigma_0),
                                                 0, kTauSizeC2Squared));
  return sigma_0 * kTauSizeConsistencyConstant * sqrt(tau_size);
}

double MLocation(const arma::vec& values, const RhoFunction& rho, const double scale,
                 const double convergence_tol, const int max_it) {
  const double scaled_conv_tol = convergence_tol * scale;
  int it = 0;
  double location = arma::median(values);

  arma::vec residuals(values.n_elem);
  arma::vec w_loc(values.n_elem);
  while (it++ < max_it) {
    residuals = values - location;
    rho.Weight(residuals, scale, &w_loc);
    const double prev_location = location;
    const double w_loc_sum = arma::accu(w_loc);

    if (w_loc_sum < convergence_tol) {
      throw ZeroWeightsException();
    }

    location = arma::accu(w_loc % values) / w_loc_sum;

    if (std::abs(prev_location - location) < scaled_conv_tol) {
      break;
    }
  }

  return location;
}

LocationScaleEstimate MLocationScale(const arma::vec& values, const Mscale& mscale,
                                     const RhoFunction& location_rho) {
  int it = 0;
  LocationScaleEstimate est {arma::median(values)};
  est.scale = robust_scale_location::InitialScaleEstimate(values - est.location, mscale.delta(), mscale.eps());

  if (est.scale < mscale.eps()) {
    est.scale = 0;
    return est;
  }

  const double convergence_tol = est.scale * mscale.eps();
  const double recip_sqrt_delta = 1. / sqrt(mscale.delta());

  arma::vec residuals(values.n_elem);
  arma::vec w_loc(values.n_elem);
  while (it++ < mscale.max_it()) {
    residuals = values - est.location;
    location_rho.Weight(residuals, est.scale, &w_loc);
    const double w_scale_mean = mscale.rho()->SumStd(residuals, est.scale) / residuals.n_elem;
    const double w_loc_sum = arma::accu(w_loc);

    if (w_loc_sum < convergence_tol) {
      throw ZeroWeightsException();
    }

    const LocationScaleEstimate prev_est = est;
    est.location = arma::accu(w_loc % values) / w_loc_sum;
    est.scale = prev_est.scale * std::sqrt(w_scale_mean) * recip_sqrt_delta;

    if (std::abs(prev_est.location - est.location) < convergence_tol &&
        std::abs(prev_est.scale - est.scale) < convergence_tol) {
      break;
    }
  }

  return est;
}

// == Mscale class implementations ====================================================
arma::vec Mscale::Derivative(const arma::vec& values) const {
  const double scale = ComputeMscale(values, InitialEstimate(values));
  if (scale < eps_) {
    return arma::vec();
  }

  const auto deriv_rho = rho_->Derivative(values, scale);
  const auto denom = sum(deriv_rho % values) / scale;
  if (denom < eps_) {
    return arma::vec(values.n_elem, arma::fill::value(R_PosInf));
  } else {
    return deriv_rho / denom;
  }
}

double Mscale::ComputeMscale(const arma::vec& values, const double init_scale) const {
  if (init_scale < kNumericZero) {
    return 0;
  }

  int iter = 0;
  double step;
  double scale = init_scale;
  // Start Newton's iterations
  do {
    step = rho_->DerivativeFixedPoint(values, scale, delta_);
    scale += scale * step;
  } while (++iter < max_it_ && std::abs(step) > eps_ && scale > kNumericZero && std::isfinite(scale));

  if (scale < kNumericZero || !std::isfinite(scale)) {
    return ComputeMscaleFallback(values, max_it_ - iter, init_scale);
  }

  return scale;
}

double Mscale::ComputeMscale(const arma::vec& values, const double init_scale) {
  if (init_scale < kNumericZero) {
    return 0;
  }

  it_ = 0;
  double step;
  double scale = init_scale;
  // Start Newton's iterations
  do {
    step = rho_->DerivativeFixedPoint(values, scale, delta_);
    scale += scale * step;
  } while (++it_ < max_it_ && std::abs(step) > eps_ && scale > kNumericZero && std::isfinite(scale));

  if (scale < kNumericZero || !std::isfinite(scale)) {
    return ComputeMscaleFallback(values, max_it_ - it_, init_scale);
  }

  return scale;
}

double Mscale::ComputeMscaleFallback(const arma::vec& values, const int max_it, double scale) const {
  const double rho_denom = 1. / (delta_ * values.n_elem);

  int iter = 0;
  double err = eps_;
  // Start iterations
  do {
    const double rho_sum = rho_->SumStd(values, scale);
    const double new_scale = scale * std::sqrt(rho_sum * rho_denom);
    err = std::abs(new_scale - scale);
    scale = new_scale;
  } while (++iter < max_it && err > eps_ * scale && std::isfinite(scale));

  if (scale < kNumericZero || !std::isfinite(scale)) {
    scale = 0;
  }

  return scale;
}

arma::mat Mscale::GradientHessian(const arma::vec& values) const {
  const double scale = this->operator()(values);
  if (scale < eps_) {
    return arma::mat(1, 1, arma::fill::value(scale));
  }
  const auto violation = rho_->SumStd(values, scale) / values.n_elem - delta_;

  arma::mat grad_hess(values.n_elem, values.n_elem + 2, arma::fill::zeros);

  // Compute the gradient and its maximum
  grad_hess.col(0) = rho_->Derivative(values, scale);
  const auto denom = arma::sum(grad_hess.col(0) % values);

  grad_hess.at(1, 2) = denom;
  grad_hess.at(2, 2) = scale;
  grad_hess.at(3, 2) = violation;

  // Compute the Hessian and its maximum
  const auto rho_2nd = rho_->SecondDerivative(values, scale);
  const auto sum_2nd = arma::sum(rho_2nd % values % values) / denom;
  grad_hess.col(1) = rho_2nd;
  double diag_offset;
  for (int i = 0; i < values.n_elem; ++i) {
    diag_offset = denom * rho_2nd[i];
    for (int k = i; k < values.n_elem; ++k) {
      grad_hess(i, k + 2) = HessianElementUnscaled(
        i, k, grad_hess.unsafe_col(0), rho_2nd, values, sum_2nd, diag_offset);

      grad_hess(i, k + 2) *= scale / (denom * denom);
      diag_offset = 0;
    }
  }

  // Final pass to get gradient right
  grad_hess.col(0) *= scale / denom;

  return grad_hess;
}

arma::vec::fixed<3> Mscale::MaxGradientHessian(const arma::vec& values) const {
  arma::vec::fixed<3> maxima(arma::fill::zeros);
  maxima[0] = this->operator()(values, InitialEstimate(values));
  if (maxima[0] < eps_) {
    return maxima;
  }
  const auto violation = rho_->SumStd(values, maxima[0]) -
    values.n_elem * delta_;

  if (violation * violation > values.n_elem * values.n_elem * eps_ * eps_) {
    return maxima;
  }

  // Compute the gradient and its maximum
  const auto rho_1st = rho_->Derivative(values, maxima[0]);
  const auto denom = sum(rho_1st % values);
  maxima[1] = (denom < eps_) ? R_PosInf :
    (arma::max(rho_1st) * maxima[0] / denom);

  // Compute the Hessian and its maximum
  const auto rho_2nd = rho_->SecondDerivative(values, maxima[0]);
  const auto sum_2nd = sum(rho_2nd % values % values) / denom;
  double diag_offset;

  for (int i = 0; i < values.n_elem; ++i) {
    diag_offset = denom * rho_2nd[i];
    for (int k = i; k < values.n_elem; ++k) {
      const auto tmp = std::abs(HessianElementUnscaled(
        i, k, rho_1st, rho_2nd, values, sum_2nd, diag_offset));

      diag_offset = 0;
      if (tmp > maxima[2]) {
        maxima[2] = tmp;
      }
    }
  }
  maxima[2] *= maxima[0] / (denom * denom);

  return maxima;
}

namespace robust_scale_location {
double InitialScaleEstimate(const vec& values, const double delta, const double eps) {
  // Try the MAD of the uncentered values.
  double mad = 0;
  try {
    mad = kMadScaleConsistencyConstant * median(abs(values));
  } catch (...) {
    mad = 0;
  }
  if (mad > eps) {
    return mad;
  } else if (static_cast<uword>((1 - delta) * values.n_elem) > values.n_elem / 2) {
    // If the MAD is also (almost) 0, but the M-scale takes into account more observations than the MAD,
    // compute the variance of the additional elements (i.e., the variance without considering the smallest
    // 50% of the observations)
    const uword lower_index = values.n_elem / 2;
    const uword upper_index = static_cast<uword>((1 - delta) * values.n_elem);
    const vec ordered_values = arma::sort(abs(values));
    const double scale = arma::var(ordered_values.rows(lower_index, upper_index));
    if (scale > eps) {
      return scale;
    }
  }
  return 0.;
}

}  // namespace robust_scale_location

}  // namespace pense
