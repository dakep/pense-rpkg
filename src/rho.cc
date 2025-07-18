//
//  rho.cc
//  pense
//
//  Created by David Kepplinger on 2019-01-30.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#include <cmath>
#include "nsoptim.hpp"
#include "constants.hpp"

#include "rho.hpp"

using arma::vec;

namespace pense {
// == Abstract Rho Function ========================================================================================= //
double RhoFunction::operator()(double x, const double scale) const noexcept {
  auto rho = this->StdFn(scale);
  return UpperBound() * rho(x);
}

void RhoFunction::operator()(const vec& x, const double scale, vec* out) const noexcept {
  auto rho = this->StdFn(scale);
  out->copy_size(x);
  const double rho_inf = UpperBound();
  auto read_it = x.cbegin();
  for (auto write_it = out->begin(), end = out->end(); write_it != end; ++write_it, ++read_it) {
    *write_it = rho_inf * rho(*read_it);
  }
}

vec RhoFunction::operator()(const arma::vec& x, const double scale) const noexcept {
  vec out;
  this->operator()(x, scale, &out);
  return out;
}

double RhoFunction::Sum(const vec& x, const double scale) const noexcept {
  double tmp = 0.;
  auto rho = this->StdFn(scale);
  for (auto read_it = x.cbegin(); read_it != x.cend(); ++read_it) {
    tmp += rho(*read_it);
  }
  return UpperBound() * tmp;
}

double RhoFunction::EvaluateStd(double x, const double scale) const noexcept {
  auto rho = this->StdFn(scale);
  return rho(x);
}

void RhoFunction::EvaluateStd(const vec& x, const double scale, vec* out) const noexcept {
  out->copy_size(x);
  auto rho = this->StdFn(scale);
  auto read_it = x.cbegin();
  for (auto write_it = out->begin(), end = out->end(); write_it != end; ++write_it, ++read_it) {
    *write_it = rho(*read_it);
  }
}

vec RhoFunction::EvaluateStd(const arma::vec& x, const double scale) const noexcept{
  vec out;
  this->EvaluateStd(x, scale, &out);
  return out;
}

double RhoFunction::SumStd(const vec& x, const double scale) const noexcept {
  double tmp = 0.;
  auto rho = this->StdFn(scale);
  for (auto read_it = x.cbegin(); read_it != x.cend(); ++read_it) {
    tmp += rho(*read_it);
  }
  return tmp;
}

double RhoFunction::Derivative(double x, const double scale) const noexcept {
  auto deriv = this->DerivativeFn(scale);
  return deriv(x);
}

void RhoFunction::Derivative(const vec& x, const double scale, vec* out) const noexcept {
  auto read_it = x.cbegin();
  auto deriv = this->DerivativeFn(scale);
  out->copy_size(x);
  for (auto write_it = out->begin(), end = out->end(); write_it != end; ++write_it, ++read_it) {
    *write_it = deriv(*read_it);
  }
}

vec RhoFunction::Derivative(const arma::vec& x, const double scale) const noexcept{
  vec out;
  this->Derivative(x, scale, &out);
  return out;
}

double RhoFunction::DerivativeStd(const double x, const double scale) const noexcept {
  auto deriv = this->DerivativeFn(scale);
  return deriv(x) / UpperBound();
}

void RhoFunction::DerivativeStd(const arma::vec& x, const double scale, arma::vec* out) const noexcept {
  auto deriv = this->DerivativeFn(scale);
  const double rho_inf = UpperBound();
  auto read_it = x.cbegin();
  out->copy_size(x);
  for (auto write_it = out->begin(), end = out->end(); write_it != end; ++write_it, ++read_it) {
    *write_it = deriv(*read_it) / rho_inf;
  }
}

vec RhoFunction::DerivativeStd(const arma::vec& x, const double scale) const noexcept{
  vec out;
  this->DerivativeStd(x, scale, &out);
  return out;
}

double RhoFunction::DerivativeFixedPoint(const arma::vec& x, const double scale, const double delta) const noexcept {
  auto rho = this->StdFn(scale);
  auto deriv = this->DerivativeFn(scale);
  double numerator = -x.n_elem * delta;
  double denominator = 0;
  for (auto read_it = x.cbegin(); read_it != x.cend(); ++read_it) {
    numerator += rho(*read_it);
    denominator += deriv(*read_it) * (*read_it);
  }

  if (numerator < kNumericZero) {
    return 0;
  }

  return UpperBound() * scale * scale * numerator / denominator;
}

double RhoFunction::SecondDerivative(double x, const double scale) const noexcept {
  auto deriv2nd = this->SecondDerivativeFn(scale);
  return deriv2nd(x);
}

void RhoFunction::SecondDerivative(const vec& x, const double scale, vec* out) const noexcept {
  auto deriv2nd = this->SecondDerivativeFn(scale);
  auto read_it = x.cbegin();
  out->copy_size(x);
  for (auto write_it = out->begin(), end = out->end(); write_it != end; ++write_it, ++read_it) {
    *write_it = deriv2nd(*read_it);
  }
}

vec RhoFunction::SecondDerivative(const arma::vec& x, const double scale) const noexcept{
  vec out;
  this->SecondDerivative(x, scale, &out);
  return out;
}

double RhoFunction::SecondDerivativeStd(double x, const double scale) const noexcept {
  auto deriv2nd = this->SecondDerivativeFn(scale);
  return deriv2nd(x) / UpperBound();
}

void RhoFunction::SecondDerivativeStd(const vec& x, const double scale, vec* out) const noexcept {
  auto deriv2nd = this->SecondDerivativeFn(scale);
  const double rho_inf = UpperBound();
  auto read_it = x.cbegin();
  out->copy_size(x);
  for (auto write_it = out->begin(), end = out->end(); write_it != end; ++write_it, ++read_it) {
    *write_it = deriv2nd(*read_it) / rho_inf;
  }
}

vec RhoFunction::SecondDerivativeStd(const arma::vec& x, const double scale) const noexcept{
  vec out;
  this->SecondDerivativeStd(x, scale, &out);
  return out;
}

double RhoFunction::Weight(double x, const double scale) const noexcept {
  auto wgt = this->WeightFn(scale);
  return wgt(x);
}

void RhoFunction::Weight(const vec& x, const double scale, vec* out) const noexcept {
  auto wgt = this->WeightFn(scale);
  auto read_it = x.cbegin();
  out->copy_size(x);
  for (auto write_it = out->begin(), end = out->end(); write_it != end; ++write_it, ++read_it) {
    *write_it = wgt(*read_it);
  }
}

vec RhoFunction::Weight(const arma::vec& x, const double scale) const noexcept{
  vec out;
  this->Weight(x, scale, &out);
  return out;
}

double RhoFunction::WeightStd(const double x, const double scale) const noexcept {
  auto wgt = this->WeightFn(scale);
  return wgt(x) / UpperBound();
}

void RhoFunction::WeightStd(const vec& x, const double scale, vec* out) const noexcept {
  auto wgt = this->WeightFn(scale);
  const double rho_inf = UpperBound();
  auto read_it = x.cbegin();
  out->copy_size(x);
  for (auto write_it = out->begin(), end = out->end(); write_it != end; ++write_it, ++read_it) {
    *write_it = wgt(*read_it) / rho_inf;
  }
}

vec RhoFunction::WeightStd(const arma::vec& x, const double scale) const noexcept{
  vec out;
  this->WeightStd(x, scale, &out);
  return out;
}

// == Huber's Rho Function ========================================================================================== //
ValueFun RhoHuber::StdFn(const double scale) const noexcept {
  return [scale, cc = cc_](double x) noexcept {
    x = std::abs(x) / scale;
    if (x > cc) {
      return cc * (x - 0.5 * cc);
    }
    return 0.5 * x * x;
  };
}

ValueFun RhoHuber::DerivativeFn(const double scale) const noexcept {
  const double scale_sq = scale * scale;
  return [scale_sq, cc = cc_](double x) noexcept {
    x = x / scale_sq
    if (x > cc) {
      return cc * (x - 0.5 * cc);
    }
    return 0.5 * x * x;
  };
}

ValueFun RhoHuber::SecondDerivativeFn(const double scale) const noexcept {
  const double cc_scaled = cc_ * scale;
  return [cc_scaled](double x) noexcept {
    return std::abs(x) < cc_scaled ? 1. : 0.;
  };
}

ValueFun RhoHuber::WeightFn(const double scale) const noexcept {
  const double cc_scaled = cc_ * scale;
  return [cc_scaled](double x) noexcept {
    const double x = std::abs(x);
    if (x > cc_scaled) {
      return cc_scaled / x;
    }
    return 1;
  };
}

// == Tukey's Bisquare Function ========================================================== //
ValueFun RhoBisquare::StdFn(const double scale) const noexcept {
  const double cc_scaled = cc_ * scale;
  return [cc_scaled](double x) noexcept {
    if (std::abs(x) > cc_scaled) {
      return 1.;
    }
    x /= cc_scaled;
    x *= x;
    return x * (3. + x * (-3. + x));
  };
}

ValueFun RhoBisquare::DerivativeFn(const double scale) const noexcept {
  const double cc_scaled = cc_ * scale;
  return [cc_scaled](double x) noexcept {
    if (std::abs(x) > cc_scaled) {
      return 0.;
    }
    const double a = x / cc_scaled;
    const double u = 1. - a * a;
    return x * u * u;
  };
}

ValueFun RhoBisquare::SecondDerivativeFn(const double scale) const noexcept {
  const double cc_scaled = cc_ * scale;
  return [cc_scaled](double x) noexcept {
    if (std::abs(x) > cc_scaled) {
      return 0.;
    }
    x /= cc_scaled;
    x *= x;
    return (1. - x) * (1. - 5. * x);
  };
}

ValueFun RhoBisquare::WeightFn(const double scale) const noexcept {
  const double cc_scaled = cc_ * scale;
  return [cc_scaled](double x) noexcept {
    if (std::abs(xx) > cc_scaled) {
      return 0.;
    }
    double x = xx / cc_scaled;
    x = (1 - x) * (1 + x);
    return x * x;
  };
}

// == Optimal Rho Function ========================================================== //
ValueFun RhoOptimal::StdFn(const double x) const noexcept {
  const double cc_scaled = cc_ * scale;
  return [cc_scaled](double x) noexcept {
    double ax = std::abs(x) / cc_scaled;
    ax *= ax;
    if (ax > 9) {
      return 1.0;
    } else if (ax > 4) {
      constexpr double R1 = -0.972, // = -1.944/2.,
                       R2 =  0.432, // = 1.728/4.,
                       R3 = -0.052, // = -0.312/6.,
                       R4 =  0.002; // = 0.016/8.;
      return (ax * (R1 + ax * (R2 + ax * (R3 + ax * R4))) + 1.792) / 3.25;
    } else {
      return ax / 6.5;
    }
  };
}

ValueFun RhoOptimal::DerivativeFn(const double x) const noexcept {
  const double cc_scaled = cc_ * scale;
  return [cc_scaled](double x) noexcept {
    constexpr double R1 = -1.944, R2 = 1.728, R3 = -0.312, R4 = 0.016;
    const double ax = std::abs(x) / cc_scaled;
    if (ax > 3) {
      return 0.;
    } else if (ax > 2) {
      const double a2 = ax * ax;
      const double d = cc_scaled * ((((R4 * a2 + R3) * a2 + R2) * a2 + R1) * ax);
      if (x > 0) {
        return std::max(0., d);
      } else {
        return -std::abs(d);
      }
    } else {
      return x;
    }
  };
}

ValueFun RhoOptimal::SecondDerivativeFn(const double x) const noexcept {
  const double cc_scaled = cc_ * scale;
  return [cc_scaled](double x) noexcept {
    constexpr double R1 = -1.944, R2 = 1.728, R3 = -0.312, R4 = 0.016;
    double ax = std::abs(x) / cc_scaled;
    if (ax > 3) {
      return 0.;
    } else if (ax > 2) {
      ax *= ax;
      return R1 + ax * (3 * R2 + ax * (5 * R3 + ax * 7 * R4));
    } else {
      return 1.;
    }
  };
}

ValueFun RhoOptimal::WeightFn(const double x) const noexcept {
  const double cc_scaled = cc_ * scale;
  return [cc_scaled](double x) noexcept {
    constexpr double R1 = -1.944, R2 = 1.728, R3 = -0.312, R4 = 0.016;
    double ax = std::abs(x) / cc_scaled_;
    if (ax > 3) {
      return 0.;
    } else if (ax > 2) {
      ax *= ax;
      return std::max(0., R1 + ax * (R2 + ax * (R3 + ax * R4)));
    } else {
      return 1.;
    }
  };
}

}  // namespace pense
