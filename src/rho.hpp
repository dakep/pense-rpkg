//
//  rho.hpp
//  pense
//
//  Created by David Kepplinger on 2019-01-30.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef RHO_HPP_
#define RHO_HPP_

#include "nsoptim.hpp"
#include <limits>
#include <functional>

namespace pense {
class RhoFunction {
 public:
  explicit RhoFunction () noexcept {}
  RhoFunction(const RhoFunction&) = default;
  RhoFunction& operator=(const RhoFunction&) = default;
  RhoFunction(RhoFunction&&) = default;
  RhoFunction& operator=(RhoFunction&&) = default;
  ~RhoFunction() {}

  //! Get the value of the rho function evaluated at x/scale.
  //! This is equivalent to `robustbase::Mpsi(x / scale, cc, psifun, deriv = -1)`.
  double    operator()(const double x, const double scale) const noexcept;
  void      operator()(const arma::vec& x, const double scale, arma::vec* out) const noexcept;
  arma::vec operator()(const arma::vec& x, const double scale) const noexcept;
  double    Sum(const arma::vec& x, const double scale) const noexcept;

  //! Get the value of the *standardized* rho function evaluated at x/scale.
  //! Standardized means that rho(inf) = 1.
  //! This is equivalent to `robustbase::Mchi(x / scale, cc, psifun, deriv = 0)`.
  double    EvaluateStd(const double x, const double scale) const noexcept;
  void      EvaluateStd(const arma::vec& x, const double scale, arma::vec* out) const noexcept;
  arma::vec EvaluateStd(const arma::vec& x, const double scale) const noexcept;
  double    SumStd(const arma::vec& x, const double scale) const noexcept;

  //! Get the derivative of the rho function evaluated at x/scale.
  //! Note that Derivative(x/scale, 1) == Derivative(x, scale) * scale
  //! This is equivalent to `robustbase::Mpsi(x / scale, cc, psifun, deriv = 0)`.
  double    Derivative(const double x, const double scale) const noexcept;
  void      Derivative(const arma::vec& x, const double scale, arma::vec* out) const noexcept;
  arma::vec Derivative(const arma::vec& x, const double scale) const noexcept;

  //! Get the derivative of the *standardized* rho function evaluated at x/scale.
  //! Note that DerivativeStd(x/scale, 1) == DerivativeStd(x, scale) / scale.
  //! This is equivalent to `robustbase::Mchi(x / scale, cc, psifun, deriv = 1)`.
  double    DerivativeStd(const double x, const double scale) const noexcept;
  void      DerivativeStd(const arma::vec& x, const double scale, arma::vec* out) const noexcept;
  arma::vec DerivativeStd(const arma::vec& x, const double scale) const noexcept;

  //! Compute the Newton step for root-finding of the M-scale equation in one pass over the vector `x`.
  //! This function computes mean(rho(x / scale; cc) - delta) / mean(rho'(x / scale; cc) * x / scale)
  double DerivativeFixedPoint(const arma::vec& x, const double scale, const double delta) const noexcept;

  //! Get the second derivative of the rho function evaluated at x/scale.
  //! Note that SecondDerivative(x/scale, 1) == SecondDerivative(x, scale)
  //! This is equivalent to `robustbase::Mpsi(x / scale, cc, psifun, deriv = 1)`.
  double    SecondDerivative(const double x, const double scale) const noexcept;
  void      SecondDerivative(const arma::vec& x, const double scale, arma::vec* out) const noexcept;
  arma::vec SecondDerivative(const arma::vec& x, const double scale) const noexcept;

  //! Get the second derivative of the *standardized* rho function evaluated at x/scale.
  //! Note that SecondDerivativeStd(x/scale, 1) == SecondDerivativeStd(x, scale) * scale^2
  //! This is equivalent to `robustbase::Mchi(x / scale, cc, psifun, deriv = 2)`.
  double    SecondDerivativeStd(const double x, const double scale) const noexcept;
  void      SecondDerivativeStd(const arma::vec& x, const double scale, arma::vec* out) const noexcept;
  arma::vec SecondDerivativeStd(const arma::vec& x, const double scale) const noexcept;

  //! Get the derivative of the rho function evaluated at (x/scale), divided by (x/scale), i.e.,
  //! Weight(x, scale) = Derivative(x/scale) / (x/scale).
  //! Note that Weight(x/scale, 1) == Weight(x, scale)
  //! This is equivalent to `robustbase::Mwgt(x / scale, cc, psifun)`.
  double    Weight(const double x, const double scale) const noexcept;
  void      Weight(const arma::vec& x, const double scale, arma::vec* out) const noexcept;
  arma::vec Weight(const arma::vec& x, const double scale) const noexcept;

  //! Get the derivative of the *standardized* rho function evaluated at (x/scale), divided by (x/scale), i.e.,
  //! WeightStd(x, scale) = DerivativeStd(x/scale) / (x/scale).
  //! Note that WeightStd(x/scale, 1) == WeightStd(x, scale)
  //! This is equivalent to
  //!   `robustbase::Mwgt(x / scale, cc, psifun) / robustbase::MrhoInf(cc, psifun)`.
  double    WeightStd(double x, const double scale) const noexcept;
  void      WeightStd(const arma::vec& x, const double scale, arma::vec* out) const noexcept;
  arma::vec WeightStd(const arma::vec& x, const double scale) const noexcept;

  //! Get the upper bound of the rho function, i.e., the limiting value of `operator()(x)` for x to infinity.
  virtual double UpperBound() const noexcept = 0;

 protected:
  using ValueFun = std::function<double(double)>;

  virtual ValueFun StdFn(const double scale) const noexcept = 0;
  virtual ValueFun DerivativeFn(const double scale) const noexcept = 0;
  virtual ValueFun SecondDerivativeFn(const double scale) const noexcept = 0;
  virtual ValueFun WeightFn(const double scale) const noexcept = 0;
};

//! Implementation of Huber's unbounded rho function defined by
//! rho(x) = 0.5 x^2 * [x < cc] + cc * (abs(x) - cc / 2) [x >= cc]
class RhoHuber : public RhoFunction {
 public:
  //! Create the Huber rho-function with cutoff `cc`
  //!
  //! @param cc cutoff value (outside the rho-function is linear)
  explicit RhoHuber(const double cc) noexcept : cc_(cc) {}

  RhoHuber(const RhoHuber&) = default;
  RhoHuber& operator=(const RhoHuber&) = default;
  RhoHuber(RhoHuber&&) = default;
  RhoHuber& operator=(RhoHuber&&) = default;

  ~RhoHuber() {}

  //! Get the value of the threshold paramter.
  //!
  //! @return threshold paramter value.
  double cc() const noexcept {
    return cc_;
  }

  double UpperBound() const noexcept override { return std::numeric_limits<double>::infinity(); }

 protected:
  ValueFun StdFn(const double scale) const noexcept override;
  ValueFun DerivativeFn(const double scale) const noexcept override;
  ValueFun SecondDerivativeFn(const double scale) const noexcept override;
  ValueFun WeightFn(const double scale) const noexcept override;

 private:
  double cc_;
};

//! Implementation of Tukey's bisquare rho function defined by
//! rho(x) = min(1, 1 - (1 - x^2/cc^2)^3)
//! All of the virtual functions are declared final and can not be overwritten again!
class RhoBisquare : public RhoFunction {
 public:
  //! Create the bisquare-rho function with threshold `cc`
  //!
  //! @param cc threshold paramter value.
  explicit RhoBisquare(const double cc) noexcept : cc_(cc) {}

  RhoBisquare(const RhoBisquare&) = default;
  RhoBisquare& operator=(const RhoBisquare&) = default;
  RhoBisquare(RhoBisquare&&) = default;
  RhoBisquare& operator=(RhoBisquare&&) = default;

  ~RhoBisquare() {}

  //! Get the value of the threshold paramter.
  //!
  //! @return threshold paramter value.
  double cc() const noexcept {
    return cc_;
  }

 protected:
  ValueFun StdFn(const double scale) const noexcept override;
  ValueFun DerivativeFn(const double scale) const noexcept override;
  ValueFun SecondDerivativeFn(const double scale) const noexcept override;
  ValueFun WeightFn(const double scale) const noexcept override;

 private:
  double cc_;
};

//! Implementation of the "Optimal" rho function given by Maronna et al. (2006, Section 5.9.1) defined by
//! rho(x) = min(1, 1 - (1 - x^2/cc^2)^3)
//! All of the virtual functions are declared final and can not be overwritten again!
class RhoOptimal : public RhoFunction {
 public:
  //! Create the optimal-rho function with threshold `cc`
  //!
  //! @param cc threshold paramter value.
  explicit RhoOptimal(const double cc) noexcept : cc_(cc) {}

  RhoOptimal(const RhoOptimal&) = default;
  RhoOptimal& operator=(const RhoOptimal&) = default;
  RhoOptimal(RhoOptimal&&) = default;
  RhoOptimal& operator=(RhoOptimal&&) = default;

  ~RhoOptimal() {}

  //! Get the value of the threshold paramter.
  //!
  //! @return threshold paramter value.
  double cc() const noexcept {
    return cc_;
  }

 protected:
  ValueFun StdFn(const double scale) const noexcept override;
  ValueFun DerivativeFn(const double scale) const noexcept override;
  ValueFun SecondDerivativeFn(const double scale) const noexcept override;
  ValueFun WeightFn(const double scale) const noexcept override;

 private:
  double cc_;
};
}  // namespace pense
#endif  // RHO_HPP_
