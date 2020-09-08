//
//  optimum.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-11-30.
//  Copyright © 2018 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_OPTIMIZER_OPTIMUM_HPP_
#define NSOPTIM_OPTIMIZER_OPTIMUM_HPP_

#include <string>
#include <limits>
#include <memory>
#include <type_traits>

#include "../container/metrics.hpp"

namespace nsoptim {

enum class OptimumStatus { kOk, kWarning, kError };

namespace optimum_internal {
using MetricsPtr = std::unique_ptr<Metrics>;

//! Wrapper around the information at an optimum point.
template <typename T, typename U, typename V>
struct Optimum {
  using LossFunction = T;
  using PenaltyFunction = U;
  using Coefficients = V;

  Optimum(const LossFunction& _loss, const PenaltyFunction& _penalty) noexcept : loss(_loss), penalty(_penalty) {}
  Optimum(const LossFunction& _loss, const PenaltyFunction& _penalty, const Coefficients& _coefs,
          const double _objf_value, MetricsPtr _metrics,
          const OptimumStatus _status, const std::string& _message) noexcept
    : loss(_loss), penalty(_penalty), coefs(_coefs), objf_value(_objf_value), metrics(std::move(_metrics)),
      status(_status), message(_message) {}
  Optimum(const Optimum& other) noexcept : loss(other.loss), penalty(other.penalty), coefs(other.coefs),
                                           objf_value(other.objf_value),
                                           metrics(other.metrics ? new Metrics(*other.metrics) : nullptr),
                                           status(other.status), message(other.message) {}

  Optimum(Optimum&& other) = default;

  //! Move-assignable operator must be explicity defined.
  Optimum& operator=(Optimum&& other) {
    loss = std::move(other.loss);
    penalty = std::move(other.penalty);
    coefs = std::move(other.coefs);
    objf_value = other.objf_value;
    status = other.status;
    message = std::move(other.message);
    metrics = std::move(other.metrics);
    return *this;
  }

  //! A copy of the loss function for which this optimum was attained.
  LossFunction loss;
  //! A copy of the penalty function for which this optimum was attained.
  PenaltyFunction penalty;
  //! The coefficients at which the objective function attains its optimum.
  Coefficients coefs;
  //! The value of the objective function at this optimum.
  double objf_value = std::numeric_limits<double>::max();
  //! Optional metrics associated with this optimum.
  MetricsPtr metrics;
  //! The status of the optimizer at the time this optimum was found.
  OptimumStatus status = OptimumStatus::kError;
  //! An optional status message of the optimizer at the time this optimum was found.
  std::string message;
};
}  // namespace optimum_internal

//! Wrapper around the information at an optimum point.
template <typename LossFunction, typename PenaltyFunction, typename Coefficients>
using Optimum = optimum_internal::Optimum<typename std::decay<LossFunction>::type,
                                          typename std::decay<PenaltyFunction>::type,
                                          typename std::decay<Coefficients>::type>;


//! Create an Optimum from the given arguments.
template <typename LossFunction, typename PenaltyFunction, typename Coefficients>
Optimum<LossFunction, PenaltyFunction, Coefficients> MakeOptimum(
    const LossFunction& loss, const PenaltyFunction& penalty, const Coefficients& coefs,
    const double objf_value, optimum_internal::MetricsPtr metrics,
    const OptimumStatus status = OptimumStatus::kOk, const std::string& message = {}) noexcept {
  return Optimum<LossFunction, PenaltyFunction, Coefficients>(loss, penalty, coefs, objf_value, std::move(metrics),
                                                              status, message);
}

template <typename LossFunction, typename PenaltyFunction, typename Coefficients>
Optimum<LossFunction, PenaltyFunction, Coefficients> MakeOptimum(
    const LossFunction& loss, const PenaltyFunction& penalty, const Coefficients& coefs,
    optimum_internal::MetricsPtr metrics,
    const OptimumStatus status = OptimumStatus::kOk, const std::string& message = {}) noexcept {
  return Optimum<LossFunction, PenaltyFunction, Coefficients>(loss, penalty, coefs, loss(coefs) + penalty(coefs),
                                                              std::move(metrics), status, message);
}

//! Create an Optimum from the given arguments.
template <typename LossFunction, typename PenaltyFunction, typename Coefficients>
Optimum<LossFunction, PenaltyFunction, Coefficients> MakeOptimum(
    const LossFunction& loss, const PenaltyFunction& penalty, const Coefficients& coefs,
    const OptimumStatus status = OptimumStatus::kOk, const std::string& message = {}) noexcept {
  return Optimum<LossFunction, PenaltyFunction, Coefficients>(loss, penalty, coefs, loss(coefs) + penalty(coefs),
                                                              nullptr, status, message);
}

//! Create an Optimum from the given arguments.
template <typename LossFunction, typename PenaltyFunction, typename Coefficients>
Optimum<LossFunction, PenaltyFunction, Coefficients> MakeOptimum(
    const LossFunction& loss, const PenaltyFunction& penalty, const Coefficients& coefs,
    const double objf_value,
    const OptimumStatus status = OptimumStatus::kOk, const std::string& message = {}) noexcept {
  return Optimum<LossFunction, PenaltyFunction, Coefficients>(loss, penalty, coefs, nullptr, objf_value, status,
                                                              message);
}

}  // namespace nsoptim

#endif  // NSOPTIM_OPTIMIZER_OPTIMUM_HPP_
