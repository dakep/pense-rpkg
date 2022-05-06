//
//  coordinate_descent.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2022-02-24.
//  Copyright © 2022 David Kepplinger. All rights reserved.
//

#ifndef PENSE_CD_PENSE_HPP_
#define PENSE_CD_PENSE_HPP_

#include <exception>
#include <forward_list>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>

#include "nsoptim.hpp"
#include "s_loss.hpp"
#include "robust_scale_location.hpp"

namespace pense {
//! Configuration options for the DAL algorithm.
struct CDPenseConfiguration {
  //! Maximum number of iterations allowed.
  int max_it;
  //! Enable/disable linesearch
  bool linesearch;
  //! Multiplier for adjusting the step size during line search.
  double linesearch_ss_multiplier;
  //! Number of step sizes to be considered for line search.
  int linesearch_ss_num;
  //! Re-compute the residuals every `reset_iter` iterations to avoid drift.
  int reset_iter;
};

namespace coorddesc {
constexpr CDPenseConfiguration kDefaultCDConfiguration = { 1000, false, 1e-6, 10, 8 };

struct SurrogateGradient {
  const double gradient;
  const double lipschitz_constant;
};

template<class Coefficients>
struct State {
  Coefficients coefs;
  arma::vec residuals;
  double mscale;
  double objf_loss;
  double objf_pen;
};

} // namespace coorddesc

//! Compute the EN regression estimate using the LARS algorithm on the
//! augmented response vector and predictor matrix.
template<class PenaltyFunction, class Coefficients>
class CDPense :
    public nsoptim::Optimizer<SLoss, PenaltyFunction, Coefficients> {
  using Base = nsoptim::Optimizer<SLoss, PenaltyFunction, Coefficients>;
  using LossFunctionPtr = std::unique_ptr<SLoss>;
  using PenaltyPtr = std::unique_ptr<PenaltyFunction>;
  using IsAdaptiveTag = typename nsoptim::traits::is_adaptive<PenaltyFunction>::type;
  using IsSparseTag = typename
    std::is_same<typename Coefficients::SlopeCoefficient, arma::sp_vec>::type;

  static_assert(nsoptim::traits::is_en_penalty<PenaltyFunction>::value,
                "PenaltyFunction must be an EN-type penalty.");

 public:
  using Optimum = typename Base::Optimum;

    //! Ininitialize the optimizer without a loss or penalty function.
  CDPense(
    const CDPenseConfiguration& config = coorddesc::kDefaultCDConfiguration) noexcept
      : config_(config) {}

  //! Ininitialize the optimizer using the given (weighted) LS loss function
  //! and penalty function.
  //! @param loss a weighted LS loss function.
  //! @param penalty penalty function.
  CDPense(const SLoss& loss,
    const PenaltyFunction& penalty,
    const CDPenseConfiguration& config = coorddesc::kDefaultCDConfiguration) noexcept
    : loss_(new SLoss(loss)),
      penalty_(new PenaltyFunction(penalty)), config_(config) {}

  //! Default copy constructor.
  //!
  //! The copied optimizer will share the identical loss and penalty
  //! functions after construction.
  CDPense(const CDPense& other) noexcept
    : loss_(other.loss_? new SLoss(*other.loss_) : nullptr),
      penalty_(other.penalty_ ? new PenaltyFunction(*other.penalty_) : nullptr),
      config_(other.config_),
      lipschitz_bounds_(other.lipschitz_bounds_),
      lipschitz_bound_intercept_(other.lipschitz_bound_intercept_),
      state_(other.state_),
      convergence_tolerance_(other.convergence_tolerance_) {}

  //! Default copy assignment.
  //!
  //! The copied optimizer will share the identical loss and penalty
  //! functions after construction.
  CDPense& operator=(const CDPense& other) = default;

  //! Default move constructor.
  CDPense(CDPense&& other) = default;

  //! Default move assignment operator.
  CDPense& operator=(CDPense&& other) = default;

  ~CDPense() = default;

  void Reset() {
    loss_.reset();
    penalty_.reset();
    state_.residuals.reset();
  }

  //! Get the convergence tolerance for the CD algorithm.
  //!
  //! @return convergence tolerance.
  double convergence_tolerance() const noexcept {
    return convergence_tolerance_;
  }

  //! Set the convergence tolerance for the CD algorithm.
  //!
  //! @param convergence_tolerance convergene tolerance for the MM algorithm.
  void convergence_tolerance(double convergence_tolerance) noexcept {
    convergence_tolerance_ = convergence_tolerance;
  }

  //! Get the current loss function.
  SLoss& loss() const {
    if (!loss_) {
      throw std::logic_error("no loss set");
    }
    return *loss_;
  }

  //! Set the new loss function.
  void loss(const SLoss& loss) noexcept {
    loss_ = std::make_unique<SLoss>(loss);
    lipschitz_bounds_.reset();
  }

  PenaltyFunction& penalty() const {
    if (!penalty_) {
      throw std::logic_error("no penalty set");
    }
    return *penalty_;
  }

  void penalty(const PenaltyFunction& penalty) noexcept {
    penalty_ = std::make_unique<PenaltyFunction>(penalty);
  }

  //! Find the minimum of the objective function, using the previous solution
  //! (or the 0-vector if no previous solution exists) as starting point.
  //!
  //! @return information about the optimum.
  Optimum Optimize() {
    return Optimize(config_.max_it);
  }

  //! Find the minimum of the objective function, using the given coefficients
  //! as starting point and at most ``max_it`` iterations.
  //!
  //! @param start where to start the optimization from.
  //! @return information about the optimum.
  Optimum Optimize(const Coefficients& start) {
    ResetState(start);
    return Optimize(config_.max_it);
  }

  //! Find the minimum of the objective function, using the given coefficients
  //! as starting point and at most ``max_it`` iterations.
  //!
  //! @param start where to start the optimization from.
  //! @param max_it maximum number of iterations.
  //! @return information about the optimum.
  Optimum Optimize(const Coefficients& start, const int max_it) {
    ResetState(start);
    return Optimize(max_it);
  }

  //! Find the minimum of the objective function, using the previous solution
  //! (or the 0-vector if no previous solution exists) as starting point and
  //! at most ``max_it`` iterations.
  //!
  //! @param max_it maximum number of iterations.
  //! @return information about the optimum.
  Optimum Optimize(const int max_it) {
    if (!loss_) {
      throw std::logic_error("no loss set");
    }
    if (!penalty_) {
      throw std::logic_error("no penalty set");
    }

    auto metrics = std::make_unique<nsoptim::Metrics>("cd-pense");

    if (state_.residuals.n_elem == 0) {
      ResetState(loss_->template ZeroCoefficients<Coefficients>());
    }

    if (lipschitz_bounds_.n_elem == 0) {
      UpdateLipschitzBounds();
    }

    int iter = 0;
    const double stepsize_start_mult = std::pow(config_.linesearch_ss_multiplier, config_.linesearch_ss_num);
    const auto& data = loss_->data();

    // if (!config_.linesearch) {
    //   config_.linesearch_ss_num = 1;
    // }

    while (iter++ < max_it) {
      double coef_change = 0;
      auto& iteration_metrics = metrics->CreateSubMetrics("cd_iteration");
      iteration_metrics.AddMetric("iter", iter);

      const double objf_before_iter = state_.objf_loss + state_.objf_pen;

      if (loss_->IncludeIntercept()) {
        // If we don't do linesearch, also compute the lipschitz constant of the surrogate WLS loss.
        auto gradlip = GradientAndSurrogateLipschitz();
        int total_mscale_iterations = 0;
        double updated_coef = state_.coefs.intercept;

        iteration_metrics.AddMetric("gradient_int", gradlip.gradient);

        if (!config_.linesearch) {
          state_.coefs.intercept -= gradlip.gradient * gradlip.lipschitz_constant;
          state_.residuals += updated_coef - state_.coefs.intercept;

          const auto eval_loss = loss_->EvaluateResiduals(state_.residuals);
          total_mscale_iterations = loss_->mscale().LastIterations();

          state_.objf_loss = eval_loss.loss;
          state_.mscale = eval_loss.scale;
          coef_change += std::abs(state_.coefs.intercept - updated_coef);

          iteration_metrics.AddMetric("ls_stepsize_int", gradlip.lipschitz_constant);
          iteration_metrics.AddMetric("ls_stepsize_int_sloss", 1 / lipschitz_bound_intercept_);
        } else {
          // Start line search with the step in the middle and increase or decrease based on the results.
          double stepsize = lipschitz_bound_intercept_ * stepsize_start_mult;
          int ls_step = -1;
          double stepsize_multiplier = 1;

          while (ls_step++ < config_.linesearch_ss_num) {
            const double try_coef = state_.coefs.intercept - gradlip.gradient / stepsize;
            state_.residuals += updated_coef - try_coef;
            const auto eval_loss = loss_->EvaluateResiduals(state_.residuals);
            total_mscale_iterations += loss_->mscale().LastIterations();

            if (ls_step == 0) {
              stepsize_multiplier = (eval_loss.loss <= state_.objf_loss) ?
                config_.linesearch_ss_multiplier : (1 / config_.linesearch_ss_multiplier);
            }

            if (stepsize_multiplier < 1) {
              if (eval_loss.loss > state_.objf_loss) {
                // Objective did not improve further. Rewind residuals and stop line search.
                state_.residuals -= updated_coef - try_coef;
                loss_->mscale().SetInitial(state_.mscale);
                break;
              }

              // Updating the objective function here ensures that we will go as long as there is more improvement.
              state_.objf_loss = eval_loss.loss;
              state_.mscale = eval_loss.scale;
            } else { // stepsize_multiplier >= 1
              if (eval_loss.loss < state_.objf_loss) {
                // We are finally at a point where the objective function improved.
                // Stop here.
                updated_coef = try_coef;
                state_.objf_loss = eval_loss.loss;
                state_.mscale = eval_loss.scale;
                break;
              } else if (stepsize > lipschitz_bound_intercept_) {
                // We are at the upper end of the step size range and haven't seen an improvement.
                // Stop here.
                ls_step = 0;
                if (std::abs(updated_coef - state_.coefs.intercept) > kNumericZero) {
                  state_.residuals += try_coef - state_.coefs.intercept;
                }
                break;
              }
            }

            updated_coef = try_coef;
            stepsize *= stepsize_multiplier;
          }

          if (ls_step > 0) {
            coef_change += std::abs(state_.coefs.intercept - updated_coef);
            state_.coefs.intercept = updated_coef;
            iteration_metrics.AddMetric("ls_stepsize_int", stepsize / lipschitz_bound_intercept_);
            iteration_metrics.AddMetric("ls_steps_int", ls_step);
          } else {
            // The coefficient value did not change.
            iteration_metrics.AddMetric("ls_stepsize_int", 0.);
            iteration_metrics.AddMetric("ls_steps_int", 0);
          }
        }

        iteration_metrics.AddMetric("mscale_iterations_int", total_mscale_iterations);
      }

      for (arma::uword j = 0; j < data.n_pred(); ++j) {
        // @TODO -- this iteration is inefficient if we have a sparse vector!
        auto& cycle_metrics = iteration_metrics.CreateSubMetrics("coordinate");
        auto gradlip = GradientAndSurrogateLipschitz(j);
        int total_mscale_iterations = 0;
        const double objf_pen_prev = state_.objf_pen - PenaltyContribution(state_.coefs.beta[j], j, IsAdaptiveTag{});
        double updated_coef = state_.coefs.beta[j];

        cycle_metrics.AddMetric("gradient", gradlip.gradient);

        if (!config_.linesearch) {
          state_.coefs.beta[j] = UpdateSlope(j, gradlip.lipschitz_constant, gradlip.gradient, IsAdaptiveTag{});

          if (std::abs(state_.coefs.beta[j] - updated_coef) > kNumericZero) {
            state_.residuals += (updated_coef - state_.coefs.beta[j]) * data.cx().col(j);
            const auto eval_loss = loss_->EvaluateResiduals(state_.residuals);
            total_mscale_iterations = loss_->mscale().LastIterations();

            state_.objf_loss = eval_loss.loss;
            state_.objf_pen = objf_pen_prev + PenaltyContribution(state_.coefs.beta[j], j, IsAdaptiveTag{});
            state_.mscale = eval_loss.scale;
            coef_change += std::abs(state_.coefs.beta[j] - updated_coef);

            cycle_metrics.AddMetric("ls_stepsize", gradlip.lipschitz_constant);
            cycle_metrics.AddMetric("ls_stepsize_sloss", 1 / lipschitz_bounds_[j]);
          }
        } else {
          double stepsize = stepsize_start_mult * lipschitz_bounds_[j];
          int ls_step = -1;
          double stepsize_multiplier = -1;

          cycle_metrics.AddMetric("index", static_cast<int>(j));
          cycle_metrics.AddMetric("gradient", gradlip.gradient);

          while (ls_step++ < config_.linesearch_ss_num) {
            const double try_coef = UpdateSlope(j, stepsize, gradlip.gradient, IsAdaptiveTag{});

            if (std::abs(try_coef - state_.coefs.beta[j]) > kNumericZero) {
              state_.residuals += (updated_coef - try_coef) * data.cx().col(j);

              const auto eval_loss = loss_->EvaluateResiduals(state_.residuals);
              const double new_objf_pen = objf_pen_prev + PenaltyContribution(try_coef, j, IsAdaptiveTag{});
              const double new_objf = eval_loss.loss + new_objf_pen;

              total_mscale_iterations += loss_->mscale().LastIterations();

              // After the first step that produces a change in the coefficient, determine
              // which way we want to go.
              if (stepsize_multiplier < 0) {
                stepsize_multiplier = (new_objf <= state_.objf_loss + state_.objf_pen) ?
                  config_.linesearch_ss_multiplier : (1 / config_.linesearch_ss_multiplier);

                if (stepsize_multiplier > 1 && ls_step > 0) {
                  // There were no changes for previous decreases in the step size,
                  // and this change leads to an increase in the objective function.
                  // Stop the line search.
                  stepsize = lipschitz_bounds_[j] + kNumericZero;
                  break;
                }
              }

              if (stepsize_multiplier < 1) {
                if (new_objf > state_.objf_loss + state_.objf_pen) {
                  // Objective did not improve further. Rewind residuals and stop line search.
                  state_.residuals -= (updated_coef - try_coef) * data.cx().col(j);
                  loss_->mscale().SetInitial(state_.mscale);
                  break;
                }

                // Updating the objective function here ensures that we will go as long as there is more improvement.
                state_.objf_loss = eval_loss.loss;
                state_.objf_pen = new_objf_pen;
                state_.mscale = eval_loss.scale;
              } else {
                if (new_objf < state_.objf_loss + state_.objf_pen) {
                  // We are finally at a point where the objective function improved.
                  // Stop here.
                  updated_coef = try_coef;
                  state_.objf_loss = eval_loss.loss;
                  state_.objf_pen = new_objf_pen;
                  state_.mscale = eval_loss.scale;
                  break;
                } else if (stepsize >= lipschitz_bounds_[j]) {
                  // We are at the upper end of the step size range and haven't seen an improvement.
                  // Stop here.
                  ls_step = 0;
                  state_.residuals += (try_coef - state_.coefs.beta[j]) * data.cx().col(j);
                  break;
                }
              }

              updated_coef = try_coef;
              stepsize *= stepsize_multiplier;
            } else if (stepsize_multiplier < 0) {
              // No change in the coefficient and there hasn't been any change in the objective
              // function yet.
              // Reduce the step size.
              stepsize *= config_.linesearch_ss_multiplier;
            } else { // stepsize_multiplier > 1 (as < 1 must lead to a change if it did before)
              break;
            }
          }

          if (ls_step > 0 && std::abs(updated_coef - state_.coefs.beta[j]) > kNumericZero) {
            coef_change += std::abs(state_.coefs.beta[j] - updated_coef);
            state_.coefs.beta[j] = updated_coef;
            cycle_metrics.AddMetric("ls_stepsize", stepsize / lipschitz_bounds_[j]);
            cycle_metrics.AddMetric("ls_steps", ls_step);
          } else {
            // The coefficient value did not change.
            cycle_metrics.AddMetric("ls_steps", 0);
            cycle_metrics.AddMetric("ls_stepsize", 0.);
          }
        }

        cycle_metrics.AddMetric("mscale_iterations", total_mscale_iterations);
      }

      const double objf_change = (state_.objf_loss + state_.objf_pen) - objf_before_iter;
      iteration_metrics.AddMetric("change", objf_change);
      iteration_metrics.AddMetric("coef_change", coef_change);

      if (objf_change * objf_change < convergence_tolerance_ * convergence_tolerance_) {
        // The objective function value did not change. Algorithm converged.
        metrics->AddMetric("iter", iter);
        return nsoptim::MakeOptimum(*loss_, *penalty_, state_.coefs, state_.residuals,
                                    std::move(metrics));
      }

      // Re-compute the residuals after every few cycles to avoid any drifts
      if (iter > 0 && iter % config_.reset_iter == 0) {
        state_.residuals = loss_->Residuals(state_.coefs);
      }
    }

    metrics->AddMetric("iter", iter);
    state_.residuals = loss_->Residuals(state_.coefs);
    return nsoptim::MakeOptimum(*loss_, *penalty_, state_.coefs, state_.residuals,
                                std::move(metrics), nsoptim::OptimumStatus::kWarning,
                                "Coordinate descent did not converge.");
  }

 private:
  void UpdateLipschitzBounds() {
    const auto data = loss_->data();
    const auto ms = loss_->mscale();
    const double eff_n = data.n_obs() * (1. - ms.delta());
    const double separation = eff_n - std::floor(eff_n);
    const double mult = std::log(separation * (1 - separation)) / std::cbrt(eff_n);
    const double u1 = std::min(80., -40. * mult) / ms.rho().cc();
    const double u2 = std::min(50., 100. * mult * mult * mult * mult) / ms.rho().cc();
    lipschitz_bounds_ = u1 * u1 * arma::square(arma::sum(data.cx())).t();
    for (arma::uword j = 0; j < data.n_pred(); ++j) {
      lipschitz_bounds_[j] += u2 * state_.mscale * std::abs(arma::accu(data.cx().col(j) * data.cx().col(j).t()));
    }

    lipschitz_bound_intercept_ = (u1 * u1 + u2 * state_.mscale) * data.n_obs() * data.n_obs();
  }

  coorddesc::SurrogateGradient GradientAndSurrogateLipschitz() {
    const arma::vec wgt = loss_->mscale().rho().Weight(state_.residuals, state_.mscale);
    const double gradient = -state_.mscale * state_.mscale * arma::dot(wgt, state_.residuals) /
      arma::dot(wgt, arma::square(state_.residuals));
    const double lipschitz = arma::mean(wgt);
    return coorddesc::SurrogateGradient { gradient, lipschitz };
  }

  coorddesc::SurrogateGradient GradientAndSurrogateLipschitz(const arma::uword j) {
    auto&& xmat = loss_->data().cx();
    const arma::vec wgt = loss_->mscale().rho().Weight(state_.residuals, state_.mscale);
    const double gradient = -state_.mscale * state_.mscale * arma::dot(wgt % xmat.col(j), state_.residuals) /
      arma::dot(wgt, arma::square(state_.residuals));
    const double lipschitz = arma::mean(wgt % arma::square(xmat.col(j)));
    return coorddesc::SurrogateGradient { gradient, lipschitz };
  }

  double Gradient() {
    arma::vec psi = loss_->mscale().rho().Derivative(state_.residuals, state_.mscale);
    return -state_.mscale * state_.mscale *
      arma::accu(psi) / (arma::dot(psi, state_.residuals));
  }

  double Gradient(const arma::uword j) {
    auto&& xmat = loss_->data().cx();
    arma::vec psi = loss_->mscale().rho().Derivative(state_.residuals, state_.mscale);
    return -state_.mscale * state_.mscale *
      arma::dot(psi, xmat.col(j)) / (arma::dot(psi, state_.residuals));
  }

  double UpdateSlope (const arma::uword j, const double stepsize, const double gradient,
                      std::false_type /* is_adaptive */) {
    const double dir = stepsize * state_.coefs.beta[j] - gradient;
    return nsoptim::SoftThreshold(dir, penalty_->lambda() * penalty_->alpha()) /
      (stepsize + penalty_->lambda() * (1 - penalty_->alpha()));
  }

  double UpdateSlope (const arma::uword j, const double stepsize, const double gradient,
                      std::true_type /* is_adaptive */) {
    const double dir = stepsize * state_.coefs.beta[j] - gradient;
    const double penalty_level = penalty_->loadings()[j] * penalty_->lambda();
    return nsoptim::SoftThreshold(dir, penalty_level * penalty_->alpha()) /
      (stepsize + penalty_level * (1 - penalty_->alpha()));
  }

  void ResetState (const Coefficients &coefs) {
    if (!loss_) {
      throw std::logic_error("no loss set");
    }
    if (!penalty_) {
      throw std::logic_error("no penalty set");
    }
    state_ = { coefs, loss_->Residuals(coefs), 0, 0, penalty_->Evaluate(coefs) };
    auto loss_eval = loss_->EvaluateResiduals(state_.residuals);
    state_.mscale = loss_eval.scale;
    state_.objf_loss = loss_eval.loss;
  }

  double PenaltyContribution(const double value, const arma::uword, std::false_type /* is_adaptive */) {
    return penalty_->lambda() * (penalty_->alpha() * std::abs(value) +
      0.5 * (1 - penalty_->alpha()) * value * value);
  }

  double PenaltyContribution(const double value, const arma::uword j, std::true_type /* is_adaptive */) {
    return penalty_->loadings()[j] * penalty_->lambda() *
      (penalty_->alpha() * std::abs(value) +
        0.5 * (1 - penalty_->alpha()) * value * value);
  }

  LossFunctionPtr loss_;
  PenaltyPtr penalty_;
  CDPenseConfiguration config_;
  arma::vec lipschitz_bounds_;
  double lipschitz_bound_intercept_;
  coorddesc::State<Coefficients> state_;
  double convergence_tolerance_ = kDefaultConvergenceTolerance;
};
} // namespace pense

#endif // PENSE_CD_PENSE_HPP_
