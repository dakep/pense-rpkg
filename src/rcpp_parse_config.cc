//
//  rcpp_parse_config.cc
//  pense
//
//  Created by David Kepplinger on 2019-05-01.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//
#include "nsoptim.hpp"

#include "rcpp_parse_config.hpp"
#include "rcpp_utils.hpp"

#include "nsoptim.hpp"

namespace {
constexpr int kAdmmMaxIt = 1000;
constexpr double kAdmmAcceleration = 1;

constexpr int kDalMaxIt = 100;
constexpr int kDalMaxInnerIt = 100;
constexpr double kDalEtaMult = 2;
constexpr double kDalEtaStartNumeratorCons = 0.01;
constexpr double kDalEtaStartNumeratorAggr = 1;
constexpr double kDalLambdaRelChangeAggr = 0.25;

constexpr int kMmMaxIt = 500;
constexpr nsoptim::MMConfiguration::TighteningType kMmTightening = nsoptim::MMConfiguration::TighteningType::kAdaptive;
constexpr int kMmTighteningSteps = 10;
}  // namespace

namespace Rcpp {
namespace traits {

nsoptim::AdmmLinearConfiguration Exporter<nsoptim::AdmmLinearConfiguration>::get() const {
  const Rcpp::List config_list = as<const Rcpp::List>(r_obj_);
  nsoptim::AdmmLinearConfiguration tmp = {
    pense::GetFallback(config_list, "max_it", kAdmmMaxIt),
    pense::GetFallback(config_list, "accelerate", kAdmmAcceleration)
  };
  return tmp;
}

nsoptim::DalEnConfiguration Exporter<nsoptim::DalEnConfiguration>::get() const {
  const Rcpp::List config_list = as<const Rcpp::List>(r_obj_);
  nsoptim::DalEnConfiguration tmp = {
      pense::GetFallback(config_list, "max_it", kDalMaxIt),
      pense::GetFallback(config_list, "max_inner_it", kDalMaxInnerIt),
      pense::GetFallback(config_list, "eta_start_numerator_conservative", kDalEtaStartNumeratorCons),
      pense::GetFallback(config_list, "eta_start_numerator_aggressive", kDalEtaStartNumeratorAggr),
      pense::GetFallback(config_list, "lambda_relchange_aggressive", kDalLambdaRelChangeAggr),
      pense::GetFallback(config_list, "eta_multiplier", kDalEtaMult)
  };
  return tmp;
}

nsoptim::MMConfiguration Exporter<nsoptim::MMConfiguration>::get() const {
  const Rcpp::List config_list = as<const Rcpp::List>(r_obj_);
  nsoptim::MMConfiguration tmp = {
      pense::GetFallback(config_list, "max_it", kMmMaxIt),
      pense::GetFallback(config_list, "tightening", kMmTightening),
      pense::GetFallback(config_list, "tightening_steps", kMmTighteningSteps)
  };
  return tmp;
}

}  // namespace traits
}  // namespace Rcpp
