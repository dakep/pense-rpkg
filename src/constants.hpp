//
//  constants.hpp
//  pense
//
//  Created by David Kepplinger on 2019-01-30.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef CONSTANTS_HPP_
#define CONSTANTS_HPP_

namespace pense {

//! Default (fallback) convergence tolerance for numeric algorithms.
constexpr double kDefaultConvergenceTolerance = 1e-6;

//! The threshold for any numeric value to be considered 0.
constexpr double kNumericZero = 1e-12;

//! Integer IDs for the supported rho-functions
enum class RhoFunctionType {
  kRhoHuber = 1,
  kRhoBisquare = 2,
  kRhoOptimal = 3
};

//! Integer IDs for supported EN algorithms
enum class EnAlgorithm {
  kLinearizedAdmm = 1,
  kVarStepAdmm = 2,
  kDal = 3,
  kRidge = 4,
  kLars = 5,
  kCoordinateDescent = 6
};

//! Integer IDs for supported EN algorithms
enum class PenseAlgorithm {
  kMm = 1,
  kAdmm = 2,
  kCoordinateDescent = 3
};

//! Integer IDs for supported EN algorithms
enum class MestEnAlgorithm {
  kMm = 1
};

//! Default tuning constant for the rho function for M-scale estimates.
constexpr double kDefaultHuberCc = 1.345;
constexpr double kDefaultBisquareCc = 2.937015;
constexpr double kDefaultOptimaleCc = 0.405;
//! Default breakdown point for the M-scale equation.
constexpr double kDefaultMscaleDelta = 0.25;

constexpr EnAlgorithm kDefaultEnAlgorithm = EnAlgorithm::kLars;
constexpr PenseAlgorithm kDefaultPenseAlgorithm = PenseAlgorithm::kMm;
constexpr MestEnAlgorithm kDefaultMestAlgorithm = MestEnAlgorithm::kMm;
constexpr bool kDefaultUseSparse = false;

}  // namespace pense

#endif  // CONSTANTS_HPP_
