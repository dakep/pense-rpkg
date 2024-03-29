//
//  armadillo_forward.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-11-30.
//  Copyright © 2018 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_ARMADILLO_FORWARD_HPP_
#define NSOPTIM_ARMADILLO_FORWARD_HPP_

#define ARMA_USE_CXX11 1
#define ARMA_DONT_USE_OPENMP 1

#ifdef __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Weverything"
#endif

#ifdef HAVE_RCPP
// For RCPP
# include <RcppArmadillo/interface/RcppArmadilloForward.h>
#else
// For stand-alone
# include <armadillo>
#endif  // HAVE_RCPP

#ifdef __clang__
#  pragma clang diagnostic pop
#endif

#endif  // NSOPTIM_ARMADILLO_FORWARD_HPP_
