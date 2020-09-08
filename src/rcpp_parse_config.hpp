//
//  rcpp_parse_config.hpp
//  pense
//
//  Created by David Kepplinger on 2019-01-30.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef RCPP_PARSE_CONFIG_HPP_
#define RCPP_PARSE_CONFIG_HPP_

#include "nsoptim_forward.hpp"

namespace Rcpp {
namespace traits {
//! Converter for an R-list to configuration options for the linearized ADMM algorithm.
template<> class Exporter< nsoptim::AdmmLinearConfiguration > {
 public:
  explicit Exporter(SEXP r_obj) noexcept : r_obj_(r_obj) {}
  nsoptim::AdmmLinearConfiguration get() const;
 private:
  SEXP r_obj_;
};

//! Converter for an R-list to configuration options for the DAL algorithm.
template<> class Exporter< nsoptim::DalEnConfiguration > {
 public:
  explicit Exporter(SEXP r_obj) noexcept : r_obj_(r_obj) {}
  nsoptim::DalEnConfiguration get() const;
 private:
  SEXP r_obj_;
};

//! Converter for an R-list to configuration options for the MM algorithm.
template<> class Exporter< nsoptim::MMConfiguration > {
 public:
  explicit Exporter(SEXP r_obj) noexcept : r_obj_(r_obj) {}
  nsoptim::MMConfiguration get() const;
 private:
  SEXP r_obj_;
};

}  // namespace traits
}  // namespace Rcpp

#endif  // RCPP_PARSE_CONFIG_HPP_
