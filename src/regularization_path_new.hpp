//
//  regularization_path.hpp
//  pense
//
//  Created by David Kepplinger on 2019-01-30.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef REGULARIZATION_PATH_NEW_HPP_
#define REGULARIZATION_PATH_NEW_HPP_

#include <memory>
#include <tuple>
#include <type_traits>

#include "nsoptim.hpp"

#include "alias.hpp"
#include "m_loss.hpp"
#include "omp_utils.hpp"

namespace pense {
namespace regpath {

//! A list of optima and associated optimizers, retaining only optima which are
//!
//! The implementation uses a single-linked list of optima, with the worst optimum at the front.
template<class Optimizer, typename... Ts>
class UniqueOptima {
  using LossFunction = typename Optimizer::LossFunction;
  using PenaltyFunction = typename Optimizer::PenaltyFunction;
  using Coefficients = typename Optimizer::Coefficients;
  using Optimum = typename Optimizer::Optimum;
  using Element = std::tuple<Optimum, Optimizer, Ts...>;

 public:
  enum class InsertResult { kGood, kBad, kDuplicate };

  //! Create a new list of optima with at most `max_size` elements.
  //!
  //! @param max_size maximum number of distinct elements contained in the list
  //! @param eps numerical tolerance for determining equivalence
  UniqueOptima(const size_t max_size, const double eps) noexcept
      : max_size_(max_size), eps_(eps), size_(0), optima_() {}

  //! Default copy constructor. Copy assignment is not possible.
  UniqueOptima(const UniqueOptima&) = default;
  UniqueOptima& operator=(const UniqueOptima&) = delete;

  //! Move constructor.
  UniqueOptima(UniqueOptima&& other) noexcept :
      max_size_(other.max_size_), eps_(other.eps_), size_(other.size_),
      optima_(std::move(other.optima_)) {
    other.size_ = 0;
  }
  UniqueOptima& operator=(UniqueOptima&&) = delete;

  //! Check if an optimum is good enough to be potentially inserted into the container.
  bool IsGoodEnough(const Optimum& optimum) {
    // Check if the optimum's objective function value is good enough.
    if (size_ >= max_size_ && optimum.objf_value > std::get<0>(optima_.front()).objf_value) {
      return false;
    }
    return true;
  }

  //! Insert an optimum and associated information into the container,
  //! if the optimum is good enough.
  template<typename T, typename... Args>
  InsertResult Insert(const Optimizer& optimizer, T&& optimum, Args&&... args) {
    // Ensure that `T` is of type `Optimum`.
    static_assert(std::is_same<typename std::decay<T>::type, Optimum>::value,
                 "Optimum is of wrong type.");

    // Check if the optimum's objective function value is good enough.
    if (size_ == max_size_ && optimum.objf_value > std::get<0>(optima_.front()).objf_value) {
      return InsertResult::kBad;
    }

    // Determine insert position.
    const auto optima_end = optima_.end();
    const auto optima_before_begin = optima_.before_begin();
    auto after_it = optima_.begin();
    auto insert_it = optima_before_begin;
    while (after_it != optima_end) {
      // Check if the two optima are equal.
      if (Equal(optimum, std::get<0>(*after_it))) {
        // Don't add.
        return InsertResult::kDuplicate;
      }

      // Check that the objective value of the given optimum is larger than the next optimum but smaller than
      // the current optimum.
      if ((optimum.objf_value > std::get<0>(*after_it).objf_value) &&
          (insert_it == optima_before_begin || std::get<0>(*insert_it).objf_value > optimum.objf_value)) {
        // Insert here.
        optima_.emplace_after(insert_it, std::forward<T>(optimum),
          std::forward<Optimizer>(optimizer), std::forward<Args>(args)...);

        // Ensure that the size keeps within the limits.
        if (++size_ > max_size_) {
          optima_.erase_after(optima_.before_begin());
          --size_;
        }
        return InsertResult::kGood;
      }
      ++after_it;
      ++insert_it;
    }

    // All other optima have larger objective function value than the given optimum. Add at the end.
    optima_.emplace_after(insert_it, std::forward<T>(optimum),
      std::forward<Optimizer>(optimizer), std::forward<Args>(args)...);
    // Ensure that the size stays within the limits.
    if (++size_ > max_size_) {
      optima_.erase_after(optima_.before_begin());
      --size_;
    }
    return InsertResult::kGood;
  }

  const alias::FwdList<Element>& Elements() const noexcept {
    return optima_;
  }

  alias::FwdList<Element>& Elements() noexcept {
    return optima_;
  }

  template<size_t I>
  alias::FwdList<typename std::tuple_element<I, Element>::type > Elements() const noexcept {
    alias::FwdList<typename std::tuple_element<I, Element>::type > out;

    for (auto&& el : optima_) {
      out.push_front(std::get<I>(el));
    }

    return out;
  }

 private:
  const size_t max_size_;
  const double eps_;
  size_t size_;
  alias::FwdList<Element> optima_;

  bool Equal(const Optimum& a, const Optimum& b) const noexcept {
    // First check if the value of the objective function is similar.
    if (std::abs(a.objf_value - b.objf_value) < eps_) {
      // Value of the objective function is similar. Check if the intercept is similar.
      const double int_diff = a.coefs.intercept - b.coefs.intercept;
      if (int_diff * int_diff < eps_) {
        // The intercept is similar. Check if the slope is also similar.
        const double beta_diff = arma::norm(a.coefs.beta - b.coefs.beta, 2);
        if (beta_diff * beta_diff < eps_) {
          // The slope is also similar. Return true.
          return true;
        }
      }
    }
    return false;
  }
};
} // namespace regpath

template<class Optimizer>
class RegularizationPath {
  using LossFunction = typename Optimizer::LossFunction;
  using PenaltyFunction = typename Optimizer::PenaltyFunction;
  using Coefficients = typename Optimizer::Coefficients;
  using PenaltyList = alias::FwdList<PenaltyFunction>;
  using Optimum = typename Optimizer::Optimum;
  using StartCoefficients = alias::FwdList<alias::FwdList<Coefficients>>;

 public:

  //! Add a starting point to be used only at the specified penalty.
  //!
  //! @param penalty penalty at which the starting point should be used.
  //! @param coefs starting point.
  void AddStartingPoint(const PenaltyFunction& penalty, const Coefficients& coefs) {}

  //! Add a starting point to be used for all penalties.
  //!
  //! @param coefs starting point.
  void AddStartingPoint(const Coefficients& coefs) {}

 private:
  void Explore()
};

} // namespace pense

#endif // REGULARIZATION_PATH_NEW_HPP_
