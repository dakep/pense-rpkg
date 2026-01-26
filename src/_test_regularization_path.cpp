#ifndef TESTTHAT_DISABLED

#include <testthat.h>

#include "nsoptim.hpp"
#include "regularization_path.hpp"

// -----------------------------------------------------------------------------
// Mock Classes to satisfy template requirements
// -----------------------------------------------------------------------------

struct MockCoefficients {
  double intercept;
  arma::vec beta;

  MockCoefficients(double i, arma::vec b) : intercept(i), beta(b) {}

  // Default constructor for containers
  MockCoefficients() : intercept(0.0), beta() {}
};

struct MockOptimum {
  MockCoefficients coefs;
  double objf_value;

  MockOptimum(MockCoefficients c, double obj) : coefs(c), objf_value(obj) {}
  MockOptimum() : coefs(), objf_value(0.0) {}
};

// Optimizer acts as a traits class
struct MockOptimizer {
  using Coefficients = MockCoefficients;
  using Optimum = MockOptimum;
};

// -----------------------------------------------------------------------------
// Tests
// -----------------------------------------------------------------------------

context("OrderedTuples Container Tests") {

  test_that("UniqueCoefficients detects duplicates correctly") {
    // UniqueCoefficients uses DuplicateCoefficients ordering (no objf comparison)
    // It simply stores unique sets of coefficients.
    using Container = pense::regpath::UniqueCoefficients<MockCoefficients>;

    double tol = 1e-5;
    Container list(pense::regpath::DuplicateCoefficients<MockCoefficients>(tol));

    arma::vec beta1 = {1.0, 2.0};
    arma::vec beta2 = {1.0, 3.0}; // Different beta

    // 1. Insert first element
    auto res1 = list.Emplace(MockCoefficients(0.5, beta1));
    expect_true(res1 == Container::InsertResult::kGood);
    expect_true(list.Size() == 1);

    // 2. Insert exact duplicate
    auto res2 = list.Emplace(MockCoefficients(0.5, beta1));
    expect_true(res2 == Container::InsertResult::kDuplicate);
    expect_true(list.Size() == 1);

    // 3. Insert approximate duplicate (within tolerance)
    // Difference is 1e-6, tolerance is 1e-5
    auto res3 = list.Emplace(MockCoefficients(0.5 + 1e-6, beta1));
    expect_true(res3 == Container::InsertResult::kDuplicate);
    expect_true(list.Size() == 1);

    // 4. Insert different element (intercept)
    auto res4 = list.Emplace(MockCoefficients(0.8, beta1));
    expect_true(res4 == Container::InsertResult::kGood);
    expect_true(list.Size() == 2);

    // 5. Insert different element (beta)
    auto res5 = list.Emplace(MockCoefficients(0.5, beta2));
    expect_true(res5 == Container::InsertResult::kGood);
    expect_true(list.Size() == 3);
  }

  test_that("UniqueOptima maintains order (Worst -> Best) and respects capacity") {
    // UniqueOptima uses OptimaOrder.
    // Logic analysis:
    // It iterates through the list. If existing_item.objf < new_item.objf (existing is better),
    // it returns kLowerObjf. The code breaks and inserts BEFORE existing.
    // This implies the list is ordered from [Worst Objective] -> [Best Objective].
    // When size > max_size, it removes the front (the Worst).

    using Container = pense::regpath::UniqueOptima<MockOptimizer>;

    size_t max_size = 3;
    double tol = 1e-5;
    pense::regpath::OptimaOrder<MockOptimizer> ordering(tol);
    Container list(max_size, ordering);

    arma::vec b = {1.0};

    // 1. Insert Value 10 (Worst)
    list.Emplace(MockOptimum(MockCoefficients(0, b), 10.0), MockOptimizer());
    expect_true(list.Size() == 1);
    expect_true(std::get<0>(list.Elements().front()).objf_value == 10.0);

    // 2. Insert Value 5 (Better)
    // Should be inserted AFTER 10. List: [10, 5]
    list.Emplace(MockOptimum(MockCoefficients(0, b), 5.0), MockOptimizer());
    expect_true(list.Size() == 2);
    expect_true(std::get<0>(list.Elements().front()).objf_value == 10.0); // Front is still worst

    // 3. Insert Value 1 (Best)
    // Should be inserted AFTER 5. List: [10, 5, 1]
    list.Emplace(MockOptimum(MockCoefficients(0, b), 1.0), MockOptimizer());
    expect_true(list.Size() == 3);

    // Verify order by iterating
    auto it = list.Elements().begin();
    expect_true(std::get<0>(*it).objf_value == 10.0); it++;
    expect_true(std::get<0>(*it).objf_value == 5.0); it++;
    expect_true(std::get<0>(*it).objf_value == 1.0);

    // 4. Capacity Test: Insert Value 3 (Intermediate)
    // Logic: 10 vs 3 -> 10>3 (Existing worse). Continue.
    //        5 vs 3  -> 5>3 (Existing worse). Continue.
    //        1 vs 3  -> 1<3 (Existing better). Break. Insert before 1.
    // Temp List: [10, 5, 3, 1]. Size 4.
    // Resize: Remove Front (10).
    // Final List: [5, 3, 1].
    list.Emplace(MockOptimum(MockCoefficients(0, b), 3.0), MockOptimizer());

    expect_true(list.Size() == 3);
    expect_true(std::get<0>(list.Elements().front()).objf_value == 5.0); // New worst

    // 5. Capacity Test: Insert Value 20 (Terrible)
    // Logic: 5 vs 20 -> 5<20 (Existing better). Break. Insert before 5.
    // Temp List: [20, 5, 3, 1].
    // Resize: Remove Front (20).
    // Final List: [5, 3, 1]. (No change effectively)
    auto res = list.Emplace(MockOptimum(MockCoefficients(0, b), 20.0), MockOptimizer());

    // The implementation might return kGood even if it effectively dropped it,
    // or it might rely on the pre-check `CompareObjf(elements_.front(), ...)`
    // logic inside Emplace implies if max_size is reached, check against front.
    // 5 vs 20 -> 5 < 20 (kLowerObjf). Returns kBad immediately.
    expect_true(res == Container::InsertResult::kBad);
    expect_true(std::get<0>(list.Elements().front()).objf_value == 5.0);
  }

  test_that("OrderedTuples handles replacements for better solutions") {
    using Container = pense::regpath::UniqueOptima<MockOptimizer>;

    double tol = 1e-5;
    pense::regpath::OptimaOrder<MockOptimizer> ordering(tol);
    Container list(10, ordering);

    arma::vec beta_a = {1.0, 1.0};
    arma::vec beta_b = {2.0, 2.0};

    // Insert Solution A with Obj 10
    list.Emplace(MockOptimum(MockCoefficients(0, beta_a), 10.0), MockOptimizer());

    // Insert Solution B with Obj 8
    list.Emplace(MockOptimum(MockCoefficients(0, beta_b), 8.0), MockOptimizer());

    // Current: [10 (A), 8 (B)]

    // 1. Try insert A again with Obj 10 (Duplicate)
    auto res1 = list.Emplace(MockOptimum(MockCoefficients(0, beta_a), 10.0), MockOptimizer());
    expect_true(res1 == Container::InsertResult::kDuplicate);
    expect_true(list.Size() == 2);

    // 2. Try insert A again with Obj 9 (Better than original A=10)
    // Logic: Compare 10 vs 9.
    // Objf comparison: 10 > 9 (SlightlyHigherObjf or HigherObjf).
    // The code checks: `equivalent(existing, new)` -> True (same coefs).
    // Then `if (objf_comparison > kEqualObjf)` -> True (10 > 9).
    // It should replace.
    auto res2 = list.Emplace(MockOptimum(MockCoefficients(0, beta_a), 9.0), MockOptimizer());

    expect_true(res2 == Container::InsertResult::kGood);
    expect_true(list.Size() == 2);
    // Verify the value of A is now 9.0 (Front of list)
    expect_true(std::get<0>(list.Elements().front()).objf_value == 9.0);

    // 3. Try insert A again with Obj 9.5 (Worse than current A=9)
    // Logic: Compare 9 vs 9.5.
    // Objf comparison: 9 < 9.5 (SlightlyLower).
    // Breaks loop (wants to insert before).
    // Wait, if it breaks loop, it inserts a duplicate?
    // The loop break condition is `objf_comparison < kSlightlyLowerObjf`.
    // If it is "SlightlyLower" or "Equal" or "SlightlyHigher", it enters the `else if`.
    // 9 vs 9.5 is kSlightlyLowerObjf (-1).
    // It does NOT break. It enters `else if`.
    // Checks `equivalent` -> True.
    // Checks `objf_comparison > kEqualObjf` (-1 > 0)? False.
    // Checks `objf_comparison < kEqualObjf` (-1 < 0)? True.
    // Breaks.
    // Inserts at `insert_it`.
    // BUT this would result in two copies of A!
    // However, looking at the header provided:
    // If equivalent is true, and new is worse (objf_comparison < kEqualObjf),
    // it executes `break`.
    // This implies it inserts the worse version alongside the better version?
    // Let's check the test expectation.
    // Usually, we don't want to insert a worse version of the same coefficients.

    // Testing specific logic in header:
    // "current element is not equivalent and has slightly better objective function. Insert here."
    // That comment is inside the `if(equivalent)` block? No, inside the else of `if(equivalent)`.
    // Wait, the code provided has:
    // } else if (objf_comparison < kEqualObjf) { break; }
    // This `else if` is inside the `if (equivalent)` block in the provided code snippet?
    // No, matching braces:
    // if (equivalent) {
    //    if (better) { replace } else { return kDuplicate }
    // } else if (objf_comparison < kEqualObjf) { break }

    // Ah, if equivalent is true:
    //    If existing is worse (comparison > 0), Replace.
    //    Else (existing is better or equal), return kDuplicate.

    // So inserting 9.5 when 9.0 exists should return kDuplicate.
    auto res3 = list.Emplace(MockOptimum(MockCoefficients(0, beta_a), 9.5), MockOptimizer());
    expect_true(res3 == Container::InsertResult::kDuplicate);
    expect_true(std::get<0>(list.Elements().front()).objf_value == 9.0);
  }
}

#endif // TESTTHAT_DISABLED
