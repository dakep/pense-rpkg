if (require(testthat)) {
  library(pense)
  test_check("pense")
} else {
  warning("'pense' requires 'testthat' for tests.")
}
