#' Compare estimates to a reference.
compare_estimates <- function (ests, ref_file, tol = 1e-6) {
  requireNamespace('testthat')
  testthat::skip_if_not_installed('jsonlite')

  requireNamespace('jsonlite')

  snapshot_file <- testthat::test_path(ref_file)
  if (!file.exists(snapshot_file)) {
    jsonlite::write_json(ests, path = snapshot_file, auto_unbox = TRUE,
                         digits = 14, pretty = TRUE)
    testthat::skip('Snapshot file did not exist and was created.')
  }

  ref <- jsonlite::read_json(snapshot_file, simplifyVector = FALSE)

  for (lai in seq_along(ref)) {
    for (lsi in seq_along(ref[[lai]])) {
      for (name in names(ref[[lai]][[lsi]])) {
        testthat::expect_success(
          testthat::expect_equal(
            drop(ests [[!!lai]] [[!!lsi]] [[!!name]]),
            unlist(ref [[!!lai]] [[!!lsi]] [[!!name]]),
            tolerance = tol))
      }
    }
  }
}
