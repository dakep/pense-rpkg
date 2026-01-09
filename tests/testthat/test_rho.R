library(pense)
library(testthat)

test_that("Rho functions", {
  skip_if_not_installed('robustbase')
  requireNamespace("robustbase")
  cc_seq <- c(0.5, 1, 1.5, 2)
  scale_seq <- c(0.5, 1, 1.5, 2)
  # scale_seq <- 1
  x <- seq(-5, 5, length.out = 1001)

  rhofun <- function (x, deriv, std, scale, cc, rho) {
    .Call(pense:::C_rho_fun, x, deriv, std, scale, list(cc = cc, rho = rho_function(rho))) |>
      drop()
  }

  for (rho in rho_function()) {
    rho_robustbase <- switch(rho,
                             mopt = 'optimal',
                             rho)

    for (cc in cc_seq) {
      for (scale in scale_seq) {
        for (deriv in 0L:2L) {
          Mchi_mult <- switch(1L + deriv, 1, scale, 1)
          expect_equal(
            rhofun(x, deriv = !!deriv, std = TRUE, scale = !!scale, cc = !!cc, rho = !!rho),
            Mchi_mult * robustbase::Mchi(x / scale, cc = cc, psi = rho_robustbase, deriv = deriv),
            tolerance = 1e-8)
          if (deriv > 0 || !identical(rho, 'huber')) {
            expect_equal(
              rhofun(x, deriv = !!deriv, std = FALSE, scale = !!scale, cc = !!cc, rho = !!rho),
              Mchi_mult * robustbase::Mpsi(x / scale, cc = cc, psi = rho_robustbase, deriv = deriv - 1L),
              tolerance = 1e-8)
          }
        }
      }
    }
  }
})

