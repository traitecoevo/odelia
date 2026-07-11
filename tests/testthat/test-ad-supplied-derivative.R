# supplied_derivative supplies the known derivative of an off-tape result to the
# reverse sweep. Differentiating g = x(a)^2, where x(a) solves x = cos(x) + a off the
# tape, must match the closed form 2*x*dx/da and finite differences of the forward
# solve -- so the known adjoint is carried through without the root-find being recorded.

# Forward-only reference: the same Newton solve of x = cos(x) + a, root squared.
g_of_a <- function(a) {
  x <- a
  for (i in 1:100) {
    f <- x - cos(x) - a
    dx <- f / (1 + sin(x))
    x <- x - dx
    if (abs(dx) < 1e-15) break
  }
  x * x
}

test_that("supplied_derivative injects the IFT adjoint", {
  testthat::skip_if(is_pkgload_dll(), "Skipping AD workflow in pkgload load_all sessions due unstable native-pointer lifecycle.")
  ensure_supplied_derivative_interface()

  for (a in c(0.3, 1.0, 2.5)) {
    res <- supplied_derivative_demo(a)

    # Value: g = x^2 reproduces the forward solve.
    expect_equal(res$g, g_of_a(a), tolerance = 1e-10)

    # Edge-propagated gradient matches the closed form ...
    expect_equal(res$dg_da, res$dg_da_analytic, tolerance = 1e-12)

    # ... and central finite differences of the forward solve.
    eps <- 1e-6
    fd <- (g_of_a(a + eps) - g_of_a(a - eps)) / (2 * eps)
    expect_equal(res$dg_da, fd, tolerance = 1e-6)
  }
})
