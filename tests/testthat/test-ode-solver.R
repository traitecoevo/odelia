# The solver over an R right-hand side (#62): OdeSolver, ode_solve() and
# domain_error(). The solver's own behaviour -- call budget, Jacobian reuse,
# domain rejection, resizing -- is proven in R-free C++ in tests/standalone/;
# these tests prove the R adapter and the front end against the same oracles.

lorenz_r_rhs <- function(t, y, p) {
  list(c(p[["sigma"]] * (y[[2]] - y[[1]]),
         p[["R"]] * y[[1]] - y[[2]] - y[[1]] * y[[3]],
         -p[["b"]] * y[[3]] + y[[1]] * y[[2]]))
}
lorenz_r_jac <- function(t, y, p) {
  matrix(c(-p[["sigma"]], p[["R"]] - y[[3]], y[[2]],
           p[["sigma"]], -1, y[[1]],
           0, -y[[1]], -p[["b"]]), 3, 3)
}
lorenz_pars <- c(sigma = 10.0, R = 28.0, b = 8.0 / 3.0)

compiled_lorenz <- function(times, method, tol = 1e-10) {
  ensure_ode_interface_loaded()
  ctrl <- odelia:::OdeControl$new()
  ctrl$set_tol_rel(tol)
  ctrl$set_tol_abs(tol)
  lz <- odelia:::LorenzSystem$new(lorenz_pars[["sigma"]], lorenz_pars[["R"]], lorenz_pars[["b"]])
  lz$set_state(c(1, 1, 1), 0.0)
  runner <- odelia:::Lorenz_Solver$new(lz$ptr, ctrl$ptr, active = FALSE, method = method)
  runner$advance_adaptive(times)
  runner$history()
}

testthat::test_that("ode_solve on Lorenz matches the compiled solver, with every stepper", {
  times <- seq(0, 2, by = 0.05)
  ref <- compiled_lorenz(times, "rodas")
  for (method in c("dopri", "rodas", "rkck")) {
    out <- ode_solve(lorenz_r_rhs, y0 = c(x = 1, y = 1, z = 1), times = times,
                     parms = lorenz_pars, method = method, rtol = 1e-10, atol = 1e-10,
                     autonomous = TRUE)
    expect_equal(dim(out), c(length(times), 4))
    expect_equal(colnames(out), c("time", "x", "y", "z"))
    expect_equal(out[, "time"], times)
    # The finite-difference Jacobian is not the AD one, so this is agreement
    # between two accurate integrators, not bit identity.
    expect_equal(out[, "x"], ref$x, tolerance = 1e-6)
    expect_equal(out[, "y"], ref$y, tolerance = 1e-6)
    expect_equal(out[, "z"], ref$z, tolerance = 1e-6)
    counts <- attr(out, "counts")
    expect_named(counts, c("n_rhs", "n_jac", "n_steps", "n_rejections"))
    expect_gt(counts$n_rhs, 0)
  }
})

testthat::test_that("ode_solve on Lorenz matches deSolve, and takes a deSolve function unchanged", {
  skip_if_not_installed("deSolve")
  times <- seq(0, 2, by = 0.05)
  ref <- deSolve::ode(y = c(1, 1, 1), times = times, func = lorenz_r_rhs,
                      parms = lorenz_pars, method = "radau", rtol = 1e-10, atol = 1e-10)
  out <- ode_solve(lorenz_r_rhs, y0 = c(1, 1, 1), times = times, parms = lorenz_pars,
                   jacfunc = lorenz_r_jac, rtol = 1e-10, atol = 1e-10, autonomous = TRUE)
  expect_equal(colnames(out), c("time", "y1", "y2", "y3"))
  expect_equal(out[, 2], as.numeric(ref[, 2]), tolerance = 1e-5)
  expect_equal(out[, 3], as.numeric(ref[, 3]), tolerance = 1e-5)
  expect_equal(out[, 4], as.numeric(ref[, 4]), tolerance = 1e-5)
})

testthat::test_that("a supplied Jacobian agrees with finite differences and is formed once per step", {
  times <- c(0, 0.5)
  fd <- ode_solve(lorenz_r_rhs, c(1, 1, 1), times, lorenz_pars, method = "rodas", autonomous = TRUE)
  an <- ode_solve(lorenz_r_rhs, c(1, 1, 1), times, lorenz_pars, method = "rodas",
                  jacfunc = lorenz_r_jac, autonomous = TRUE)
  expect_equal(an[2, -1], fd[2, -1], tolerance = 1e-5)
  ca <- attr(an, "counts")
  cf <- attr(fd, "counts")
  expect_equal(ca$n_jac, ca$n_steps)
  # The budget from tests/standalone: one evaluation to seed, six per attempt,
  # and for the finite-difference Jacobian n = 3 more per accepted step.
  expect_equal(ca$n_rhs, 1 + 6 * (ca$n_steps + ca$n_rejections))
  expect_equal(cf$n_rhs, 1 + 6 * (cf$n_steps + cf$n_rejections) + 3 * cf$n_steps)
  # Not declared autonomous: one more per accepted step for df/dt.
  na <- ode_solve(lorenz_r_rhs, c(1, 1, 1), times, lorenz_pars, method = "rodas",
                  jacfunc = lorenz_r_jac, autonomous = FALSE)
  cn <- attr(na, "counts")
  expect_equal(cn$n_rhs, 1 + 6 * (cn$n_steps + cn$n_rejections) + cn$n_steps)
})

testthat::test_that("RODAS through R takes far fewer steps than RKCK on stiff Van der Pol", {
  eps <- 1e-4
  vdp <- function(t, y, eps) list(c(y[[2]], ((1 - y[[1]]^2) * y[[2]] - y[[1]]) / eps))
  times <- seq(0, 2, by = 0.2)
  rodas <- ode_solve(vdp, c(2, 0), times, eps, method = "rodas", autonomous = TRUE)
  rkck <- ode_solve(vdp, c(2, 0), times, eps, method = "rkck", autonomous = TRUE)
  expect_true(all(is.finite(rodas)))
  expect_equal(rodas[, 2], rkck[, 2], tolerance = 1e-3)
  expect_lt(attr(rodas, "counts")$n_steps, attr(rkck, "counts")$n_steps / 5)
})

testthat::test_that("dense output under dopri agrees with landing on every time, for far fewer evaluations", {
  times <- seq(0, 2, by = 0.001)
  dense <- ode_solve(lorenz_r_rhs, c(1, 1, 1), times, lorenz_pars, method = "dopri",
                     rtol = 1e-8, atol = 1e-8, autonomous = TRUE, dense = TRUE)
  landed <- ode_solve(lorenz_r_rhs, c(1, 1, 1), times, lorenz_pars, method = "dopri",
                      rtol = 1e-8, atol = 1e-8, autonomous = TRUE, dense = FALSE)
  expect_equal(dense[, "time"], times)
  expect_equal(dense[, -1], landed[, -1], tolerance = 1e-5)
  expect_lt(attr(dense, "counts")$n_rhs, attr(landed, "counts")$n_rhs / 2)
  # A quartic is reproduced exactly, whatever the steps: the dense output has
  # the method's order.
  grid <- seq(0, 2, by = 0.01)
  quartic <- ode_solve(function(t, y, p) 4 * t^3, 0, grid, method = "dopri")
  expect_equal(quartic[, 2], grid^4, tolerance = 1e-10)
  expect_lt(attr(quartic, "counts")$n_steps, 30)
  # Under the other steppers the interpolant is cubic Hermite: exact for a
  # cubic, and one order short otherwise.
  cubic <- ode_solve(function(t, y, p) 3 * t^2, 0, grid, method = "rkck")
  expect_equal(cubic[, 2], grid^3, tolerance = 1e-10)
})

testthat::test_that("dopri through R matches deSolve's ode45, the same method", {
  skip_if_not_installed("deSolve")
  times <- seq(0, 2, by = 0.05)
  ref <- deSolve::ode(y = c(1, 1, 1), times = times, func = lorenz_r_rhs,
                      parms = lorenz_pars, method = "ode45", rtol = 1e-10, atol = 1e-10)
  out <- ode_solve(lorenz_r_rhs, c(1, 1, 1), times, lorenz_pars, method = "dopri",
                   rtol = 1e-10, atol = 1e-10, autonomous = TRUE)
  expect_equal(out[, 2], as.numeric(ref[, 2]), tolerance = 1e-6)
  expect_equal(out[, 3], as.numeric(ref[, 3]), tolerance = 1e-6)
  expect_equal(out[, 4], as.numeric(ref[, 4]), tolerance = 1e-6)
})

testthat::test_that("advance_collect returns the states at the times asked for, in one call", {
  s <- OdeSolver$new(function(t, y) -y, c(a = 1, b = 2), autonomous = TRUE)
  out <- s$advance_collect(c(0, 0.5, 1))
  expect_equal(dim(out), c(3, 3))
  expect_equal(out[, 1], c(0, 0.5, 1))
  expect_equal(out[, 2], exp(-c(0, 0.5, 1)), tolerance = 1e-6)
  expect_equal(out[, 3], 2 * exp(-c(0, 0.5, 1)), tolerance = 1e-6)
  expect_equal(s$time(), 1)
  expect_error(s$advance_collect(c(0, 2)), "same as current time")
  expect_error(s$advance_collect(numeric(0)), "at least length 1")
  expect_error(s$advance_collect(c(1, 3, 2)), "strictly increasing")
  landed <- s$advance_collect(c(1, 1.5, 2), dense = FALSE)
  expect_equal(landed[, 2], exp(-c(1, 1.5, 2)), tolerance = 1e-6)
  expect_equal(s$times()[length(s$times())], 2)
  expect_true(1.5 %in% s$times())
})

testthat::test_that("a non-autonomous right-hand side is integrated correctly", {
  times <- seq(0, 3, by = 0.5)
  out <- ode_solve(function(t, y, p) cos(t), y0 = 0, times = times, rtol = 1e-8, atol = 1e-8)
  expect_equal(out[, 2], sin(times), tolerance = 1e-6)
  for (method in c("rkck", "rodas")) {
    out2 <- ode_solve(function(t, y, p) cos(t), y0 = 0, times = times,
                      method = method, rtol = 1e-8, atol = 1e-8, dense = FALSE)
    expect_equal(out2[, 2], sin(times), tolerance = 1e-6)
  }
})

testthat::test_that("an error in a callback surfaces with its own message and poisons the solver until set_state", {
  expect_error(OdeSolver$new(function(t, y) stop("boom"), 1), "boom")
  s <- OdeSolver$new(function(t, y) if (t > 0.3) stop("boom later") else -y, 1)
  expect_error(s$advance_adaptive(c(0, 1)), "boom later")
  expect_error(s$step(), "set_state")
  s$set_state(1, 0)
  expect_silent(s$step(0.1))
  expect_error(OdeSolver$new(function(t, y) c(1, 2), 1), "returned 2 rates, expected 1")
  expect_error(OdeSolver$new(function(t, y) -y, 1, jac = function(t, y) matrix(0, 2, 2),
                             method = "rodas")$step(), "2 x 2 matrix, expected 1 x 1")
  expect_error(OdeSolver$new("not a function", 1), "rhs must be a function")
  expect_error(OdeSolver$new(function(t, y) -y, 1, method = "euler"), "Unknown method")
  expect_s3_class(OdeSolver$new(function(t, y) -y, 1, method = "ode45"), "OdeSolver")
})

testthat::test_that("domain_error() rejects the step rather than ending the solve", {
  logistic <- function(t, y) {
    if (y[1] < 0 || y[1] > 1) domain_error(paste("y =", y[1], "is outside [0, 1]"))
    50 * y[1] * (1 - y[1])
  }
  ctrl <- OdeControl$new()
  ctrl$set_controls(1e-2, 1e-2, 1, 0, 1e-8, 10, 1)
  for (method in c("rkck", "rodas")) {
    s <- OdeSolver$new(logistic, 0.5, control = ctrl, method = method,
                       jac = function(t, y) matrix(50 * (1 - 2 * y[1])))
    expect_silent(s$advance_adaptive(c(0, 1)))
    expect_true(s$state() >= 0 && s$state() <= 1)
    expect_equal(s$time(), 1)
    expect_gt(s$counts()$n_rejections, 0)
  }
  # And through a validity predicate instead.
  s <- OdeSolver$new(function(t, y) 50 * y * (1 - y), 0.5, control = ctrl, method = "rkck",
                     state_valid = function(t, y) y[1] >= 0 && y[1] <= 1)
  s$advance_adaptive(c(0, 1))
  expect_true(s$state() >= 0 && s$state() <= 1)
  expect_gt(s$counts()$n_rejections, 0)
  # A domain the exact flow leaves cannot be rescued, and says so.
  ramp <- function(t, y) { if (y[1] > 1) domain_error("over the ceiling"); 1 }
  s <- OdeSolver$new(ramp, 0.9, method = "rkck")
  expect_error(s$advance_adaptive(c(0, 1)), "over the ceiling")
})

testthat::test_that("the state can be re-seeded at a new length, and the step size read and set", {
  decay <- function(t, y) -y
  ctrl <- OdeControl$new()
  # Loose enough that a step of 0.125 passes the error estimate: at the default
  # 1e-8 a fourth-order method rejects it even on a linear system.
  ctrl$set_controls(1e-6, 1e-6, 1, 0, 1e-12, 10, 1e-2)
  s <- OdeSolver$new(decay, c(1, 2, 3), control = ctrl, autonomous = TRUE)
  s$step()
  expect_length(s$state(), 3)
  t1 <- s$time()
  s$set_state(c(0.5, 0.25), t1)
  expect_length(s$state(), 2)
  expect_equal(s$time(), t1)
  expect_equal(s$step_size(), 1e-2)
  expect_equal(s$times(), t1)
  s$set_step_size(0.125)
  s$step()
  expect_equal(s$time() - t1, 0.125)
  expect_equal(s$state(), c(0.5, 0.25) * exp(-0.125), tolerance = 1e-6)
  s$set_state(1, 0)
  while (s$time() < 0.5) s$step(0.5)
  expect_equal(s$time(), 0.5)
  expect_error(s$step(0.5), "already at time_max")
  expect_error(s$set_step_size(-1), "positive")
})

testthat::test_that("after every accepted step the last evaluation was at the solver's state", {
  seen <- NULL
  recording <- function(t, y) { seen <<- list(t = t, y = y); -y }
  for (method in c("rodas", "rkck")) {
    s <- OdeSolver$new(recording, c(1, 2), method = method, autonomous = TRUE)
    expect_identical(seen, list(t = 0, y = c(1, 2)))
    for (i in 1:5) {
      s$step()
      expect_identical(seen$t, s$time())
      expect_identical(seen$y, s$state())
    }
    s$set_state(c(3, 4, 5), 2)
    expect_identical(seen, list(t = 2, y = c(3, 4, 5)))
    expect_identical(s$rates(), -c(3, 4, 5))
  }
})
