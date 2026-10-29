# A Solver steps its own copy of the System it was built from. Parameters set
# on that System afterwards do not reach the solver; the solver's own
# set_params() does.

lorenz_at <- function(p) {
  lz <- LorenzSystem$new(sigma = p[1], R = p[2], b = p[3])
  lz$set_initial_state(c(1, 1, 1), t0 = 0)
  lz
}

run_to_one <- function(solver) {
  solver$advance_adaptive(c(0, 1))
  solver$state()
}

p_true <- c(10, 28, 8 / 3)
p_other <- c(12, 30, 3)

testthat::test_that("parameters set on the System after the solver is built do not reach it", {
  ctrl <- OdeControl$new()
  lz <- lorenz_at(p_true)
  s <- Lorenz_Solver$new(lz$ptr, ctrl$ptr)
  testthat::expect_identical(s$pars(), p_true)

  lz$set_params(p_other)
  testthat::expect_identical(lz$pars(), p_other)
  testthat::expect_identical(s$pars(), p_true)

  reference <- Lorenz_Solver$new(lorenz_at(p_true)$ptr, ctrl$ptr)
  testthat::expect_identical(run_to_one(s), run_to_one(reference))
})

testthat::test_that("parameters set on the solver give the run of a solver built at them", {
  ctrl <- OdeControl$new()
  s <- Lorenz_Solver$new(lorenz_at(p_true)$ptr, ctrl$ptr)
  s$set_params(p_other)
  testthat::expect_identical(s$pars(), p_other)

  reference <- Lorenz_Solver$new(lorenz_at(p_other)$ptr, ctrl$ptr)
  testthat::expect_identical(run_to_one(s), run_to_one(reference))
  testthat::expect_false(isTRUE(all.equal(
    run_to_one(Lorenz_Solver$new(lorenz_at(p_true)$ptr, ctrl$ptr)), s$state())))
})

testthat::test_that("set_params after a run, with reset, starts the new run from the initial state", {
  ctrl <- OdeControl$new()
  s <- Lorenz_Solver$new(lorenz_at(p_true)$ptr, ctrl$ptr)
  run_to_one(s)
  s$set_params(p_other)
  s$reset()
  testthat::expect_identical(s$time(), 0)
  testthat::expect_identical(s$state(), c(1, 1, 1))

  reference <- Lorenz_Solver$new(lorenz_at(p_other)$ptr, ctrl$ptr)
  testthat::expect_identical(run_to_one(s), run_to_one(reference))
  testthat::expect_identical(s$times(), reference$times())
})

testthat::test_that("set_params refuses the wrong number of parameters", {
  s <- Lorenz_Solver$new(lorenz_at(p_true)$ptr, OdeControl$new()$ptr)
  testthat::expect_error(s$set_params(c(1, 2)), "one entry per parameter")
  testthat::expect_identical(s$pars(), p_true)
})
