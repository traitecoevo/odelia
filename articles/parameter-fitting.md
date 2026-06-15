# Parameter fitting with automatic differentiation

A distinctive feature of `odelia` is that ODE systems are templated on
their scalar type, so the solver can run with an *automatic
differentiation* (AD) number type as well as ordinary `double`. This
lets the solver return the exact gradient of a loss function — the
mismatch between the solution and a set of target observations — with
respect to the system’s parameters. Those gradients can then be handed
to a gradient-based optimiser to calibrate the model.

This vignette demonstrates the workflow by recovering known Lorenz
parameters from a reference trajectory.

``` r

library(odelia)
```

## Generate a reference trajectory

First we solve the system with known “true” parameters to produce the
data we will later try to recover.

``` r

true_pars <- c(sigma = 10.0, R = 28.0, b = 8.0 / 3.0)

lz <- LorenzSystem$new(sigma = true_pars[1], R = true_pars[2], b = true_pars[3])
lz$set_initial_state(c(1, 1, 1), t0 = 0)

ctrl <- OdeControl$new()

runner <- Lorenz_Solver$new(lz$ptr, ctrl$ptr)
runner$advance_adaptive(c(0, 1))

times <- runner$times()
hist <- runner$history()

# The internal solver times that coincide with stored output rows
obs_index <- which(times %in% hist$time)
target_vals <- as.matrix(hist[, c("x", "y", "z")])
```

## Set up an AD-enabled solver

We move the system parameters away from the truth (this is our starting
guess), then build a solver with `active = TRUE` to enable AD, and
register the calibration target.

``` r

initial_guess <- c(sigma = 12.0, R = 30.0, b = 3.0)
lz$set_params(initial_guess)

ad_runner <- Lorenz_Solver$new(lz$ptr, ctrl$ptr, active = TRUE)
ad_runner$set_target(times, target_vals, obs_index)
```

The `$fit()` method solves the system for a given parameter vector and
returns both the loss and its exact gradient:

``` r

res0 <- ad_runner$fit(params = initial_guess)
res0$loss
#> [1] 3.480593
res0$gradient
#> [1] 1.0865743 0.4861987 1.9423226
```

## Optimise

Because each `$fit()` call returns both loss and gradient, we cache the
most recent result so the optimiser does not have to re-run the solver
to obtain the gradient separately.

``` r

memoise_fit <- function(runner) {
  cache_p <- NULL
  cache_res <- NULL

  fn <- function(p) {
    cache_res <<- runner$fit(params = p)
    cache_p <<- p
    cache_res$loss
  }
  gr <- function(p) {
    if (!identical(p, cache_p)) fn(p) # populate cache if stale
    cache_res$gradient
  }
  list(obj = fn, grad = gr)
}

helper <- memoise_fit(ad_runner)

res <- optim(
  par = initial_guess,
  fn = helper$obj,
  gr = helper$grad,
  method = "L-BFGS-B",
  control = list(maxit = 100)
)

res$convergence # 0 indicates successful convergence
#> [1] 0
res$value       # final loss, should be ~0
#> [1] 2.697385e-13
```

## Check the recovered parameters

``` r

rbind(true = true_pars, recovered = res$par)
#>           sigma  R        b
#> true         10 28 2.666667
#> recovered    10 28 2.666667
```

The optimiser recovers the true parameters to high precision, confirming
that the gradients produced by AD are correct and usable for
calibration.
