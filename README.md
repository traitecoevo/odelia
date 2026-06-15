# odelia: ODE solver with automatic differentiation, in C++ header files

<!-- badges: start -->
[![R-CMD-check](https://github.com/traitecoevo/odelia/workflows/R-CMD-check/badge.svg)](https://github.com/traitecoevo/odelia/master)
<!-- badges: end -->

`odelia` is an ODE solver implemented in C++ header files, using an adaptive-step
Runge-Kutta 4-5 method, with an interface to R via Rcpp. The solver runs entirely
in compiled code, so it is fast, and ODE systems can be templated on their scalar
type to support **automatic differentiation (AD)** — letting you compute exact
gradients of a solution with respect to its parameters for use in optimisation and
calibration.

The core solver was first developed by Rich FitzJohn as part of the
[plant package](https://github.com/traitecoevo/plant/). This package spins that
code out so it can be used more widely.

## Features

- Adaptive-step **RK4-5** integrator running entirely in C++ (~30-90x faster than
  the equivalent solved via `deSolve`; see the Lorenz example).
- **Automatic differentiation** of ODE solutions w.r.t. parameters and initial
  conditions, via the vendored [XAD](https://github.com/auto-differentiation/xad)
  library — enabling gradient-based parameter fitting.
- **External drivers**: time-varying forcing variables, smoothly interpolated with
  cubic splines and queried by the system at each step.
- Header-only C++ core that other Rcpp packages can link against.
- Friendly **R6** wrappers around the C++ objects.

## Installation

`odelia` compiles C++ from source, so you need a working C++ toolchain
(Rtools on Windows, Xcode command-line tools on macOS, a recent g++/clang on
Linux). Then:

```r
# install.packages("remotes")
remotes::install_github("traitecoevo/odelia")
```

## Quick start

`odelia` ships with one ready-to-run system, the classic
[Lorenz system](https://en.wikipedia.org/wiki/Lorenz_system), so you can try the
solver out of the box. **It is bundled purely as a demonstration** — to model
your own problem you write a small C++ system class (and a thin Rcpp/R6
interface) that the same solver then drives. See
[Building your own model](https://traitecoevo.github.io/odelia/articles/leaf-thermal.html)
for a complete worked example.

The example below solves the bundled Lorenz system:

```r
library(odelia)

# Define the system (parameters) and its initial state
lz <- LorenzSystem$new(sigma = 10, R = 28, b = 8 / 3)
lz$set_state(c(1, 1, 1))   # autonomous system, so no time needed

# Solver control settings (tolerances, step sizes) use sensible defaults
ctrl <- OdeControl$new()

# Build a solver (a "runner") and advance with adaptive stepping
runner <- Lorenz_Solver$new(lz$ptr, ctrl$ptr)
runner$advance_adaptive(seq(0, 100, by = 0.01))

# Collect output as a tibble
out <- runner$history()
```

For more, see the vignettes (also rendered on the
[package website](https://traitecoevo.github.io/odelia/)):

- `vignette("odelia")` — getting started with the core abstractions
  ([source](vignettes/odelia.Rmd)).
- `vignette("parameter-fitting")` — recovering parameters with automatic
  differentiation ([source](vignettes/parameter-fitting.Rmd)).

And the worked examples:

- [Building your own model with external drivers](https://traitecoevo.github.io/odelia/articles/leaf-thermal.html)
  — a leaf-thermal model showing how to define your own ODE system in C++ and
  drive it with time-varying forcing
  ([source](vignettes/articles/leaf-thermal.Rmd)).
- [Lorenz system](examples/lorenz/readme.qmd) — also includes a speed comparison
  against `deSolve`.

## Parameter fitting with automatic differentiation

Because the solver can be templated on an AD scalar type, you can recover exact
gradients of a loss (the mismatch between the solution and a set of target
observations) with respect to the system parameters, and hand them to a
gradient-based optimiser such as `optim()`:

```r
# An AD-enabled runner exposes a $fit() method returning loss and gradient
ad_runner <- Lorenz_Solver$new(lz$ptr, ctrl$ptr, active = TRUE)
ad_runner$set_target(times, target_vals, obs_index)

res <- ad_runner$fit(params = c(sigma = 12, R = 30, b = 3))
res$loss      # scalar mismatch with the target trajectory
res$gradient  # exact gradient w.r.t. each parameter
```

A complete optimisation workflow (recovering known Lorenz parameters) is walked
through in `vignette("parameter-fitting")`.

## Vocabulary

`odelia` is organised around a few core abstractions:

- **System** — your ODE model. A C++ class that holds parameters and state and
  knows how to compute its rates (right-hand side) `dy/dt`. Templated on its scalar
  type so it works with both `double` and AD types.
- **Stepper** — the numerical integration scheme (adaptive RK4-5) that takes one
  step of the system, estimating the error to choose the next step size.
- **Solver / runner** — drives the stepper forward over a requested set of times,
  applying step-size control and collecting the solution history.
- **Drivers** — external, time-varying forcing variables (e.g. air temperature)
  that the system queries during integration. Supplied as time series and
  interpolated with cubic splines.
- **Control** (`OdeControl`) — the solver's tuning knobs: absolute and relative
  tolerances, state/derivative scaling, and minimum/maximum/initial step sizes.

## License

This package is released under the
[GNU Affero General Public License v3 (AGPL-3)](LICENSE.md). The vendored XAD
library retains its own license; see [inst/include/XAD/LICENSE.md](inst/include/XAD/LICENSE.md).

## Contributing

Contributions are welcome. By submitting a pull request or code to this repository,
you agree to the terms of the [Contributor License Agreement](CLA.md).

### ODE System Structure

An ODE system in odelia consists of:

- **A system header** (`inst/include/odelia/examples/*.hpp`) — C++ class defining the ODE. Must implement `ode_size()`, `set_ode_state()`, `ode_state()`, and `ode_rates()`. Template on scalar type for AD support.

- **An Rcpp interface** (`src/*_interface.cpp`) — `[[Rcpp::export]]` functions exposing the system to R.

- **An R wrapper** (`R/*-interface.R`, optional) — R6 classes providing a friendlier API around the external pointers.

- **A demo script** (`examples/*/demo.R`) — Runnable demonstration.

### Testing variants

Use the Makefile targets below depending on the level of test coverage you want.

```bash
make test
```

Runs package tests after local compile. Fast default for day-to-day work.

```bash
make test-local
```

Runs `testthat::test_local()` from the source tree.

```bash
make test-installed
```

Installs with tests and runs `testthat::test_package()` against the installed package copy.
