# odelia: ODE solver with automatic differentiation, in C++ header files

<!-- badges: start -->
[![R-CMD-check](https://github.com/traitecoevo/odelia/workflows/R-CMD-check/badge.svg)](https://github.com/traitecoevo/odelia/actions)
<!-- badges: end -->

`odelia` is an ODE solver implemented in C++ header files, with an interface to R
via Rcpp. Three adaptive steppers (Cash-Karp RK4(5), Dormand-Prince 5(4), and
RODAS4(3) for stiff systems) run entirely in compiled code, so it is fast. ODE
systems are templated on their scalar type, so a run can be recorded and
differentiated in **reverse mode**: one sweep over the recorded run gives the
exact derivative of a solution with respect to every parameter and initial
condition at once, for optimisation and calibration.

The core solver was first developed by Rich FitzJohn as part of the
[plant package](https://github.com/traitecoevo/plant/). This package spins that
code out so it can be used more widely.

## Features

- Three adaptive steppers running entirely in C++: explicit Cash-Karp RK4(5),
  Dormand-Prince 5(4) with dense output, and the implicit RODAS4(3) for stiff
  systems. A compiled system solves about a hundred times faster than the same
  right-hand side written in R (Lorenz: 1 ms against 25-60 ms under `deSolve`).
- **Reverse-mode automatic differentiation** of a recorded run, via the vendored
  [XAD](https://github.com/auto-differentiation/xad) library: one sweep, every
  parameter and initial condition. Forward mode supplies the implicit stepper's
  Jacobian. See the article
  [Reverse mode](https://traitecoevo.github.io/odelia/articles/reverse-mode.html).
- **External drivers**: time-varying forcing variables, interpolated with a cubic
  Hermite spline and queried by the system at each step.
- A C++ core (header-only but for the XAD tape runtime) that other packages link
  against, with no R needed to compile or test it.
- Friendly **R6** wrappers around the C++ objects.
- A right-hand side **written in R** can be solved by the same steppers, through
  `ode_solve()` (shaped like `deSolve::ode()`) or the step-at-a-time `OdeSolver`,
  with genuinely adaptive Runge-Kutta stepping, dense output of the method's
  order (Dormand-Prince 5(4)), the accepted steps exposed, and a count of what
  the solve cost. On a fine output grid that is about twice as fast as
  `deSolve::ode45` with the same R function.

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

The example below solves the bundled Lorenz system. (A right-hand side written
in R can be solved without any C++ at all; see `ode_solve()` and
`vignette("odelia")`.)

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
- [Reverse mode: one solve, every parameter](https://traitecoevo.github.io/odelia/articles/reverse-mode.html)
  — the gradient machinery from C++: what a run records, how the sweep walks it,
  what a System must provide ([source](vignettes/articles/reverse-mode.Rmd)).

## From C++

A System is a class templated on its scalar type that holds parameters and
state and knows its rates. The solver drives it:

```cpp
// [[Rcpp::plugins(cpp20)]]   (or CXX_STD = CXX20 in a package's Makevars)
#include <odelia/ode_solver.hpp>
#include <examples/lorenz_system.hpp>   // the shipped example System

LorenzSystem<double> system(10.0, 28.0, 8.0 / 3.0);
odelia::ode::Solver<LorenzSystem<double>> s(system, odelia::ode::OdeControl());
s.advance_adaptive({0.0, 10.0});
std::vector<double> y = s.state();
```

A gradient of that run is three more lines: keep the record, seed the output
wanted, and sweep.

```cpp
s.set_keep_states(true);                       // before the run
s.advance_adaptive({0.0, 10.0});
auto lambda = odelia::ode::adjoint_rows::one_row({1.0, 0.0, 0.0});   // d x(10)
odelia::ode::adjoint_rows dp(1, 3);            // one row per seed, zeroed
s.solve_adjoint(lambda, dp);                   // dp[0][j] = d x(10) / d parameter j
```

A package that uses the headers adds `odelia` to `LinkingTo` and `Imports`,
compiles as C++20 with `-DXAD_NO_THREADLOCAL -DXAD_USE_STRONG_INLINE`, and on
Windows links against odelia's DLL; [ARCHITECTURE.md](ARCHITECTURE.md) has the
details and a map of the headers. `inst/include/examples/lorenz_system.hpp` is
the template for a System of your own, with its members marked by which
contract they serve; the article
[Building your own model](https://traitecoevo.github.io/odelia/articles/leaf-thermal.html)
walks through one with external drivers.

## Parameter fitting with automatic differentiation

The shipped Lorenz and leaf-thermal runners expose a `$fit()` method that
returns the least-squares loss against a target trajectory and its exact
gradient with respect to the parameters (and initial conditions), by one
reverse sweep over a replay of the reference run. Hand it to a gradient-based
optimiser such as `optim()`:

```r
fit_runner <- Lorenz_Solver$new(lz$ptr, ctrl$ptr)
fit_runner$set_target(times, target_vals, obs_index)

res <- fit_runner$fit(params = c(sigma = 12, R = 30, b = 3))
res$loss      # scalar mismatch with the target trajectory
res$gradient  # exact gradient w.r.t. each parameter
```

A complete optimisation workflow (recovering known Lorenz parameters) is walked
through in `vignette("parameter-fitting")`. A runner for your own System gets
`$fit()` by instantiating `Solver_fit_impl` from `solver_interface.hpp`, as the
leaf-thermal example does.

## Vocabulary

`odelia` is organised around a few core abstractions:

- **System** — your ODE model. A C++ class that holds parameters and state and
  knows how to compute its rates (right-hand side) `dy/dt`. Templated on its scalar
  type so it works with both `double` and AD types.
- **Stepper** — the numerical integration scheme (Cash-Karp RK4(5), Dormand-Prince
  5(4) or RODAS4(3)) that takes one step of the system, estimating the error to
  choose the next step size.
- **Solver** — `Solver<System>` in C++ drives the stepper forward over a requested
  set of times, applying step-size control, and sweeps a recorded run backward
  for its gradient. In R, `Lorenz_Solver` wraps that for the shipped system and
  `OdeSolver` steps a right-hand side written in R.
- **Drivers** — external, time-varying forcing variables (e.g. air temperature)
  that the system queries during integration. Supplied as time series and
  interpolated with a cubic Hermite spline.
- **Control** (`OdeControl`) — the solver's tuning knobs: absolute and relative
  tolerances, state/derivative scaling, minimum/maximum/initial step sizes, and
  the step-size rule.
- **Recording** — what a run keeps when asked (`set_keep_states`): one row per
  accepted step with its time, step size, state, and whatever the step solved
  for inside a stage. A sweep reads it; a replay reproduces the run from it.
- **Sweep** — the backward pass over a recording (`solve_adjoint`). It is seeded
  with the derivative wanted of the final state, one **seed** per output, and
  returns one **row** of derivatives per seed against the parameters and the
  initial state; `adjoint_rows` is the batch.
- **Insertion** — a scheduled growth of the state vector mid-run. The sweep
  transposes the System's own widening map across it.
- **Supplied derivative** — a value put on the tape with rows obtained by other
  means (`record_with_derivatives`, `implicit_value`), for a quantity a submodel
  solved for by iteration.

## License

This package is released under the
[GNU Affero General Public License v3 (AGPL-3)](LICENSE.md). The vendored XAD
library retains its own license; see [inst/include/XAD/LICENSE.md](inst/include/XAD/LICENSE.md).

## Contributing

Contributions are welcome. By submitting a pull request or code to this repository,
you agree to the terms of the [Contributor License Agreement](CLA.md).

> **Linking from another package?** If you `LinkingTo: odelia` and instantiate
> `Solver`, read [ARCHITECTURE.md](ARCHITECTURE.md) — it documents how the
> compiled XAD `Tape` symbols are resolved per platform, what your package must
> do, and the invariants that must not be broken.

### ODE System Structure

An ODE system in odelia consists of:

- **A system header** (`inst/include/examples/lorenz_system.hpp`, `inst/examples/leaf_thermal/src/leaf_thermal_system.hpp`) — a C++ class templated on its scalar type. To be solved it implements `ode_size()`, `set_ode_state()`, `ode_state()` and `ode_rates()`; to be swept it adds `rebind_from<U>()`, `ad_parameters()`, `for_each_active(f)` and `set_recorded_state(y, time)` (the `Sweepable` concept in `ode_interface.hpp`).

- **An Rcpp interface** (`src/*_interface.cpp`, `inst/examples/leaf_thermal/src/*_interface.cpp`) — `[[Rcpp::export]]` functions exposing the system to R, built on the `Solver_*_impl` templates in `solver_interface.hpp`.

- **An R wrapper** (`R/*-interface.R`, optional) — R6 classes providing a friendlier API around the external pointers.

- **A worked example** — the article for the leaf-thermal model, built by pkgdown.

### Testing

```bash
make test
```

Installs the package with its tests and runs the full suite against the installed copy. The AD and DLL-lifecycle tests need an installed package.

```bash
make test-local
```

The fast development loop: `testthat::test_local()` via `load_all()`, skipping those tests.

```bash
make test-cpp
```

Builds and runs the solver core as plain C++ with no R on the include path (`tests/standalone/`), and compiles every core header on its own. Anything compiled against the headers outside the package's own `src/` needs C++20 and `-DXAD_NO_THREADLOCAL -DXAD_USE_STRONG_INLINE`, matching `src/Makevars`.

## Plant family

`odelia` is part of the **plant family** of packages in the
[`traitecoevo`](https://github.com/traitecoevo) org, built around the
[`plant`](https://github.com/traitecoevo/plant) forest model. Docs hub:
<https://traitecoevo.github.io/overstorey/>.

**Contributing:** please skim the family
[issue guide](https://github.com/traitecoevo/plant-meta/blob/main/governance/issue-guide.md)
before filing — issues across the family are triaged on
[board #5](https://github.com/orgs/traitecoevo/projects/5), and cross-package context lives in
[`plant-meta`](https://github.com/traitecoevo/plant-meta).
