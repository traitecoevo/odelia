# odelia: ODE solver using Runge-Kutta 4-5 adaptive step, C++ header files

<!-- badges: start -->
[![R-CMD-check](https://github.com/traitecoevo/odelia/workflows/R-CMD-check/badge.svg)](https://github.com/traitecoevo/odelia/master)
<!-- badges: end -->

The odelia package for R is an ODE solver in C++ header files, using RK4-5. There's potential to link to R via Rcpp. Code was first developed by Rich FitzJohn as part of the [plant package](https://github.com/traitecoevo/plant/). I'm spinning that code out into a package, as I want to use it elsewhere.

## Vocabulary
- system
- stepper
- solver
- runner

## License

This package is released under the [GNU Affero General Public License v3 (AGPL-3)](LICENSE).

## Contributing

Contributions are welcome. By submitting a pull request or code to this repository, you agree to the terms of the [Contributor License Agreement](CLA.md).

### ODE System Structure

An ODE system in odelia consists of:

- **A system header** (`inst/include/odelia/examples/*.hpp`) — C++ class defining the ODE. Must implement `ode_size()`, `set_ode_state()`, `ode_state()`, and `ode_rates()`. Template on scalar type for AD support.

- **An Rcpp interface** (`src/*_interface.cpp`) — `[[Rcpp::export]]` functions exposing the system to R.

- **An R wrapper** (`R/*-interface.R`, optional) — R6 classes providing a friendlier API around the external pointers.

- **A demo script** (`examples/*/demo.R`) — Runnable demonstration.
