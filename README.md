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

## Examples

Example systems are in `examples/`:

| Example | Description |
|---------|-------------|
| `lorenz` | Classic Lorenz attractor ODE |
| `parabola` | AD gradient optimization |

### Development

```bash
make
```

```r
devtools::load_all()
```

`make` regenerates Rcpp exports and compiles C++. `load_all()` compiles code and loads it into R, but doesn't regenerate functions with `[[Rcpp::export]]` annotations - always run `make` if you're editing Rcpp interfaces.

### ODE System Structure

An ODE system in odelia consists of:

- **A system header** (`inst/include/odelia/examples/*.hpp`) — C++ class defining the ODE. Must implement `ode_size()`, `set_ode_state()`, `ode_state()`, and `ode_rates()`. Template on scalar type for AD support.

- **An Rcpp interface** (`src/*_interface.cpp`) — `[[Rcpp::export]]` functions exposing the system to R.

- **An R wrapper** (`R/*-interface.R`, optional) — R6 classes providing a friendlier API around the external pointers.

- **A demo script** (`examples/*/demo.R`) — Runnable demonstration.
