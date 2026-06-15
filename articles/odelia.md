# Getting started with odelia

This vignette walks through the minimal setup for solving an ODE system
with `odelia`, using the classic [Lorenz
system](https://en.wikipedia.org/wiki/Lorenz_system) as the worked
example. The Lorenz system is built into the package, so everything
below runs against
[`library(odelia)`](https://github.com/traitecoevo/odelia) directly — no
C++ compilation needed.

``` r

library(odelia)
```

## The core abstractions

A simulation in `odelia` is assembled from a few pieces:

- a **system** — the ODE model itself (parameters, state, and a rule for
  the rates `dy/dt`);
- a **control** object (`OdeControl`) — the solver’s tolerances and step
  sizes;
- a **solver** (or *runner*) — which drives the adaptive RK4-5 stepper
  forward over the times you request and collects the solution history.

## Define the system

The Lorenz system has three state variables and three parameters
(`sigma`, `R`, `b`). We create it and set an initial state. The Lorenz
system is autonomous (its rates do not depend on time directly), so no
time is required.

``` r

lz <- LorenzSystem$new(sigma = 10, R = 28, b = 8 / 3)

# Inspect stored parameters
lz$pars()
#> [1] 10.000000 28.000000  2.666667

# Set the initial state (x, y, z)
lz$set_state(c(1, 1, 1))

# We can query the current state and rates directly
lz$state()
#> [1] 1 1 1
lz$rates()
#> [1]  0.000000 26.000000 -1.666667
```

## Solve the system

Create a control object (the defaults are sensible) and a solver, then
advance with adaptive time stepping over the times we want output at.

``` r

ctrl <- OdeControl$new()

runner <- Lorenz_Solver$new(lz$ptr, ctrl$ptr)
runner$advance_adaptive(seq(0, 100, by = 0.01))

out <- runner$history()
head(out)
#> # A tibble: 6 × 7
#>    time     x     y     z  dxdt  dydt    dzdt
#>   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>   <dbl>
#> 1  0     1     1    1      0     26   -1.67  
#> 2  0.01  1.01  1.26 0.985  2.47  26.1 -1.35  
#> 3  0.02  1.05  1.52 0.973  4.75  26.8 -0.997 
#> 4  0.03  1.11  1.80 0.965  6.91  28.1 -0.583 
#> 5  0.04  1.19  2.09 0.962  9.02  30.0 -0.0858
#> 6  0.05  1.29  2.40 0.964 11.1   32.4  0.520
```

[`history()`](https://rdrr.io/r/utils/savehistory.html) returns a tibble
with one row per requested output time. Because the stepper is adaptive,
the solver may take more internal steps than there are output rows — the
internal evaluation times are available via `runner$times()`.

## Plot the attractor

``` r

plot(out$x, out$z,
  type = "l",
  xlab = "x", ylab = "z",
  main = "Lorenz attractor (via odelia)"
)
```

![](odelia_files/figure-html/plot-1.png)

## Controlling the solver

`OdeControl` exposes the solver’s tuning knobs — absolute/relative
tolerances, state and derivative scaling, and minimum/maximum/initial
step sizes. You can inspect the current settings or tighten them:

``` r

ctrl$get_controls()
#> $tol_abs
#> [1] 1e-08
#> 
#> $tol_rel
#> [1] 1e-08
#> 
#> $a_y
#> [1] 1
#> 
#> $a_dydt
#> [1] 0
#> 
#> $step_size_min
#> [1] 1e-08
#> 
#> $step_size_max
#> [1] 10
#> 
#> $step_size_initial
#> [1] 1e-06

# Tighten tolerances for a more accurate (but slower) solve
ctrl$set_tol_abs(1e-8)
ctrl$set_tol_rel(1e-8)
```

## Where next

- For models with **external time-varying forcing**, and how to define
  your own ODE system in C++, see [Building your own model with external
  drivers](https://traitecoevo.github.io/odelia/articles/leaf-thermal.html).
- To recover parameters from data using exact gradients, see
  [`vignette("parameter-fitting", package = "odelia")`](https://traitecoevo.github.io/odelia/articles/parameter-fitting.md).
