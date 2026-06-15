# Package index

## Systems and solvers

Define an ODE system and advance it through time.

- [`LorenzSystem`](https://traitecoevo.github.io/odelia/reference/LorenzSystem.md)
  : Lorenz System R6 Class
- [`Lorenz_Solver`](https://traitecoevo.github.io/odelia/reference/Lorenz_Solver.md)
  : Lorenz Solver R6 Class

## Solver control

Tolerances and step-size settings for the adaptive stepper.

- [`OdeControl`](https://traitecoevo.github.io/odelia/reference/OdeControl.md)
  : OdeControl R6 Class

## Drivers

External, time-varying forcing variables.

- [`Drivers`](https://traitecoevo.github.io/odelia/reference/Drivers.md)
  : Drivers R6 Class

## Package and utilities

- [`odelia-package`](https://traitecoevo.github.io/odelia/reference/odelia-package.md)
  [`odelia`](https://traitecoevo.github.io/odelia/reference/odelia-package.md)
  : ODE Solver with Automatic Differentiation, in C++ Header Files
- [`odelia_load_dll()`](https://traitecoevo.github.io/odelia/reference/odelia_load_dll.md)
  : Load odelia shared library with configurable symbol visibility

## Low-level wrappers

Thin wrappers over the underlying C++ functions. Most users should
prefer the R6 classes above.

- [`low_level_wrappers`](https://traitecoevo.github.io/odelia/reference/low_level_wrappers.md)
  [`Drivers_new`](https://traitecoevo.github.io/odelia/reference/low_level_wrappers.md)
  [`Drivers_set_constant`](https://traitecoevo.github.io/odelia/reference/low_level_wrappers.md)
  [`Drivers_set_variable`](https://traitecoevo.github.io/odelia/reference/low_level_wrappers.md)
  [`Drivers_set_extrapolate`](https://traitecoevo.github.io/odelia/reference/low_level_wrappers.md)
  [`Drivers_evaluate`](https://traitecoevo.github.io/odelia/reference/low_level_wrappers.md)
  [`Drivers_evaluate_range`](https://traitecoevo.github.io/odelia/reference/low_level_wrappers.md)
  [`Drivers_get_names`](https://traitecoevo.github.io/odelia/reference/low_level_wrappers.md)
  [`Drivers_clear`](https://traitecoevo.github.io/odelia/reference/low_level_wrappers.md)
  [`OdeControl_new`](https://traitecoevo.github.io/odelia/reference/low_level_wrappers.md)
  [`OdeControl_get_controls`](https://traitecoevo.github.io/odelia/reference/low_level_wrappers.md)
  [`OdeControl_set_controls`](https://traitecoevo.github.io/odelia/reference/low_level_wrappers.md)
  [`OdeControl_set_tol_abs`](https://traitecoevo.github.io/odelia/reference/low_level_wrappers.md)
  [`OdeControl_set_tol_rel`](https://traitecoevo.github.io/odelia/reference/low_level_wrappers.md)
  [`OdeControl_set_a_y`](https://traitecoevo.github.io/odelia/reference/low_level_wrappers.md)
  [`OdeControl_set_a_dydt`](https://traitecoevo.github.io/odelia/reference/low_level_wrappers.md)
  [`OdeControl_set_step_size_min`](https://traitecoevo.github.io/odelia/reference/low_level_wrappers.md)
  [`OdeControl_set_step_size_max`](https://traitecoevo.github.io/odelia/reference/low_level_wrappers.md)
  [`OdeControl_set_step_size_initial`](https://traitecoevo.github.io/odelia/reference/low_level_wrappers.md)
  : Low-level Drivers and OdeControl wrappers
