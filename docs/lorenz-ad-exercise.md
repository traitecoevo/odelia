# ODE Adjoint Optimization Implementation Plan

**Goal:** Build a general-purpose AD interface for `odelia` that enables gradient-based optimization of ODE models, demonstrated on the Lorenz system.

**Design Invariants (Must Hold After Implementation):**
- The `odelia` solver core (`Solver`, `SolverInternal`, `Step`) must remain unchanged, AD support is additive.
- ODE systems must be templated on scalar type `T` to work with both `double` and XAD's `active_type`.
- The R interface must expose loss+gradient computation through R6 methods, matching the existing `LorenzSystem` pattern.
- Only `advance_fixed` may be used during AD recording (see Design Constraints).

**Architecture:** Template the system and solver infrastructure on scalar type. Wrap simulation in an AD tape (XAD), register inputs, run forward, compute loss, backpropagate. Expose to R via Rcpp.

**Tech Stack:** C++17, XAD (adjoint AD library), Rcpp, R6 classes.

---

## Background: Automatic Differentiation

### The Problem
You want to minimize a loss function $L(\theta)$ where $\theta$ are ODE parameters (initial conditions, rate constants). Optimizers need the gradient $\nabla_\theta L$, but computing this by hand is tedious because $L$ depends on the ODE solution, which depends on $\theta$ in a complex way.

### The Solution: AD
**Automatic Differentiation (AD)** computes exact gradients automatically. The idea:
1. Replace `double` with a special type (`active_type`) that remembers its history
2. Every arithmetic operation (`+`, `-`, `*`, `/`) is logged to a **tape**
3. After computing the final result, "rewind" the tape using the chain rule to get derivatives

This is **not** finite differences (which is approximate) or symbolic differentiation (which produces expressions). AD gives you exact numeric gradients at ~2-4x the cost of the forward computation.

### Forward vs Adjoint Mode
- **Forward mode**: Compute $\partial y / \partial x$ while computing $y$. Efficient when few inputs, many outputs.
- **Adjoint (reverse) mode**: Record forward pass, then propagate derivatives backward. Efficient when many inputs, few outputs (like a scalar loss).

We use **adjoint mode** (XAD's `adj<double>`) because we have many parameters but a single scalar loss.

### The Tape
Think of the tape as a log of operations:
```
Step 1: a = x + y        → record (a, x, y, "+")
Step 2: b = a * z        → record (b, a, z, "*")
Step 3: loss = b - target → record (loss, b, target, "-")
```

After the forward pass, the tape is "rewound":
- Start with $\partial L / \partial L = 1$
- At step 3: $\partial L / \partial b = 1$
- At step 2: $\partial L / \partial a = z$, $\partial L / \partial z = a$
- At step 1: $\partial L / \partial x = \partial L / \partial a = z$, same for $y$

XAD does all of this automatically, you just write normal arithmetic code.

---

## Design Constraints

### Why `advance_fixed` only?

The adaptive solver (`advance_adaptive`) decides step sizes at runtime: "is the error too large? If so, shrink the step and retry."

**Problem for AD:** The tape records exactly what happened during the forward pass. If the adaptive solver took 47 steps with specific sizes, the tape assumes those 47 steps are "the computation." If you ran it with a different number of steps, the tape would assume the gradient depended on a different computation. But the gradient should reflect how the loss changes with parameters, not how the step count changes. If changing $\theta$ slightly would cause the solver to take a different path (different number of steps), the gradient becomes discontinuous or wrong.

**Constraint:** Use `advance_fixed` with pre-specified time steps. Every run takes exactly the same steps, no branching, no path ambiguity.

### Why template on `T`?

XAD works by replacing `double` with `active_type`, a class that records every arithmetic operation. For `dy = sigma * (y1 - y0)` to be recorded:
- `sigma`, `y0`, `y1`, `dy` must all be `active_type`
- The `*` and `-` operators must be XAD's overloaded versions

By templating `LorenzSystem<T>`, the same source code works for:
- `T = double`: Normal simulation
- `T = xad::adj<double>::active_type`: AD-enabled simulation

---

## Task 1: Understand the Existing Interface

**Intent:** Familiarize yourself with how `odelia` and the Lorenz example work before modifying anything.

**Steps:**
1. Read `inst/examples/lorenz/readme.qmd` and run the examples.
2. Study `inst/examples/lorenz/src/lorenz_system.hpp`, understand the ODE interface methods (`ode_size`, `set_ode_state`, `ode_state`, `ode_rates`).
3. Study `inst/include/odelia/ode_solver.hpp`, understand how `Solver` uses the system.

**Deliverable:** Be able to explain how `Solver::advance_fixed` updates the system state.

---

## Task 2: Template the Lorenz System

**Intent:** Make `LorenzSystem` work with generic scalar types so XAD's `active_type` can flow through.

**Files:**
- Modify: `inst/examples/lorenz/src/lorenz_system.hpp`

**Constraints:**
- Template the class as `LorenzSystem<T = double>`
- All state variables (`y0`, `y1`, `y2`), rates, and parameters must be type `T`
- Iterator-based methods (`set_ode_state`, `ode_state`, `ode_rates`) must remain templated on iterator type
- Define `system_traits<LorenzSystem<T>>` specialization 

**Gap for you to solve:** The `record_step()` method returns `vector<double>` for history storage. How should this work when `T` is an AD type? (Hint: look at `odelia::util::value`.)

---

## Task 3: Create the Simulation Wrapper

**Intent:** Factor out the "run ODE and get result" logic so it can be reused for both forward simulation and AD.

**Files:**
- Create: `inst/include/odelia/ode_simulate.hpp`

**Constraints:**
- Implement `simulate_trajectory<System>(sys, ic, times)` returning `vector<state_type>`
- Must use `advance_fixed` only
- Must work for any `System` that satisfies the ODE interface

**Gap for you to solve:** How do you extract the trajectory from `Solver::get_history()`? The history stores `System` objects, you need to call `ode_state()` on each.

---

## Task 4: Create the AD Wrapper

**Intent:** Build the core AD logic that wraps simulation in an XAD tape.

**Files:**
- Create: `inst/include/odelia/ode_ad.hpp`

**Constraints:**
- Implement `compute_ode_gradient<SystemTemplate, LossFn>(ic, params, times, loss_fn)`
- Must: create tape, register inputs, call `simulate_trajectory`, compute loss, backpropagate, extract gradients

**XAD API you'll need:**
```cpp
using ad = xad::adj<double>;
using T = ad::active_type;
ad::tape_type tape;
tape.registerInput(x);  // Mark x as an input
tape.newRecording();
// ... computation ...
tape.registerOutput(y);
xad::derivative(y) = 1.0;  // Seed
tape.computeAdjoints();
double grad_x = xad::derivative(x);
double val_y = xad::value(y);
```

**Gap for you to solve:** The system constructor takes parameters. How do you pass them as `T` so they can be differentiated?

---

## Task 5: Define Loss Functors

**Intent:** Create reusable loss functions that work with the AD wrapper.

**Files:**
- Create: `inst/examples/lorenz/src/lorenz_fit.hpp`

**Constraints:**
- Define `TerminalLoss { target }` with `operator()(trajectory) -> T`
- Define `TrajectoryLoss { data }` with `operator()(trajectory) -> T`
- Both must be templated on the trajectory's scalar type

**Gap for you to solve:** How do you handle the case where `trajectory.size() != data.size()`? e.g. if you only have observations at some timepoints.

---

## Task 6: Create Rcpp Bindings

**Intent:** Expose the AD functionality to R.

**Files:**
- Modify: `inst/examples/lorenz/src/lorenz_interface.cpp`

**Constraints:**
- Export `lorenz_gradient(ic, params, target, times)`
- Return `List::create(Named("value") = ..., Named("gradient") = ...)`

**Gap for you to solve:** How do you convert an R `NumericMatrix` to `vector<vector<double>>`?

---

## Task 7: Extend R6 Interface

**Intent:** Add the `compute_loss` method to the R6 class so users can call `lz$compute_loss(...)`.

**Files:**
- Modify: `inst/examples/lorenz/R/lorenz_interface.R`

**Constraints:**
- Add `compute_loss(ic, data, times)` method
- Must call the C++ function and return the result

---

## Task 8: Create Demo and Verify

**Intent:** Demonstrate the full workflow and verify correctness.

**Files:**
- Create: `inst/examples/lorenz/demo_fit.R`

**Constraints:**
- Generate synthetic data from known IC/params
- Recover IC using `nlm()` or `optim()`
- Recover params using the same pattern

**Verification criteria:**
- IC estimation error < 1.0 from true
- Param estimation error < 1.0 from true
- Gradient norm at solution < 1e-6

---

## References
- `inst/examples/lorenz/readme.qmd`: Existing Lorenz demo
- XAD documentation: https://auto-differentiation.github.io/XAD/
