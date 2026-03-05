# ODE Adjoint Optimization Implementation Plan

**Goal:** Build a general-purpose AD interface for `odelia` that enables gradient-based optimization of ODE models, demonstrated on the Lorenz system.

**Design Invariants (Must Hold After Implementation):**
- The `odelia` solver core (`Solver`, `SolverInternal`, `Step`) must remain usable for both `double` and AD types—AD support is additive, achieved through templating.
- ODE systems must be templated on scalar type `T` to work with both `double` and XAD's `active_type`.
- The R interface must expose loss+gradient computation through R6 methods, matching the existing `LorenzSystem` pattern.
- Only `advance_fixed` may be used during AD recording (see Design Constraints).
- History collection must be disabled (`collect = false`) during AD to avoid polluting the tape.

**Architecture:** Template the system and solver infrastructure on scalar type. Wrap simulation in an AD tape (XAD), register inputs, run forward, compute loss at observation times, backpropagate. Expose to R via Rcpp.

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

XAD does all of this automatically when you write normal arithmetic code.

### Tape Pollution: Why Copies Matter
When the tape is recording, **copying an `active_type` value creates a new tape entry**. This is because XAD tracks data flow: when you copy `a` to `b`, the tape must record that `b` depends on `a`.

This has implications for odelia's history collection. The solver stores history as `std::vector<System>`, which copies the system at each step. With AD types, every copy pollutes the tape with spurious entries. **Solution:** Disable history collection (`collect = false`) during AD runs.

---

## Design Constraints

### Why `advance_fixed` only?

The adaptive solver (`advance_adaptive`) decides step sizes at runtime: "is the error too large? If so, shrink the step and retry."

**Problem for AD:** The tape records exactly what happened during the forward pass. If the adaptive solver took 47 steps with specific sizes, the tape assumes those 47 steps are "the computation." Changing parameters slightly might cause the solver to take different steps, making the gradient discontinuous or wrong.

**Constraint:** Use `advance_fixed` with pre-specified time steps. Every run takes exactly the same steps, no branching, no path ambiguity.

**Important:** The AD pass must use the **exact same time steps** as the forward pass. Different step sizes give different numerical results, so the gradient would be computing derivatives for a different numerical path.

### Why template on `T`?

XAD works by replacing `double` with `active_type`, a class that records every arithmetic operation. For `dy = sigma * (y1 - y0)` to be recorded:
- `sigma`, `y0`, `y1`, `dy` must all be `active_type`
- The `*` and `-` operators must be XAD's overloaded versions

By templating `LorenzSystem<T>`, the same source code works for:
- `T = double`: Normal simulation
- `T = xad::adj<double>::active_type`: AD-enabled simulation

---

## Task 1: Trace the Solver and Discover Templating Points

**Intent:** Understand how state flows through the solver and identify exactly where `double` is hardcoded.

**Guiding Trace:**

Follow this path through the code, noting where types are defined:

1. **Entry point:** `Solver::advance_fixed(times)` in `ode_solver.hpp`
   - What type is `history`? Where does it come from?
   - What calls does it make?

2. **Internal stepping:** `SolverInternal::step_to(system, time)` in `ode_solver_internal.hpp`
   - Find the member variables `y`, `yerr`, `dydt_in`, `dydt_out` (lines 50-53)
   - What is their type? Where is `state_type` defined?

3. **RK45 stepping:** `Step::step(...)` in `ode_step.hpp`
   - Find the intermediate storage `k1`, `k2`, ..., `ytmp` (line 35)
   - What is their type?

4. **System interface:** `ode::derivs(system, y, dydt, t, index)` in `ode_interface.hpp`
   - Find the global `typedef` for `state_type` (line 13)
   - Find the iterator typedefs (lines 14-15)

5. **System implementation:** `LorenzSystem` in `examples/lorenz_system.hpp`
   - Find all `double` declarations for state, rates, and parameters

**Expected Findings:**

| Location | Current | Problem |
|:--|:--|:--|
| `ode_interface.hpp:13` | `typedef std::vector<double> state_type` | Hardcoded double |
| `ode_interface.hpp:14-15` | Iterator typedefs | Derived from hardcoded state_type |
| `ode_step.hpp:35` | `state_type k1, k2, ...` | Uses global typedef |
| `ode_solver_internal.hpp:50-53` | `state_type y, yerr, ...` | Uses global typedef |
| `lorenz_system.hpp:78-81` | `double time, sigma, R, b, y0-y2` | All hardcoded |

**What can stay `double`?**
- RK45 coefficients in `ode_step.hpp` (lines 38-52): `const double ah[]`, etc. These are mathematical constants, not parameters we differentiate. `active_type * double` is valid.
- `OdeControl` methods: Only used by `advance_adaptive`, which we bypass.
- Time variables: Lorenz is autonomous (rates don't depend on time), and we don't differentiate w.r.t. time.

**Outcome:** A clear mental map of what needs templating and why.

---

## Task 2: Template the System and Define Traits

**Intent:** Make `LorenzSystem` work with generic scalar types and establish the traits pattern for propagating type information.

### The Problem: How Does the Solver Know What Types to Use?

Currently, `Step<System>` and `SolverInternal<System>` use a global `typedef std::vector<double> state_type`. When we template `LorenzSystem<T>`, we need these solver components to use `std::vector<T>` instead.

But `Step` is templated on `System`, not on `T` directly. It only sees `LorenzSystem<active_type>` as a type. How does it extract `active_type` from that?

**Answer: Traits.** A traits class provides a compile-time mapping from a `System` type to its associated types:

```cpp
// Primary template - default for non-templated systems
template<typename System>
struct system_traits {
    using value_type = double;
    using state_type = std::vector<value_type>;
};

Add this to `ode_interface.hpp`. After adding traits, the global typedef becomes obsolete for templated code:

```cpp
// AFTER: Can remove. Solver components should use system_traits<System>::state_type instead
typedef std::vector<double> state_type;
```

For `LorenzSystem<T>`, we need a **specialization** that extracts `T`:

```cpp
// Specialization - extracts T from LorenzSystem<T>
template<typename T>
struct system_traits<LorenzSystem<T>> {
    using value_type = T;
    using state_type = std::vector<T>;
};
```

**Why the specialization is required:** The primary template always returns `double`. Only the specialization can "look inside" `LorenzSystem<T>` to extract `T`. Without it, `Step<LorenzSystem<active_type>>` would incorrectly use `std::vector<double>`.

Now solver components can deduce the correct types:

```cpp
template<class System>
class Step {
    using value_type = typename system_traits<System>::value_type;
    using state_type = typename system_traits<System>::state_type;
    
    state_type k1, k2, ...;  // Correctly typed as vector<T>
};
```

The existing functions (like `ode::derivs`) are already templated on `System`, so they'll work automatically once the solver components use trait-based types.

### Files to Modify
- `inst/include/examples/lorenz_system.hpp` (template the class)
- `inst/include/odelia/ode_interface.hpp` (add traits, remove or deprecate global typedef)

### Requirements
- Template the class as `LorenzSystem<T = double>`
- All state variables (`y0`, `y1`, `y2`), rates (`dy0dt`, ...), and parameters must be type `T`
- Define `system_traits<LorenzSystem<T>>` specialization in `ode_interface.hpp`
- The iterator-based methods (`set_ode_state`, `ode_state`, `ode_rates`) should remain templated on iterator type

### Verification (Test Code)

Write a simple test to confirm the templating works. This is throwaway verification code, not production code:

```cpp
// Test: Explicit double template should work identically to before
LorenzSystem<double> sys(10, 28, 8.0/3.0);
```

---

## Task 3: Template the Solver Infrastructure

**Intent:** Propagate the scalar type `T` through `Step`, `SolverInternal`, and related components using the traits defined in Task 2.

### Files to Modify
- `inst/include/odelia/ode_interface.hpp`
- `inst/include/odelia/ode_step.hpp`
- `inst/include/odelia/ode_solver_internal.hpp`
- `inst/include/odelia/ode_solver.hpp`

### Key Changes

1. **ode_interface.hpp**: Replace the global `typedef` with trait-based types, or make the interface functions deduce types from their arguments.
2. **ode_step.hpp**: `Step<System>` should deduce `state_type` from `system_traits<System>`.
3. **ode_solver_internal.hpp**: Member variables `y`, `yerr`, `dydt_in`, `dydt_out` should use `system_traits<System>::state_type`.
4. **ode_solver.hpp**: `history` storage type needs consideration—but since we'll use `collect = false` for AD, this is less critical.

### Verification: Active Types Without Tape (Test Code)

Before introducing the tape, write a test to verify that AD types work correctly in forward mode. This is throwaway verification code:

```cpp
using ad = xad::adj<double>;
using T = ad::active_type;

// No tape active - active_type behaves like expensive double
LorenzSystem<double> sys_double(10, 28, 8.0/3.0);
LorenzSystem<T> sys_active(10, 28, 8.0/3.0);

// Set same initial conditions on both
// Run advance_fixed with same times on both
// Compare final states - should be IDENTICAL (bitwise)
// If they match, templating is correct. Discard test.
```

When no tape is active, `active_type` values don't record operations—they just carry values. This comparison verifies the templating is correct before adding tape complexity.

---

## Task 4: Implement the Gradient Computation

**Intent:** Build the AD integration that computes loss and gradients with respect to initial conditions.

### R Interface Design

The design uses an `active` flag to create an AD-capable solver:

```r
# Forward simulation (unchanged pattern)
lz <- LorenzSystem$new(sigma = 10, R = 28, b = 8/3)
lz$set_state(c(1, 1, 1))
ctrl <- OdeControl_new()
runner <- Lorenz_Solver$new(lz$ptr, ctrl$ptr)
runner$advance_fixed(times)
target <- runner$history()

# AD solver - explicit flag creates active-typed solver
runner_ad <- Lorenz_Solver$new(lz$ptr, ctrl$ptr, active = TRUE)
runner_ad$set_target(times, target, obs_indices)
result <- runner_ad$fit(ic = c(2, 2, 2))
# result$loss, result$gradient
```

`fit(ic)` forces `collect = FALSE` (to avoid polluting tape with copies, but you can set `collect = TRUE` afterwards. This means you can visualize the optimized trajectory using the same solver:

```r
# After optimization completes...
runner_ad$set_state(optimal_ic, 0.0)
runner_ad$set_collect(TRUE)
runner_ad$advance_fixed(times)
plot(runner_ad$history())
```

### C++ Architecture

**Rcpp layer** - dispatches based on `active` flag:

```cpp
// [[Rcpp::export]]
SEXP Solver_new(SEXP system_xp, SEXP control_xp, bool active = false) {
    auto sys = get_system(system_xp);
    auto ctrl = get_control(control_xp);
    
    if (active) {
        auto params = sys->pars();
        using ActiveSys = LorenzSystem<active_type>;
        auto* sys_active = new ActiveSys(params[0], params[1], params[2]);
        return XPtr<Solver<ActiveSys>>(new Solver<ActiveSys>(*sys_active, *ctrl));
    } else {
        return XPtr<Solver<System>>(new Solver<System>(*sys, *ctrl));
    }
}
```

**Solver class** - stores fit configuration:

```cpp
template <typename System>
class Solver {
    // Fit configuration (stored for repeated calls)
    std::vector<double> fit_times_;
    std::vector<size_t> obs_indices_;
    std::vector<std::vector<double>> targets_;
    
public:
    void set_target(times, targets, obs_indices);
    
    // TODO: define a suitable type for the return value
    std::pair<T, std::vector<T>> fit(const std::vector<double>& ic); // calls compute gradient
};
```
**Note on `times`**: Unlike `advance_fixed(times)` which takes times as an argument, `set_target` stores times for repeated `fit()` calls. In a future phase, this pattern enables replacing `times` with a step schedule recorded from adaptive solving.

**ode_fit.hpp** - AD logic (illustrative pseudocode):

```cpp
template<typename System>
FitResult compute_gradient(
    Solver<System>& solver,
    const std::vector<double>& ic
) {
    using ad = xad::adj<double>;
    using T = ad::active_type;
    
    // 1. Fresh tape (required each call)
    ad::tape_type tape;
    
    // 2. Register IC as inputs
    std::vector<T> ic_active(ic.begin(), ic.end());
    for (auto& x : ic_active) tape.registerInput(x);
    tape.newRecording();
    
    // 3. Reset solver with new IC
    solver.system().set_state(ic_active, 0.0);
    solver.reset();
    solver.set_collect(false);  // Don't pollute tape
    
    // 4. Step through times, accumulate loss at observation points
    T loss = 0.0;
    size_t obs_idx = 0;
    for (size_t i = 1; i < solver.fit_times().size(); ++i) {
        solver.advance_fixed(solver.fit_times()[i]);
        
        // if obs_idx  in this step 
        // accumulate the loss
            
    }
    
    // 5. Backpropagate
    tape.registerOutput(loss);
    xad::derivative(loss) = 1.0;
    tape.computeAdjoints();
    
    // 6. Extract results
    return {xad::value(loss), extract_gradients(ic_active)};
}
```

### Files to Create/Modify
- Create: `inst/include/odelia/ode_fit.hpp` (reusable AD logic)
- Modify: `inst/include/odelia/ode_solver.hpp` (add `set_target`, `fit` methods)
- Modify: `src/lorenz_interface.cpp` (Rcpp bindings with `active` flag)

---

## Task 5: R Interface and Verification

**Intent:** Expose the gradient computation to R and verify correctness using optimization.

### Rcpp Bindings

The key changes to `lorenz_interface.cpp`:

1. **Modify `Solver_new`** to accept `active` flag (see Task 4)
2. **Add `Solver_set_target`** to configure fit parameters
3. **Add `Solver_fit`** that takes only IC:

```cpp
// [[Rcpp::export]]
void Solver_set_target(SEXP solver_xp, 
                       Rcpp::NumericVector times,
                       Rcpp::NumericMatrix target,
                       Rcpp::IntegerVector obs_indices) {
    // Store in solver for repeated fit() calls
}

// [[Rcpp::export]]
Rcpp::List Solver_fit(SEXP solver_xp, Rcpp::NumericVector ic) {
    // Call compute_gradient with stored target config
    return List::create(
        Named("loss") = result.loss,
        Named("gradient") = wrap(result.gradient)
    );
}
```

### R6 Extension

Add to the `Lorenz_Solver` R6 class:

```r
initialize = function(system_ptr, control_ptr, active = FALSE) {
    private$ptr <- Solver_new(system_ptr, control_ptr, active)
    private$active <- active
},

set_target = function(times, target, obs_indices) {
    Solver_set_target(private$ptr, times, target, obs_indices)
    invisible(self)
},

fit = function(ic) {
    if (!private$active) stop("fit() requires active = TRUE solver")
    Solver_fit(private$ptr, ic)
}
```

### Verification

1. Compare double and active_types are numerically equivalent
2. Compare AD gradient vs finite difference; AD gradient should match finite differences to ~1e-6 relative error.
3. Run the fit and minimise the loss
4. Minimising the loss recovers the true ICs.


## Future Work: Adaptive Solver Integration

This exercise uses `advance_fixed` with predetermined time steps. For real applications, you often want adaptive stepping in the forward pass (for efficiency and accuracy) while maintaining valid gradients.

The solution is **schedule recording and replay**:
1. **Forward pass (double):** Run `advance_adaptive`, record the exact step schedule (times and step sizes)
2. **AD pass (active_type):** Replay the recorded schedule with `advance_fixed`, ensuring identical numerical path

This ensures the gradient reflects the actual computation that produced the forward solution.

---

## References
- `inst/examples/lorenz/readme.qmd`: Existing Lorenz demo
- XAD documentation: https://auto-differentiation.github.io/XAD/
