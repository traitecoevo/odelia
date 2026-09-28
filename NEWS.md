## odelia 0.6.0

**A recorded run can now be differentiated in reverse mode, so one solve yields the
derivative with respect to every parameter at once.** Forward mode answers for one
parameter per solve; reverse mode answers for all of them, but it must first record
every operation the solve performed, and a run of several thousand adaptive steps
performs far more of them than memory will hold. The solver stores its state at each
accepted step and replays one step's arithmetic at a time while the derivative is
taken, so what is held is bounded by a single step rather than by the run. Memory
then grows with the number of steps, which is cheap, rather than with the
arithmetic, which is not.

A value a submodel solved for rather than computed is lifted onto the record
carrying a derivative obtained by other means — `implicit_value` and
`record_with_derivatives` — which keeps an iterative solver's own iterations out of
the answer. What *selects* rather than *moves* (a step size, a knot position, an
arm) stays a plain `double` and is replayed, because differentiating a selector
manufactures a discontinuity the model does not have.

**A forward pass can now REPLAY a recording instead of only taking one.** The
store/load channel already picked its direction by the constness of what a walk
handed over, and `step_adjoint` already handed a const row — but every forward
entry point funnelled through one place that zeroed its scratch and passed it
mutable, so nothing on the way forward could load. `Step::step`'s row parameter is
templated, so a const row loads and a mutable one stores with no change to the
body; `advance_recorded()` gains an overload taking a **recording** rather than a
program, pairing each step with its own row so the two cannot be crossed. No new
concept and no new name. `method='rodas'` refuses rather than silently re-deriving:
the Rosenbrock stepper keeps no per-stage row.

What this is for: a pass that must take the run's answers rather than its own. The
first consumer is plant's invasion run, where an invader stands in a resident's
recorded field — exogenous to it, so its derivative is zero rather than severed.

**⚠️ The recorded row is six long, not five, and that changes `step_record`.** Five
of the six are a step's stages; the sixth is the evaluation at the state the step
ends at, which first-same-as-last hands the next step as its own k1. A sweep
re-derives that one at the state it was handed and still reads only the first five.
A forward replay cannot — re-deriving is the thing it replays to avoid — and a step
whose k1 was re-derived is wrong at first order in `h`.

The forward-mode scalar lives in `tangent.hpp`, apart from the reverse-mode
machinery in `adjoint.hpp`, so a consumer wanting a directional derivative and no
record is never handed vocabulary for one. `tangent.hpp` rejects a forward scalar
nested above a reverse one at compile time: three kernels costing 31 statements flat
cost 566 nested.

`ode_fit.hpp` is removed; `compute_gradient` has no drop-in, and `vector_jacobian_product` plus `sweep.hpp` replace it in C++. The R fitting interface is kept on the new machinery: `$set_target(times, target, obs_indices)` and `$fit(ic, params)` return the same least-squares loss and gradient as before, now by one reverse sweep over a replay of `times`, for both the Lorenz and the leaf thermal examples. What goes is the `active` argument on the `Solver_*` bindings and the separate AD solver it made, since every solver can now fit. A System needs `set_recorded_state()` and `for_each_active()` to be swept. The System contract becomes C++20 concepts (`HasOdeTime`, `SolvesForValues`, `ChecksState`, `Rebindable`), and `rebind()` becomes `rebind_from()`.

A **minor** bump, and a breaking one for a System written against 0.5.0's traits or for a caller that passed `active = TRUE`.
## odelia 0.5.1

**The spline reads as fast as 0.4.0's again, with the same numbers.** On a graded knot grid, 0.5.0 found a query's span by binary search, where 0.4.0 guessed from the mean spacing and stepped from there (#21). plant's adaptive light field is graded and read in height order, so the search made its FF16 runs 17% slower than on 0.4.0. `hermite_spline` now uses the guess-and-step lookup again, which returns exactly the span the search did: 2.55 million reads on random graded grids, knots and their neighbouring doubles included, are bit-identical. The front end's unchecked `operator()` also skips the initialisation check, as 0.4.0's did.

Against plant on 0.4.0, interleaved: FF16 full lifetime 0.106 to 0.102 s, K93 0.036 to 0.035 s, and TF24 6 years 2.35 to 2.31 s.

## odelia 0.5.0

**The spline's backend becomes a cubic Hermite that can take a slope at each knot; the front end does not change, and neither do its numbers.** A model that integrates over crowns knows the derivative at each knot as well as the value, and the global natural spline had no way to accept it. `spline.hpp` now holds `hermite_spline<S>`, a local cubic Hermite built from a value and a slope per knot. `interpolator.hpp` keeps every call the family makes — `init(x, y)`, `eval` with its out-of-domain refusal, `deriv`, `min`/`max`, `set_extrapolate`, `r_eval` — on `hermite_interpolator<S>` (was `basic_interpolator<S>`), which derives from the backend so a caller that has slopes reaches `init(x, y, m)`, `value_and_slope`, `set_nodes`/`set_data` and `refine`.

A caller with values alone gets the natural cubic spline's own knot slopes (`natural_slopes`), and a Hermite read through those slopes is the natural spline on every span, so `init(x, y)` reads the curve it read before: `test-drivers.R` and `test-spline.R` are unchanged, and phylloptim 0.8.1 compiled against this release makes the same number of solver evaluations at the same speed (82.3 per solve, 3.0–3.1 µs, interleaved against 0.4.0), with outputs equal to rounding (worst 2e-10 absolute on its golden grid, from evaluating the same cubic in a different order). With supplied slopes the Hermite reproduces a cubic exactly and a read reaches four knots rather than all of them, which is what keeps an adjoint through it O(1). `monotone_slopes` (Fritsch–Carlson) is there for a caller that needs an interpolant that never leaves its data's range.

`value_with_slope<T>` pairs a value with its slope so the two cannot be handed over separately and paired wrongly; `util::to_passive` strips every derivative layer from a scalar.

A **minor** bump: breaking only for a consumer that named `basic_interpolator<S>` or `spline::basic_spline<S>` directly. Split out of the reverse-mode change (#59).

## odelia 0.4.0

**Step rejection now works when the integration is pinned to fixed times (plant#642).** `advance_fixed()` steps exactly to a caller-supplied set of times — how a replay reproduces a trajectory recorded earlier. It called the stepper bare, so both of the ways #55 gave a system to refuse a state were unreachable from it: a `util::DomainError` from a stage killed the solve, and a declared `ode_state_valid()` was never consulted at all.

The endpoints of a pinned step are the caller's and cannot be moved, so #55's answer — take a smaller step — becomes *take several smaller ones to the same endpoint*. On a refusal the sub-step shrinks through `control.reject_step()`, the same rule and the same `step_size_min` floor as the adaptive path, and the walk continues to the original time. That time is still hit exactly, so the times a caller records are unchanged. A system that raises no objection is stepped exactly as before: one RKCK step per interval, six stage evaluations, asserted in `tests/standalone/`.

An unreachable domain still fails, and now says where it gave up and why, rather than surfacing as whatever the offending stage happened to throw.

This was not a rare corner. `plant`'s mutant replay pins the stepper to a resident's recorded times, and its TF24 model reports an empty carbon pool this way as a matter of routine — ~480 rejections in a resident run that goes on to complete normally — so a replay was near-certain to meet one and die. Invasion-fitness analysis was impossible for that model, not merely slow.

One caveat for systems that cache per-stage data through `cache(system, rk_step)`: the stage indices restart at 0 on each sub-step, so a subdivided interval leaves the system holding the last sub-step's stages rather than stages spanning the whole interval. A consumer recording such a cache for later replay gets a coarser record of a subdivided step than of a plain one. **⚠️ 0.5.0 deletes `cache(system, rk_step)` and this caveat's machinery with it** — the per-stage record is now `step_record::solved`, addressed by the walk rather than by a cursor in the System. The hazard it names did not go away with the spelling: a replay driven by `step_by` takes the recorded size and has no subdivision path at all, so a system that refuses a state mid-replay fails rather than shrinking.

A **minor** bump, so downstreams can pin against the capability (`odelia (>= 0.4.0)`). Systems that never refuse a state are unaffected.

## odelia 0.3.1

**Removes `util::to_string_g()`, added in 0.3.0 an hour earlier as a duplicate of
`util::format_double()`.** Both formatted a double for an error message with six
significant figures — `"%g"` and `"%.6g"` are the same format string, since `%g`'s
default precision is 6 — and they were byte-identical on every value tried. The one
call site, in the invariant-rejection failure message, now uses `format_double`; its
output is unchanged.

The duplicate arose because #55 was developed on a branch stacked below the 0.2.2 work
that introduced `format_double`. Worth recording as a hazard rather than a review miss:
two helpers with *different names* in different parts of the same file merge with no
textual conflict, so neither the rebase nor the diff had anything to show.

A **patch** bump, although the header core lost a public symbol — the situation that
earned 0.2.0 a minor one. The difference is what a version number can usefully say.
0.2.0's removals had been reachable across released versions, so the bump warned of a
break a consumer could actually hit. `to_string_g` never left this repository:
**0.3.0 was never tagged** — `v0.2.1` remains the only tag and the only release — and
both consumers are still on `odelia (>= 0.2.x)` with `@v0.2.1` remotes. Nothing can
have compiled against it, so a minor bump would announce an incompatibility that has
no possible victim, and spend the signal for nothing.

Downstream floors are deliberately **not** raised. Neither `plant` nor `phylloptim`
uses `ode_state_valid()` yet, and per the precedent set for `leaf` in 0.2.1 — checked
rather than aligned for symmetry — raising a floor forces an upgrade for a change the
consumer does not use. They should move to `odelia (>= 0.3.0)` when plant#609 or
plant#599 actually adopts the domain check; 0.3.0 is the version that introduced it,
and this release does not change it.

The 0.3.0 entry below has been corrected accordingly; it advertised a function that
no longer exists.

## odelia 0.3.0

**Invariant-aware step rejection (#55).** A system may now declare the domain its
state lives in, and the adaptive stepper will refuse to commit a step that leaves it.
Two ways to say so, both opt-in:

- an optional `bool ode_state_valid(const state_type&) const` on the system, checked
  after each completed step;
- `util::stop_domain(msg)`, which throws the new `util::DomainError`, from anywhere in
  a stage.

Either one turns the step into a *rejection* — shrink and retry — rather than a
committed out-of-domain state or, in the throwing case, a dead solve. Previously a
throw from a stage ended the whole integration even though the pre-step state was
still on the stack one frame up.

Only `DomainError` is caught. `util::stop()` and everything else still propagate,
which is the point: absorbing them would turn a programming error into step-shrinking
until "Cannot achieve the desired accuracy", a diagnostic that points at the solver
instead of at the bug.

This matters for bounded quantities that a finite RK step can overshoot even when the
exact flow cannot — a carbon pool at zero, soil water at saturation, a probability at
one. It is a discretisation guard, **not** a way to fix a model whose exact flow leaves
its domain: that case shrinks to the minimum step and raises, now naming the reason and
the location rather than blaming accuracy.

Version bumped so downstreams can pin against the capability
(`odelia (>= 0.3.0)`); systems declaring neither hook are unaffected, and a Lorenz
trajectory over 4127 steps is bit-identical across the change.

Numbers in that message are rendered with `util::format_double()`, so a step size at
its floor reads `1e-08` rather than the `"0.000000"` `std::to_string` would give.


## odelia 0.2.2

**An out-of-domain interpolator lookup now says which point, how far out, and what
the domain was.** `Interpolator::eval()` had all three in hand — `u`, `min()`
and `max()` — and reported none of them:

```
Extrapolation disabled and evaluation point outside of interpolated domain.
```

That sentence is the same whichever spline threw it, so a consumer holding several
of them learns nothing about which one, and nothing about whether the point missed
the near end or the far end. It cost real time downstream: localising
[traitecoevo/plant#576](https://github.com/traitecoevo/plant/issues/576) meant
instrumenting four call sites by hand to discover which spline was being asked and
at what value, and the answer — the **lower** end, not past the far end as everyone
had assumed — inverted the fix. Now:

```
Extrapolation disabled and evaluation point outside of interpolated domain:
u = -0.0023 lies 0.0023 beyond the lower end of [0, 6.8918].
```

Which spline, and which caller, is the one thing this layer cannot know; consumers
that build several should catch and add it. A patch bump so downstreams can pin
against the message.

Behaviour is otherwise unchanged. In particular the guard is still written
`u < min() || u > max()` rather than the negation of an in-range test, because every
comparison against NaN is false and a non-finite `u` must keep falling through to the
spline — plant relies on that.



## odelia 0.2.1

A patch bump for one reason: **#46 has no version number, and a downstream needs
one.** `d8235d1` ("Let a system compute its rates when the solver reads them", #46)
landed *after* `3bdfcf7` bumped the version to 0.2.0, and the version has not moved
since — so `0.2.0` names two different header sets, one with #46 and one without, and
no `>= ` requirement can tell them apart.

That is not academic. `traitecoevo/plant`'s `develop` **does not compile** against the
released 0.2.0:

```
odelia/ode_interface.hpp:212:3: error: 'this' argument to member function 'ode_rates'
  has type 'const plant::Patch<plant::FF16_Strategy, plant::FF16_Environment>',
  but function is not marked const
```

plant #585 made `Patch::ode_rates` non-const; `r_ode_rates(const T& obj)` and
`ode_solver_internal.hpp:155` both call it on a `const&`. #46 is the fix. Eight errors,
four templated `<Strategy, Environment>` pairs × two call sites, and they surface
*inside these headers*, which points nowhere near the cause. Confirmed pre-existing by
syntax-checking plant's `origin/develop` unmodified.

This is the same job 0.2.0 was bumped for, and the 0.2.0 entry below says so in as many
words: it exists "to give downstream packages something to pin against ... so a build
against an older odelia fails at dependency resolution with a clear message rather than
at compile time". #46 needed the same courtesy and did not get it.

**No header changes.** Only `DESCRIPTION`. Downstream floors after this:

- `plant` -> `odelia (>= 0.2.1)`, because it links the ODE solver and needs #46.
- `leaf` stays at `odelia (>= 0.2.0)`. Checked rather than aligned for symmetry: leaf
  includes exactly one odelia header, `odelia/interpolator.hpp`, and never touches
  `ode_rates` or the solver. Raising its floor would force an upgrade for a fix in a
  header it does not include.

Closes #48.

## odelia 0.2.0

A minor-version bump rather than a patch, because the header core lost public
symbols: `util::index`, `util::index_vector()` and the `base_1_to_0` /
`base_0_to_1` helpers are gone, and `util::stop` / `util::warning` no longer call
into Rcpp. Nothing in the family used them — `plant` has its own `plant::util`
equivalents — but a consumer that did would fail to compile, which is exactly what
a version number is for. It also gives downstream packages something to pin
against: `leaf` now requires `odelia (>= 0.2.0)`, so a build against an older
odelia fails at dependency resolution with a clear message rather than at compile
time with `RcppCommon.h: No such file or directory`.

* The **header-only solver core is now free of R** (#43). `ode_util.hpp` no longer includes `RcppCommon.h`, so everything reachable through `interpolator.hpp`, `ode_control.hpp` and `ode_solver.hpp` compiles and runs as plain C++ with no R installation — see `tests/standalone/`, which integrates the Lorenz system on a runner that has no R on it. `util::stop()` now throws `std::runtime_error` instead of calling `Rcpp::stop()`; Rcpp converts that into an ordinary R error with the same message at the package boundary, so R-level behaviour is unchanged apart from the condition's class vector, which gains `std::runtime_error` in place of `Rcpp::exception`. `util::warning()` writes to `std::cerr` rather than raising an R warning; it had no callers. The unused `util::index` struct, its undefined `Rcpp::as`/`wrap` specializations, `util::index_vector()` and the `base_1_to_0`/`base_0_to_1` helpers are removed — `plant` has its own. R remains where it belongs, in `src/` and in `solver_interface.hpp` / `rcpp_interface_helpers.hpp`.

## odelia 0.1.0

Odelia is a new package, arising out of https://github.com/traitecoevo/plant/. In that project, Rich FitzJohn built a custom ODE solver, using a Runge-Kutta 4-5 method, in C++. I'm spinning that code out into a package, as I want to use it elsewhere.

* New implicit, adaptive-step **RODAS4(3)** Rosenbrock stepper for stiff systems, selectable via `method = "rodas"` when constructing a solver (#35). It reuses the existing adaptive step-size controller and obtains an exact Jacobian by forward-mode automatic differentiation; systems opt in by providing a `template<class U> rebind<U> rebind_from()` method (the same double->AD lift the gradient driver uses). The explicit RKCK 4(5) method (`method = "rkck"`) remains the default.

* `odelia` now loads its shared library with global symbol visibility in `.onLoad`, so packages that `LinkingTo: odelia` and instantiate `Solver` can resolve the compiled XAD runtime symbols at load time without per-package linker hacks (#26).

