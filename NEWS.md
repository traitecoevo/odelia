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

* New implicit, adaptive-step **RODAS4(3)** Rosenbrock stepper for stiff systems, selectable via `method = "rodas"` when constructing a solver (#35). It reuses the existing adaptive step-size controller and obtains an exact Jacobian by forward-mode automatic differentiation; systems opt in by providing a `template<class U> System<U> rebind()` method. The explicit RKCK 4(5) method (`method = "rkck"`) remains the default.

* `odelia` now loads its shared library with global symbol visibility in `.onLoad`, so packages that `LinkingTo: odelia` and instantiate `Solver` can resolve the compiled XAD runtime symbols at load time without per-package linker hacks (#26).

