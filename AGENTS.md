# odelia — developer guide for agents

`odelia` is an **ODE solver with automatic differentiation, implemented
in C++ header files**, with an R interface via Rcpp. Three adaptive
steppers (Cash–Karp RK4(5), Dormand–Prince 5(4), RODAS4(3) for stiff
systems) run entirely in compiled code. ODE systems are templated on
their scalar type, so a recorded run can be differentiated in **reverse
mode** (one sweep gives the derivative with respect to every parameter
and initial condition) and the implicit stepper’s Jacobian is taken in
forward mode, both through the vendored
[XAD](https://github.com/auto-differentiation/xad) library. Time-varying
external drivers are interpolated with a cubic Hermite spline. A
right-hand side written in R can be solved by the same steppers.

The solver core was first written by Rich FitzJohn inside
[`plant`](https://github.com/traitecoevo/plant); `odelia` spins it out
as a reusable library (header-only but for the XAD tape runtime, see
ARCHITECTURE.md) that other Rcpp packages link against. The
next-generation `plant` core links against it.

**Start here.** `ARCHITECTURE.md` has the map of the headers in
dependency order. `vignettes/articles/reverse-mode.Rmd` is the design of
the gradient machinery and the reading order through the headers. The
System contract is the concepts in
`inst/include/odelia/ode_interface.hpp`.

## Layout

- `inst/include/odelia/` — the **C++ core** (the solver; this is the
  reusable artifact); see the map in `ARCHITECTURE.md`.
  `inst/include/examples/lorenz_system.hpp` is the shipped example
  System; `inst/examples/leaf_thermal/` the one the article builds.
- `inst/include/XAD/` — the vendored autodiff library, with local
  patches recorded in `tools/`. `src/Tape.cpp` is the one compiled copy
  of its tape runtime.
- `src/` — Rcpp glue compiled into the package. `r_system.h` /
  `r_system_interface.cpp` are the R adapter over the callback system
  and the `OdeSolver` exports.
- `R/` — friendly **R6** wrappers around the C++ objects.
- `vignettes/` — `odelia.Rmd` (getting started),
  `parameter-fitting.Rmd`; `articles/` holds the website-only pieces
  that compile C++ (`leaf-thermal`, `reverse-mode`).
- `tests/testthat/` — the R suite; `tests/standalone/` — the R-free C++
  suite, which also holds the smallest complete sweepable System
  (`Grow`).

## Build & test (Makefile)

- `make compile` — compile C++ after C++-only changes.
- `make Rcpp` / `make roxygen` — regenerate Rcpp exports / roxygen docs
  (don’t hand-edit generated files: `R/RcppExports.R`,
  `src/RcppExports.cpp`, `NAMESPACE`, `man/`).
- `make test` — install with tests and run the full suite against the
  installed package (the AD and DLL-lifecycle tests need that).
  `make test-local` — the fast `load_all` loop, which skips those.
  `make check` — `R CMD check`.
- `make test-cpp` — build and run the core as plain C++ with no R on the
  include path (`tests/standalone/`), and compile every core header on
  its own; fast, and the guard that keeps the headers R-free.
- Everything is C++20. Anything compiled against the headers outside
  `src/` (a `sourceCpp` snippet, a consumer) needs
  `-DXAD_NO_THREADLOCAL -DXAD_USE_STRONG_INLINE` to match
  `src/Makevars`; the tests take them from `odelia_cppflags()` in
  `tests/testthat/helper-load-odelia.R`.

## Gotchas

- It compiles C++ from source — a working toolchain is required, and a
  header change can break dependents at **compile time**, not just
  runtime.
- The header core is a cross-boundary artifact: changing a solver
  signature ripples to anything that `LinkingTo` it (notably the
  next-gen `plant`). Treat such changes as `cross-package` / `breaking`.
- **The header core must stay free of R.** Everything in
  `inst/include/odelia/` bar `solver_interface.hpp` and
  `rcpp_interface_helpers.hpp` compiles with no R installed —
  `util::stop` throws, it does not call `Rcpp::stop`. Consumers depend
  on this (`leaf` runs its C++ tests without R). Keep Rcpp in `src/` and
  the two interface headers, and remember that these headers no longer
  get `<cassert>`, `<string>` and friends for free via R — include what
  you use. See ARCHITECTURE.md.

## Code & comment style

New code should be indistinguishable from the existing header core (Rich
FitzJohn’s): terse, template-heavy, `const` by default, 2-space indent.
AD code is **glue around the vendored XAD facilities**
(`computeJacobian`, `CheckpointCallback`, the tape drivers) — invoke
them, don’t re-implement them. Modify the type that already exists
rather than adding a parallel one, and prefer a mechanism that scales (a
System hands back its fields; no per-index switch) over a special case.

Comments say what the code **is** and what must hold — the invariant,
the reason behind a non-obvious choice — and nothing else. The bar the
AD surface is held to:

- **State the thing, not its history.** No issue/PR numbers, no “was
  renamed from…”, no “the old X did Y”. A stable external anchor (a
  paper, `#472`, a GSL routine) is fine; process references drift the
  moment the code moves.
- **Present odelia’s design as its own fact.** Don’t explain it via
  plant, “the spike”, or how we got here.
- **Plain and direct — no metaphor, no flourish.** Name things for what
  they are; avoid decorative nouns (`contract`, `oracle`, `surface`) and
  cute metaphors (`frozen`/`mutant`/`live`/`comb`). If a name needs a
  metaphor to make sense, rename it.
- **Be sparing.** The code carries most of the meaning; a comment earns
  its place by helping the reader over a genuine hump. Don’t narrate a
  counter for a paragraph.
- **Generic machinery is background.** In a concrete System (Lorenz,
  leaf, canopy) the members required by the AD contract should read as
  ordinary code, not as the point of the file — the physics is the
  point. Give an example a real applied domain, not an abstract
  stand-in.

The contract a System implements is stated as concepts in
`inst/include/odelia/ode_interface.hpp` – `HasOdeTime`, `Rebindable`,
`SolvesForValues`, `ChecksState` for solving, and `Sweepable` for a
reverse sweep, which `Solver::solve_adjoint` asserts – so a System that
does not satisfy one fails to compile naming the requirement it missed.
Read those rather than any prose account: a prose copy of a
compiler-checked contract drifts. The two members a widening System adds
(`apply_insertion`) and a System solving inside a stage adds
(`solved_values`) are documented beside the concepts. Every header in
`inst/include/odelia/` meets the comment standard above; keep it that
way, and state an invariant once, at its owning definition, referring to
it by name elsewhere. Don’t hand-edit generated files
(`R/RcppExports.R`, `src/RcppExports.cpp`, `NAMESPACE`, `man/`).

## Plant family

`odelia` is part of the **plant family** in the
[`traitecoevo`](https://github.com/traitecoevo) org — a hub-and-spoke
set of packages built around the
[`plant`](https://github.com/traitecoevo/plant) size- and
trait-structured forest model.

- **Docs hub** — family user guides & theory:
  <https://traitecoevo.github.io/overstorey/>
- **Cross-package orientation** — how the family fits together (who
  depends on whom, source-of-truth rules, cross-repo gotchas) lives in
  [`plant-meta`](https://github.com/traitecoevo/plant-meta); start with
  its
  [`AGENTS.md`](https://github.com/traitecoevo/plant-meta/blob/main/AGENTS.md).
  Keep family-wide concerns there, not here.
- **Issues & board** — follow the [issue
  guide](https://github.com/traitecoevo/plant-meta/blob/main/governance/issue-guide.md);
  work is tracked on [board
  \#5](https://github.com/orgs/traitecoevo/projects/5) (new issues
  auto-add with no Status = the triage queue). Labels: `bug` / `task` /
  `epic` plus `blocked`, `needs-info`, `cross-package`, `breaking`,
  `question`.
- **Commit messages** — the repo squash-merges, so a PR’s title and body
  are copied verbatim into permanent history. Keep them short and
  durable, and put the working detail (measurements, alternatives
  rejected, what you tried first) in the first PR comment instead — see
  [`commit-messages.md`](https://github.com/traitecoevo/plant-meta/blob/main/governance/commit-messages.md).
