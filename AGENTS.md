# odelia — developer guide for agents

`odelia` is an **ODE solver with automatic differentiation, implemented
in C++ header files**, with an R interface via Rcpp. It uses an
adaptive-step Runge–Kutta 4–5 integrator that runs entirely in compiled
code, and templates ODE systems on their scalar type so solutions can be
differentiated w.r.t. parameters and initial conditions (forward-mode AD
via the vendored [XAD](https://github.com/auto-differentiation/xad)
library). It also supports time-varying external drivers (cubic-spline
interpolated).

The solver core was first written by Rich FitzJohn inside
[`plant`](https://github.com/traitecoevo/plant); `odelia` spins it out
as a reusable, header-only library that other Rcpp packages can link
against. The next-generation `plant` core links against it.

## Layout

- `inst/include/` — the **header-only C++ core** (the solver; this is
  the reusable artifact).
- `src/` — Rcpp glue compiled into the package.
- `R/` — friendly **R6** wrappers around the C++ objects.
- `tools/`, `vendor` (XAD) — the vendored autodiff library.
- `ARCHITECTURE.md` — read this for the C++ design; `vignettes/` for
  worked examples (e.g. the Lorenz benchmark).

## Build & test (Makefile)

- `make compile` — compile C++ after C++-only changes.
- `make Rcpp` / `make roxygen` — regenerate Rcpp exports / roxygen docs
  (don’t hand-edit generated files: `R/RcppExports.R`,
  `src/RcppExports.cpp`, `NAMESPACE`, `man/`).
- `make test` — run the test suite (`testthat`). `make check` —
  `R CMD check`.
- `make test-cpp` — build and run the core as plain C++ with no R on the
  include path (`tests/standalone/`); fast, and the guard that keeps the
  headers R-free.

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
