# `odelia` — developer guide for agents

`odelia` is an **ODE solver with automatic differentiation, implemented in C++ header
files**, with an R interface via Rcpp. It uses an adaptive-step Runge–Kutta 4–5 integrator
that runs entirely in compiled code, and templates ODE systems on their scalar type so
solutions can be differentiated w.r.t. parameters and initial conditions (forward-mode AD
via the vendored [XAD](https://github.com/auto-differentiation/xad) library). It also
supports time-varying external drivers (cubic-spline interpolated).

The solver core was first written by Rich FitzJohn inside [`plant`](https://github.com/traitecoevo/plant);
`odelia` spins it out as a reusable, header-only library that other Rcpp packages can link
against. The next-generation `plant` core links against it.

## Layout

- `inst/include/` — the **header-only C++ core** (the solver; this is the reusable artifact).
- `src/` — Rcpp glue compiled into the package.
- `R/` — friendly **R6** wrappers around the C++ objects.
- `tools/`, `vendor` (XAD) — the vendored autodiff library.
- `ARCHITECTURE.md` — read this for the C++ design; `vignettes/` for worked examples
  (e.g. the Lorenz benchmark).

## Build & test (Makefile)

- `make compile` — compile C++ after C++-only changes.
- `make Rcpp` / `make roxygen` — regenerate Rcpp exports / roxygen docs (don't hand-edit
  generated files: `R/RcppExports.R`, `src/RcppExports.cpp`, `NAMESPACE`, `man/`).
- `make test` — run the test suite (`testthat`). `make check` — `R CMD check`.

## Gotchas

- It compiles C++ from source — a working toolchain is required, and a header change can
  break dependents at **compile time**, not just runtime.
- The header core is a cross-boundary artifact: changing a solver signature ripples to
  anything that `LinkingTo` it (notably the next-gen `plant`). Treat such changes as
  `cross-package` / `breaking`.

## Code & comment style

New code should be indistinguishable from the existing header core (Rich FitzJohn's):
terse, template-heavy, `const` by default, 2-space indent. AD code is **glue around the
vendored XAD facilities** (`computeJacobian`, `CheckpointCallback`, the tape drivers) —
invoke them, don't re-implement them. Modify the type that already exists rather than
adding a parallel one, and prefer a mechanism that scales (a System hands back its fields;
no per-index switch) over a special case.

Comments say what the code **is** and what must hold — the invariant, the reason behind a
non-obvious choice — and nothing else. The bar the AD surface is held to:

- **State the thing, not its history.** No issue/PR numbers, no "was renamed from…", no
  "the old X did Y". A stable external anchor (a paper, `#472`, a GSL routine) is fine;
  process references drift the moment the code moves.
- **Present odelia's design as its own fact.** Don't explain it via plant, "the spike", or
  how we got here.
- **Plain and direct — no metaphor, no flourish.** Name things for what they are; avoid
  decorative nouns (`contract`, `oracle`, `surface`) and cute metaphors
  (`frozen`/`mutant`/`live`/`comb`). If a name needs a metaphor to make sense, rename it.
- **Be sparing.** The code carries most of the meaning; a comment earns its place by
  helping the reader over a genuine hump. Don't narrate a counter for a paragraph.
- **Generic machinery is background.** In a concrete System (Lorenz, leaf, canopy) the
  members required by the AD contract should read as ordinary code, not as the point of
  the file — the physics is the point. Give an example a real applied domain, not an
  abstract stand-in.

`AUTODIFF.md` is the reference for the AD surface a System implements; `ARCHITECTURE.md`
for the XAD `Tape` link. Don't hand-edit generated files (`R/RcppExports.R`,
`src/RcppExports.cpp`, `NAMESPACE`, `man/`).

## Plant family

`odelia` is part of the **plant family** in the [`traitecoevo`](https://github.com/traitecoevo)
org — a hub-and-spoke set of packages built around the
[`plant`](https://github.com/traitecoevo/plant) size- and trait-structured forest model.

- **Docs hub** — family user guides & theory: <https://traitecoevo.github.io/overstorey/>
- **Cross-package orientation** — how the family fits together (who depends on whom,
  source-of-truth rules, cross-repo gotchas) lives in
  [`plant-meta`](https://github.com/traitecoevo/plant-meta); start with its
  [`AGENTS.md`](https://github.com/traitecoevo/plant-meta/blob/main/AGENTS.md). Keep
  family-wide concerns there, not here.
- **Issues & board** — follow the
  [issue guide](https://github.com/traitecoevo/plant-meta/blob/main/governance/issue-guide.md);
  work is tracked on [board #5](https://github.com/orgs/traitecoevo/projects/5) (new issues
  auto-add with no Status = the triage queue). Labels: `bug` / `task` / `epic` plus `blocked`,
  `needs-info`, `cross-package`, `breaking`, `question`.
