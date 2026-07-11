# Automatic differentiation in odelia: exact gradients of an ODE run

**Scope:** `traitecoevo/odelia` — the reverse-mode automatic-differentiation API. A
new developer with only this code should come away understanding what a gradient is made
of, who owns what, the two axes every workflow is a point on, the contract a System
implements to be differentiable, and the one way to get an adaptive-component gradient
wrong. Self-contained.

**Companion:** [`ARCHITECTURE.md`](ARCHITECTURE.md) — how the XAD `Tape` runtime is
compiled once and linked across the DLL boundary. Read it before touching the build;
this document is about the API, not the link.

---

## The one-sentence design

odelia differentiates **any reduction of an ODE run** — on one scalar-templated
`Solver`, with reverse-mode XAD compiled once — and only `double` ever crosses back to R.

Reverse mode (a *tape* records each operation, then sweeps backward to get every input's
derivative from one output) is the right engine here because there are many inputs
(parameters, initial conditions) and few outputs: one backward sweep yields the whole
gradient. XAD is the vendored AD library; odelia compiles its `Tape` once and is glue
around it, never a second AD engine.

Everything else follows: the System is templated on its scalar so the same code simulates
(`double`) and differentiates (`active`); the thing differentiated is a caller's
functional, not a built-in loss; adaptive numerics are recorded once and replayed on fixed
nodes so the tape carries no adaptive branching; and the active machinery is born, used,
and destroyed inside one C++ call.

---

## The mental model: two orthogonal axes

Most of the apparent complexity dissolves once you see that a gradient workflow is a choice
on **two independent axes**:

- **Replay** — which adaptive constructions the *system* must record on the double pass and
  replay *fixed* on the active pass. A property of the system, not the question. A plain
  ODE needs only the step schedule (L1). A system whose rates read an adaptively-built
  background also records that background's node positions (L2), and optionally its values
  (L3) for a variant that holds the background fixed.
- **Functional** — what scalar you reduce the replayed run to: an emergent summary of the
  state, or a likelihood over measured data. Orthogonal to replay.

Three representative points:

| Application | What you differentiate | Replay | Functional |
|---|---|---|---|
| **Sensitivity** | d(output)/d(param or IC) of a plain ODE | L1 | any reduction |
| **Adaptive-component gradient** | the same, when a rate reads an adaptively-built background | L1·L2 (+L3 to freeze it) | any reduction |
| **Calibration / inference** | fit params to measured data | whatever the system needs | a likelihood |

The first two differ only on the replay axis (how complex the system is); the third differs
only on the functional axis. Calibration is not a different replay — it is the system's own
replay with a likelihood functional on top.

---

## Ownership

R holds one object — the ordinary (double) `Solver`. Everything active hangs off it and
never surfaces:

```
R  ── holds ──▶  d : Solver<System<double>>          the double solver
                 │
                 ├─ System<double>                    immutable after the adaptive pass
                 ├─ recording                          produced once, read per call:
                 │    ├─ L1  recorded_steps()          the discovered schedule   (Solver state)
                 │    ├─ L2  positions_history[step]   adaptive node positions   (System state)
                 │    └─ L3  values_history[step][stage] recorded background values (System state)
                 └─ active_solver : Solver<System<active>>
                      built once, reused; owns its tape; re-seeded and re-fed the
                      recording every call; holds no semantics between calls
```

Three ownership rules carry the whole design:

- **The Solver owns the schedule (L1); the System owns its background (L2/L3).** The stepper
  steps any system either adaptively (discovering node positions) or on a fixed grid
  (replaying them). What a rate reads is the System's own business; the stepper only signals
  cadence through the hooks.
- **The double solver owns the recording, and it is immutable after the adaptive pass.** Its
  only remaining job is to feed replays. The recording is keyed to that run's ICs and
  params; change those and you must re-record. The active solver never records — it is
  handed the recording, per call.
- **The Solver owns the AD scratch.** The active solver and its tape live on the double
  `Solver` object (not an R handle), so a C++ caller holding the solver as a plain member
  gets the reuse for free. Reuse is pure speed — it never changes a number. The recording,
  read per call, is what carries semantics.

## Control flow of a gradient

```
1. d.advance_adaptive({0, T})          discover the schedule, record L1/L2/L3.
                                        d is immutable hereafter.
2. compute_gradient(d, targets, schedule, functional):
     active = d.active_solver()         lift-or-reuse
     feed active the recording
     tape on; per output row:
         seed targets → active.system   (ad_parameters / ad_initial_state)
         active.reset()
         active.advance_fixed(schedule)  ◀── the DRIVER owns the replay
         functional(active)              ◀── a PURE REDUCTION: reads state, returns scalar(s)
     sweep adjoints → gradient
     tape off → return doubles
```

The functional never drives the solver and never carries the schedule. `recorded_steps()`
is the single source of the replay grid, so it can't go inconsistent, and the "forgot to
record" guard is one check at the single place the replay happens (`schedule.empty()`).

---

## The System contract

A System is an ODE right-hand side plus the state it integrates. To simulate, it provides
the ODE interface (`ode_size`, `set_ode_state`, `ode_rates`, `ode_state`). To be
**differentiable**, it adds:

```cpp
using value_type = S;                                  // the scalar it carries
template <class S2> using rebind = System<…, S2>;      // the double -> active mould
template <class S2> rebind<S2> rebind_from() const;    // copy config (values only) into S2
std::vector<S*> ad_parameters();                       // handles to the active parameters
std::vector<S*> ad_initial_state();                    // handles to the active initial state
```

That is the whole contract for a bare ODE. `rebind_from` copies configuration values only,
so the active system starts free of tape identity; the driver then seeds a chosen subset of
inputs by writing through the handles `ad_parameters()` / `ad_initial_state()` return (in a
fixed order, so a 30-parameter System is one `return {&a, &b, …}` rather than a switch).
Nothing forces a System to be differentiable — one without these still simulates.
(`LorenzSystem` is the worked example.)

## Differentiation targets and the drivers

```cpp
compute_jacobian(solver, targets, schedule, functional);  // m x n, row i = d(out_i)/d(in)
compute_gradient(solver, targets, schedule, functional);  // the one-row case
```

**`targets`** names the inputs directly — no opaque flat-slot space:

```cpp
struct DifferentiationTargets {
  std::vector<int>    params;   // param indices to seed active
  std::vector<int>    ics;      // initial-state indices to seed active
  std::vector<double> values;   // seed values, params-then-ics
};
```

Seed only `params`, only `ics`, or both. The driver writes each through the handles
`ad_parameters()` / `ad_initial_state()` return.

> **Column order is a contract.** Jacobian column *j* is `d(output)/d(input_j)` for the
> *j*-th seeded input, in `values` order (params then ics). A caller that resolves names →
> indices must seed in the same order it reads columns, or the columns transpose silently.

## Functionals: what you differentiate

A functional is **any callable** that reads an already-replayed `Solver` and returns the
scalar(s) to differentiate, plus a `codomain()` reporting how many outputs it returns. It
does not drive the solver.

The whole final state — `m = ode_size` outputs:

```cpp
struct final_state {
  std::size_t codomain() const { return m; }
  template <class Solver>
  std::vector<typename Solver::value_type> operator()(Solver& s) const {
    return s.state();
  }
};
```

A scalar summary — the summed final state:

```cpp
struct sum_final_state {
  std::size_t codomain() const { return 1; }
  template <class Solver>
  typename Solver::value_type operator()(Solver& s) const {
    typename Solver::value_type total(0);
    for (auto const& x : s.state()) total += x;
    return total;
  }
};
```

Anything you can compute from the replayed state is a valid functional: one component, a
weighted sum, a nonlinear metric, a residual against data. You write a struct with
`codomain()` and `operator()(solver)`; odelia differentiates whatever it returns and never
learns what the scalars mean. `least_squares` is just a prebuilt one — a scalar loss that
holds measured data and scores the trajectory against it. Reading the count off the
functional (rather than passing it) means the count and the outputs cannot silently
disagree, and it saves XAD an extra forward callback — a full wasted replay — just to size
the output. The record-once / row-sweep is the vendored `xad::computeJacobian`'s; odelia
supplies only the forward callback.

---

## Record → replay: any adaptive component

A bare ODE needs nothing more than the above. But if a rate reads a background built by an
**adaptive** construction — one that makes parameter-dependent decisions about where to
place nodes — differentiating through those decisions corrupts the tape. On a gradient pass
the adaptive run has already happened, so the placement is already known: record it once,
replay pinned to it.

The **interpolator** is the running example. On the double pass it refines its knots
adaptively; on the active pass it rebuilds on the recorded knots with active values. But
nothing here is interpolator-specific: a quadrature's subdivision, or any adaptive refiner,
records through the identical hooks. The System records whatever *positions* its adaptive
machinery chose; odelia never learns a node is a knot. A System opts in with the
`Replayable` hooks:

| Hook | Fires | Job |
|---|---|---|
| `record_stage(k)` | per RK stage, record pass | stash node positions / this stage's background value |
| `record_ode_step()` | per accepted step, record pass | commit the step's recording |
| `replay_step()` | per step, active pass | load this step's recorded slice |
| `has_recorded_field()` | query | is the L3 background cache populated? |

Two cadences, and they are the whole model:

- **Positions are recorded per step** (L2). The only job is to keep the adaptive branching
  off the tape; one representative position set per accepted step does it, and it is
  representative by construction — the step controller only accepts a step whose state (hence
  the background's structure) varies within tolerance across it.
- **Values are recorded per stage** (L3). A background held fixed must read back the *exact*
  value each RK stage consumed — a per-step value would be wrong, not merely coarse.

**Whether a background's derivative flows is a data question.** `has_recorded_field()` asks
one thing: are recorded background values present (the L3 cache)?

- **empty** → the System recomputes the background on the recorded positions with the active
  scalar, so its feedback into the rates is differentiated (the derivative flows through it).
- **populated** → the System reads the recorded values as `double`, off the tape, so that
  background's derivative is zero by construction — you have held it fixed.

The caller chooses per replay whether to populate the cache. Reach for the hold-fixed
(populated) variant when a background is a shared or coupling quantity you deliberately want
to hold constant — its contribution to the derivative is then zero, not merely small.

> **The one way to get it wrong.** If you want a background's feedback in the gradient, you
> must leave L3 *empty* so it is recomputed with the active scalar on the recorded positions.
> Populating it — or otherwise reading the recorded values — drops that background's
> contribution to zero *silently*: the run succeeds and returns a plausible number that is
> missing the feedback cross term. The positions are fixed so the tape stays clean; the
> values are recomputed so the gradient flows through *what the nodes hold*, never through
> *where they sit*.

---

## Injecting a known derivative

**When you need it.** Sometimes the forward pass computes a value the tape never recorded —
an iterative solver, a root-find, an optimizer result. Differentiating *through* the
iteration is wasteful (it tapes every iterate) or impossible (the iteration count is
data-dependent), but you know the result's derivative analytically — typically from the
implicit function theorem. `SuppliedDerivative` lets you register the off-tape value as a
fresh input and hand the reverse sweep its known partials, so the internal solve is never
recorded.

**Example.** A root-find gives `x` solving `x = cos(x) + a`, and you differentiate
`g = x²` w.r.t. `a`. Nothing about the Newton loop is on the tape. By the implicit function
theorem, `dx/da = 1 / (1 + sin(x))`. Register `x` as a leaf carrying that partial; then
`g = x*x` records normally and the reverse sweep returns `dg/da = 2x · dx/da`. Built on
`xad::CheckpointCallback`, it is a free function called from *within* the forward pass —
where the off-tape value exists — not a pre-declared target. odelia sees only inputs, an
output value, and partials.

---

## The R boundary

Only `double` crosses to R. R holds the double `Solver`; a gradient call builds the active
solver internally, differentiates, and returns doubles — R never holds an active type, so
the boundary can only fail loudly, never reinterpret a handle as the wrong scalar. The
active solver and its tape are cached on the double `Solver` and reused across calls, so an
optimiser loop (asking for value and gradient each iteration) amortizes them.

---

## Building an AD-compatible System

**A bare ODE (no adaptive sub-numerics).** Implement the ODE interface plus the four
contract members — `value_type`, `rebind` / `rebind_from`, `ad_parameters`,
`ad_initial_state`. That's it: `compute_gradient` / `compute_jacobian` work, doubles in
and out, and a downstream package gets `Solver_gradient` / `Solver_jacobian` for free.

The complete member set:

```cpp
template <class S = double>
class MySystem {
public:
  using value_type = S;

  // Simulate: the ODE interface.
  std::size_t ode_size() const;                          // number of state variables
  double      ode_time() const;                          // current time
  template <class It> It set_ode_state(It y, double t);  // load state at time t; recompute rates
  template <class It> It ode_state(It out) const;        // write the current state
  template <class It> It ode_rates(It out) const;        // write the current dy/dt
  void reset();                                          // return to the initial state

  // Differentiate: copy onto another scalar, and expose the inputs to seed active.
  template <class S2> using rebind = MySystem<S2>;
  template <class S2> rebind<S2> rebind_from() const;    // xad::value the config into an S2 copy
  std::vector<S*> ad_parameters();                       // addresses of the active parameters
  std::vector<S*> ad_initial_state();                    // addresses of the active initial state
};
```

`LorenzSystem` is the worked example to copy from.

**Rates that read an adaptively-refined background.** Add the `Replayable` hooks. Record the
node **positions** your refiner chose (per accepted step) and the background **value** per
RK stage; on replay, rebuild on the recorded positions with the active scalar. This is
switched on by *doing reverse-mode AD*, not by any user flag — miss it and gradients are
silently wrong wherever the adaptive component bites.

```cpp
void record_stage(int stage);     // per RK stage, record pass: stash a value
void record_ode_step();           // per accepted step, record pass: commit the step
void replay_step();               // per step, replay pass: load the recorded slice
bool has_recorded_field() const;  // are recorded values present to reuse?
```

`CanopySystem` is the worked example.

**A variant that holds a background fixed.** If you also want a cheap run that reads some
recomputable quantity as a constant (its derivative zero) — a coupling field, held-fixed
state variables — record it per stage and have the variant entry populate the L3 cache;
`has_recorded_field()` reports it. L3 is "read a recorded `double` from a container"; the
quantity need not be an interpolated field.

odelia stays agnostic to all of it — it records "some positions" and "some values" and never
learns a node is a height or a value a field.
