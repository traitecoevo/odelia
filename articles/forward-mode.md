# Forward mode: few parameters, many outputs

Fitting a small ODE model to a time series needs `dy(t)/dθ` at every
observation time, for a handful of parameters. [Reverse
mode](https://traitecoevo.github.io/odelia/articles/reverse-mode.md)
gives the derivative of one output with respect to every input from one
sweep; forward mode gives the derivative of every output with respect to
one input from one pass. With few inputs and many outputs, forward mode
is the cheaper of the two. It is also the only way `odelia`
differentiates a run of the implicit stepper, RODAS, which cannot be
swept.

`odelia` does it by walking a run’s recorded steps a second time, with
the System on a scalar that carries one lane per direction. Nothing a
lane carries reaches a step’s size or whether a step is taken, so the
derivative is that of the scheme as the run stepped it. The sweep keeps
the same convention, so over the same run forward and reverse give the
same numbers to rounding.

> This is a website-only article rather than a packaged vignette,
> because it compiles C++ at build time.

## The whole of it

``` r

library(odelia)
Sys.setenv(PKG_CPPFLAGS = paste0(
  "-I", shQuote(system.file("include", package = "odelia")),
  " -DXAD_NO_THREADLOCAL -DXAD_USE_STRONG_INLINE"))
```

``` cpp
// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <odelia/ode_solver.hpp>
#include <examples/lorenz_system.hpp>

using namespace odelia::ode;

// d (x, y, z)(t) / d (sigma, R, b, x(0)) at each of `times`, from one solve and
// one forward pass, and the same derivative of x(t_end) from a reverse sweep.
// [[Rcpp::export]]
Rcpp::List lorenz_forward(std::vector<double> times, std::string method) {
  LorenzSystem<double> system(10.0, 28.0, 8.0 / 3.0);
  const Method m = method == "rodas" ? Method::rodas : Method::dopri;
  Solver<LorenzSystem<double>> s(system, OdeControl(), m);
  s.set_keep_states(true);
  s.advance_adaptive(times);                    // a step lands on every time

  const std::vector<std::size_t> rows = s.rows_at(times);   // the row of each time
  tangent_directions wrt{{0, 1, 2}, {0}};       // sigma, R, b; and x(0)
  std::vector<tangent_row> d = s.sweep_tangent<4>(wrt, rows);

  Rcpp::NumericMatrix dx(rows.size(), 4);       // d x(t) / d each input
  for (std::size_t k = 0; k < rows.size(); ++k) {
    for (std::size_t j = 0; j < 4; ++j) dx(k, j) = d[k].at(0, j);
  }
  Rcpp::NumericVector swept(4, NA_REAL);
  if (m != Method::rodas) {
    adjoint_batch lambda = adjoint_batch::one({1.0, 0.0, 0.0});
    adjoint_batch dp(1, 3);
    s.sweep_adjoint(lambda, dp);
    swept = Rcpp::NumericVector::create(dp[0][0], dp[0][1], dp[0][2], lambda[0][0]);
  }
  return Rcpp::List::create(Rcpp::Named("dx") = dx, Rcpp::Named("swept") = swept);
}
```

Against an ordinary solve, three things are new: the record is kept, the
directions are named, and `sweep_tangent<W>` is called. `W` is how many
directions one pass carries; with more directions than lanes it walks
again for the rest. `tangent_directions` names parameters by their place
in the System’s `ad_parameters()`, then components of the state the run
started from. Each answer row holds the state at a recorded row and
`at(i, j)`, the derivative of component `i` along direction `j`.

``` r

r <- lorenz_forward(c(0, 0.5, 1), "dopri")
dimnames(r$dx) <- list(c("t=0", "t=0.5", "t=1"), c("sigma", "R", "b", "x0"))
r$dx
#>            sigma          R         b        x0
#> t=0    0.0000000  0.0000000  0.000000  1.000000
#> t=0.5 -1.8984006 -1.0452473  1.502801 -4.769701
#> t=1   -0.1681545  0.3099964 -1.717643  0.459014
r$swept               # the last row again, from one reverse sweep
#> [1] -0.1681545  0.3099964 -1.7176430  0.4590140
r$dx[3, ] - r$swept
#>         sigma             R             b            x0 
#> -1.471046e-15 -1.498801e-15 -1.243450e-14 -4.940492e-15
```

At `t = 0` the only nonzero derivative is that of `x` with respect to
its own starting value. At the last row, forward and reverse agree to
rounding.

## Through RODAS

The implicit stepper keeps no per-stage record, so it cannot be swept.
Forward mode runs through it unchanged: each RODAS step is taken at the
lane scalar, and the step’s Jacobian, which RODAS takes by forward-mode
AD on the System, is then a tangent of a tangent that carries the lanes
too. The derivative is the RODAS step’s own, and it converges with the
solution as the tolerance tightens.

``` r

rodas <- lorenz_forward(c(0, 0.5, 1), "rodas")$dx
signif(rodas[3, ] - r$dx[3, ], 2)    # RODAS against Dormand-Prince, each at its own steps
#>    sigma        R        b       x0 
#> -2.1e-08 -1.0e-08 -4.3e-07 -1.4e-07
```

The two differ by about the integration error of the two runs, as two
solutions do. A gradient-based fit needs the derivative of the loss it
evaluates, which is what each gives for its own method.

## What it costs

A forward pass costs about one solve, plus a share for each lane: every
operation computes the value once and then each lane. A reverse sweep
re-runs each step at the adjoint scalar inside a recording and walks it
backwards, at a cost that does not depend on the number of parameters.
The table below times both on a forced, damped oscillator with `P`
forcing coefficients and 200 output times, as multiples of the plain
solve. Forward mode also needs that plain solve first, to fix the steps,
so add one to its numbers.

``` cpp
// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <chrono>
#include <odelia/ode_solver.hpp>

using namespace odelia::ode;

// x'' + c x' + k x = sum_j theta_j sin(j t), with P forcing coefficients.
template <typename T = double>
struct Forced {
  using value_type = T;
  std::vector<T> theta;
  T c = T(0.3), k = T(2.0), x = T(1.0), v = T(0.0), dx = T(0.0), dv = T(0.0);
  double time = 0.0;
  explicit Forced(std::size_t P = 1) : theta(P, T(0.1)) { compute_rates(); }
  template <class S> Forced<S> rebind() const {
    Forced<S> o(theta.size());
    for (std::size_t j = 0; j < theta.size(); ++j) o.theta[j] = S(xad::value(theta[j]));
    std::vector<S> y{S(xad::value(x)), S(xad::value(v))};
    o.set_ode_state(y.begin(), time);
    return o;
  }
  std::size_t ode_size() const { return 2; }
  double ode_time() const { return time; }
  template <class It> It set_ode_state(It it, double t) {
    x = *it++; v = *it++; time = t; compute_rates(); return it;
  }
  void set_recorded_state(const std::vector<T>& y, double t) { set_ode_state(y.begin(), t); }
  void compute_rates() {
    T f = T(0.0);
    for (std::size_t j = 0; j < theta.size(); ++j) f += theta[j] * std::sin((j + 1.0) * time);
    dx = v;
    dv = -c * v - k * x + f;
  }
  template <class It> It ode_state(It it) const { *it++ = x; *it++ = v; return it; }
  template <class It> It ode_rates(It it) const { *it++ = dx; *it++ = dv; return it; }
  std::vector<T*> ad_parameters() {
    std::vector<T*> p;
    for (T& t : theta) p.push_back(&t);
    return p;
  }
  template <class F> void for_each_active(F&& f) {
    for (T& t : theta) f(t);
    f(c); f(k); f(x); f(v); f(dx); f(dv);
  }
};

template <class F> double fastest(F f) {
  double best = 1e300;
  for (int r = 0; r < 7; ++r) {
    const auto t0 = std::chrono::steady_clock::now();
    f();
    best = std::min(best, std::chrono::duration<double>(
                              std::chrono::steady_clock::now() - t0).count());
  }
  return best;
}

// Reverse (one sweep seeded with a sum of squares' adjoints, interval by
// interval) and forward at 8 lanes, each over the plain solve.
// [[Rcpp::export]]
Rcpp::NumericVector forced_costs(int P) {
  std::vector<double> times;
  for (int i = 0; i <= 200; ++i) times.push_back(0.1 * i);
  OdeControl control;
  control.set_tol_abs(1e-8);
  control.set_tol_rel(1e-8);
  Forced<double> system(P);
  const double plain = fastest([&] {
    Solver<Forced<double>> r(system, control, Method::dopri);
    r.set_collect(false);
    r.set_keep_states(true);
    r.advance_adaptive(times);
  });
  Solver<Forced<double>> s(system, control, Method::dopri);
  s.set_collect(false);
  s.set_keep_states(true);
  s.advance_adaptive(times);
  const auto rec = s.recording();
  const std::vector<std::size_t> rows = s.rows_at(times);
  std::vector<adjoint_batch> at;                // 2 y(t_i): a sum of squares' adjoints
  for (const std::size_t k : rows) {
    at.push_back(adjoint_batch::one({2 * rec[k].state[0], 2 * rec[k].state[1]}));
  }
  const double reverse = fastest([&] {
    adjoint_batch lambda, dp(1, P);
    s.sweep_adjoint_at(rows, at, lambda, dp);
  });
  tangent_directions wrt;
  for (int j = 0; j < P; ++j) wrt.parameters.push_back(j);
  const double forward = fastest([&] { s.sweep_tangent<8>(wrt, rows); });
  return Rcpp::NumericVector::create(Rcpp::Named("reverse") = reverse / plain,
                                     Rcpp::Named("forward") = forward / plain);
}
```

``` r

P <- c(1, 2, 4, 8, 16, 32, 64)
costs <- t(vapply(P, forced_costs, numeric(2)))
data.frame(P = P, reverse = round(costs[, "reverse"], 1),
           forward = round(costs[, "forward"], 1))
#>    P reverse forward
#> 1  1     3.8     2.2
#> 2  2     3.8     2.1
#> 3  4     3.3     1.7
#> 4  8     2.6     1.4
#> 5 16     2.3     2.5
#> 6 32     2.1     4.4
#> 7 64     1.7     7.6
```

On a small System the sweep costs several plain solves whatever `P` is.
Forward mode starts near one, grows with `P`, and passes the sweep at a
few tens of parameters. A larger state, or rates that cost more against
the bookkeeping, moves where they meet, so time both on the model in
hand. Choose the lane width at or just above the number of directions:
idle lanes still cost, and a width below it means further passes.

## From R

A compiled solver exposes both:

``` r

lz <- LorenzSystem$new(10, 28, 8 / 3)
s <- Lorenz_Solver$new(lz$ptr, OdeControl$new()$ptr, method = "rodas")
out <- s$sensitivities(c(0, 0.5, 1), parameters = c("sigma", "R"), initial = "x")
out$dy[3, , ]
#>        sigma         R         x0
#> x -0.1681545 0.3099964 0.45901384
#> y -0.5792739 0.3131371 1.01430284
#> z -0.2646244 0.4461349 0.07709175
```

`$sensitivities(times, parameters, initial)` solves from where the
solver stands, without moving it, and returns `y` and `dy`, an array
over time, state component and input.
`$loss_and_gradient(mode = "forward")` takes a least-squares loss’s
gradient the same way over the target’s schedule; it is the default for
a `"rodas"` solver, where `mode = "reverse"` is refused.

## What a System needs

A System walked forward needs less than one swept in reverse. The
concept `TangentSweepable` in `ode_interface.hpp`, which `sweep_tangent`
asserts, asks for `rebind<U>()` to the lane scalar and
`set_recorded_state(y, time)`, which is where the pass seats the
directions it starts from. A direction on a parameter also needs
`ad_parameters()`, and the rates must read each parameter live, as for a
sweep. A System that widens is walked across its insertions as the run
took them.

Under RODAS the System also declares
`static constexpr bool nests_tangents = true;` (`NestsTangents`), unless
it supplies its own `ode_jacobian()`: the step’s Jacobian is then a
tangent of a tangent, and whether a System’s rates compile at that
scalar shows only when they are instantiated, so it is declared rather
than detected. The Lorenz example declares it.

Two limits:

- A System that solves for values inside a stage (`SolvesForValues`) is
  refused at compile time. A forward pass would solve them again at the
  lane scalar instead of loading what the run solved for, and
  differentiate a function the run never computed.
- Under RODAS, a System’s own `ode_jacobian()` hook is called at the
  lane scalar. One written over the System’s `value_type`, as
  `fd_jacobian()` is, carries the lanes. One that computes in `double`
  drops them, and the pass then differentiates a step whose Jacobian
  ignores the inputs. Leave the Jacobian to AD, or write the hook over
  `value_type`.
