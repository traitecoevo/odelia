/* Multirate (MRI-GARK) interface for the two-rate demonstrator.
 *
 * Exposes the model-agnostic MRI macro stepper on the TwoRateSystem example so
 * the native multirate capability can be checked from R: a macro-stepped run, a
 * single-rate adaptive reference to check it against, and a reverse-mode gradient
 * through the stepper by record->replay, validated against a frozen-schedule
 * finite difference.
 */

#include <Rcpp.h>
#include <XAD/XAD.hpp>
#include <odelia/ode_solver.hpp>
#include <odelia/mri.hpp>
#include <examples/two_rate_system.hpp>
#include <examples/drainage_system.hpp>

using namespace Rcpp;
using namespace odelia;
using namespace odelia::ode;

static MRICoupling pick_table(const std::string& name) {
  if (name == "forward_euler") return mri_forward_euler();
  if (name == "heun")          return mri_heun();
  if (name == "kutta3")        return mri_kutta3();
  if (name == "erk33a")        return mri_erk33a();
  Rcpp::stop("unknown MRI table: " + name);
}

static OdeControl make_control(double tol) {
  return OdeControl(tol, tol, 1.0, 0.0, 1e-10, 1e10, 1e-4);
}

template <class State>
static Rcpp::NumericMatrix to_matrix(const std::vector<State>& h) {
  const int nt = (int)h.size(), nd = h.empty() ? 0 : (int)h[0].size();
  Rcpp::NumericMatrix out(nt, nd);
  for (int i = 0; i < nt; ++i)
    for (int j = 0; j < nd; ++j) out(i, j) = xad::value(h[i][j]);
  return out;
}

// [[Rcpp::export]]
Rcpp::List two_rate_mri(double k, int n_slow, std::string table,
                        Rcpp::NumericVector macro_times, double tol) {
  TwoRateSystem<double> sys(k, n_slow);
  MRICoupling M = pick_table(table);
  OdeControl control = make_control(tol);
  std::vector<double> times(macro_times.begin(), macro_times.end());
  MRISchedule sched;
  auto hist = mri_advance(sys, M, control, times, sched, /*replay=*/false);
  return Rcpp::List::create(_["states"] = to_matrix(hist),
                            _["n_fast"] = (double)mri_fast_steps(sched),
                            _["order"]  = M.order);
}

// Drive the two-rate system through the Solver `method` surface: method="mri"
// selects the multirate stepper exactly as "rkck"/"rodas" select theirs, so this
// is the path the SCM uses. Fixed macro grid (the forcing-kink grid).
// [[Rcpp::export]]
Rcpp::List two_rate_solver(double k, int n_slow, std::string method,
                           Rcpp::NumericVector times, double tol) {
  TwoRateSystem<double> sys(k, n_slow);
  OdeControl control = make_control(tol);
  Method m = method == "mri"   ? Method::mri
           : method == "rodas" ? Method::rodas
           : method == "imex"  ? Method::imex
                               : Method::rkck;
  Solver<TwoRateSystem<double>> solver(sys, control, m);
  std::vector<double> t(times.begin(), times.end());
  solver.advance_fixed(t);
  auto hist = solver.get_history();
  const int nt = (int)hist.size(), nd = (int)sys.ode_size();
  Rcpp::NumericMatrix out(nt, nd);
  for (int i = 0; i < nt; ++i) {
    std::vector<double> s(nd);
    hist[i].ode_state(s.begin());
    for (int j = 0; j < nd; ++j) out(i, j) = s[j];
  }
  return Rcpp::List::create(_["states"] = out);
}

// [[Rcpp::export]]
Rcpp::List two_rate_reference(double k, int n_slow, Rcpp::NumericVector times, double tol) {
  TwoRateSystem<double> sys(k, n_slow);
  OdeControl control = make_control(tol);
  Solver<TwoRateSystem<double>> solver(sys, control);
  std::vector<double> t(times.begin(), times.end());
  solver.advance_adaptive(t);

  auto hist = solver.get_history();
  const int nt = (int)hist.size(), nd = (int)sys.ode_size();
  Rcpp::NumericMatrix out(nt, nd);
  for (int i = 0; i < nt; ++i) {
    std::vector<double> s(nd);
    hist[i].ode_state(s.begin());
    for (int j = 0; j < nd; ++j) out(i, j) = s[j];
  }
  const long nsteps = (long)solver.times().size() - 1;
  return Rcpp::List::create(_["states"] = out, _["nsteps"] = (double)nsteps);
}

// Macro-stepped drainage run, with the inner sub-cycle either adaptive (unsplit)
// or operator-split (exact drainage flow + ROS34PW2 residual). Same model both
// ways, so the sub-step counts are directly comparable.
// [[Rcpp::export]]
Rcpp::List drainage_mri(double c, int n_fast, int n_slow, bool split,
                        std::string table, Rcpp::NumericVector macro_times, double tol) {
  DrainageSystem<double> sys(c, n_fast, n_slow);
  MRICoupling M = pick_table(table);
  OdeControl control = make_control(tol);
  std::vector<double> times(macro_times.begin(), macro_times.end());
  MRISchedule sched;
  auto hist = split
    ? mri_advance(sys, M, control, times, sched, false, SplitSubcycle{})
    : mri_advance(sys, M, control, times, sched, false, AdaptiveSubcycle{});
  return Rcpp::List::create(_["states"] = to_matrix(hist),
                            _["n_fast"] = (double)mri_fast_steps(sched));
}

// [[Rcpp::export]]
Rcpp::List drainage_reference(double c, int n_fast, int n_slow,
                              Rcpp::NumericVector times, double tol) {
  DrainageSystem<double> sys(c, n_fast, n_slow);
  OdeControl control = make_control(tol);
  Solver<DrainageSystem<double>> solver(sys, control);
  std::vector<double> t(times.begin(), times.end());
  solver.advance_adaptive(t);
  auto hist = solver.get_history();
  const int nt = (int)hist.size(), nd = (int)sys.ode_size();
  Rcpp::NumericMatrix out(nt, nd);
  for (int i = 0; i < nt; ++i) {
    std::vector<double> s(nd);
    hist[i].ode_state(s.begin());
    for (int j = 0; j < nd; ++j) out(i, j) = s[j];
  }
  return Rcpp::List::create(_["states"] = out, _["nsteps"] = (double)(solver.times().size() - 1));
}

// J = sum of the final slow states -- a smooth scalar of the whole trajectory.
template <class T>
static T functional(const std::vector<T>& final_state, size_t n_slow_) {
  T s(0.0);
  for (size_t i = 0; i < n_slow_; ++i) s += final_state[i];   // slow block is [0, n_slow)
  return s;
}

// Gradient of J w.r.t. [k, initial state] through the MRI stepper by
// record->replay, checked against a frozen-schedule central finite difference.
// [[Rcpp::export]]
Rcpp::List two_rate_gradient(double k, int n_slow, std::string table,
                             Rcpp::NumericVector macro_times, double tol, double eps_fd) {
  using ad = xad::adj<double>;
  using ad_type = ad::active_type;
  MRICoupling M = pick_table(table);
  OdeControl control = make_control(tol);
  std::vector<double> times(macro_times.begin(), macro_times.end());

  // Record the sub-cycle schedule with a double adaptive pass; capture the
  // baseline initial state the perturbations start from.
  TwoRateSystem<double> sys_d(k, n_slow);
  const size_t slow = sys_d.slow_size(), nd = sys_d.ode_size();
  std::vector<double> x0(nd);
  sys_d.ode_state(x0.begin());
  MRISchedule sched;
  mri_advance(sys_d, M, control, times, sched, /*replay=*/false);

  // Active replay under the tape: the tape is the discrete adjoint of the run.
  ad::tape_type tape(false);
  tape.activate();
  TwoRateSystem<ad_type> sys_a(k, n_slow);
  std::vector<double> kv{k};
  std::vector<ad_type*> inputs = sys_a.set_params(tape, kv.begin());
  auto ic = sys_a.set_initial_state(tape, x0.begin(), times.front());
  inputs.insert(inputs.end(), ic.begin(), ic.end());
  tape.newRecording();
  sys_a.reset();
  auto hist_a = mri_advance(sys_a, M, control, times, sched, /*replay=*/true);
  ad_type J = functional(hist_a.back(), slow);
  tape.registerOutput(J);
  xad::derivative(J) = 1.0;
  tape.computeAdjoints();
  std::vector<double> grad_adj(inputs.size());
  for (size_t i = 0; i < inputs.size(); ++i) grad_adj[i] = xad::derivative(*inputs[i]);
  const double value = xad::value(J);
  tape.deactivate();

  // Frozen-schedule central finite difference over the same inputs [k, state].
  auto replay_J = [&](TwoRateSystem<double>& s) {
    auto h = mri_advance(s, M, control, times, sched, /*replay=*/true);
    return functional(h.back(), slow);
  };
  std::vector<double> grad_fd(inputs.size());
  {
    TwoRateSystem<double> sp(k + eps_fd, n_slow), sm(k - eps_fd, n_slow);
    grad_fd[0] = (replay_J(sp) - replay_J(sm)) / (2 * eps_fd);
  }
  for (size_t i = 0; i < nd; ++i) {
    std::vector<double> xp = x0, xm = x0;
    xp[i] += eps_fd; xm[i] -= eps_fd;
    TwoRateSystem<double> sp(k, n_slow), sm(k, n_slow);
    sp.set_initial_state(xp.begin(), times.front()); sp.reset();
    sm.set_initial_state(xm.begin(), times.front()); sm.reset();
    grad_fd[1 + i] = (replay_J(sp) - replay_J(sm)) / (2 * eps_fd);
  }

  return Rcpp::List::create(_["value"] = value,
                            _["grad_adjoint"] = Rcpp::wrap(grad_adj),
                            _["grad_fd"] = Rcpp::wrap(grad_fd));
}

// Gradient of J w.r.t. [c, initial state] through the operator-split inner by
// record->replay: a double split pass records the sub-cycle schedule, an active
// split pass replays it under the tape (the ROS Jacobian is passive, the exact
// drainage flow and the stages tape), checked against a frozen-schedule central
// difference. Demonstrates reverse mode through the exact flow + ROS34PW2 split.
// [[Rcpp::export]]
Rcpp::List drainage_gradient_split(double c, int n_fast, int n_slow, std::string table,
                                   Rcpp::NumericVector macro_times, double tol, double eps_fd) {
  using ad = xad::adj<double>;
  using ad_type = ad::active_type;
  MRICoupling M = pick_table(table);
  OdeControl control = make_control(tol);
  std::vector<double> times(macro_times.begin(), macro_times.end());

  DrainageSystem<double> sys_d(c, n_fast, n_slow);
  const size_t slow = sys_d.slow_size(), nd = sys_d.ode_size();
  std::vector<double> x0(nd);
  sys_d.ode_state(x0.begin());
  MRISchedule sched;
  mri_advance(sys_d, M, control, times, sched, /*replay=*/false, SplitSubcycle{});

  ad::tape_type tape(false);
  tape.activate();
  DrainageSystem<ad_type> sys_a(c, n_fast, n_slow);
  std::vector<double> cv{c};
  std::vector<ad_type*> inputs = sys_a.set_params(tape, cv.begin());
  auto ic = sys_a.set_initial_state(tape, x0.begin(), times.front());
  inputs.insert(inputs.end(), ic.begin(), ic.end());
  tape.newRecording();
  sys_a.reset();
  auto hist = mri_advance(sys_a, M, control, times, sched, /*replay=*/true, SplitSubcycle{});
  ad_type J = functional(hist.back(), slow);
  tape.registerOutput(J);
  xad::derivative(J) = 1.0;
  tape.computeAdjoints();
  std::vector<double> grad_adj(inputs.size());
  for (size_t i = 0; i < inputs.size(); ++i) grad_adj[i] = xad::derivative(*inputs[i]);
  const double value = xad::value(J);
  tape.deactivate();

  auto replay_J = [&](DrainageSystem<double>& s) {
    auto h = mri_advance(s, M, control, times, sched, /*replay=*/true, SplitSubcycle{});
    return functional(h.back(), slow);
  };
  std::vector<double> grad_fd(inputs.size());
  {
    DrainageSystem<double> sp(c + eps_fd, n_fast, n_slow), sm(c - eps_fd, n_fast, n_slow);
    grad_fd[0] = (replay_J(sp) - replay_J(sm)) / (2 * eps_fd);
  }
  for (size_t i = 0; i < nd; ++i) {
    std::vector<double> xp = x0, xm = x0;
    xp[i] += eps_fd; xm[i] -= eps_fd;
    DrainageSystem<double> sp(c, n_fast, n_slow), sm(c, n_fast, n_slow);
    sp.set_initial_state(xp.begin(), times.front()); sp.reset();
    sm.set_initial_state(xm.begin(), times.front()); sm.reset();
    grad_fd[1 + i] = (replay_J(sp) - replay_J(sm)) / (2 * eps_fd);
  }

  return Rcpp::List::create(_["value"] = value,
                            _["grad_adjoint"] = Rcpp::wrap(grad_adj),
                            _["grad_fd"] = Rcpp::wrap(grad_fd));
}
