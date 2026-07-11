/* CanopySystem interface -- the record -> replay demonstrator, exercising the
 * Replayable hooks and a fixed-node interpolator on the AD path.
 *
 * Names are Canopy_-prefixed: this TU links into the same odelia.so as the Lorenz
 * interface, so the exported symbols must not collide.
 */

#include <Rcpp.h>
#include <XAD/XAD.hpp>
#include <odelia/ode_solver.hpp>
#include <odelia/gradient.hpp>
#include <odelia/solver_interface.hpp>
#include <examples/canopy_system.hpp>

using namespace Rcpp;
using namespace odelia;

typedef CanopySystem<double> SystemType;
typedef CanopySystem<xad::adj<double>::active_type> ActiveSystemType;

// [[Rcpp::export]]
SEXP Canopy_new(double gain, double y0 = 1.0) {
  return Rcpp::XPtr<SystemType>(new SystemType(gain, y0), true);
}

// A fully adaptive (double) run to Tmax; the finite-difference reference the replay
// gradient is checked against. Returns the final state y(Tmax).
// [[Rcpp::export]]
double Canopy_adaptive_final(SEXP system_xp, SEXP control_xp, double Tmax) {
  Rcpp::XPtr<SystemType> sys(system_xp);
  Rcpp::XPtr<ode::OdeControl> ctrl(control_xp);
  ode::Solver<SystemType> solver(*sys, *ctrl);
  solver.advance_adaptive({0.0, Tmax});
  return solver.state()[0];
}

// The functional to differentiate: the final canopy state. Like any functional it is
// a pure reduction -- the driver hands it a positioned solver and it returns a value.
struct canopy_final {
  std::size_t codomain() const { return 1; }
  template<typename Solver>
  typename Solver::value_type operator()(Solver& solver) const {
    return solver.state()[0];
  }
};

// Record one adaptive double pass, then replay it on the active System and
// differentiate the final state w.r.t. `gain`. `reuse_light` picks the replay:
// FALSE recomputes the light with the active scalar on the recorded node positions
// (the canopy re-shades), TRUE reuses the recorded light values as constants.
// Returns value, gradient, and the recorded step count.
// [[Rcpp::export]]
Rcpp::List Canopy_record_replay_gradient(SEXP system_xp, SEXP control_xp,
                                         double Tmax, bool reuse_light = false) {
  Rcpp::XPtr<SystemType> sys(system_xp);
  Rcpp::XPtr<ode::OdeControl> ctrl(control_xp);

  // Record: one adaptive double pass on its own double Solver.
  SystemType dbl(*sys);
  dbl.start_recording();
  ode::Solver<SystemType> dbl_solver(dbl, *ctrl);
  dbl_solver.advance_adaptive({0.0, Tmax});

  const std::vector<double> times = dbl_solver.times();
  auto& rec_sys = dbl_solver.get_system_ref();
  const double gain = rec_sys.pars();

  // Replay: lift to the active System, hand over the recording, differentiate.
  ActiveSystemType act = rec_sys.rebind_from<xad::adj<double>::active_type>();
  act.set_recording(rec_sys.recorded_positions(), rec_sys.recorded_values(), reuse_light);
  ode::Solver<ActiveSystemType> active_solver(act, *ctrl);
  active_solver.set_schedule(times);   // the recorded L1 schedule to replay

  ode::DifferentiationTargets ind;
  ind.params.push_back(0);     // gain
  ind.values.push_back(gain);

  auto [value, gradient] =
      ode::compute_gradient(active_solver, ind, canopy_final{});

  return Rcpp::List::create(
      Rcpp::Named("value") = value,
      Rcpp::Named("gradient") = Rcpp::wrap(gradient),
      Rcpp::Named("n_steps") = static_cast<int>(times.size() - 1));
}

//-------------------------------------------------------------------------
// A persistent double Solver held across calls. The one-shot above rebuilds the
// active solver every call; here the active solver and its tape are cached on the
// double Solver and reused, while the recording is read fresh on each gradient call.

// [[Rcpp::export]]
SEXP Canopy_Solver_new(SEXP system_xp, SEXP control_xp) {
  Rcpp::XPtr<SystemType> sys(system_xp);
  Rcpp::XPtr<ode::OdeControl> ctrl(control_xp);
  return Rcpp::XPtr<ode::Solver<SystemType>>(
      new ode::Solver<SystemType>(*sys, *ctrl), true);
}

// One adaptive double pass, recording node positions and light into the solver's
// System. Returns the resolved step count.
// [[Rcpp::export]]
int Canopy_record(SEXP solver_xp, double Tmax) {
  Rcpp::XPtr<ode::Solver<SystemType>> d(solver_xp);
  d->reset();
  d->get_system_ref().start_recording();
  d->advance_adaptive({0.0, Tmax});
  return static_cast<int>(d->recorded_steps().size() - 1);
}

// Differentiate the replayed final state w.r.t. `gain`, reusing the cached active
// solver. The recording is read from the double System on this call, so a fresh
// Canopy_record() is always picked up rather than replayed against a stale schedule.
// [[Rcpp::export]]
Rcpp::List Canopy_replay_gradient(SEXP solver_xp, bool reuse_light = false) {
  Rcpp::XPtr<ode::Solver<SystemType>> d(solver_xp);
  if (!d->has_recording()) {
    util::stop("Canopy_replay_gradient: call Canopy_record() first");
  }
  auto& rec = d->get_system_ref();
  const std::vector<double> times = d->recorded_steps();

  auto& active = solver::active_solver<SystemType, ActiveSystemType>(*d);
  active.set_schedule(times);   // the recorded L1 schedule to replay
  active.get_system_ref().set_recording(rec.recorded_positions(), rec.recorded_values(), reuse_light);

  ode::DifferentiationTargets ind;
  ind.params.push_back(0);     // gain
  ind.values.push_back(rec.pars());

  auto [value, gradient] = ode::compute_gradient(active, ind, canopy_final{});
  return Rcpp::List::create(
      Rcpp::Named("value") = value,
      Rcpp::Named("gradient") = Rcpp::wrap(gradient),
      Rcpp::Named("n_steps") = static_cast<int>(times.size() - 1));
}
