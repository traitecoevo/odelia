#ifndef CANOPY_SYSTEM_HPP_
#define CANOPY_SYSTEM_HPP_

#include <odelia/ode_solver.hpp>
#include <odelia/interpolator.hpp>
#include <XAD/XAD.hpp>
#include <vector>
#include <cmath>

using namespace odelia;

// A one-state canopy that relaxes toward the light it captures -- the demonstrator
// for record -> replay, which Lorenz and leaf_thermal don't exercise.
//
//   dy/dt = -turnover * y + L(y)
//
// y is the canopy state (say, leaf area). L is the light captured, read from a
// vertical light profile over relative depth x in [0, 1]:
//
//   light_profile(x) = gain * exp(-extinction * x)         (incident light attenuated
//                                                            with depth, Beer-Lambert)
//                    + y * exp(-shade_conc * (x - 0.5)^2)  (a leaf layer at mid-canopy
//                                                            that grows with the state)
//
// The profile is built by adaptively refining a spline over depth, then read at
// ref_depth. Because the leaf-layer term moves with y, the refinement places its
// nodes differently every step -- awkward to differentiate through directly, cheap to
// record once and replay. The ODE is linear in (gain, y), so dy(T)/dgain has a clean
// finite-difference reference.
//
// On the adaptive pass the System records, per step, the node positions the
// refinement chose (and, if asked, the light at each RK stage). A later pass replays
// on those fixed positions: either recomputing the light with the active scalar (its
// derivative flows -- the canopy re-shades in response to the parameter), or reusing
// the recorded light values as fixed doubles (the derivative through the light is
// then zero).
template <typename T = double>
class CanopySystem {
public:
  using value_type = T;

  CanopySystem(T gain_, double y0_ = 1.0,
               double turnover_ = 1.0, double extinction_ = 1.0, double shade_conc_ = 40.0,
               double ref_depth_ = 0.5, double tol_ = 1e-4, int max_depth_ = 12)
    : gain(gain_), y0_init(y0_), turnover(turnover_), extinction(extinction_),
      shade_conc(shade_conc_), ref_depth(ref_depth_), tol(tol_), max_depth(max_depth_),
      t0(0.0) {
    reset();
  }

  // Copy the configuration onto another scalar so the driver can build the active
  // version; the recording is handed in per call (set_recording), not carried here.
  template <class S2> using rebind = CanopySystem<S2>;

  template <class S2>
  rebind<S2> rebind_from() const {
    return rebind<S2>(xad::value(gain), y0_init, turnover, extinction, shade_conc,
                      ref_depth, tol, max_depth);
  }

  // ---- ODE interface -------------------------------------------------------
  size_t ode_size() const { return 1; }
  double ode_time() const { return time; }
  double ode_t0() const { return t0; }

  // Adaptive pass, and replay that recomputes the light: set the state and (re)build
  // the light profile.
  template <typename Iterator>
  Iterator set_ode_state(Iterator it, double time_) {
    y = *it++;
    time = time_;
    light = compute_light();
    compute_rates();
    return it;
  }

  // Replay that reuses recorded light: read this RK stage's recorded value as a plain
  // double (off the tape, so its derivative is zero) instead of rebuilding the
  // profile. derivs routes here when has_recorded_field() is true.
  template <typename Iterator>
  Iterator set_ode_state(Iterator it, int stage) {
    y = *it++;
    bg_light = stage_light.at(stage);
    compute_rates();
    return it;
  }

  // Turnover minus captured light. On a reuse replay the light is the recorded double;
  // otherwise it is the value read from the (active) profile.
  void compute_rates() {
    dydt = -turnover * y;
    if (has_recorded_field()) dydt += bg_light;
    else                      dydt += light;
  }

  template <typename Iterator>
  Iterator set_initial_state(Iterator it, double t0_ = 0.0) {
    t0 = t0_;
    y0_init = xad::value(*it++);
    return it;
  }

  template <typename Iterator>
  Iterator ode_state(Iterator it) const { *it++ = y; return it; }

  template <typename Iterator>
  Iterator ode_rates(Iterator it) const { *it++ = dydt; return it; }

  void reset() {
    y = T(y0_init);
    time = t0;
    step = 0;
    if (has_recorded_field()) {
      // light_history[step][stage]: take this first step's row of per-RK-stage values,
      stage_light = light_history.front();
      bg_light = stage_light.front();     // and the value the first stage will read.
    } else {
      if (replaying()) node_positions = positions_history.front();
      light = compute_light();
    }
    compute_rates();
  }

  // The single differentiable input; this canopy has no seedable initial state.
  std::vector<T*> ad_parameters()    { return {&gain}; }
  std::vector<T*> ad_initial_state() { return {}; }

  // ---- Record / replay hooks (Replayable) ---------------------------------
  // record_* run only while recording; replay_step only on replay. Record the node
  // positions once per accepted step, and the light at each of the six RK stages.
  void record_stage(int stage) {
    if (!recording) return;
    if (stage == 0) {
      pending_positions = node_positions;   // hold this step's positions until it commits
      stage_light.assign(6, 0.0);
    }
    if (stage >= 0 && stage < 6) {
      stage_light[stage] = xad::value(light);
    }
  }

  void record_ode_step() {
    if (!recording) return;
    positions_history.push_back(pending_positions);
    light_history.push_back(stage_light);
  }

  // Load this step's recorded positions (and light, if recorded), then advance the
  // cursor. A no-op on the recording pass. Replay runs on the recorded schedule, so
  // the cursor walks [0, size) without clamping.
  void replay_step() {
    if (!replaying()) return;
    node_positions = positions_history.at(step);
    if (has_recorded_field()) stage_light = light_history.at(step);
    ++step;
  }

  // ---- Record / replay channel (System -> System) -------------------------
  // The double Solver produces a recording; the active replay reads it per call.
  // Positions and light values are plain doubles, so they cross the double->active
  // boundary directly.
  void start_recording() {
    recording = true;
    positions_history.clear();
    light_history.clear();
  }
  const std::vector<std::vector<double>>& recorded_positions() const { return positions_history; }
  const std::vector<std::vector<double>>& recorded_values()    const { return light_history; }

  // Hand a recording to the active System for one replay. `reuse_light` = reuse the
  // recorded light values as constants (derivative through the light is zero);
  // otherwise the light is recomputed with the active scalar on the recorded
  // positions. has_recorded_field() then reports which path applies.
  void set_recording(std::vector<std::vector<double>> positions,
                     std::vector<std::vector<double>> values, bool reuse_light) {
    positions_history = std::move(positions);
    light_history = reuse_light ? std::move(values) : std::vector<std::vector<double>>{};
    recording = false;
    step = 0;
  }
  bool has_recording() const { return !positions_history.empty(); }

  double pars() const { return xad::value(gain); }

  // Whether recorded light values are present to reuse: the query derivs reads to
  // choose the replay path. Guarded against the recording pass, which is still
  // building the values.
  bool has_recorded_field() const { return !recording && !light_history.empty(); }

private:
  bool replaying() const { return !recording && has_recording(); }

  // Build the light profile as a spline over depth and read it at ref_depth. On the
  // adaptive pass the nodes are refined and stored; on replay they are rebuilt with
  // the active scalar on the recorded positions.
  T compute_light() {
    interpolator::basic_interpolator<T> interp;
    if (replaying()) {
      std::vector<T> vals;
      vals.reserve(node_positions.size());
      for (double x : node_positions) vals.push_back(light_profile(x));
      interp.init(node_positions, vals);
    } else {
      interp.construct([this](double x) { return light_profile(x); }, 0.0, 1.0,
                       tol, 0.0, 5, static_cast<std::size_t>(max_depth));
      if (recording) node_positions = interp.get_x();
    }
    return interp.eval(ref_depth);
  }

  T light_profile(double x) const {
    using std::exp;
    const double g = x - 0.5;
    return gain * exp(-extinction * x) + y * exp(-shade_conc * g * g);
  }

  T gain;
  double y0_init;
  double turnover, extinction, shade_conc, ref_depth, tol;
  int max_depth;
  double t0;

  T y, dydt, light;
  double bg_light = 0.0;   // on a reuse replay, the recorded light for the current RK stage
  double time;

  bool recording = false;
  std::vector<std::vector<double>> positions_history;  // node positions per ODE step
  std::vector<std::vector<double>> light_history;      // light per [ODE step][RK stage]

  std::size_t step = 0;
  std::vector<double> node_positions;   // positions built (recording) or loaded (replay)
  std::vector<double> stage_light;      // this step's light at each of the six RK stages
  std::vector<double> pending_positions;
};

#endif
