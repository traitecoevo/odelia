// -*-c++-*-
#ifndef ODELIA_ODE_SUBSTEPS_HPP_
#define ODELIA_ODE_SUBSTEPS_HPP_

// A step that integrates a System's subsystem (HasSubsystem) in substeps of its
// own, for an explicit Runge-Kutta stepper with a stage at t + h: the stages and
// the end of such a step, the choice of its substeps, and what it leaves behind
// for the solver's decision on the next step. The two explicit steppers derive
// from this, each handing over its tableau:
//
//   static constexpr int stage_at_end;       the stage i, in 1..5, at t + h
//   static constexpr int stage_rates;        rate vectors a step carries, k[0..)
//   static constexpr bool error_reads_end_rates;
//                                            whether step_error reads the rates
//                                            at the end state, in the last slot
//   static constexpr double stability_limit; -h times a real rate the solution
//                                            is stable down to
//   static constexpr double stage_fraction(int i);
//                                            stage i's fraction of the step
//   static double end_weight(int i);         stage i's weight in step_end
//   stage_state(i, y, k, h, out), step_end(y, k, h, out), step_error(k, h, out),
//   stage_time(i, time, h), order()
//
// The subsystem needs its inputs at every moment of the step, and the rest of
// the System computes them only at the stages. So the step takes the inputs
// along a line, and integrates the subsystem twice: a predictor under the inputs
// extrapolated along their slope over the step before, whose values at the stage
// fractions the stages read, and a corrector under the line from the inputs at
// the start to those the stage at t + h computed, whose end is the step's. The
// error on the subsystem is the larger of the gap between the two and the
// stages' inputs' defect against the corrector's line, so it shrinks with the
// square of the step where the stepper's own shrinks with the fifth power: a
// plain step is cheaper where the subsystem is not stiff.

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <utility>
#include <vector>
#include <odelia/ode_control.hpp>
#include <odelia/ode_interface.hpp>
#include <odelia/ode_util.hpp>

namespace odelia {
namespace ode {

template <class Stepper, class System>
class SubsystemSubsteps {
public:
  using value_type = typename System::value_type;
  using state_type = std::vector<value_type>;

  // How stiff the System's subsystem is at the end of the last step: its rates'
  // change over its state's change, from the stage at t + h to the end. Two
  // equal states measure nothing, and the last measure stands.
  double subsystem_stiffness = 0.0;
  // The last step's own error estimate on the subsystem, whether or not it
  // substepped, and the subsystem's inputs at its stage at t + h.
  state_type subsystem_own_error, inputs_at_end_stage;

  // The substeps of a step's predictor: the stepper on the subsystem from `y`
  // under inputs extrapolated from `u0` along `input_slope`, controlled at
  // substep_tol and no smaller than `step_size_min`, landing on every stage
  // fraction and every knot of the subsystem's forcing.
  subsystem_substeps choose_substeps(const System& system, double time,
                                     double step_size, const state_type& y,
                                     const state_type& dydt,
                                     const state_type& u0,
                                     std::vector<double> input_slope,
                                     double step_size_min) const
    requires HasSubsystem<System>;

protected:
  // Stage i of a step being evaluated at `at`: where it is the stage at t + h,
  // the subsystem's part of it and the inputs there are kept for the step's end.
  void note_stage(int i, const System& system, const state_type& at)
    requires HasSubsystem<System>;
  // After the evaluation at the step's end: the stiffness, the step's own error
  // on the subsystem's rows of `yerr`, and those rows replaced by the coupling
  // error `subsystem_error` where the step substepped (empty where it did not).
  void record_subsystem(const System& system, const state_type& y,
                        const std::vector<state_type>& k,
                        const state_type& dydt_out, state_type& yerr,
                        const std::vector<double>& subsystem_error)
    requires HasSubsystem<System>;

  // The stages and the end for a step that substeps the subsystem from its
  // inputs `u0`: the rest from the tableau, the subsystem from the predictor and
  // the corrector. `rates_at(i, stage, out)` evaluates stage i. Returns, at
  // double, the subsystem's error from the coupling, and nothing at any other
  // scalar.
  template <class Sys, class S, class RatesAt>
  std::vector<double> take_substepped_step(const Sys& sys, double time, double h,
                                           const std::vector<S>& y0,
                                           std::vector<std::vector<S>>& k,
                                           std::vector<S>& y_end,
                                           const subsystem_substeps& substeps,
                                           const std::vector<S>& u0,
                                           RatesAt&& rates_at) const;

private:
  const Stepper& self() const { return static_cast<const Stepper&>(*this); }

  // The subsystem from y0's, over the substeps ending at `ends`, its first rates
  // the step's own, under `rates(u, s, out)`; `at_end(u, s)` sees it where each
  // substep ends.
  template <class S, class Rates, class AtEnd>
  void integrate_substeps(std::size_t first, const std::vector<double>& ends,
                          double h, const std::vector<S>& y0,
                          const std::vector<std::vector<S>>& k, Rates&& rates,
                          std::vector<S>& own, AtEnd&& at_end) const;
  // One substep from fraction `from` to `to` of a step of h, from `own` whose
  // rates r[0] holds, into `out`; r[1..5] take its stage rates.
  template <class S, class Rates>
  void take_substep(const std::vector<S>& own, std::vector<std::vector<S>>& r,
                    double from, double to, double h, Rates&& rates,
                    std::vector<S>& out) const;
  // The inputs extrapolated over a step of h from u0 along `input_slope`, held
  // at zero rather than changing sign.
  template <class S>
  static std::vector<S> extrapolated_inputs(const std::vector<S>& u0,
                                            const std::vector<double>& input_slope,
                                            double h);
  // The subsystem's rates at fraction u of a step of h from `time`, under inputs
  // on the line from `from` to `to`. `from` outlives the call; `to` is kept.
  template <class Sys, class S>
  static auto rates_on_line(const Sys& sys, double time, double h,
                            const std::vector<S>& from, std::vector<S> to) {
    return [&sys, time, h, &from, to = std::move(to),
            line = std::vector<S>(from.size())](
             double u, const std::vector<S>& s, std::vector<S>& out) mutable {
      for (std::size_t q = 0; q < line.size(); ++q) {
        line[q] = from[q] + (to[q] - from[q]) * u;
      }
      sys.subsystem_rates(time + u * h, s, line, out);
    };
  }
  // The stages' fractions of the step, ascending: where a predictor's substeps
  // must land.
  static std::array<double, 5> stage_fractions() {
    std::array<double, 5> out{};
    for (int i = 1; i <= 5; ++i) {
      out[static_cast<std::size_t>(i - 1)] = Stepper::stage_fraction(i);
    }
    std::sort(out.begin(), out.end());
    return out;
  }

  // The substeps' relative tolerance; the absolute one is 1e-4 of it.
  static constexpr double substep_tol = 1e-9;
  // The subsystem's part of the stage at t + h.
  state_type subsystem_at_end_stage;
};

template <class Stepper, class System>
void SubsystemSubsteps<Stepper, System>::note_stage(int i, const System& system,
                                                    const state_type& at)
  requires HasSubsystem<System> {
  if (i == Stepper::stage_at_end) {
    const auto [first, n] = system.subsystem();
    const auto from = at.begin() + static_cast<std::ptrdiff_t>(first);
    subsystem_at_end_stage.assign(from, from + static_cast<std::ptrdiff_t>(n));
    system.subsystem_inputs(inputs_at_end_stage);
  }
}

template <class Stepper, class System>
void SubsystemSubsteps<Stepper, System>::record_subsystem(
    const System& system, const state_type& y, const std::vector<state_type>& k,
    const state_type& dydt_out, state_type& yerr,
    const std::vector<double>& subsystem_error)
  requires HasSubsystem<System> {
  const auto [first, n] = system.subsystem();
  double rates = 0.0, states = 0.0;
  for (std::size_t q = 0; q < n; ++q) {
    const double rate = util::to_passive(dydt_out[first + q]) -
      util::to_passive(k[Stepper::stage_at_end][first + q]);
    const double state = util::to_passive(y[first + q]) -
      util::to_passive(subsystem_at_end_stage[q]);
    rates += rate * rate;
    states += state * state;
  }
  if (states > 0.0) {
    subsystem_stiffness = std::sqrt(rates / states);
  }
  const auto from = yerr.begin() + static_cast<std::ptrdiff_t>(first);
  subsystem_own_error.assign(from, from + static_cast<std::ptrdiff_t>(n));
  if (!subsystem_error.empty()) {
    std::copy(subsystem_error.begin(), subsystem_error.end(), from);
  }
}

template <class Stepper, class System>
subsystem_substeps SubsystemSubsteps<Stepper, System>::choose_substeps(
    const System& system, double time, double step_size, const state_type& y,
    const state_type& dydt, const state_type& u0, std::vector<double> input_slope,
    double step_size_min) const
  requires HasSubsystem<System> {
  const double h = step_size;
  const auto [first, n] = system.subsystem();
  util::check_length(u0.size(), n);
  const auto at = static_cast<std::ptrdiff_t>(first);
  const auto last = at + static_cast<std::ptrdiff_t>(n);
  auto rates = rates_on_line(system, time, h, u0,
                             extrapolated_inputs(u0, input_slope, h));
  std::vector<double> s(y.begin() + at, y.begin() + last), next(n), error(n);
  std::vector<std::vector<double>> r(Stepper::stage_rates, std::vector<double>(n));
  std::copy(dydt.begin() + at, dydt.begin() + last, r[0].begin());
  // GSL's rule, in fractions of this step and down to the solver's smallest step.
  OdeControl control(1e-4 * substep_tol, substep_tol, 1.0, 0.0, step_size_min / h,
                     1.0, 1.0);
  control.set_controller(Controller::gsl);
  const std::array<double, 5> fractions = stage_fractions();
  std::vector<double> stops(fractions.begin(), fractions.end());
  for (const double knot : system.subsystem_knots(time, time + h)) {
    const double fraction = (knot - time) / h;
    if (fraction > 0.0 && fraction < 1.0) {
      stops.push_back(fraction);
    }
  }
  std::sort(stops.begin(), stops.end());
  stops.erase(std::unique(stops.begin(), stops.end()), stops.end());
  // The rates at a substep's end, where the error reads them, go in the last
  // slot and carry over as the next substep's first.
  constexpr std::size_t end_slot = Stepper::stage_rates - 1;
  subsystem_substeps chosen{std::move(input_slope), {}};
  double from = 0.0, substep = 1.0;
  for (const double stop : stops) {
    while (from < stop) {
      const double end = std::min(from + substep, stop);
      take_substep(s, r, from, end, h, rates, next);
      if constexpr (Stepper::error_reads_end_rates) {
        rates(end, next, r[end_slot]);
      }
      self().step_error(r, (end - from) * h, error);
      // Compared with the substep asked for: `end - from` rounds away from it.
      const double asked = substep;
      substep = control.adjust_step_size(n, self().order(), end - from, next,
                                         error, r[0]);
      if (control.step_size_shrank()) {
        if (!(substep < asked)) {
          throw util::DomainError(
            "The subsystem fails its error test at the smallest substep");
        }
        continue;
      }
      chosen.ends.push_back(end);
      from = end;
      std::swap(s, next);
      if (from < 1.0) {
        if constexpr (Stepper::error_reads_end_rates) {
          std::swap(r[0], r[end_slot]);
        } else {
          rates(from, s, r[0]);
        }
      }
    }
  }
  return chosen;
}

template <class Stepper, class System>
template <class S>
std::vector<S> SubsystemSubsteps<Stepper, System>::extrapolated_inputs(
    const std::vector<S>& u0, const std::vector<double>& input_slope, double h) {
  util::check_length(input_slope.size(), u0.size());
  std::vector<S> out(u0.size());
  for (std::size_t q = 0; q < u0.size(); ++q) {
    out[q] = u0[q] + h * input_slope[q];
    if ((util::to_passive(out[q]) < 0.0) != (util::to_passive(u0[q]) < 0.0)) {
      out[q] = S(0.0);
    }
  }
  return out;
}

template <class Stepper, class System>
template <class Sys, class S, class RatesAt>
std::vector<double> SubsystemSubsteps<Stepper, System>::take_substepped_step(
    const Sys& sys, double time, double h, const std::vector<S>& y0,
    std::vector<std::vector<S>>& k, std::vector<S>& y_end,
    const subsystem_substeps& substeps, const std::vector<S>& u0,
    RatesAt&& rates_at) const {
  const state_range subsystem = sys.subsystem();
  const auto first = static_cast<std::ptrdiff_t>(subsystem.first);
  util::check_length(u0.size(), subsystem.size);
  if (substeps.ends.empty() || substeps.ends.back() != 1.0) {
    util::stop("Substeps must end at the step's end");
  }
  const std::array<double, 5> fractions = stage_fractions();

  // Predictor: the subsystem over the substeps, its inputs extrapolated from u0
  // along the recorded slope, kept at each stage fraction.
  std::array<std::vector<S>, 5> predicted;
  std::size_t kept = 0;
  std::vector<S> s(subsystem.size);
  integrate_substeps(subsystem.first, substeps.ends, h, y0, k,
                     rates_on_line(sys, time, h, u0,
                                   extrapolated_inputs(u0, substeps.input_slope, h)),
                     s, [&](double u, const std::vector<S>& at) {
                       if (kept < predicted.size() && u == fractions[kept]) {
                         predicted[kept++] = at;
                       }
                     });
  if (kept != predicted.size()) {
    util::stop("Substeps must land on every stage fraction");
  }
  auto predicted_at = [&](double u) -> const std::vector<S>& {
    const auto at = std::find(fractions.begin(), fractions.end(), u);
    return predicted.at(static_cast<std::size_t>(at - fractions.begin()));
  };

  // Stages: the rest from the tableau, the subsystem from the predictor. Each
  // evaluation leaves that stage's inputs on the System.
  std::vector<S> stage(y0.size());
  std::array<std::vector<S>, 5> stage_inputs;
  for (int i = 1; i < 6; ++i) {
    self().stage_state(i, y0, k, h, stage);
    const std::vector<S>& p = predicted_at(Stepper::stage_fraction(i));
    std::copy(p.begin(), p.end(), stage.begin() + first);
    rates_at(i, stage, k[i]);
    sys.subsystem_inputs(stage_inputs[static_cast<std::size_t>(i - 1)]);
  }

  // Corrector: the subsystem again over the same substeps, its inputs on the
  // line from u0 to those at the stage at t + h. Its end is the step's.
  const std::vector<S>& inputs_at_end =
    stage_inputs[static_cast<std::size_t>(Stepper::stage_at_end - 1)];
  std::vector<S> end(subsystem.size);
  integrate_substeps(subsystem.first, substeps.ends, h, y0, k,
                     rates_on_line(sys, time, h, u0, inputs_at_end), end,
                     [](double, const std::vector<S>&) {});
  self().step_end(y0, k, h, y_end);
  std::copy(end.begin(), end.end(), y_end.begin() + first);

  // Error, at double: how far the corrector ends from the predictor, or how far
  // the stages' inputs stray from the corrector's line, weighed as step_end
  // weighs the stages, whichever is larger.
  std::vector<double> error;
  if constexpr (std::same_as<S, double>) {
    const std::vector<S>& predicted_end = predicted_at(1.0);
    auto corrector_line = rates_on_line(sys, time, h, u0, inputs_at_end);
    std::vector<double> stage_rates(subsystem.size), line_rates(subsystem.size);
    error.assign(subsystem.size, 0.0);
    for (int i = 1; i < 6; ++i) {
      const double weight = Stepper::end_weight(i);
      if (weight == 0.0) {
        continue;
      }
      const double fraction = Stepper::stage_fraction(i);
      const std::vector<S>& p = predicted_at(fraction);
      sys.subsystem_rates(self().stage_time(i, time, h), p,
                          stage_inputs[static_cast<std::size_t>(i - 1)],
                          stage_rates);
      corrector_line(fraction, p, line_rates);
      for (std::size_t q = 0; q < subsystem.size; ++q) {
        error[q] += h * weight * std::abs(stage_rates[q] - line_rates[q]);
      }
    }
    for (std::size_t q = 0; q < subsystem.size; ++q) {
      error[q] = std::max(std::abs(end[q] - predicted_end[q]), error[q]);
    }
  }
  return error;
}

template <class Stepper, class System>
template <class S, class Rates, class AtEnd>
void SubsystemSubsteps<Stepper, System>::integrate_substeps(
    std::size_t first, const std::vector<double>& ends, double h,
    const std::vector<S>& y0, const std::vector<std::vector<S>>& k,
    Rates&& rates, std::vector<S>& own, AtEnd&& at_end) const {
  const std::size_t width = own.size();
  const auto at = static_cast<std::ptrdiff_t>(first);
  const auto last = at + static_cast<std::ptrdiff_t>(width);
  std::vector<std::vector<S>> r(Stepper::stage_rates, std::vector<S>(width));
  std::copy(y0.begin() + at, y0.begin() + last, own.begin());
  double from = 0.0;
  for (const double to : ends) {
    if (!(to > from)) {
      util::stop("Substeps must end in ascending order inside the step");
    }
    // At the step's start the subsystem's rates are the step's own first ones.
    if (from == 0.0) {
      std::copy(k[0].begin() + at, k[0].begin() + last, r[0].begin());
    } else {
      rates(from, own, r[0]);
    }
    take_substep(own, r, from, to, h, rates, own);
    at_end(to, own);
    from = to;
  }
}

template <class Stepper, class System>
template <class S, class Rates>
void SubsystemSubsteps<Stepper, System>::take_substep(
    const std::vector<S>& own, std::vector<std::vector<S>>& r, double from,
    double to, double h, Rates&& rates, std::vector<S>& out) const {
  const double substep = (to - from) * h;
  std::vector<S> stage(own.size());
  for (int i = 1; i < 6; ++i) {
    self().stage_state(i, own, r, substep, stage);
    rates(from + Stepper::stage_fraction(i) * (to - from), stage, r[i]);
  }
  self().step_end(own, r, substep, out);
}

} // namespace ode
} // namespace odelia

#endif
