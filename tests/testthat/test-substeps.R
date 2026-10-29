# Tests for a step that substeps a System's subsystem: where a stiff subsystem
# reads the rest only through a few inputs, the rest takes long steps while the
# subsystem takes substeps of its own, which replays and the sweep take as
# recorded.

# A population x drawing on a two-layer store w, which a seasonal source fills
# and drainage K w^2 carries layer to layer: stiff in w where K is large. The
# uptake c x w is what the rest of the state hands the store. With
# `WithSubsystem` false it is the same System with no subsystem. `Pad` adds that
# many copies of the population, each with x's own dynamics and nothing the
# store reads: more components in the rest of the state, and no other change.
store_system <- '
  template <typename T, bool WithSubsystem = true, int Pad = 0>
  class Store {
  public:
    using value_type = T;
    Store(T a_ = T(0), T m_ = T(0), T K_ = T(0), T s0_ = T(0), T c_ = T(0))
      : a(a_), m(m_), K(K_), s0(s0_), c(c_) { reset(); }

    template <typename, bool, int> friend class Store;

    template <class S2>
    void assign_from(const Store<S2, WithSubsystem, Pad>& src) {
      a = T(odelia::util::to_passive(src.a));
      m = T(odelia::util::to_passive(src.m));
      K = T(odelia::util::to_passive(src.K));
      s0 = T(odelia::util::to_passive(src.s0));
      c = T(odelia::util::to_passive(src.c));
      x = T(odelia::util::to_passive(src.x));
      w1 = T(odelia::util::to_passive(src.w1));
      w2 = T(odelia::util::to_passive(src.w2));
      for (int j = 0; j < Pad; ++j) {
        xs[j] = T(odelia::util::to_passive(src.xs[j]));
      }
      time = src.time;
      fail_after = src.fail_after;
      pulses = src.pulses;
      compute_rates();
    }

    template <class S2>
    Store<S2, WithSubsystem, Pad> rebind_from() const {
      Store<S2, WithSubsystem, Pad> out;
      out.assign_from(*this);
      return out;
    }

    std::vector<T*> ad_parameters() { return {&a, &m, &K, &s0, &c}; }

    template <class F>
    void for_each_active(F&& f) {
      f(a); f(m); f(K); f(s0); f(c);
      f(x); f(w1); f(w2);
      f(dx); f(dw1); f(dw2); f(u1); f(u2);
      for (T& v : xs) f(v);
      for (T& v : dxs) f(v);
    }

    size_t ode_size() const { return 3 + Pad; }
    double ode_time() const { return time; }
    void reset() {
      x = 1.0; w1 = 0.5; w2 = 0.5; time = 0.0;
      xs.assign(Pad, T(1.0));
      dxs.assign(Pad, T(0.0));
      compute_rates();
    }

    template <typename It> It set_ode_state(It it, double time_) {
      time = time_;
      x = *it++;
      w1 = *it++;
      w2 = *it++;
      for (T& v : xs) v = *it++;
      compute_rates();
      return it;
    }
    template <class V> void set_recorded_state(const V& y, double time_) {
      set_ode_state(y.begin(), time_);
    }
    template <typename It> It ode_state(It it) const {
      *it++ = x; *it++ = w1; *it++ = w2;
      for (const T& v : xs) *it++ = v;
      return it;
    }
    template <typename It> It ode_rates(It it) const {
      *it++ = dx; *it++ = dw1; *it++ = dw2;
      for (const T& v : dxs) *it++ = v;
      return it;
    }

    // The store, from index 1.
    odelia::ode::state_range subsystem() const requires WithSubsystem {
      return {1, 2};
    }
    std::vector<double> subsystem_knots(double from, double to) const
      requires WithSubsystem {
      return pulses != nullptr ? pulses->knots(from, to) : std::vector<double>{};
    }
    void subsystem_inputs(std::vector<T>& u) const requires WithSubsystem {
      u = {u1, u2};
    }
    template <class U>
    void subsystem_rates(double t, const std::vector<U>& w,
                         const std::vector<U>& u,
                         std::vector<U>& out) const requires WithSubsystem {
      store_rates(t, w, out);
      out[0] -= u[0];
      out[1] -= u[1];
      if (t > fail_after) {
        out[0] = U(std::nan(""));
      }
    }

    // Past this time the store has NaN rates.
    double fail_after = 1e300;
    // Pulses of inflow on top of the seasonal one.
    const odelia::drivers::Function* pulses = nullptr;

  private:
    template <class U> U at(const T& v) const {
      if constexpr (std::is_same_v<U, T>) {
        return v;
      } else {
        return odelia::util::to_passive(v);
      }
    }
    double pulse_at(double t) const {
      return pulses != nullptr ? std::max(0.0, pulses->evaluate(t)) : 0.0;
    }
    template <class U>
    void store_rates(double t, const std::vector<U>& w, std::vector<U>& out) const {
      const U k = at<U>(K);
      const U inflow = at<U>(s0) * (1.0 + 0.5 * std::sin(2.0 * M_PI * t)) + pulse_at(t);
      out.resize(2);
      out[0] = inflow - k * w[0] * w[0];
      out[1] = k * w[0] * w[0] - k * w[1] * w[1];
    }
    void compute_rates() {
      std::vector<T> w{w1, w2}, f(2);
      store_rates(time, w, f);
      u1 = c * x * w1;
      u2 = c * x * w2;
      dx = a * x * (w1 + w2) - m * x;
      dw1 = f[0] - u1;
      dw2 = f[1] - u2;
      for (int j = 0; j < Pad; ++j) {
        dxs[j] = a * xs[j] * (w1 + w2) - m * xs[j];
      }
    }

    T a, m, K, s0, c;
    T x = 0, w1 = 0, w2 = 0, dx = 0, dw1 = 0, dw2 = 0, u1 = 0, u2 = 0;
    std::vector<T> xs, dxs;
    double time = 0;
  };
'

compile_substep_interface <- function() {
  ensure_ode_interface_loaded()

  odelia_so <- .odelia_test_cache$odelia_so
  # Without contraction, because the tests below hold a replay at another
  # scalar to the run bit for bit (see compile_step_adjoint_interface).
  withr::local_envvar(
    PKG_CPPFLAGS = odelia_cppflags(odelia_include_dir()),
    PKG_CXXFLAGS = "-ffp-contract=off",
    PKG_LIBS = shQuote(normalizePath(odelia_so, winslash = "/", mustWork = TRUE))
  )
  Rcpp::sourceCpp(code = paste0('
    // [[Rcpp::plugins(cpp20)]]
    // Built with -ffp-contract=off (see compile_substep_interface).
    #include <Rcpp.h>
    #include <algorithm>
    #include <cmath>
    #include <limits>
    #include <vector>
    #include <odelia/drivers.hpp>
    #include <odelia/ode_solver.hpp>
    #include <odelia/tangent.hpp>
    ', store_system, '

    using namespace odelia;
    using StoreA = Store<double, true>;
    using StoreP = Store<double, false>;
    // The store among forty more populations.
    using StoreMany = Store<double, true, 40>;

    static_assert(ode::HasSubsystem<StoreA>);
    static_assert(!ode::HasSubsystem<StoreP>);

    ode::Method method_of(const std::string& name) {
      return name == "dopri" ? ode::Method::dopri :
        name == "rodas" ? ode::Method::rodas : ode::Method::rkck;
    }

    // The step-size rule by name, gsl unless asked otherwise.
    ode::OdeControl control_at(double tol, const std::string& rule = "gsl") {
      ode::OdeControl ctl;
      ctl.set_tol_abs(tol);
      ctl.set_tol_rel(tol);
      ctl.set_controller(rule == "hairer" ? ode::Controller::hairer
                                          : ode::Controller::gsl);
      return ctl;
    }

    template <class S = StoreA>
    S store(const std::vector<double>& p) {
      S s(p[0], p[1], p[2], p[3], p[4]);
      return s;
    }

    std::vector<ode::instruction> program_of(const std::vector<double>& pars,
                                             const std::vector<double>& y0,
                                             double t_end, double tol,
                                             const std::string& method,
                                             const std::string& rule) {
      ode::Solver<StoreA> solver(store(pars), control_at(tol, rule),
                                 method_of(method));
      solver.set_state(y0, 0.0);
      solver.advance_adaptive({0.0, t_end});
      return solver.schedule();
    }

    // A run of the store, with the rest of the state `Pad` wide: where it
    // lands, its steps, which of them substepped and how, and what a replay of
    // its recording lands at. `fractions` are the stage fractions every
    // substepped step must land a substep on.
    template <class S>
    Rcpp::List run_store(double tol, double t_end, const std::vector<double>& pars,
                         const std::vector<double>& y0, double fail_after,
                         const std::string& method, const std::string& rule,
                         bool substep_subsystem,
                         const std::vector<double>& fractions) {
      S s = store<S>(pars);
      s.fail_after = fail_after;
      ode::OdeControl ctl = control_at(tol, rule);
      ctl.substep_subsystem = substep_subsystem;
      ode::Solver<S> solver(s, ctl, method_of(method));
      solver.set_keep_states(true);
      solver.set_state(y0, 0.0);
      solver.advance_adaptive({0.0, t_end});
      const std::vector<ode::instruction> program = solver.schedule();
      int substepped = 0, substeps = 0;
      bool fractions_landed = true, ascending = true, sized = true;
      for (const ode::instruction& row : program) {
        if (!row.subsystem) {
          continue;
        }
        const ode::subsystem_substeps& taken = *row.subsystem;
        ++substepped;
        substeps += taken.ends.size();
        sized = sized && taken.input_slope.size() == 2;
        ascending = ascending && std::is_sorted(taken.ends.begin(),
                                                taken.ends.end());
        for (double fraction : fractions) {
          fractions_landed = fractions_landed &&
            std::find(taken.ends.begin(), taken.ends.end(), fraction) !=
              taken.ends.end();
        }
      }
      // A replay of the recording, its rows loaded, repeats the run; the
      // implicit stepper keeps no rows, and replays the program.
      ode::Solver<S> replay(store<S>(pars), ctl, method_of(method));
      replay.set_state(y0, 0.0);
      if (method == "rodas") {
        replay.advance_recorded(program);
      } else {
        replay.advance_recorded(solver.recording());
      }
      return Rcpp::List::create(
        Rcpp::_["y"] = solver.state(),
        Rcpp::_["steps"] = (int) program.size() - 1,
        Rcpp::_["substepped"] = substepped, Rcpp::_["substeps"] = substeps,
        Rcpp::_["fractions_landed"] = fractions_landed,
        Rcpp::_["ascending"] = ascending,
        Rcpp::_["sized"] = sized,
        Rcpp::_["recording_replayed"] = replay.state());
    }

    // [[Rcpp::export]]
    Rcpp::List store_run_(double tol, double t_end,
                          std::vector<double> pars, std::vector<double> y0,
                          double fail_after, std::string method,
                          std::string rule, bool substep,
                          std::vector<double> fractions) {
      return run_store<StoreA>(tol, t_end, pars, y0, fail_after, method, rule,
                               substep, fractions);
    }

    // The error of a step as the rule reads it (OdeControl::error_ratio).
    // [[Rcpp::export]]
    double control_error_ratio(std::string rule, std::vector<double> y,
                               std::vector<double> yerr,
                               std::vector<double> dydt, double h, double tol) {
      ode::OdeControl ctl = control_at(tol, rule);
      return ctl.error_ratio(y.size(), h, y, yerr, dydt);
    }

    // [[Rcpp::export]]
    Rcpp::List store_many_run(double tol, double t_end,
                              std::vector<double> pars, std::vector<double> y0,
                              std::string rule) {
      return run_store<StoreMany>(tol, t_end, pars, y0, 1e300, "rkck", rule,
                                  true, {});
    }

    // A run under pulses of inflow through (x, y), stopping at `stops`: where it
    // lands, the times its steps end, whether each step substepped the store, and
    // the times its substeps end.
    // [[Rcpp::export]]
    Rcpp::List store_pulsed_run(double tol, std::vector<double> stops,
                                std::vector<double> pars, std::vector<double> y0,
                                std::vector<double> x, std::vector<double> y,
                                std::string method = "rkck",
                                std::string rule = "gsl") {
      const drivers::Function pulses(x, y, drivers::Slopes::monotone);
      StoreA s = store(pars);
      s.pulses = &pulses;
      ode::Solver<StoreA> solver(s, control_at(tol, rule), method_of(method));
      solver.set_state(y0, 0.0);
      solver.advance_adaptive(stops);
      const std::vector<ode::instruction> program = solver.schedule();
      std::vector<double> times, substep_ends;
      std::vector<bool> substepped;
      for (size_t k = 0; k < program.size(); ++k) {
        times.push_back(program[k].time);
        if (k > 0) {
          substepped.push_back(program[k].subsystem.has_value());
        }
        if (k > 0 && program[k].subsystem) {
          for (double end : program[k].subsystem->ends) {
            substep_ends.push_back(program[k - 1].time + end * program[k].step_size);
          }
        }
      }
      return Rcpp::List::create(Rcpp::_["y"] = solver.state(),
                                Rcpp::_["steps"] = (int) program.size() - 1,
                                Rcpp::_["times"] = times,
                                Rcpp::_["substepped"] = substepped,
                                Rcpp::_["substep_ends"] = substep_ends);
    }

    // The program of a pulsed run, replayed with its substeps dropped, which a step
    // spanning a knot refuses.
    // [[Rcpp::export]]
    void store_pulsed_replay_without_substeps(double tol,
                                              std::vector<double> pars,
                                              std::vector<double> y0,
                                              std::vector<double> x,
                                              std::vector<double> y,
                                              std::string rule = "gsl") {
      const drivers::Function pulses(x, y, drivers::Slopes::monotone);
      StoreA s = store(pars);
      s.pulses = &pulses;
      ode::Solver<StoreA> run(s, control_at(tol, rule));
      run.set_state(y0, 0.0);
      run.advance_adaptive({0.0, 1.0});
      std::vector<ode::instruction> program = run.schedule();
      for (ode::instruction& row : program) {
        row.subsystem.reset();
      }
      ode::Solver<StoreA> replay(s, control_at(tol, rule));
      replay.set_state(y0, 0.0);
      replay.advance_recorded(program);
    }

    // The same System with no subsystem, stepped plainly.
    // [[Rcpp::export]]
    Rcpp::List store_plain_run(double tol, double t_end,
                               std::vector<double> pars,
                               std::vector<double> y0,
                               std::string method = "rkck",
                               std::string rule = "gsl") {
      ode::Solver<StoreP> solver(store<StoreP>(pars), control_at(tol, rule),
                                 method_of(method));
      solver.set_state(y0, 0.0);
      solver.advance_adaptive({0.0, t_end});
      return Rcpp::List::create(Rcpp::_["y"] = solver.state(),
                                Rcpp::_["steps"] = (int) solver.times().size() - 1);
    }

    // The program a run at `program_pars` took, replayed at `pars` from `y0`:
    // where it lands, and whether its own rows took the substeps the run took.
    // `method` is the replay\'s; the run\'s is `program_method`.
    // [[Rcpp::export]]
    Rcpp::List store_replay(std::vector<double> program_pars,
                            std::vector<double> pars, std::vector<double> y0,
                            std::vector<double> y0_program, double t_end,
                            double tol, std::string method,
                            std::string program_method = "rkck",
                            std::string rule = "gsl") {
      const std::vector<ode::instruction> program =
        program_of(program_pars, y0_program, t_end, tol, program_method, rule);
      ode::Solver<StoreA> replay(store(pars), control_at(tol, rule),
                                 method_of(method));
      replay.set_state(y0, 0.0);
      replay.advance_recorded(program);
      const std::vector<ode::instruction> taken = replay.schedule();
      bool same_steps = taken.size() == program.size();
      for (size_t k = 0; same_steps && k < taken.size(); ++k) {
        same_steps = taken[k].subsystem == program[k].subsystem;
      }
      return Rcpp::List::create(Rcpp::_["y"] = replay.state(),
                                Rcpp::_["same_steps"] = same_steps);
    }

    // The program a run took, replayed at the tangent scalar along `direction`
    // in the parameters: where it lands, and its derivative there.
    // [[Rcpp::export]]
    Rcpp::List store_tangent_replay(std::vector<double> pars, std::vector<double> y0,
                                    double t_end, double tol,
                                    std::vector<double> direction,
                                    std::string method = "rkck",
                                    std::string rule = "gsl") {
      const std::vector<ode::instruction> program =
        program_of(pars, y0, t_end, tol, method, rule);
      using tangent = ode::tangent_scalar<>;
      auto active = store(pars).rebind_from<tangent>();
      std::vector<tangent*> p = active.ad_parameters();
      for (size_t j = 0; j < p.size(); ++j) {
        ode::seed_direction(*p[j], direction[j]);
      }
      ode::Solver<Store<tangent, true>> replay(active, control_at(tol, rule),
                                               method_of(method));
      replay.set_collect(false);
      replay.set_state(y0, 0.0);
      replay.advance_recorded(program);
      std::vector<double> y, dy;
      for (const tangent& v : replay.state()) {
        y.push_back(xad::value(v));
        dy.push_back(ode::derivative_along(v));
      }
      return Rcpp::List::create(Rcpp::_["y"] = y, Rcpp::_["dy"] = dy);
    }

    // The times a run took, replayed as a grid, which a step that would
    // substep the store refuses.
    // [[Rcpp::export]]
    Rcpp::List store_grid_replay(std::vector<double> pars, std::vector<double> y0,
                                 double t_end, double tol,
                                 std::string rule = "gsl") {
      std::vector<ode::instruction> grid =
        program_of(pars, y0, t_end, tol, "rkck", rule);
      for (ode::instruction& row : grid) {
        row = {row.time, std::numeric_limits<double>::quiet_NaN()};
      }
      ode::Solver<StoreA> replay(store(pars), control_at(tol, rule));
      replay.set_state(y0, 0.0);
      replay.advance_recorded(grid);
      return Rcpp::List::create(Rcpp::_["y"] = replay.state());
    }

    // [[Rcpp::export]]
    Rcpp::List store_adjoint(std::vector<double> pars, std::vector<double> y0,
                             double t_end, std::vector<double> lambda_end,
                             double tol, std::string method = "rkck",
                             std::string rule = "gsl") {
      ode::Solver<StoreA> run(store(pars), control_at(tol, rule),
                              method_of(method));
      run.set_keep_states(true);
      run.set_state(y0, 0.0);
      run.advance_adaptive({0.0, t_end});
      ode::adjoint_rows lambda = ode::adjoint_rows::one_row(lambda_end);
      ode::adjoint_rows rows(1, pars.size());
      run.clear_recorded_rates();
      run.solve_adjoint(lambda, rows);
      return Rcpp::List::create(
        Rcpp::_["lambda"] = lambda.to_rows()[0],
        Rcpp::_["parameter_adjoint"] = rows.to_rows()[0],
        Rcpp::_["n_steps"] = (int) run.recording().size() - 1,
        Rcpp::_["recorded_rates"] = (int) run.recorded_rates());
    }
  '), verbose = FALSE)
}

store_run <- function(tol, t_end, pars, y0, fail_after = 1e300, method = "rkck",
                      rule = "gsl", substep = TRUE, fractions = numeric(0)) {
  store_run_(tol, t_end, pars, y0, fail_after, method, rule, substep, fractions)
}

store_pars <- c(1.0, 0.8, 1e4, 200, 0.3)
store_y <- c(1.0, 0.5, 0.5)
# Drainage K = 10: a store inside the stability limit, whose coupling error
# exceeds the stepper's own on it, so the stepper's error never sets the step.
slow <- store_pars * c(1, 1, 1e-3, 1, 1)
# The stages' fractions of a step, where every substepped step lands a substep.
rkck_fractions <- c(1 / 5, 0.3, 3 / 5, 7 / 8, 1)
dopri_fractions <- c(1 / 5, 3 / 10, 4 / 5, 8 / 9, 1)
# Both step-size rules: hairer is the default, gsl the one before it.
rules <- c("gsl", "hairer")

for (rule in rules) testthat::test_that(sprintf("[%s] where the store's error never sets the step, every step is the plain System's", rule), {
  compile_substep_interface()
  off <- store_run(1e-6, 0.5, slow, store_y, rule = rule)
  plain <- store_plain_run(1e-6, 0.5, slow, store_y, rule = rule)
  expect_identical(off$substepped, 0L)
  expect_identical(off$steps, plain$steps)
  expect_identical(off$y, plain$y)
})

for (rule in rules) testthat::test_that(sprintf("[%s] substepped, the stiff store lets the population take long steps, as accurately", rule), {
  compile_substep_interface()
  reference <- store_plain_run(1e-12, 2.0, store_pars, store_y)$y
  ck <- store_plain_run(1e-6, 2.0, store_pars, store_y, rule = rule)
  substepped <- store_run(1e-6, 2.0, store_pars, store_y, rule = rule)
  # Cash-Karp is held to its stability limit by the store's drainage; a step
  # past it takes substeps of the store instead.
  expect_gt(ck$steps, 1000L)
  expect_gt(substepped$substepped, substepped$steps / 2)
  expect_lt(substepped$steps, ck$steps / 10)
  expect_lt(max(abs(substepped$y - reference)), 1e-5)
  # It converges on the same answer as its tolerance tightens.
  tight <- store_run(1e-9, 2.0, store_pars, store_y, rule = rule)
  expect_lt(max(abs(tight$y - reference)), 1e-8)
})

for (rule in rules) testthat::test_that(sprintf("[%s] a store the population takes nothing from keeps substepping", rule), {
  compile_substep_interface()
  # With no uptake the predictor and the corrector end at the same state, which
  # measures no stiffness; the last measure stands.
  no_uptake <- store_pars * c(1, 1, 1, 1, 0)
  ck <- store_plain_run(1e-4, 2.0, no_uptake, store_y, rule = rule)
  substepped <- store_run(1e-4, 2.0, no_uptake, store_y, rule = rule)
  expect_lt(substepped$steps, ck$steps / 10)
})

for (rule in rules) testthat::test_that(sprintf("[%s] each row records its input slope and substeps, landing on every stage fraction", rule), {
  compile_substep_interface()
  r <- store_run(1e-6, 2.0, store_pars, store_y, rule = rule,
                 fractions = rkck_fractions)
  expect_true(r$sized)
  expect_true(r$ascending)
  expect_true(r$fractions_landed)
  expect_gt(r$substeps, 5L * r$substepped)
  expect_identical(r$recording_replayed, r$y)
})

testthat::test_that("a store whose rates fail at every substep stops the run, saying why", {
  compile_substep_interface()
  expect_error(store_run(1e-6, 2.0, store_pars, store_y, fail_after = 1.0),
               "subsystem fails its error test at the smallest substep")
})

for (rule in rules) testthat::test_that(sprintf("[%s] a replay takes the run's substeps, at the run's parameters to the bit", rule), {
  compile_substep_interface()
  run <- store_run(1e-6, 2.0, store_pars, store_y, rule = rule)
  same <- store_replay(store_pars, store_pars, store_y, store_y, 2.0, 1e-6,
                       "rkck", rule = rule)
  expect_identical(same$y, run$y)
  expect_true(same$same_steps)
  moved <- store_pars * c(1, 1, 1.01, 1, 1)
  away <- store_replay(store_pars, moved, store_y, store_y, 2.0, 1e-6, "rkck",
                       rule = rule)
  expect_true(away$same_steps)
  expect_false(identical(away$y, run$y))
  expect_error(store_replay(store_pars, store_pars, store_y, store_y, 2.0, 1e-6,
                            "rodas", rule = rule),
               "taken only by method='rkck' or method='dopri'")
})

for (rule in rules) testthat::test_that(sprintf("[%s] a grid of a run's times is refused where a step would substep the store", rule), {
  compile_substep_interface()
  expect_error(store_grid_replay(store_pars, store_y, 0.5, 1e-6, rule),
               "replay the run's program")
  # Where the store's error never sets the step no step substeps, and the grid
  # is stepped.
  off <- store_run(1e-6, 0.5, slow, store_y, rule = rule)
  expect_equal(store_grid_replay(slow, store_y, 0.5, 1e-6, rule)$y, off$y,
               tolerance = 1e-12)
})

for (rule in rules) testthat::test_that(sprintf("[%s] a tangent replay takes the run's substeps, its value the run's to the bit", rule), {
  compile_substep_interface()
  direction <- c(0, 0, 1, 0, 0)
  r <- store_tangent_replay(store_pars, store_y, 2.0, 1e-6, direction,
                            rule = rule)
  run <- store_run(1e-6, 2.0, store_pars, store_y, rule = rule)
  expect_identical(r$y, run$y)
  e <- 1e-4 * store_pars[3]
  replay <- function(pars) {
    store_replay(store_pars, pars, store_y, store_y, 2.0, 1e-6, "rkck",
                 rule = rule)$y
  }
  expected <- (replay(store_pars + e * direction) -
               replay(store_pars - e * direction)) / (2 * e)
  expect_equal(r$dy, expected, tolerance = 1e-5)
})

# The sweep's adjoints against central differences of the replayed program, for
# a run under `method` and `rule`.
expect_sweep_matches_differences <- function(method, rule) {
  t_end <- 2.0
  tol <- 1e-6
  lambda_end <- c(0.4, -1.3, 0.7)
  r <- store_adjoint(store_pars, store_y, t_end, lambda_end, tol, method, rule)
  expect_identical(r$recorded_rates, 6L * r$n_steps)

  eps <- 1e-4
  replay <- function(pars, y0) {
    store_replay(store_pars, pars, y0, store_y, t_end, tol, method, method,
                 rule)$y
  }
  expected_state <- sapply(seq_along(store_y), function(j) {
    e <- eps * store_y[j]
    up <- store_y; up[j] <- up[j] + e
    down <- store_y; down[j] <- down[j] - e
    sum((replay(store_pars, up) - replay(store_pars, down)) / (2 * e) * lambda_end)
  })
  expected_pars <- sapply(seq_along(store_pars), function(j) {
    e <- eps * store_pars[j]
    up <- store_pars; up[j] <- up[j] + e
    down <- store_pars; down[j] <- down[j] - e
    sum((replay(up, store_y) - replay(down, store_y)) / (2 * e) * lambda_end)
  })
  expect_equal(r$lambda / expected_state, rep(1, 3), tolerance = 1e-5)
  expect_equal(r$parameter_adjoint / expected_pars, rep(1, 5), tolerance = 1e-5)
}

for (rule in rules) testthat::test_that(sprintf("[%s] the sweep through substepped steps matches central differences of the replayed program", rule), {
  compile_substep_interface()
  expect_sweep_matches_differences("rkck", rule)
})

# Pulses of inflow two hundredths wide, every quarter of a unit, through a
# monotone cubic; its knots are those beside a non-zero value, where the cubic
# changes.
pulse_x <- seq(0, 2.5, by = 0.01)
pulse_y <- ifelse(seq_along(pulse_x) %% 25 == 5, 300, 0)
pulse_knots <- local({
  nz <- pulse_y != 0
  n <- length(nz)
  pulse_x[nz | c(FALSE, nz[-n]) | c(nz[-1], FALSE)]
})

for (rule in rules) testthat::test_that(sprintf("[%s] under Dormand-Prince the store substeps, replays and sweeps as under Cash-Karp", rule), {
  compile_substep_interface()
  reference <- store_plain_run(1e-12, 2.0, store_pars, store_y)$y
  dp <- store_plain_run(1e-6, 2.0, store_pars, store_y, "dopri", rule)
  substepped <- store_run(1e-6, 2.0, store_pars, store_y, method = "dopri",
                          rule = rule, fractions = dopri_fractions)
  expect_gt(dp$steps, 1000L)
  expect_gt(substepped$substepped, substepped$steps / 2)
  expect_lt(substepped$steps, dp$steps / 10)
  expect_lt(max(abs(substepped$y - reference)), 1e-5)
  expect_true(substepped$fractions_landed)
  expect_true(substepped$ascending)
  expect_identical(substepped$recording_replayed, substepped$y)
  # A replay of the program takes the run's substeps, to the bit at the run's
  # parameters; a tangent replay too, with the derivative of what it took.
  same <- store_replay(store_pars, store_pars, store_y, store_y, 2.0, 1e-6,
                       "dopri", "dopri", rule)
  expect_identical(same$y, substepped$y)
  expect_true(same$same_steps)
  direction <- c(0, 0, 1, 0, 0)
  tangent <- store_tangent_replay(store_pars, store_y, 2.0, 1e-6, direction,
                                  "dopri", rule)
  expect_identical(tangent$y, substepped$y)
  e <- 1e-4 * store_pars[3]
  replay <- function(pars) {
    store_replay(store_pars, pars, store_y, store_y, 2.0, 1e-6, "dopri",
                 "dopri", rule)$y
  }
  expected <- (replay(store_pars + e * direction) -
               replay(store_pars - e * direction)) / (2 * e)
  expect_equal(tangent$dy, expected, tolerance = 1e-5)
  expect_sweep_matches_differences("dopri", rule)
  # A step spanning a knot of the inflow substeps under Dormand-Prince too.
  knots <- pulse_knots[pulse_knots < 1.0]
  free <- store_pulsed_run(1e-6, c(0, 1.0), slow, store_y, pulse_x, pulse_y,
                           "dopri", rule)
  t <- free$times
  spans <- vapply(seq_along(free$substepped), function(k)
    any(knots > t[k] & knots < t[k + 1]), TRUE)
  expect_true(any(spans))
  expect_true(all(free$substepped[spans]))
})

for (rule in rules) testthat::test_that(sprintf("[%s] under the implicit stepper the store is taken with the rest, and no step substeps", rule), {
  compile_substep_interface()
  reference <- store_plain_run(1e-12, 0.5, store_pars, store_y)$y
  r <- store_run(1e-6, 0.5, store_pars, store_y, method = "rodas", rule = rule)
  expect_identical(r$substepped, 0L)
  expect_lt(max(abs(r$y - reference) / (1 + abs(reference))), 1e-5)
})

testthat::test_that("the control reads a step's error in its rule's norm", {
  compile_substep_interface()
  y <- c(1, 2, 3)
  yerr <- c(1e-6, 2e-6, 0)
  dydt <- c(0, 0, 0)
  ratios <- abs(yerr) / (1e-6 * abs(y) + 1e-6)
  expect_equal(control_error_ratio("gsl", y, yerr, dydt, 0.1, 1e-6), max(ratios))
  expect_equal(control_error_ratio("hairer", y, yerr, dydt, 0.1, 1e-6),
               sqrt(mean(ratios^2)))
  expect_identical(control_error_ratio("hairer", y, c(1e-6, NaN, 0), dydt, 0.1,
                                       1e-6), Inf)
})

testthat::test_that("the store's error sets the step in the control's norm", {
  compile_substep_interface()
  # K = 100: inside the stability limit, so the store substeps only where its
  # own error set the step before.
  mild <- store_pars * c(1, 1, 1e-2, 1, 1)
  many_y <- c(store_y, rep(1, 40))
  alone <- store_run(1e-6, 2.0, mild, store_y, rule = "gsl")
  expect_gt(alone$substepped, 0L)
  expect_identical(store_run(1e-6, 2.0, mild, store_y, rule = "gsl",
                             substep = FALSE)$substepped, 0L)
  # Under the max norm forty copies of the population change nothing: the same
  # steps, substeps and numbers as the store alone.
  many <- store_many_run(1e-6, 2.0, mild, many_y, "gsl")
  expect_identical(many$steps, alone$steps)
  expect_identical(many$substepped, alone$substepped)
  expect_identical(many$y[1:3], alone$y)
  # Under the RMS norm the copies dilute the error and the steps grow; the
  # store's error sets the step where it is the larger part of the whole, and
  # the run stays within tolerance.
  alone_rms <- store_run(1e-6, 2.0, mild, store_y, rule = "hairer")
  many_rms <- store_many_run(1e-6, 2.0, mild, many_y, "hairer")
  expect_lt(many_rms$steps, alone_rms$steps)
  expect_gt(many_rms$substepped, 0L)
  reference <- store_plain_run(1e-12, 2.0, mild, store_y)$y
  expect_lt(max(abs(many_rms$y[1:3] - reference) / (1 + abs(reference))), 1e-5)
})

for (rule in rules) testthat::test_that(sprintf("[%s] a step spanning the inflow's knots substeps the store, landing a substep on each", rule), {
  compile_substep_interface()
  t_end <- 1.0
  knots <- pulse_knots[pulse_knots < t_end]
  # The slow store holds each pulse as a rise in its level.
  reference <- store_pulsed_run(1e-11, c(0, knots, t_end), slow, store_y,
                                pulse_x, pulse_y)$y
  stopped <- store_pulsed_run(1e-6, c(0, knots, t_end), slow, store_y,
                              pulse_x, pulse_y, "rkck", rule)
  free <- store_pulsed_run(1e-6, c(0, t_end), slow, store_y,
                           pulse_x, pulse_y, "rkck", rule)
  # Stopped at every knot, no step spans one and none substeps the store.
  expect_false(any(stopped$substepped))
  t <- free$times
  spans <- vapply(seq_along(free$substepped), function(k)
    any(knots > t[k] & knots < t[k + 1]), TRUE)
  expect_true(any(spans))
  expect_true(all(free$substepped[spans]))
  landed <- vapply(knots, function(k) any(abs(c(t, free$substep_ends) - k) < 1e-12), TRUE)
  expect_true(all(landed))
  # Within the tolerance in the solver's norm, |error| / (1 + |y|).
  error <- function(y) max(abs(y - reference) / (1 + abs(reference)))
  expect_lt(error(stopped$y), 1e-6)
  expect_lt(error(free$y), 1e-6)
  # A replay that lost those substeps is refused at the first step across a knot.
  expect_error(store_pulsed_replay_without_substeps(1e-6, slow, store_y,
                                                    pulse_x, pulse_y, rule),
               "carries no substeps")
})
