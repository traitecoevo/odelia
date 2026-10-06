// -*-c++-*-
#ifndef ODELIA_ODE_CALLBACK_SYSTEM_HPP_
#define ODELIA_ODE_CALLBACK_SYSTEM_HPP_

// A System whose right-hand side is a function handed in at run time (#62),
// rather than a class compiled against odelia. This is how an R closure, a
// Python callable or a plain C++ lambda is stepped: the solver sees an ordinary
// System, and the language binding is a few lines that wrap its callable in a
// std::function. Nothing here knows about R -- this header is part of the
// R-free core and is exercised as plain C++ in tests/standalone/.
//
// What the callable must do: write f(t, y) into dydt, which arrives sized n.
// Optionally a Jacobian callable fills J row-major, J[row * n + col] =
// d f_row / d y_col, with f(t, y) handed in as `dydt` for a finite-difference
// implementation to use; when there is none, fd_jacobian() (ode_jacobian.hpp)
// forms it against the right-hand side at n evaluations. Optionally a
// validity predicate over (t, y) maps to ode_state_valid() (#55). A callable
// that wants the current step rejected throws util::DomainError
// (util::stop_domain); anything else it throws propagates and ends the solve,
// which is the core's rule for telling a bug from a domain refusal.
//
// Evaluation is lazy: set_ode_state() records (t, y) and nothing more, and the
// right-hand side is called when the rates are first read. Three places in the
// solver set the state without wanting rates -- the restore after a caught
// DomainError, and the two set_state_from_system() calls inside
// Solver::set_state() -- and for a right-hand side that is itself an expensive
// computation (regnans's is a demographic equilibrium solve) an evaluation
// nobody reads is the whole cost of a step wasted. The count of right-hand-side
// calls is kept, with the count of Jacobian formations, because for such a
// system those two numbers are what a solve costs.
//
// The last call made during an accepted step is at the accepted (t, y): both
// steppers end with the derivative at the new point and the validity check
// makes no call. A consumer that keeps extra results from its last evaluation
// may rely on that, and should still compare the solver's state with the one it
// was called at.

#include <cstddef>
#include <functional>
#include <string>
#include <vector>
#include <odelia/ode_util.hpp>
#include <odelia/ode_interface.hpp>
#include <odelia/ode_jacobian.hpp>

namespace odelia {
namespace ode {

class CallbackSystem {
public:
  using value_type = double;
  using state_type = std::vector<double>;
  using rhs_type =
      std::function<void(double t, const state_type& y, state_type& dydt)>;
  using jac_type = std::function<void(double t, const state_type& y,
                                      const state_type& dydt, state_type& J)>;
  using valid_type = std::function<bool(double t, const state_type& y)>;

  CallbackSystem() = default;

  // y0 and t0 are the initial state the solver starts from and that reset()
  // returns to; the callable is sized by y0.
  // jac_fd_step and jac_fd_floor are fd_jacobian()'s rel_step and y_floor;
  // see there for why the floor is 1e-5 and not the absolute tolerance.
  CallbackSystem(rhs_type rhs, state_type y0, double t0 = 0.0,
                 jac_type jac = jac_type(), valid_type valid = valid_type(),
                 bool autonomous = false, double jac_fd_step = 1e-6,
                 double jac_fd_floor = 1e-5)
      : rhs_(std::move(rhs)), jac_(std::move(jac)), valid_(std::move(valid)),
        autonomous_(autonomous), jac_fd_step_(jac_fd_step), jac_fd_floor_(jac_fd_floor),
        y_(std::move(y0)), t_(t0), y0_(y_), t0_(t0) {
    if (!rhs_) {
      util::stop("CallbackSystem needs a right-hand side");
    }
    if (y_.empty()) {
      util::stop("CallbackSystem needs a state of at least one element");
    }
    if (!(jac_fd_step_ > 0.0)) {
      util::stop("jac_fd_step must be positive");
    }
    if (!(jac_fd_floor_ > 0.0)) {
      util::stop("jac_fd_floor must be positive");
    }
    dydt_.assign(y_.size(), 0.0);
    dirty_ = true;
  }

  // --- The System interface ------------------------------------------------

  size_t ode_size() const { return y_.size(); }
  double ode_time() const { return t_; }

  template <typename Iterator>
  Iterator set_ode_state(Iterator it, double t) {
    for (size_t i = 0; i < y_.size(); ++i) {
      y_[i] = *it++;
    }
    t_ = t;
    dirty_ = true;
    return it;
  }

  template <typename Iterator>
  Iterator ode_state(Iterator it) const {
    for (size_t i = 0; i < y_.size(); ++i) {
      *it++ = y_[i];
    }
    return it;
  }

  template <typename Iterator>
  Iterator ode_rates(Iterator it) {
    if (dirty_) {
      evaluate();
    }
    for (size_t i = 0; i < dydt_.size(); ++i) {
      *it++ = dydt_[i];
    }
    return it;
  }

  template <typename Iterator>
  Iterator ode_aux(Iterator it) const { return it; }
  size_t aux_size() const { return 0; }

  bool ode_autonomous() const { return autonomous_; }

  bool ode_state_valid(const state_type& y) const {
    return valid_ ? valid_(t_, y) : true;
  }

  void ode_jacobian(const state_type& y, double t, const state_type& dydt,
                    state_type& J) {
    ++n_jac;
    const size_t n = y.size();
    if (jac_) {
      J.assign(n * n, 0.0);
      jac_(t, y, dydt, J);
      if (J.size() != n * n) {
        util::stop("Jacobian callback returned " + std::to_string(J.size()) +
                   " elements, expected " + std::to_string(n * n));
      }
    } else {
      fd_jacobian(*this, y, t, dydt, J, jac_fd_step_, jac_fd_floor_);
    }
  }

  // --- What the Solver wrapper and a driving consumer need -----------------

  // Back to the initial state. The counters are not touched: they describe the
  // callable's life, not the trajectory's.
  void reset() {
    y_ = y0_;
    t_ = t0_;
    dydt_.assign(y_.size(), 0.0);
    dirty_ = true;
  }

  // A different number of unknowns from here on: the consumer's state has
  // changed shape (regnans drops or splits residents between steps). The state
  // itself is set afterwards through set_ode_state(); what this does is make
  // the sizes agree so that Solver::set_state() accepts the new vector.
  void resize(size_t n) {
    if (n == 0) {
      util::stop("CallbackSystem needs a state of at least one element");
    }
    y_.assign(n, 0.0);
    dydt_.assign(n, 0.0);
    y0_ = y_;
    dirty_ = true;
  }

  bool autonomous() const { return autonomous_; }
  double jac_fd_step() const { return jac_fd_step_; }
  double jac_fd_floor() const { return jac_fd_floor_; }
  bool has_jacobian_callback() const { return static_cast<bool>(jac_); }

  // Right-hand-side evaluations and Jacobian formations since construction.
  // Finite-difference Jacobians count their evaluations in n_rhs as well.
  size_t n_rhs = 0;
  size_t n_jac = 0;

private:
  void evaluate() {
    ++n_rhs;
    const size_t n = y_.size();
    dydt_.resize(n);
    rhs_(t_, y_, dydt_);
    if (dydt_.size() != n) {
      util::stop("right-hand side returned " + std::to_string(dydt_.size()) +
                 " rates, expected " + std::to_string(n));
    }
    dirty_ = false;
  }

  rhs_type rhs_;
  jac_type jac_;
  valid_type valid_;
  bool autonomous_ = false;
  double jac_fd_step_ = 1e-6;
  double jac_fd_floor_ = 1e-5;

  state_type y_;
  double t_ = 0.0;
  state_type dydt_;
  bool dirty_ = true;

  state_type y0_;
  double t0_ = 0.0;
};

} // namespace ode
} // namespace odelia

#endif
