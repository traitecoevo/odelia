// -*-c++-*-
#ifndef ODELIA_ODE_JACOBIAN_HPP_
#define ODELIA_ODE_JACOBIAN_HPP_

// Jacobian J = d(dydt)/dy for the implicit (Rosenbrock) stepper, from one of two
// sources, in order of preference:
//
//   1. The System's own hook (#62),
//        void ode_jacobian(const state_type& y, double t,
//                          const state_type& dydt, state_type& J);
//      written row-major, J[row * n + col] = d f_row / d y_col. Whatever the
//      system knows is accepted here: an analytic matrix, a closure handed in
//      from R, or finite differences through fd_jacobian() below. `dydt` is
//      f(t, y), already in hand at the start of a step, so a finite-difference
//      implementation costs n evaluations rather than n + 1.
//
//   2. Forward-mode (tangent) automatic differentiation on an active "twin" of
//      the System whose scalar type is FReal<value_type>. Forward mode is used
//      (not adjoint) because: for a square N->N Jacobian both cost N sweeps, but
//      forward mode needs no tape (no recording, no allocation, no interaction
//      with the single thread-local active-tape pointer). It therefore composes
//      cleanly as FReal<AReal<double>> when the solver itself is being
//      differentiated by an outer adjoint fit -- the tangent layer never contends
//      with the outer tape. Obtaining the twin requires the System to expose
//        template <class U> System<U> rebind() const;
//      which returns a copy of itself with the scalar type swapped to U
//      (parameters carried over via xad::value + U(...)).
//
// A system declaring neither cannot use the implicit stepper: `supported` says so
// at compile time and the stepper raises a clear error. When both exist the hook
// wins -- whoever wrote it meant it.
//
// The AD route also yields the parameter Jacobian d(dydt)/dtheta, used by the
// implicit-function-theorem steady-state sensitivity in ode_steady_state.hpp:
// the same forward sweep on the same twin, differing only in where the unit
// tangent seed is placed. It needs the twin (a hook knows nothing about
// parameters), so it is gated on `params_supported` below.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>
#include <XAD/XAD.hpp>
#include <odelia/ode_interface.hpp>

namespace odelia {
namespace ode {

// Detect `template<class U> ... rebind()` on a System, probed at the System's own
// scalar type (every system can at least rebind to itself).
template <typename S, typename = void>
struct has_rebind : std::false_type {};

template <typename S>
struct has_rebind<
    S, std::void_t<decltype(std::declval<const S>()
                                .template rebind<typename S::value_type>())>>
    : std::true_type {};

// Detect the `ode_jacobian(y, t, dydt, J)` hook (#62). Same shape of probe as
// has_state_check in ode_interface.hpp: a system that omits the member is
// unaffected, and nothing is called on its behalf.
template <typename S>
class has_jacobian {
  typedef char true_type;
  typedef long false_type;
  template <typename C> static true_type test(decltype(&C::ode_jacobian));
  template <typename C> static false_type test(...);
public:
  enum { value = sizeof(test<S>(0)) == sizeof(true_type) };
};

// Detect a `std::vector<scalar*> ode_parameters()` hook: pointers to the
// system's differentiable parameters, in a fixed order, used to seed parameter
// tangents for the forward-mode parameter Jacobian df/dtheta. A system that
// omits it simply cannot have its parameter sensitivity taken (gated below).
template <typename S, typename = void>
struct has_ode_parameters : std::false_type {};

template <typename S>
struct has_ode_parameters<
    S, std::void_t<decltype(std::declval<S&>().ode_parameters())>>
    : std::true_type {};

// The System type rebound to scalar U, i.e. decltype(system.rebind<U>()). When
// the System has no rebind() the type is not evaluated (a harmless placeholder
// is used instead), so that Jacobian<System> can still be *class*-instantiated
// for systems that will never use the AD route -- the actual use is gated on
// `ad_supported` below.
template <typename S, typename U, bool = has_rebind<S>::value>
struct rebound_system {
  using type = decltype(std::declval<const S>().template rebind<U>());
};
template <typename S, typename U>
struct rebound_system<S, U, false> {
  using type = S;
};

// Jacobian helper. Owns the active twin's scratch buffers so that repeated
// evaluations (once per accepted step) reuse storage.
template <typename System>
class Jacobian {
public:
  using value_type = typename System::value_type;
  // Tangent scalar: one forward-mode layer on top of the solver's scalar type.
  using tangent_type = typename xad::fwd<value_type>::active_type;
  using twin_type = typename rebound_system<System, tangent_type>::type;

  // Whether the forward-AD route is instantiable and usable for this System.
  // Requires (a) a rebind() hook and (b) that the tangent twin can be built from
  // the current scalar type. (b) is currently false when value_type is itself an
  // active AD type (nested tangent-over-adjoint, e.g. FReal<AReal<double>>, is
  // not yet wired up -- see issue #36).
  static constexpr bool ad_supported =
      has_rebind<System>::value &&
      std::is_constructible<tangent_type, value_type>::value;

  // Whether a Jacobian can be had at all: the system's own hook, or the AD
  // route. Callers gate on this, so Jacobian can be class-instantiated even for
  // systems that never use the implicit stepper.
  static constexpr bool supported = has_jacobian<System>::value || ad_supported;

  // Whether the parameter Jacobian df/dtheta is additionally available: needs
  // the AD route (a twin to differentiate on; a hook says nothing about
  // parameters) plus an `ode_parameters()` hook on that twin exposing pointers
  // to the differentiable parameters to seed.
  static constexpr bool params_supported =
      ad_supported && has_ode_parameters<twin_type>::value;

  void resize(size_t size_) {
    size = size_;
    v.resize(size);
    dydt_ad.resize(size);
  }

  // Compute J = d f / d y at (y, t), written row-major into `J` (size n*n),
  // J[row * n + col] = d f_row / d y_col. `dydt` is f(t, y), handed through to
  // a hook that can use it. Parameters are held fixed (on the AD route they are
  // seeded with zero tangent), so J is the state Jacobian only.
  //
  // The system is taken mutable because a hook may need to evaluate it (finite
  // differences do), and the stepper already holds it that way.
  void compute(System& system, const std::vector<value_type>& y, double t,
               const std::vector<value_type>& dydt,
               std::vector<value_type>& J) {
    if constexpr (has_jacobian<System>::value) {
      J.assign(size * size, value_type(0.0));
      system.ode_jacobian(y, t, dydt, J);
    } else {
      // Refresh the twin from the live system each call so current parameters
      // are reflected (cheap: a small value copy). The twin's scalar is the
      // tangent type; its parameters carry zero derivative.
      twin_type twin = system.template rebind<tangent_type>();

      for (size_t j = 0; j < size; ++j) {
        v[j] = tangent_type(y[j]);
      }

      J.assign(size * size, value_type(0.0));
      for (size_t col = 0; col < size; ++col) {
        xad::derivative(v[col]) = 1.0;
        ode::derivs(twin, v, dydt_ad, t);
        for (size_t row = 0; row < size; ++row) {
          J[row * size + col] = xad::derivative(dydt_ad[row]);
        }
        xad::derivative(v[col]) = 0.0;
      }
    }
  }

  // Compute the parameter Jacobian Jp = d f / d theta at (y, t), written
  // row-major into `Jp` (size n * n_params), Jp[row * n_params + col] =
  // d f_row / d theta_col. `n_params` is set to the number of parameters the
  // System exposes via ode_parameters().
  //
  // Same forward-mode sweep as compute()'s AD route, but the tangent seed is
  // placed on a *parameter* of the twin rather than a state component: state
  // carries zero derivative, one parameter carries unit derivative per column,
  // so the output tangent is exactly that parameter's column of df/dtheta. This
  // reuses the twin and buffers -- seeding parameters instead of routing through
  // the reverse-mode set_params-on-tape path keeps the whole sweep tape-free, so
  // it never contends with an outer adjoint tape and composes under one later.
  void compute_params(const System& system, const std::vector<value_type>& y,
                      double t, std::vector<value_type>& Jp,
                      size_t& n_params) {
    static_assert(params_supported,
                  "compute_params requires rebind() and an ode_parameters() "
                  "hook on the System twin");
    twin_type twin = system.template rebind<tangent_type>();

    for (size_t j = 0; j < size; ++j) {
      v[j] = tangent_type(y[j]);
    }

    // Pointers into the twin's own parameter storage; valid for the lifetime of
    // `twin`. Seeding a tangent here propagates through compute_rates().
    std::vector<tangent_type*> params = twin.ode_parameters();
    n_params = params.size();
    Jp.assign(size * n_params, value_type(0.0));

    for (size_t col = 0; col < n_params; ++col) {
      xad::derivative(*params[col]) = 1.0;
      ode::derivs(twin, v, dydt_ad, t);
      for (size_t row = 0; row < size; ++row) {
        Jp[row * n_params + col] = xad::derivative(dydt_ad[row]);
      }
      xad::derivative(*params[col]) = 0.0;
    }
  }

private:
  size_t size = 0;
  std::vector<tangent_type> v;
  std::vector<tangent_type> dydt_ad;
};

// Forward-difference Jacobian of the right-hand side, for a system that has no
// rebind() to differentiate through: the one-line body of an ode_jacobian() hook
// on such a system. One evaluation per column, at y + h_j e_j with
// h_j = rel_step * max(|y_j|, y_floor), against the `dydt` already known at y,
// written row-major like Jacobian::compute(). The system is left on the last
// perturbed point; the stepper sets it again before anything reads it.
//
// rel_step is a trade between truncation and round-off and the right value
// depends on how clean the right-hand side is: 1e-6 suits an exactly evaluated
// function, while one that is itself an iterative solve (a tolerance inside it)
// wants a larger step, 1e-5 or so, to stay above its noise floor.
//
// y_floor is the size below which a component is treated as "small" and the
// perturbation stops shrinking with it. The default 1e-5 is rodas.f's. It has
// to sit between two failures: a floor of 1 perturbs a component of size 1e-5
// by a tenth of itself (on Robertson's kinetics RODAS then takes 19000 steps
// where 430 do, measured), while a floor of the solve's absolute tolerance
// makes the perturbation rel_step * tol_abs, which at tol_abs = 1e-10 is 1e-16
// and vanishes in the subtraction (J = 0 for a unit derivative, measured).
// Lower it only for a state whose components legitimately live below 1e-5,
// keeping rel_step * y_floor well above eps * |f|.
//
// A right-hand side that refuses y + h e_j (util::DomainError: a component
// sitting on the upper edge of its domain) is asked at y - h e_j instead, so
// that a bounded state on its boundary does not make every retry fail the
// same way; the Jacobian does not depend on the step size, so a smaller step
// could never have helped. If both sides refuse, the refusal propagates and
// the stepper treats it as it treats any domain refusal.
template <typename System>
void fd_jacobian(System& system,
                 const std::vector<typename System::value_type>& y, double t,
                 const std::vector<typename System::value_type>& dydt,
                 std::vector<typename System::value_type>& J,
                 double rel_step = 1e-6, double y_floor = 1e-5) {
  using value_type = typename System::value_type;
  const size_t n = y.size();
  J.assign(n * n, value_type(0.0));
  std::vector<value_type> yj(y);
  std::vector<value_type> fj(n);
  for (size_t col = 0; col < n; ++col) {
    double h = rel_step * std::max(std::abs(util::to_passive(y[col])), y_floor);
    yj[col] = y[col] + value_type(h);
    try {
      ode::derivs(system, yj, fj, t);
    } catch (const util::DomainError&) {
      h = -h;
      yj[col] = y[col] + value_type(h);
      ode::derivs(system, yj, fj, t);
    }
    for (size_t row = 0; row < n; ++row) {
      J[row * n + col] = (fj[row] - dydt[row]) / value_type(h);
    }
    yj[col] = y[col];
  }
}

// Finite-difference partial derivative of the RHS with respect to time,
// d f / d t at (y, t), against `dydt` = f(t, y) already in hand. The System
// stores time as a plain double (not the scalar type), so this term cannot be
// seeded through an AD twin; a one-sided difference is used. It is (near) zero
// for autonomous systems, and a system that declares ode_autonomous() is not
// asked for it at all (see ode_step_rodas.hpp). Uses value_type arithmetic
// throughout, so it tapes correctly under an outer adjoint fit.
template <typename System>
void dfdt_fd(System& system, const std::vector<typename System::value_type>& y,
             double t, const std::vector<typename System::value_type>& dydt,
             std::vector<typename System::value_type>& out) {
  using value_type = typename System::value_type;
  const size_t n = y.size();
  out.resize(n);
  std::vector<value_type> f1(n);

  // Scale the perturbation to the magnitude of t (with a floor for t near 0).
  const double dt = 1e-7 * (std::abs(t) + 1.0);

  ode::derivs(system, y, f1, t + dt);
  for (size_t i = 0; i < n; ++i) {
    out[i] = (f1[i] - dydt[i]) / dt;
  }
}

} // namespace ode
} // namespace odelia

#endif
