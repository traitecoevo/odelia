// -*-c++-*-
#ifndef ODELIA_ODE_JACOBIAN_HPP_
#define ODELIA_ODE_JACOBIAN_HPP_

// Exact Jacobians of the ODE right-hand side via forward-mode (tangent)
// automatic differentiation: the state Jacobian J = d(dydt)/dy (used by the
// implicit Rosenbrock stepper) and the parameter Jacobian d(dydt)/dtheta (used
// by the implicit-function-theorem steady-state sensitivity in
// ode_steady_state.hpp). Both are the same forward sweep on the same active
// twin, differing only in where the unit tangent seed is placed.
//
// The RHS is differentiated on an active "twin" of the System whose scalar type
// is the tangent type FReal<value_type>. Forward mode is used (not adjoint)
// because: for a square N->N Jacobian both cost N sweeps, but forward mode needs
// no tape (no recording, no allocation, no interaction with the single
// thread-local active-tape pointer). It therefore composes cleanly as
// FReal<AReal<double>> when the solver itself is being differentiated by an outer
// adjoint fit -- the tangent layer never contends with the outer tape.
//
// Obtaining the twin requires the System to expose
//     template <class U> System<U> rebind() const;
// which returns a copy of itself with the scalar type swapped to U (parameters
// carried over via xad::value + U(...)). This is the one concept extension the
// implicit stepper adds; a clear error fires below if it is missing.

#include <vector>
#include <cstddef>
#include <type_traits>
#include <utility>
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
// for systems that will never use the implicit stepper -- the actual use is
// gated on `supported` below.
template <typename S, typename U, bool = has_rebind<S>::value>
struct rebound_system {
  using type = decltype(std::declval<const S>().template rebind<U>());
};
template <typename S, typename U>
struct rebound_system<S, U, false> {
  using type = S;
};

// Forward-mode AD Jacobian helper. Owns the active twin and scratch buffers so
// that repeated evaluations (once per accepted step) reuse storage.
template <typename System>
class Jacobian {
public:
  using value_type = typename System::value_type;
  // Tangent scalar: one forward-mode layer on top of the solver's scalar type.
  using tangent_type = typename xad::fwd<value_type>::active_type;
  using twin_type = typename rebound_system<System, tangent_type>::type;

  // Whether the forward-AD Jacobian is instantiable and usable for this System.
  // Requires (a) a rebind() hook and (b) that the tangent twin can be built from
  // the current scalar type. (b) is currently false when value_type is itself an
  // active AD type (nested tangent-over-adjoint, e.g. FReal<AReal<double>>, is
  // not yet wired up -- see issue #35). Callers gate on this, so Jacobian can be
  // class-instantiated even for systems that never use the implicit stepper.
  static constexpr bool supported =
      has_rebind<System>::value &&
      std::is_constructible<tangent_type, value_type>::value;

  // Whether the parameter Jacobian df/dtheta is additionally available: needs
  // `supported` (a twin to differentiate on) plus an `ode_parameters()` hook on
  // that twin exposing pointers to the differentiable parameters to seed.
  static constexpr bool params_supported =
      supported && has_ode_parameters<twin_type>::value;

  void resize(size_t size_) {
    size = size_;
    v.resize(size);
    dydt_ad.resize(size);
  }

  // Compute J = d f / d y at (y, t), written row-major into `J` (size n*n),
  // J[row * n + col] = d f_row / d y_col. Parameters are held fixed (they are
  // seeded with zero tangent), so J is the state Jacobian only.
  void compute(const System& system, const std::vector<value_type>& y,
               double t, std::vector<value_type>& J) {
    // Refresh the twin from the live system each call so current parameters are
    // reflected (cheap: a small value copy). The twin's scalar is the tangent
    // type; its parameters carry zero derivative.
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

  // Compute the parameter Jacobian Jp = d f / d theta at (y, t), written
  // row-major into `Jp` (size n * n_params), Jp[row * n_params + col] =
  // d f_row / d theta_col. `n_params` is set to the number of parameters the
  // System exposes via ode_parameters().
  //
  // Same forward-mode sweep as compute(), but the tangent seed is placed on a
  // *parameter* of the twin rather than a state component: state carries zero
  // derivative, one parameter carries unit derivative per column, so the output
  // tangent is exactly that parameter's column of df/dtheta. This reuses the
  // twin and buffers -- seeding parameters instead of routing through the
  // reverse-mode set_params-on-tape path keeps the whole sweep tape-free, so it
  // never contends with an outer adjoint tape and composes under one later.
  void compute_params(const System& system, const std::vector<value_type>& y,
                      double t, std::vector<value_type>& Jp,
                      size_t& n_params) {
    static_assert(params_supported,
                  "compute_params requires an ode_parameters() hook on the "
                  "System twin");
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

// Finite-difference partial derivative of the RHS with respect to time,
// d f / d t at (y, t). The System stores time as a plain double (not the scalar
// type), so this term cannot be seeded through the twin; a one-sided difference
// is used. It is (near) zero for autonomous systems. Uses value_type arithmetic
// throughout, so it tapes correctly under an outer adjoint fit.
template <typename System>
void dfdt_fd(System& system, const std::vector<typename System::value_type>& y,
             double t, std::vector<typename System::value_type>& out) {
  using value_type = typename System::value_type;
  const size_t n = y.size();
  out.resize(n);
  std::vector<value_type> f0(n), f1(n);

  // Scale the perturbation to the magnitude of t (with a floor for t near 0).
  const double dt = 1e-7 * (std::abs(t) + 1.0);

  ode::derivs(system, y, f0, t);
  ode::derivs(system, y, f1, t + dt);
  for (size_t i = 0; i < n; ++i) {
    out[i] = (f1[i] - f0[i]) / dt;
  }
}

} // namespace ode
} // namespace odelia

#endif
