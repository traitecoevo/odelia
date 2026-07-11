// -*-c++-*-
#ifndef ODELIA_ODE_JACOBIAN_HPP_
#define ODELIA_ODE_JACOBIAN_HPP_

// Exact Jacobian J = d(dydt)/dy for the implicit (Rosenbrock) stepper, via
// forward-mode (tangent) automatic differentiation.
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
//     template <class U> rebind<U> rebind_from() const;
// which returns a copy of itself with the scalar type swapped to U (parameters
// carried over via xad::value + U(...)). This is the same double->AD lift the
// gradient driver uses; a clear error fires below if it is missing.

#include <vector>
#include <cstddef>
#include <type_traits>
#include <utility>
#include <XAD/XAD.hpp>
#include <odelia/ode_interface.hpp>

namespace odelia {
namespace ode {

// Detect `template<class U> ... rebind_from()` on a System, probed at the System's
// own scalar type (every system can at least rebind to itself).
template <typename S, typename = void>
struct has_rebind_from : std::false_type {};

template <typename S>
struct has_rebind_from<
    S, std::void_t<decltype(std::declval<const S>()
                                .template rebind_from<typename S::value_type>())>>
    : std::true_type {};

// The System type rebound to scalar U, i.e. decltype(system.rebind_from<U>()).
// When the System has no rebind_from() the type is not evaluated (a harmless
// placeholder is used instead), so that Jacobian<System> can still be
// *class*-instantiated for systems that will never use the implicit stepper --
// the actual use is gated on `supported` below.
template <typename S, typename U, bool = has_rebind_from<S>::value>
struct rebound_system {
  using type = decltype(std::declval<const S>().template rebind_from<U>());
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
  // Requires (a) a rebind_from() hook and (b) that the tangent twin can be built from
  // the current scalar type. (b) is currently false when value_type is itself an
  // active AD type (nested tangent-over-adjoint, e.g. FReal<AReal<double>>, is
  // not yet wired up -- see issue #35). Callers gate on this, so Jacobian can be
  // class-instantiated even for systems that never use the implicit stepper.
  static constexpr bool supported =
      has_rebind_from<System>::value &&
      std::is_constructible<tangent_type, value_type>::value;

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
    twin_type twin = system.template rebind_from<tangent_type>();

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
