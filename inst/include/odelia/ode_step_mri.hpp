// -*-c++-*-
#ifndef ODELIA_ODE_STEP_MRI_HPP_
#define ODELIA_ODE_STEP_MRI_HPP_

// MRI-GARK multirate step, matching the `Step` (RKCK) interface so SolverInternal
// can drive it through the same loop, selectable with `method = "mri"`. One step
// advances the whole state by step_size: an explicit Runge-Kutta on the slow
// block, the fast block sub-cycled by the existing adaptive solver over each leg.
// The macro grid should be the forcing-kink grid, so drive it with fixed steps
// (advance_fixed / step_to); adaptivity lives in the fast sub-cycle.
//
// A System opts in by declaring its fast/slow split on top of the plain ODE
// interface (fast_size / slow_size / coupling_size, slow_rates / fast_rates /
// aggregate; see mri.hpp), with the state laid out fast-block first. Forward
// (double) only; for gradients use the record->replay path.
//
// step() is defined out of line in ode_step_mri_impl.hpp (included by
// ode_solver_internal.hpp after SolverInternal is defined): its body sub-cycles
// the fast block through SolverInternal, so it cannot be defined before it.

#include <vector>
#include <cstddef>
#include <type_traits>
#include <odelia/ode_control.hpp>

namespace odelia {
namespace ode {

// A System is multirate when it declares the fast/slow partition hooks.
template <class, class = void>
struct is_multirate : std::false_type {};
template <class S>
struct is_multirate<S, std::void_t<decltype(std::declval<S&>().fast_size()),
                                   decltype(std::declval<S&>().coupling_size())>>
  : std::true_type {};

template <class System>
class MriStep {
public:
  using value_type = typename System::value_type;
  using state_type = std::vector<value_type>;

  // Available when the System declares the partition and the scalar is the
  // passive double (the record->replay path handles active types).
  static constexpr bool supported =
    is_multirate<System>::value && std::is_same<value_type, double>::value;

  static const bool can_use_dydt_in = false;
  static const bool first_same_as_last = false;

  void resize(size_t size_) { size = size_; }
  size_t order() const { return 3; }   // MRI-GARK-ERK33a

  // One macro step over [time, time+step_size]. Defined in ode_step_mri_impl.hpp.
  void step(System& system, double time, double step_size,
            state_type& y, state_type& yerr,
            const state_type& dydt_in, state_type& dydt_out);

  OdeControl inner_control{1e-6, 1e-6, 1.0, 0.0, 1e-10, 1e10, 1e-4};

private:
  size_t size = 0;
};

}
}

#endif
