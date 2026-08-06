// -*-c++-*-
#ifndef ODELIA_ODE_STEP_MRI_IMPL_HPP_
#define ODELIA_ODE_STEP_MRI_IMPL_HPP_

// Out-of-line definition of MriStep::step. Split from ode_step_mri.hpp so it can
// be included after SolverInternal is defined: the macro step sub-cycles the fast
// block through SolverInternal (via mri_macro_step in mri.hpp), which therefore
// must already exist. Included at the bottom of ode_solver_internal.hpp.

#include <odelia/ode_step_mri.hpp>
#include <odelia/mri.hpp>

namespace odelia {
namespace ode {

template <class System>
void MriStep<System>::step(System& system, double time, double step_size,
                           state_type& y, state_type& yerr,
                           const state_type& /*dydt_in*/, state_type& dydt_out) {
  if constexpr (supported) {
    const size_t ns = system.slow_size();
    // State is laid out [slow | fast], fast block last.
    state_type x(y.begin(), y.begin() + ns), u(y.begin() + ns, y.end());
    MRICoupling coupling = mri_erk33a();
    MRISchedule sched;
    // Select the fast-block inner: the exact-flow split (Lever 1) when the
    // System opts in at runtime, else the adaptive black-box inner. Both are
    // instantiated only for a System that provides mri_split() (which implies
    // the analytic_flow / residual_rhs hooks the split needs).
    auto run = [&](auto subcycle) {
      mri_macro_step(system, coupling, inner_control, time, step_size,
                     x, u, sched, /*replay=*/false, subcycle);
    };
    if constexpr (mri_wants_split<System>::value) {
      if (system.mri_split()) run(SplitSubcycle{});
      else                    run(AdaptiveSubcycle{});
    } else {
      run(AdaptiveSubcycle{});
    }
    std::copy(x.begin(), x.end(), y.begin());
    std::copy(u.begin(), u.end(), y.begin() + ns);
    // The macro grid is the fixed forcing-kink grid, so there is no macro-level
    // error to report; the fast sub-cycle is error-controlled internally.
    std::fill(yerr.begin(), yerr.end(), value_type(0.0));
    // Rates at the new state -- also syncs the System onto y for history.
    ode::derivs(system, y, dydt_out, time + step_size);
  }
}

}
}

#endif
