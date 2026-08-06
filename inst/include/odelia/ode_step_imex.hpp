// -*-c++-*-
#ifndef ODELIA_ODE_STEP_IMEX_HPP_
#define ODELIA_ODE_STEP_IMEX_HPP_

// IMEX single-step stepper: RODAS4(3) Rosenbrock machinery driven by a Jacobian
// that is populated only on the System's fast (soil) block, by finite
// differences through the full coupled RHS (BlockFdJacobian). The result is
// implicit on the L<=5 soil states -- including the O(N) state-dependent root
// uptake a(theta), which is where the stiffness lives -- and explicit on the
// cohorts, in a single global step (no sub-cycling, unlike MRI).
//
// It reuses RODAS wholesale via the Jacobian-policy template, so there is no new
// integrator to validate: same order, same embedded error pair, same adaptive
// controller. It differs from RODAS in exactly one respect -- the Jacobian
// source -- and that difference is what lets it run on the real coupled Patch,
// which has no rebind<U>() hook and so cannot supply RODAS's forward-AD Jacobian.

#include <odelia/ode_step_rodas.hpp>
#include <odelia/ode_jacobian.hpp>

namespace odelia {
namespace ode {

template <class System>
using ImexStep = RodasStep<System, BlockFdJacobian<System>>;

} // namespace ode
} // namespace odelia

#endif
