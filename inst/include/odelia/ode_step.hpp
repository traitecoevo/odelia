// -*-c++-*-
#ifndef ODELIA_ODE_STEP_HPP_
#define ODELIA_ODE_STEP_HPP_

// The Cash-Karp stepper moved to ode_step_rkck.hpp when it stopped being the
// only one (ode_step_rodas.hpp, ode_step_dopri.hpp). This name is kept because
// plant includes it directly (plant.h); include ode_step_rkck.hpp in new code.

#include <odelia/ode_step_rkck.hpp>

#endif
