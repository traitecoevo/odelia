#ifndef ODELIA_ODE_FIT_HPP_
#define ODELIA_ODE_FIT_HPP_

#include <XAD/XAD.hpp>
#include <odelia/ode_solver.hpp>
#include <vector>

namespace odelia {
namespace ode {

// Result from gradient computation
struct FitResult {
    double loss;
    std::vector<double> gradient;
};

// Compute loss and gradient w.r.t. initial conditions
// System must be templated on xad::adj<double>::active_type
template <typename System>
FitResult compute_gradient_ic(
    System& system,
    const std::vector<double>& ic,
    const std::vector<double>& times,
    const std::vector<std::vector<double>>& targets,
    const std::vector<size_t>& obs_indices,
    OdeControl& control
) {
    using ad = xad::adj<double>;
    using ADType = ad::active_type;
    
    // 1. Create fresh tape
    ad::tape_type tape;
    
    // 2. Convert IC to AD types and register as inputs
    std::vector<ADType> ic_ad(ic.begin(), ic.end());
    for (auto& x : ic_ad) {
        tape.registerInput(x);
    }
    tape.newRecording();
    
    // 3. Set initial state
    system.set_ode_state(ic_ad.begin(), times[0]);
    
    // 4. Create solver (no history collection to avoid tape pollution)
    Solver<System> solver(system, control);
    solver.set_collect(false);
    
// 5. Run ODE and accumulate loss at observation points
    ADType loss = ADType(0.0);
    size_t target_idx = 0;
    
    // First check if initial time is an observation point
    if (target_idx < obs_indices.size() && obs_indices[target_idx] == 0) {
        auto state = solver.state();
        for (size_t j = 0; j < targets[target_idx].size() && j < state.size(); ++j) {
            ADType diff = state[j] - ADType(targets[target_idx][j]);
            loss += diff * diff;
        }
        target_idx++;
    }
    
    // Now step through time and check for observations
    for (size_t i = 1; i < times.size(); ++i) {
        // Step to next time
        std::vector<double> step_times = {times[i-1], times[i]};
        solver.advance_fixed(step_times);
        
        // Check if current time index is an observation point
        while (target_idx < obs_indices.size() && obs_indices[target_idx] == i) {
            auto state = solver.state();
            for (size_t j = 0; j < targets[target_idx].size() && j < state.size(); ++j) {
                ADType diff = state[j] - ADType(targets[target_idx][j]);
                loss += diff * diff;
            }
            target_idx++;
        }
    }
    
    // 6. Backpropagate
    tape.registerOutput(loss);
    xad::derivative(loss) = 1.0;
    tape.computeAdjoints();
    
    // 7. Extract results
    FitResult result;
    result.loss = xad::value(loss);
    result.gradient.resize(ic.size());
    for (size_t i = 0; i < ic.size(); ++i) {
        result.gradient[i] = xad::derivative(ic_ad[i]);
    }
    
    return result;
}

// Compute loss and gradient w.r.t. parameters
// Parameters are part of the System, so we need to recreate it
template <typename System, typename... Pars>
FitResult compute_gradient_params(
    const std::vector<double>& params,
    const std::vector<double>& ic,
    const std::vector<double>& times,
    const std::vector<std::vector<double>>& targets,
    const std::vector<size_t>& obs_indices,
    OdeControl& control
) {
    using ad = xad::adj<double>;
    using ADType = ad::active_type;
    
    // 1. Create fresh tape
    ad::tape_type tape;
    
    // 2. Convert params to AD types and register as inputs
    std::vector<ADType> params_ad(params.begin(), params.end());
    for (auto& p : params_ad) {
        tape.registerInput(p);
    }
    tape.newRecording();
    
    // 3. Create system with AD parameters
    System system(params_ad[0], params_ad[1], params_ad[2]);
    
    // 4. Fixed IC (convert to AD but not differentiated)
    std::vector<ADType> ic_ad(ic.begin(), ic.end());
    system.set_ode_state(ic_ad.begin(), times[0]);
    
    // 5. Create solver
    Solver<System> solver(system, control);
    solver.set_collect(false);
    
// 6. Run ODE and accumulate loss
    ADType loss = ADType(0.0);
    size_t target_idx = 0;
    
    // First check if initial time is an observation point
    if (target_idx < obs_indices.size() && obs_indices[target_idx] == 0) {
        auto state = solver.state();
        for (size_t j = 0; j < targets[target_idx].size() && j < state.size(); ++j) {
            ADType diff = state[j] - ADType(targets[target_idx][j]);
            loss += diff * diff;
        }
        target_idx++;
    }
    
    // Now step through time and check for observations
    for (size_t i = 1; i < times.size(); ++i) {
        std::vector<double> step_times = {times[i-1], times[i]};
        solver.advance_fixed(step_times);
        
        // Check if current time index is an observation point
        while (target_idx < obs_indices.size() && obs_indices[target_idx] == i) {
            auto state = solver.state();
            for (size_t j = 0; j < targets[target_idx].size() && j < state.size(); ++j) {
                ADType diff = state[j] - ADType(targets[target_idx][j]);
                loss += diff * diff;
            }
            target_idx++;
        }
    }
    
    // 7. Backpropagate
    tape.registerOutput(loss);
    xad::derivative(loss) = 1.0;
    tape.computeAdjoints();
    
    // 8. Extract results
    FitResult result;
    result.loss = xad::value(loss);
    result.gradient.resize(params.size());
    for (size_t i = 0; i < params.size(); ++i) {
        result.gradient[i] = xad::derivative(params_ad[i]);
    }
    
    return result;
}

} // namespace ode
} // namespace odelia

#endif