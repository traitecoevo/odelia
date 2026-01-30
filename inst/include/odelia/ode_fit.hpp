#ifndef ODELIA_ODE_FIT_HPP_
#define ODELIA_ODE_FIT_HPP_

#include <XAD/XAD.hpp>
#include <odelia/ode_solver.hpp>
#include <optional>
#include <vector>

namespace odelia {
namespace ode {

// Loss function: sum of squared errors
template<typename T>
T sum_of_squares(const std::vector<std::vector<T>>& obs, 
                 const std::vector<std::vector<double>>& targets) {
    T loss(0.0);
    for (size_t i = 0; i < obs.size(); ++i) {
        for (size_t j = 0; j < obs[i].size(); ++j) {
            T diff = obs[i][j] - targets[i][j]; 
            loss += diff * diff;
        }
    }
    return loss;
}

// Compute gradient w.r.t. initial conditions and/or parameters
// Returns pair of (loss, gradient)
template<typename Solver>
std::pair<double, std::vector<double>> compute_gradient(
    Solver& solver,
    std::optional<std::vector<double>> ic = std::nullopt,
    std::optional<std::vector<double>> params = std::nullopt
) {
    using ad = xad::adj<double>;
    using ad_type = ad::active_type;
    
    // At least one must be provided
    if (!ic && !params) {
        util::stop("Must provide at least one of 'ic' or 'params'");
    }
    
    // Create tape
    ad::tape_type tape;
    
    // Vector to track all inputs
    std::vector<ad_type> inputs;
    
    // Get reference to system
    auto& system = solver.get_system_ref();
    
    // Register parameters if provided
    if (params) {
        std::vector<ad_type> params_ad(params->begin(), params->end());
        for (auto& p : params_ad) {
            tape.registerInput(p);
            inputs.push_back(p);
        }
    }
    
    // Register initial conditions if provided
    if (ic) {
        std::vector<ad_type> ic_ad(ic->begin(), ic->end());
        for (auto& x : ic_ad) {
            tape.registerInput(x);
            inputs.push_back(x);
        }
    }
    
    // Start recording
    tape.newRecording();
    
    // NOW set params/IC 
    if (params) {
        system.set_params(inputs.begin());
    }
    
    if (ic) {
        auto ic_start = params ? inputs.begin() + params->size() : inputs.begin();
        system.set_ode_state(ic_start, solver.fit_times()[0]);
    } else {
    // If IC not provided, use current solver state but ensure time matches fit_times[0]
    auto current_state = solver.state();
    std::vector<ad_type> state_ad(current_state.begin(), current_state.end());
    system.set_ode_state(state_ad.begin(), solver.fit_times()[0]);
}
    // Reset solver to sync internal state with system
    solver.reset();
    
    // Forward pass - returns states only at observation times
    auto obs = solver.advance_target();
    
    // Compute loss
    ad_type loss = sum_of_squares(obs, solver.targets());
    
    // Backpropagate
    tape.registerOutput(loss);
    xad::derivative(loss) = 1.0;
    tape.computeAdjoints();
    
    // Extract gradients
    std::vector<double> gradient(inputs.size());
    for (size_t i = 0; i < inputs.size(); ++i) {
        gradient[i] = xad::derivative(inputs[i]);
    }
    
    return {xad::value(loss), gradient};
}

} // namespace ode
} // namespace odelia

#endif
