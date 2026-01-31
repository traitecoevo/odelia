#ifndef ODELIA_ODE_FIT_HPP_
#define ODELIA_ODE_FIT_HPP_

#include <XAD/XAD.hpp>
#include <odelia/ode_solver.hpp>
#include <optional>
#include <vector>

namespace odelia {
namespace ode {

// Loss function
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
    
    // Collect input pointers for gradient extraction
    std::vector<ad_type*> inputs;
    
    auto& system = solver.get_system_ref();
    
    // Set params and register
    if (params) {
        auto refs = system.set_params(tape, params->begin());
        inputs.insert(inputs.end(), refs.begin(), refs.end());
    }
    
    // Set initial conditions and register
    if (ic) {
        auto refs = system.set_ode_state(tape, ic->begin(), solver.fit_times()[0]);
        inputs.insert(inputs.end(), refs.begin(), refs.end());
    } else {
        system.compute_rates();
    }
    
    solver.reset();
    
    tape.newRecording();
    
    // Forward pass - returns states at observation times
    auto obs = solver.advance_target();
    
    // Compute loss
    ad_type loss = sum_of_squares(obs, solver.targets());
    
    // Backpropagate
    tape.registerOutput(loss);
    xad::derivative(loss) = 1.0;
    tape.computeAdjoints();
    
    // Extract gradients from registered inputs
    std::vector<double> gradient(inputs.size());
    for (size_t i = 0; i < inputs.size(); ++i) {
        gradient[i] = xad::derivative(*inputs[i]);
    }
    
    return {xad::value(loss), gradient};
}

} // namespace ode
} // namespace odelia

#endif
