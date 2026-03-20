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

    // Reuse tape across calls to avoid invalidating slots
    if (!solver.tape) {
        solver.tape = new ad::tape_type(false); // create without activating
    }
    solver.tape->activate();

    try {      
        // Collect input pointers for gradient extraction
        std::vector<ad_type*> inputs;

        // Configure the System
        auto& system = solver.get_system_ref();

        if (params) {
            auto refs = system.set_params(*solver.tape, params->begin());
            inputs.insert(inputs.end(), refs.begin(), refs.end());
        }

        if (ic) {
            auto refs = system.set_initial_state(*solver.tape, ic->begin(), solver.fit_times()[0]);
            inputs.insert(inputs.end(), refs.begin(), refs.end());
        }

        // Start recording operations on AD tape
        solver.tape->newRecording();

        // Reset system to propagate initial conditions into ODE buffers
        solver.reset();
        
        // Forward pass - returns states at observation times
        auto obs = solver.advance_target();
        
        // Compute loss
        ad_type loss = sum_of_squares(obs, solver.targets());
        
        // Backpropagate
        solver.tape->registerOutput(loss);
        xad::derivative(loss) = 1.0;
        solver.tape->computeAdjoints();
        
        // Extract gradients from registered inputs
        std::vector<double> gradient(inputs.size());
        for (size_t i = 0; i < inputs.size(); ++i) {
            gradient[i] = xad::derivative(*inputs[i]);
        }
        
        // Stop recording
        solver.tape->deactivate();
        return {xad::value(loss), gradient};
        
    } catch (const std::exception& e) {
        solver.tape->deactivate();
        throw;
    }
}

} // namespace ode
} // namespace odelia

#endif
