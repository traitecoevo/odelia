#ifndef ODELIA_SOLVER_WRAPPER_HPP_
#define ODELIA_SOLVER_WRAPPER_HPP_

#include <odelia/ode_solver.hpp>
#include <XAD/XAD.hpp>
#include <memory>
#include <vector>

namespace odelia {
namespace ode {

// Abstract base class for type-erased solver
class SolverBase {
public:
    virtual ~SolverBase() {}
    
    // Virtual methods for all Solver operations
    virtual void reset() = 0;
    virtual double time() const = 0;
    virtual std::vector<double> state() const = 0;
    virtual std::vector<double> times() const = 0;
    virtual void set_state(const std::vector<double>& y, double time) = 0;
    virtual void advance_adaptive(const std::vector<double>& times) = 0;
    virtual void advance_fixed(const std::vector<double>& times) = 0;
    virtual void step() = 0;
    virtual bool get_collect() const = 0;
    virtual void set_collect(bool x) = 0;
    virtual std::size_t get_history_size() const = 0;
    virtual std::vector<double> get_history_step(std::size_t i) const = 0;
    
    // Fit-related methods
    virtual void set_target(const std::vector<double>& times,
                           const std::vector<std::vector<double>>& targets,
                           const std::vector<size_t>& obs_indices) = 0;
    virtual const std::vector<double>& fit_times() const = 0;
    virtual const std::vector<std::vector<double>>& targets() const = 0;
    virtual const std::vector<size_t>& obs_indices() const = 0;
    
    // For getting system parameters
    virtual std::vector<double> get_system_pars() const = 0;
    
    // Check if this is an active (AD-enabled) solver
    virtual bool is_active() const = 0;
};
// Helper to convert AD state to double
template <typename T>
inline std::vector<double> to_double_vector(const std::vector<T>& v) {
    std::vector<double> result(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        result[i] = xad::value(v[i]);
    }
    return result;
}

// Specialization for double (no conversion needed)
template <>
inline std::vector<double> to_double_vector<double>(const std::vector<double>& v) {
    return v;
}
// Concrete implementation that wraps a Solver<System>
template <typename System>
class SolverWrapper : public SolverBase {
private:
    std::unique_ptr<Solver<System>> solver_;
    ode::OdeControl control_;
    
public:
SolverWrapper(System* sys, const ode::OdeControl& ctrl) 
    : control_(ctrl), 
      solver_(new Solver<System>(*sys, control_)) {
    delete sys; // Transfer ownership to Solver
}
    
    void reset() override { solver_->reset(); }
    double time() const override { return solver_->time(); }
    std::vector<double> state() const override { 
    return to_double_vector(solver_->state()); 
}
    std::vector<double> times() const override { return solver_->times(); }
    
    void set_state(const std::vector<double>& y, double time) override {
        solver_->set_state(y, time);
    }
    
    void advance_adaptive(const std::vector<double>& times) override {
        solver_->advance_adaptive(times);
    }
    
    void advance_fixed(const std::vector<double>& times) override {
        solver_->advance_fixed(times);
    }
    
    void step() override { solver_->step(); }
    bool get_collect() const override { return solver_->get_collect(); }
    void set_collect(bool x) override { solver_->set_collect(x); }
    std::size_t get_history_size() const override { return solver_->get_history_size(); }
    
    std::vector<double> get_history_step(std::size_t i) const override {
        auto sys = solver_->get_history_step(i);
        return sys.record_step();
    }
    
    void set_target(const std::vector<double>& times,
                   const std::vector<std::vector<double>>& targets,
                   const std::vector<size_t>& obs_indices) override {
        solver_->set_target(times, targets, obs_indices);
    }
    
    const std::vector<double>& fit_times() const override {
        return solver_->fit_times();
    }
    
    const std::vector<std::vector<double>>& targets() const override {
        return solver_->targets();
    }
    
    const std::vector<size_t>& obs_indices() const override {
        return solver_->obs_indices();
    }
    
    std::vector<double> get_system_pars() const override {
        return solver_->get_system().pars();
    }
    
    // Access the underlying solver (for fit functions)
    Solver<System>* get_solver() { return solver_.get(); }
    
    bool is_active() const override { return false; }

    // Cast to specific wrapper type (for runtime dispatch)
    template<typename S>
    SolverWrapper<S>* as() { 
        return dynamic_cast<SolverWrapper<S>*>(this); 
    }
};

} // namespace ode
} // namespace odelia

#endif
