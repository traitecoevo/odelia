// Must NOT compile. sweep_tangent refuses, by name, a System that solves for
// inner values: a forward pass would solve them again rather than load the
// run's (TangentSweepable). The Makefile's `refusals` target checks the message.
#include <odelia/ode_solver.hpp>
#include <cstddef>
#include <vector>

struct inner { double x = 0.0; };

template <typename T = double>
struct Solving {
  using value_type = T;
  using solved_values = inner;
  T y = T(0.5), dydt = T(0.0);
  double time = 0.0;
  template <class S> Solving<S> rebind() const { return Solving<S>(); }
  std::size_t ode_size() const { return 1; }
  template <class It> It set_ode_state(It it, double t) {
    y = *it++; time = t; dydt = -y; return it;
  }
  void set_recorded_state(const std::vector<T>& s, double t) {
    set_ode_state(s.begin(), t);
  }
  template <class It> It ode_state(It it) const { *it++ = y; return it; }
  template <class It> It ode_rates(It it) const { *it++ = dydt; return it; }
  void store_solved(inner&, std::size_t) {}
  void load_solved(const inner&, std::size_t) {}
  void end_solved() {}
};

int main() {
  odelia::ode::Solver<Solving<>> s(Solving<>{}, odelia::ode::OdeControl());
  s.set_keep_states(true);
  s.sweep_tangent<1>(odelia::ode::tangent_directions{{}, {0}});
}
