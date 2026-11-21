// A reverse sweep allocates nothing per step after its first: its buffers are
// sized once per segment and reused. Under an allocator that returns every
// freed block to the operating system, a sweep that allocated per step would
// cost more in the allocator than in the arithmetic.

#include <odelia/ode_solver.hpp>
#include <odelia/adjoint.hpp>
#include <examples/lorenz_system.hpp>

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <string>
#include <vector>

// Every allocation through operator new is counted. XAD's tape storage goes
// around it (AllocHelper::aligned_alloc), and keeps its chunks across a clear,
// so a tape that regrew every step would not show here.
static std::size_t allocations = 0;

void* operator new(std::size_t n) {
  ++allocations;
  if (void* p = std::malloc(n == 0 ? 1 : n)) {
    return p;
  }
  throw std::bad_alloc();
}
void operator delete(void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }

namespace {

using namespace odelia::ode;

int failures = 0;

void check(bool ok, const std::string& what) {
  std::printf("  %s %s\n", ok ? "ok  " : "FAIL", what.c_str());
  if (!ok) {
    ++failures;
  }
}

// Allocations a sweep of rows [first, last] makes, for one output.
std::size_t sweep_allocations(Solver<LorenzSystem<double>>& s, std::size_t first,
                              std::size_t last) {
  adjoint_batch lambda = adjoint_batch::one({1.0, 0.0, 0.0});
  adjoint_batch dp(1, 3);
  const std::size_t before = allocations;
  s.sweep_adjoint(lambda, dp, first, last);
  return allocations - before;
}

void test_a_sweep_allocates_nothing_per_step(Method m) {
  std::printf("a sweep allocates nothing per step under %s\n",
              m == Method::dopri ? "dopri" : "rkck");
  LorenzSystem<double> system(10.0, 28.0, 8.0 / 3.0);
  Solver<LorenzSystem<double>> s(system, OdeControl(), m);
  s.set_collect(false);
  s.set_keep_states(true);
  s.advance_adaptive(std::vector<double>{0.0, 10.0});
  const std::size_t last = s.recording().size() - 1;
  check(last > 200, "the run took more than 200 steps");

  // Each sweep builds its own tape and active System, a fixed cost; comparing
  // two lengths cancels it and leaves what a step adds.
  const std::size_t whole = sweep_allocations(s, 0, last);
  const std::size_t half = sweep_allocations(s, last / 2, last);
  const std::size_t steps_apart = last - last / 2;
  std::printf("       (%zu steps: %zu allocations; %zu steps: %zu)\n", last,
              whole, steps_apart, half);
  check(whole == half,
        "a sweep over " + std::to_string(last) + " steps allocates as much as "
        "one over " + std::to_string(steps_apart) + ": nothing per step");
}

}  // namespace

int main() {
  std::printf("odelia reverse sweep, allocations per step\n");
  test_a_sweep_allocates_nothing_per_step(Method::rkck);
  test_a_sweep_allocates_nothing_per_step(Method::dopri);
  if (failures > 0) {
    std::printf("%d check(s) failed\n", failures);
    return 1;
  }
  std::printf("all passed\n");
  return 0;
}
