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

// sweep_adjoint_at over a run's first `n` output times, one adjoint at each.
std::size_t sweep_at_allocations(Solver<LorenzSystem<double>>& s,
                                 const std::vector<std::size_t>& rows,
                                 std::size_t n) {
  const std::vector<std::size_t> first(rows.begin(), rows.begin() + n);
  const std::vector<adjoint_batch> at(n, adjoint_batch::one({1.0, 0.5, 0.0}));
  adjoint_batch lambda(1, 3);
  adjoint_batch dp(1, 3);
  const std::size_t before = allocations;
  s.sweep_adjoint_at(first, at, lambda, dp);
  return allocations - before;
}

// sweep_adjoint_at sweeps one interval per output time, and a time series has
// hundreds. Its intervals share one tape and one active System, as the segments
// of one sweep do (#103): a tape per interval reserved and freed a large block
// each time, which under CmdStan's tbbmalloc proxy cost 17 times the sweep.
//
// The tape's block goes around operator new (see above), so this counts what is
// held beside it: the active System and its buffers, which an interval that built
// its own workspace would build again. A sweep allocates nothing per step, so a
// whole sweep's count is its fixed cost; an interval must cost less than that.
void test_intervals_share_one_workspace(Method m) {
  std::printf("sweep_adjoint_at's intervals share a tape and an active System "
              "under %s\n", m == Method::dopri ? "dopri" : "rkck");
  LorenzSystem<double> system(10.0, 28.0, 8.0 / 3.0);
  Solver<LorenzSystem<double>> s(system, OdeControl(), m);
  s.set_collect(false);
  s.set_keep_states(true);
  std::vector<double> times;
  for (int i = 0; i <= 400; ++i) {
    times.push_back(0.025 * i);
  }
  s.advance_adaptive(times);
  const std::vector<std::size_t> rows =
      s.rows_at(std::vector<double>(times.begin() + 1, times.end()));

  const std::size_t setup = sweep_allocations(s, 0, s.recording().size() - 1);
  const std::size_t many = sweep_at_allocations(s, rows, 400);
  const std::size_t few = sweep_at_allocations(s, rows, 200);
  const double per_interval = static_cast<double>(many - few) / 200.0;
  std::printf("       (a sweep's setup: %zu allocations; per interval: %.1f)\n",
              setup, per_interval);
  check(per_interval < static_cast<double>(setup),
        "an interval allocates less than a sweep's setup: the intervals do not "
        "each build a workspace");
}

}  // namespace

int main() {
  std::printf("odelia reverse sweep, allocations per step\n");
  test_a_sweep_allocates_nothing_per_step(Method::rkck);
  test_a_sweep_allocates_nothing_per_step(Method::dopri);
  test_intervals_share_one_workspace(Method::rkck);
  test_intervals_share_one_workspace(Method::dopri);
  if (failures > 0) {
    std::printf("%d check(s) failed\n", failures);
    return 1;
  }
  std::printf("all passed\n");
  return 0;
}
