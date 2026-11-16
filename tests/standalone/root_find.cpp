// The bracketed root-finders in root_find.hpp: what each returns, when each
// refuses, and that a search allocates nothing. That they visit the points
// Boost's do is tests/testthat/test-root-find.R's, which has Boost to hand.

#include <odelia/root_find.hpp>

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <string>

// Every allocation in this program is counted, so a test can require that a
// search makes none.
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

namespace root = odelia::root;

int failures = 0;

void check(bool ok, const std::string& what) {
  std::printf("  %s %s\n", ok ? "ok  " : "FAIL", what.c_str());
  if (!ok) {
    ++failures;
  }
}

// What a call threw, or "none".
template <class Call>
std::string thrown(Call call) {
  try {
    call();
  } catch (const odelia::util::RootFindError&) {
    return "RootFindError";
  } catch (const odelia::util::BracketError&) {
    return "BracketError";
  } catch (...) {
    return "other";
  }
  return "none";
}

// A target that counts its evaluations.
struct counted {
  double (*f)(double);
  int calls = 0;
  double operator()(double x) {
    ++calls;
    return f(x);
  }
};

double sqrt2_target(double x) { return x * x - 2.0; }
double kepler_target(double x) { return x - 0.5 * std::sin(x) - 1.0; }
double line_target(double x) { return x - 0.3; }
// A step: no method can do better than halve the bracket around it.
double step_target(double x) { return x < 0.3 ? -1.0 : 1.0; }

void test_each_form_meets_the_tolerance() {
  std::printf("each form meets the tolerance\n");
  const root::tolerance tol{1e-12, 1e-12};
  const double root2 = std::sqrt(2.0);

  const double x3 = root::toms748(sqrt2_target, 1.0, 2.0, tol);
  check(std::abs(x3 - root2) < 1e-12, "toms748 evaluating the ends finds sqrt(2)");

  const double x5 = root::toms748(sqrt2_target, 1.0, 2.0, -1.0, 2.0, tol);
  check(x5 == x3, "toms748 handed the ends' values finds the same midpoint");

  const auto br = root::toms748_bracket(sqrt2_target, 1.0, 2.0, -1.0, 2.0, tol);
  check(tol.met(br.first, br.second) && br.first <= root2 && root2 <= br.second,
        "its final bracket meets the tolerance and holds the root");
  check((br.first + br.second) / 2.0 == x5,
        "the midpoint form returns that bracket's midpoint");

  const double xb = root::bisect(sqrt2_target, 1.0, 2.0, tol);
  check(std::abs(xb - root2) < 1e-12, "bisect finds sqrt(2)");

  counted c748{kepler_target};
  counted cbis{kepler_target};
  root::toms748(c748, 0.0, 3.0, tol);
  root::bisect(cbis, 0.0, 3.0, tol);
  check(c748.calls * 3 < cbis.calls,
        "toms748 needs " + std::to_string(c748.calls) +
            " evaluations on a smooth target where bisect needs " +
            std::to_string(cbis.calls));
}

void test_an_exact_zero_at_an_end_is_the_root() {
  std::printf("an exact zero at an end is the root\n");
  const root::tolerance tol{1e-12, 1e-12};
  counted f{line_target};
  const auto at_a = root::toms748_bracket(f, 0.3, 1.0, 0.0, 0.7, tol);
  check(at_a.first == 0.3 && at_a.second == 0.3 && f.calls == 0,
        "toms748_bracket returns [a, a] for f(a) = 0 without evaluating f");
  const auto at_b = root::toms748_bracket(f, -1.0, 0.3, -1.3, 0.0, tol);
  check(at_b.first == 0.3 && at_b.second == 0.3,
        "and [b, b] for f(b) = 0");
  const auto at_a_none =
      root::toms748_bracket(f, 0.3, 1.0, 0.0, 0.7, {1e-12, 1e-12, 0});
  check(at_a_none.first == 0.3 && at_a_none.second == 0.3,
        "with no evaluations to spend, as with any");
  check(root::toms748(line_target, 0.3, 1.0, tol) == 0.3,
        "toms748 evaluating the ends returns that end");
  check(root::bisect(line_target, 0.3, 1.0, tol) == 0.3 &&
            root::bisect(line_target, -1.0, 0.3, tol) == 0.3,
        "bisect returns the end at either end");
  check(root::bisect(line_target, 0.3, 0.3, tol) == 0.3,
        "bisect accepts a zero-width bracket on the root, as Boost does");
}

void test_a_bracket_that_brackets_nothing_is_refused() {
  std::printf("a bracket that brackets nothing is refused\n");
  const root::tolerance tol{1e-12, 1e-12};
  check(thrown([&] { root::toms748(sqrt2_target, 2.0, 3.0, tol); }) ==
            "BracketError",
        "toms748 refuses f of one sign at both ends");
  check(thrown([&] {
          root::toms748_bracket(sqrt2_target, 2.0, 3.0, 2.0, 7.0, tol);
        }) == "BracketError",
        "toms748_bracket refuses the same");
  check(thrown([&] { root::bisect(sqrt2_target, 2.0, 3.0, tol); }) ==
            "BracketError",
        "bisect refuses the same");
  check(thrown([&] { root::toms748(sqrt2_target, 2.0, 1.0, tol); }) ==
            "BracketError",
        "toms748 refuses ends out of order");
  check(thrown([&] { root::bisect(sqrt2_target, 2.0, 1.0, tol); }) ==
            "BracketError",
        "bisect refuses ends out of order");
  const root::tolerance none{1e-12, 1e-12, 0};
  check(thrown([&] {
          root::toms748_bracket(sqrt2_target, 2.0, 3.0, 2.0, 7.0, none);
        }) == "BracketError",
        "and with no evaluations to spend, as with any");
  const root::tolerance two{1e-12, 1e-12, 2};
  check(thrown([&] { root::toms748(sqrt2_target, 1.0, 1.0, two); }) ==
            "BracketError",
        "toms748 evaluating the ends refuses a zero-width bracket at any cap");
  bool is_domain_error = false;
  try {
    root::toms748(sqrt2_target, 2.0, 3.0, tol);
  } catch (const std::domain_error&) {
    is_domain_error = true;
  } catch (...) {
  }
  check(is_domain_error, "the refusal is a std::domain_error, as Boost's is");
}

void test_a_search_that_cannot_converge_is_refused() {
  std::printf("a search that cannot converge within its cap is refused\n");
  const root::tolerance tight{1e-15, 1e-15, 20};
  check(thrown([&] { root::bisect(step_target, 0.0, 1.0, tight); }) ==
            "RootFindError",
        "bisect of a step needs ~50 halvings to 1e-15, and 20 are refused");
  check(thrown([&] { root::toms748(step_target, 0.0, 1.0, tight); }) ==
            "RootFindError",
        "toms748 on the same step is refused");
  check(thrown([&] {
          root::toms748_bracket(step_target, 0.0, 1.0, -1.0, 1.0, tight);
        }) == "RootFindError",
        "and so is its bracket form");
  const root::tolerance two{1e-12, 1e-12, 2};
  check(thrown([&] { root::toms748(sqrt2_target, 1.0, 2.0, two); }) ==
            "RootFindError",
        "a cap of 2 spends itself on the ends, and [a, b] is refused");
  const root::tolerance loose{1.5, 0.0, 2};
  check(root::toms748(sqrt2_target, 1.0, 2.0, loose) == 1.5,
        "unless [a, b] already meets the tolerance");
  bool is_runtime_error = false;
  try {
    root::bisect(step_target, 0.0, 1.0, tight);
  } catch (const std::runtime_error& e) {
    is_runtime_error = std::string(e.what()).find("still wider") !=
                       std::string::npos;
  }
  check(is_runtime_error, "the refusal is a std::runtime_error naming why");
}

// The cap counts evaluations, and a search that meets the tolerance with its
// last one has converged: whether it converged is read off the bracket, not
// off the count.
void test_converging_on_the_last_iteration_is_accepted() {
  std::printf("converging on the last allowed iteration is accepted\n");
  const double t = 1e-12;
  {
    counted f{line_target};
    root::bisect(f, 0.0, 1.0, {t, t});
    const std::size_t n = static_cast<std::size_t>(f.calls);
    // Bisection's cap counts its two end evaluations and one more.
    const root::tolerance exact{t, t, n + 1};
    const root::tolerance short_one{t, t, n};
    check(thrown([&] { root::bisect(line_target, 0.0, 1.0, exact); }) == "none",
          "bisect converging on its last iteration (" + std::to_string(n) +
              " evaluations) returns");
    check(thrown([&] { root::bisect(line_target, 0.0, 1.0, short_one); }) ==
              "RootFindError",
          "one fewer is refused");
  }
  {
    counted f{kepler_target};
    root::toms748_bracket(f, 0.0, 3.0, kepler_target(0.0), kepler_target(3.0),
                          {t, t});
    const std::size_t n = static_cast<std::size_t>(f.calls);
    const root::tolerance exact{t, t, n};
    const root::tolerance short_one{t, t, n - 1};
    check(thrown([&] {
            root::toms748_bracket(kepler_target, 0.0, 3.0, kepler_target(0.0),
                                  kepler_target(3.0), exact);
          }) == "none",
          "toms748 converging on its last iteration (" + std::to_string(n) +
              " evaluations) returns");
    check(thrown([&] {
            root::toms748_bracket(kepler_target, 0.0, 3.0, kepler_target(0.0),
                                  kepler_target(3.0), short_one);
          }) == "RootFindError",
          "one fewer is refused");
  }
}

void test_a_search_allocates_nothing() {
  std::printf("a search allocates nothing\n");
  const root::tolerance tol{1e-12, 1e-12};
  // A target holding state by reference, as a step's dense output would.
  double scale = 0.5;
  auto target = [&scale](double x) { return x - scale * std::sin(x) - 1.0; };
  const std::size_t before = allocations;
  double sum = 0.0;
  for (int i = 0; i < 100; ++i) {
    sum += root::toms748(target, 0.0, 3.0, target(0.0), target(3.0), tol);
    sum += root::toms748_bracket(target, 0.0, 3.0, target(0.0), target(3.0),
                                 tol).first;
    sum += root::toms748(target, 0.0, 3.0, tol);
    sum += root::bisect(target, 0.0, 3.0, tol);
  }
  check(allocations == before && std::isfinite(sum),
        "400 converging searches made " + std::to_string(allocations - before) +
            " allocations");
}

}  // namespace

int main() {
  std::printf("odelia bracketed root-finding, standalone\n");
  test_each_form_meets_the_tolerance();
  test_an_exact_zero_at_an_end_is_the_root();
  test_a_bracket_that_brackets_nothing_is_refused();
  test_a_search_that_cannot_converge_is_refused();
  test_converging_on_the_last_iteration_is_accepted();
  test_a_search_allocates_nothing();
  if (failures > 0) {
    std::printf("%d FAILED\n", failures);
    return 1;
  }
  std::printf("all passed\n");
  return 0;
}
