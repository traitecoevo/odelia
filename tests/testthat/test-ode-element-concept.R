# Tests for the range helpers that walk a container of elements, and for the
# constraint that decides which iterators they accept.
# The interface header is header-only, so this needs the include path and no link
# against the odelia shared library.

compile_ode_element_interface <- function() {
  include_dir <- odelia_include_dir()
  withr::local_envvar(PKG_CPPFLAGS = odelia_cppflags(include_dir))
  Rcpp::sourceCpp(code = '
    // [[Rcpp::plugins(cpp20)]]
    #include <Rcpp.h>
    #include <vector>
    #include <odelia/ode_interface.hpp>

    // One state, one rate, one aux, moved through std::vector<value_type> iterators.
    struct Point {
      using value_type = double;
      double x = 0.0;

      std::vector<double>::const_iterator
      set_ode_state(std::vector<double>::const_iterator it) {
        x = *it++;
        return it;
      }
      std::vector<double>::iterator ode_state(std::vector<double>::iterator it) const {
        *it++ = x;
        return it;
      }
      std::vector<double>::iterator ode_rates(std::vector<double>::iterator it) const {
        *it++ = 2.0 * x;
        return it;
      }
      std::vector<double>::iterator ode_aux(std::vector<double>::iterator it) const {
        *it++ = x + 1.0;
        return it;
      }
      std::vector<double>::iterator set_ode_aux(std::vector<double>::iterator it) {
        return ++it;
      }
      size_t ode_size() const { return 1; }
      size_t aux_size() const { return 1; }
    };

    // Moves its state through int iterators rather than double ones, which is what
    // a double-typed signature does to an element carrying an active scalar.
    struct WrongPoint {
      using value_type = double;
      std::vector<int>::const_iterator
      set_ode_state(std::vector<int>::const_iterator it) { return it; }
      std::vector<int>::iterator ode_state(std::vector<int>::iterator it) const { return it; }
      std::vector<int>::iterator ode_rates(std::vector<int>::iterator it) const { return it; }
      std::vector<int>::iterator ode_aux(std::vector<int>::iterator it) const { return it; }
    };

    // What the helper accepts, asserted at the call rather than on the element.
    static_assert(requires(std::vector<Point> v, std::vector<double>::iterator it) {
      odelia::ode::ode_state(v.begin(), v.end(), it); });

    // Per member, not per element: WrongPoint moves its state through int iterators
    // and has no set_ode_aux at all, and each of those rules it out of exactly the
    // helpers that ask for what it lacks. The mismatches themselves are compile
    // failures, asserted as such below.
    static_assert(requires(std::vector<WrongPoint> v, std::vector<int>::iterator it) {
      odelia::ode::ode_rates(v.begin(), v.end(), it); });

    // [[Rcpp::export]]
    Rcpp::List element_range_round_trip(std::vector<double> state) {
      std::vector<Point> pts(state.size());
      const size_t n = odelia::ode::ode_size(pts.begin(), pts.end());
      odelia::ode::set_ode_state(pts.begin(), pts.end(), state.cbegin());

      std::vector<double> got(n), rates(n),
          aux(odelia::ode::aux_size(pts.begin(), pts.end()));
      auto state_end = odelia::ode::ode_state(pts.begin(), pts.end(), got.begin());
      odelia::ode::ode_rates(pts.begin(), pts.end(), rates.begin());
      auto aux_end = odelia::ode::ode_aux(pts.begin(), pts.end(), aux.begin());
      auto set_aux_end = odelia::ode::set_ode_aux(pts.begin(), pts.end(), aux.begin());

      return Rcpp::List::create(
          Rcpp::_["state"] = got,
          Rcpp::_["rates"] = rates,
          Rcpp::_["aux"] = aux,
          Rcpp::_["state_advance"] = (int)(state_end - got.begin()),
          Rcpp::_["aux_advance"] = (int)(aux_end - aux.begin()),
          Rcpp::_["set_aux_advance"] = (int)(set_aux_end - aux.begin()));
    }', verbose = FALSE)
}

testthat::test_that("the range helpers move each element's state through its own iterator", {
  compile_ode_element_interface()

  s <- c(1.0, 2.0, 3.0)
  r <- element_range_round_trip(s)

  expect_equal(r$state, s)
  expect_equal(r$rates, 2 * s)
  expect_equal(r$aux, s + 1)
})

testthat::test_that("every helper advances the iterator by the size it reports", {
  compile_ode_element_interface()

  s <- c(1.0, 2.0, 3.0, 4.0)
  r <- element_range_round_trip(s)

  expect_equal(r$state_advance, length(s))
  expect_equal(r$aux_advance, length(s))
  expect_equal(r$set_aux_advance, length(s))
})

# A program that calls one helper on one element with one iterator. Compiling it is
# the assertion: the constraint either admits the call or rejects it.
compile_helper_call <- function(element, member, iterator) {
  include_dir <- odelia_include_dir()
  withr::local_envvar(PKG_CPPFLAGS = odelia_cppflags(include_dir))
  Rcpp::sourceCpp(code = sprintf('
    // [[Rcpp::plugins(cpp20)]]
    #include <Rcpp.h>
    #include <vector>
    #include <odelia/ode_interface.hpp>

    struct Double {
      using value_type = double;
      double x = 0.0;
      std::vector<double>::const_iterator
      set_ode_state(std::vector<double>::const_iterator it) { x = *it++; return it; }
      std::vector<double>::iterator ode_state(std::vector<double>::iterator it) const {
        *it++ = x; return it; }
      std::vector<double>::iterator ode_rates(std::vector<double>::iterator it) const {
        *it++ = x; return it; }
      std::vector<double>::iterator ode_aux(std::vector<double>::iterator it) const {
        *it++ = x; return it; }
      size_t ode_size() const { return 1; }
      size_t aux_size() const { return 1; }
    };

    // Its state moves through int iterators, which is what a double-typed signature
    // does to an element carrying an active scalar.
    struct Int {
      using value_type = double;
      std::vector<int>::const_iterator
      set_ode_state(std::vector<int>::const_iterator it) { return it; }
      std::vector<int>::iterator ode_state(std::vector<int>::iterator it) const { return it; }
      std::vector<int>::iterator ode_rates(std::vector<int>::iterator it) const { return it; }
      std::vector<int>::iterator ode_aux(std::vector<int>::iterator it) const { return it; }
      size_t ode_size() const { return 1; }
      size_t aux_size() const { return 1; }
    };

    // [[Rcpp::export]]
    void helper_call() {
      std::vector<%s> v(1);
      std::vector<%s> buffer(1);
      odelia::ode::%s(v.begin(), v.end(), buffer.begin());
    }', element, iterator, member), verbose = FALSE)
}

testthat::test_that("a helper admits exactly the iterator its element moves state through", {
  # The element and the iterator agree, so the call compiles.
  expect_no_error(compile_helper_call("Double", "ode_state", "double"))

  # They do not, and the rejection is the constraint on the call: an element-level
  # requirement cannot see this, because Double is well formed on its own.
  expect_error(compile_helper_call("Double", "ode_state", "int"))
  expect_error(compile_helper_call("Int", "ode_state", "double"))

  # Per member: Int has no set_ode_aux, which rules it out of that helper only.
  expect_error(compile_helper_call("Int", "set_ode_aux", "int"))
  expect_no_error(compile_helper_call("Int", "ode_rates", "int"))
})
