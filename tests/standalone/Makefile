# Build and run the R-free core check as plain C++ (#43), and the check that
# every copy of a stage time rounds it the same way (stage_times.cpp).
#
#   make          build and run
#   make build    build only
#   make clean
#
# The include path is deliberately just ../../inst/include -- no R, no Rcpp, no
# BH. If this needs any of those added back in order to build, the solver core
# has stopped being usable outside an R package and something needs fixing.
#
# src/Tape.cpp is compiled in because the XAD tape runtime is defined in exactly
# one object file by design (ARCHITECTURE.md); anything that instantiates Solver
# has to link it. It is R-free too.

# `?=`, so a CXX exported into the environment wins. That is deliberate, and it is
# also how these probes stop running: an environment naming a compiler that is not
# installed fails with `No such file or directory` from make rather than from the
# compiler, which reads as a broken Makefile. `make CXX=g++ <probe>` overrides it,
# and g++ is what the packages themselves build with.
CXX     ?= c++
# The core uses concepts, so this must stay at C++20 or later: below it every
# `concept` in the headers reads as `'concept' does not name a type` and the
# hundred errors after that are cascade.
CXXSTD  ?= -std=c++20
OPT     ?= -O0
WARN    ?= -Wall -Wextra -Wno-unused-parameter

INCLUDES = -I../../inst/include
HEADERS  = $(wildcard ../../inst/include/odelia/*.hpp)

.PHONY: all build selfcontained layering clean

all: build selfcontained layering
	./r_free
	./stage_times

build: r_free stage_times

# Every core header must also compile on its own -- no leaning on whichever
# sibling header happened to be included first to have pulled in <vector> or
# <cassert>. The two R interface headers are exempt; they are meant to need Rcpp.
selfcontained:
	@for h in ../../inst/include/odelia/*.hpp; do \
	  b=$$(basename $$h); \
	  case $$b in solver_interface.hpp|rcpp_interface_helpers.hpp) continue;; esac; \
	  printf '#include <odelia/%s>\n' $$b \
	    | $(CXX) $(CXXSTD) $(INCLUDES) -fsyntax-only -x c++ - \
	    || { echo "not self-contained: $$b"; exit 1; }; \
	  echo "  self-contained: $$b"; \
	done

# implicit_node.hpp serves a consumer that records partials on a tape and takes
# no solver (phylloptim), so an edit to the System contract must not reach it.
# Compiling on its own does not show that -- it would with the contract included
# too -- so this checks the contract's include guard is still unset after it.
layering:
	@printf '#include <odelia/implicit_node.hpp>\n#ifdef ODELIA_ODE_INTERFACE_HPP_\n#error implicit_node.hpp reaches ode_interface.hpp\n#endif\n' \
	  | $(CXX) $(CXXSTD) $(INCLUDES) -fsyntax-only -x c++ - \
	  || { echo "implicit_node.hpp includes the System contract"; exit 1; }
	@echo "  implicit_node.hpp: no System contract"

r_free: r_free.o Tape.o
	$(CXX) $(CXXSTD) -o $@ $^

r_free.o: r_free.cpp $(HEADERS)
	$(CXX) $(CXXSTD) $(OPT) $(WARN) $(INCLUDES) -c -o $@ $<

# Optimised, and under the compiler's own contraction default: a fused
# multiply-add the optimiser chooses in one copy of a stage time and not in
# another is what this checks the steppers leave it no room for.
stage_times: stage_times.o Tape.o
	$(CXX) $(CXXSTD) -o $@ $^

stage_times.o: stage_times.cpp $(HEADERS)
	$(CXX) $(CXXSTD) -O2 $(WARN) $(INCLUDES) -c -o $@ $<

# Vendored third-party code: build it quietly, we do not police its warnings.
Tape.o: ../../src/Tape.cpp
	$(CXX) $(CXXSTD) $(OPT) -w $(INCLUDES) -c -o $@ $<

clean:
	rm -f r_free r_free.o stage_times stage_times.o Tape.o
