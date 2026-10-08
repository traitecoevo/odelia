# Build and run the R-free core check as plain C++ (#43).
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

.PHONY: all build selfcontained clean

all: build selfcontained
	./r_free

build: r_free

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

r_free: r_free.o Tape.o
	$(CXX) $(CXXSTD) -o $@ $^

r_free.o: r_free.cpp $(HEADERS)
	$(CXX) $(CXXSTD) $(OPT) $(WARN) $(INCLUDES) -c -o $@ $<

# Vendored third-party code: build it quietly, we do not police its warnings.
Tape.o: ../../src/Tape.cpp
	$(CXX) $(CXXSTD) $(OPT) -w $(INCLUDES) -c -o $@ $<

clean:
	rm -f r_free r_free.o Tape.o
