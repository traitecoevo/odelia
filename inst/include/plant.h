// -*-c++-*-
#ifndef _PLANT_H_
#define _PLANT_H_

#include <plant/ode_control.h>
#include <plant/ode_step.h>
#include <plant/ode_solver.h>
#include <plant/ode_runner.h>

// Purely for testing
#include <plant/lorenz.h>

// Include this early on.  It can be either after classes have been
// declared (but before Rcpp has been loaded) or first.  This file will
// attempt to provide declarations for the classes and namespaces that
// you use, but this might be fragile.
#include <plant/RcppR6_pre.hpp>

// Anything after this point is OK to include Rcpp.h.  This is
// probably where the meat of the included material goes if your
// classes directly use Rcpp types.  Otherwise you can just declare
// them earlier up.

#include <Rcpp.h>
#include <plant/ode_r.h>

// This line can safely be the last line in the file, but may go any
// point after RcppR6_pre.hpp is included.
#include <plant/RcppR6_post.hpp>
#include <plant/util_post_rcpp.h>

#endif
