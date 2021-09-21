#ifndef __SRC_LIB_SOLVER_HPP__
#define __SRC_LIB_SOLVER_HPP__

#include "basis.hpp"
#include "eventhandler.hpp"
#include "factor.hpp"
#include "instance.hpp"
#include "runtime.hpp"

struct Solver {
  Solver(Runtime& rt);

  void solve(const Vector& x0, const Vector& ra, Basis& b0);

  void solve();

 private:
  Runtime& runtime;

  void loginformation(Runtime& rt, Basis& basis, NewCholeskyFactor& factor);
};

#endif
