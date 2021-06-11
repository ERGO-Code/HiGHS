#ifndef __SRC_LIB_SOLVER_HPP__
#define __SRC_LIB_SOLVER_HPP__

#include "basis.hpp"
#include "eventhandler.hpp"
#include "runtime.hpp"

#include "redhes/factor.hpp"
#include "nullspace.hpp"

#include "instance.hpp"

struct Solver {
   Solver(Runtime& rt);
   
   void solve(const Vector& x0, const Vector& ra, Basis& b0);

   void solve();

private:
   Runtime& runtime;

   void loginformation(Runtime& rt, Basis& basis, Nullspace& ns, NewCholeskyFactor& factor);
};

#endif