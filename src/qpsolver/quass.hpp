#ifndef __SRC_LIB_QUASS_HPP__
#define __SRC_LIB_QUASS_HPP__

#include "qpsolver/basis.hpp"
#include "qpsolver/eventhandler.hpp"
#include "qpsolver/factor.hpp"
#include "qpsolver/instance.hpp"
#include "qpsolver/runtime.hpp"

struct Quass {
  Quass(Runtime& rt);

  void solve(const Vector& x0, const Vector& ra, Basis& b0);

  void solve();

 private:
  Runtime& runtime;

  void loginformation(Runtime& rt, Basis& basis, CholeskyFactor& factor);
};

#endif
