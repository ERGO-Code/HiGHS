#ifndef __SRC_LIB_REDUCEDCOSTS_HPP__
#define __SRC_LIB_REDUCEDCOSTS_HPP__

#include "qpsolver/basis.hpp"
#include "qpsolver/gradient.hpp"
#include "qpsolver/runtime.hpp"
#include "qpsolver/vector.hpp"

class ReducedCosts {
  Basis& basis;

  Gradient& gradient;

  Vector reducedcosts;
  bool uptodate;

 public:
  ReducedCosts(Runtime& rt, Basis& bas, Gradient& grad)
      : basis(bas),
        gradient(grad),
        reducedcosts(Vector(rt.instance.num_var)),
        uptodate(false) {}

  void recompute() {
    basis.ftran(gradient.getGradient(), reducedcosts);
    uptodate = true;
  }

  Vector& getReducedCosts() {
    if (!uptodate) {
      recompute();
    }
    return reducedcosts;
  }

  void update() { uptodate = false; }
};

#endif
