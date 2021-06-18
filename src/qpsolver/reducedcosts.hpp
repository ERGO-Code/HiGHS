#ifndef __SRC_LIB_REDUCEDCOSTS_HPP__
#define __SRC_LIB_REDUCEDCOSTS_HPP__

#include "basis.hpp"
#include "vector.hpp"
#include "runtime.hpp"
#include "gradient.hpp"

class ReducedCosts {
   Runtime& runtime;
   Basis& basis;

   Gradient& gradient;

   Vector reducedcosts;
   bool uptodate;

   void recompute() {
      basis.ftran(gradient.getGradient(), reducedcosts);
      uptodate = true;
   }

public:
   ReducedCosts(Runtime& rt, Basis& bas, Gradient& grad) : runtime(rt), basis(bas), gradient(grad), reducedcosts(Vector(rt.instance.num_var)), uptodate(false) {
   
   }

   Vector& getReducedCosts() {
      if (!uptodate) {
         recompute();
      }
      return reducedcosts;
   }

   void update() {
      uptodate = false;
   }
};

#endif
