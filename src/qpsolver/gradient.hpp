#ifndef __SRC_LIB_GRADIENT_HPP__
#define __SRC_LIB_GRADIENT_HPP__

#include "vector.hpp"
#include "runtime.hpp"

class Gradient {
   Runtime& runtime;

   Vector gradient;
   bool uptodate;
   int numupdates = 0;

   Vector buffer_temp;

   void recompute() {
      runtime.instance.Q.vec_mat(runtime.primal, gradient);
      gradient += runtime.instance.c;
      uptodate = true;
      numupdates = 0;
   }

public:
   Gradient(Runtime& rt) : runtime(rt),  gradient(Vector(rt.instance.num_var)), uptodate(false), buffer_temp(rt.instance.num_var) {
   
   }

   Vector& getGradient() {
      if (!uptodate || numupdates >= runtime.settings.gradientrecomputefrequency) {
         recompute();
      }
      return gradient;
   }

   void update(Vector& p, double stepsize) {
      runtime.instance.Q.mat_vec(p, buffer_temp);
      gradient.saxpy(stepsize, buffer_temp);
      numupdates++;
   }
};

#endif
