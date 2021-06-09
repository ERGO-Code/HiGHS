#ifndef __SRC_LIB_REDUCEDGRADIENT_HPP__
#define __SRC_LIB_REDUCEDGRADIENT_HPP__

#include "vector.hpp"
#include "runtime.hpp"

#include "nullspace.hpp"

class ReducedGradient {
   Vector rg;
   bool uptodate = false;
   Gradient& gradient;
   Nullspace& nullspace;

   void recompute() {
      rg.dim = nullspace.getNullspace().mat.num_col;
      nullspace.Ztprod(gradient.getGradient(), rg);
      uptodate = true;
   }

public:
   ReducedGradient(Runtime& rt, Nullspace& ns, Gradient& grad) : rg(rt.instance.num_var), gradient(grad) , nullspace(ns){
      
   }

   Vector& get() {
      if (!uptodate) {
         recompute();
      }
      return rg;
   }

   void reduce(NullspaceReductionResult& nrr) {
      if (!uptodate) {
         return;
      }
      // Vector r(rg.dim-1); 
      // for (int col=0; col<nrr.maxabsd; col++) {
      //    r.index[col] = col;
      //    r.value[col] = -nrr.d[col] / nrr.d[nrr.maxabsd];
      // }
      // for (int col=nrr.maxabsd+1; col<rg.dim; col++) {
      //    r.index[col-1] = col-1;
      //    r.value[col-1] = -nrr.d[col] / nrr.d[nrr.maxabsd];
      // }
      // r.num_nz = rg.dim-1;

      for (int i=0; i<nrr.d.num_nz; i++) {
         int idx = nrr.d.index[i];
         if (idx== nrr.maxabsd) {
            continue;
         }
         rg.value[idx] -= rg.value[nrr.maxabsd] * nrr.d.value[idx] / nrr.d.value[nrr.maxabsd];
      }

      rg.resparsify();

      uptodate = true;
   }

   void expand(const Vector& yp) {
      if (!uptodate) {
         return;
      }

      double newval = yp * gradient.getGradient();
      rg.value.push_back(newval);
      rg.index.push_back(0);
      rg.index[rg.num_nz++] = rg.dim++;

      uptodate = true;
   }

   void update(double alpha, bool minor) {
      if (!uptodate) {
         return;
      }
      if (minor) {
         for (int i=0; i<rg.num_nz; i++) {
            rg.value[rg.index[i]] *= (1.0-alpha);
         }
         uptodate = true;
      } else {
         uptodate = false;
      }
   }

};

#endif
