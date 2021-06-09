#ifndef __SRC_LIB_PRICING_DEVEXPRICING_HPP__
#define __SRC_LIB_PRICING_DEVEXPRICING_HPP__

#include "pricing.hpp"
#include "../runtime.hpp"
#include "../basis.hpp"
#include "../reducedcosts.hpp"

// 42726, 78965776.391299, 559, 104.321553, 0.000669, 7937

class DevexPricing : public Pricing {
private:
   Runtime& runtime;
   Basis& basis; 
   ReducedCosts& redcosts;

   std::vector<double> weights;

   int chooseconstrainttodrop(const Vector& lambda) {
      auto activeconstraintidx = basis.getactive();
      auto constraintindexinbasisfactor = basis.getindexinfactor(); 
      
      int minidx = -1;
      double maxabslambda = 0.0;
      for (int i = 0; i < activeconstraintidx.size(); i++) {
         int indexinbasis = constraintindexinbasisfactor[activeconstraintidx[i]];
         if (indexinbasis == -1) {
            printf("error\n");
         }
         assert(indexinbasis != -1);

         double val = lambda.value[indexinbasis] * lambda.value[indexinbasis] / weights[indexinbasis];
         if (val > maxabslambda && fabs(lambda.value[indexinbasis]) > runtime.settings.lambda_zero_threshold) {
            if (basis.getstatus(activeconstraintidx[i]) == BasisStatus::ActiveAtLower &&
               -lambda.value[indexinbasis] > 0) {
               minidx = activeconstraintidx[i];
               maxabslambda = val;
            } else if (basis.getstatus(activeconstraintidx[i]) == BasisStatus::ActiveAtUpper &&
                        lambda.value[indexinbasis] > 0) {
               minidx = activeconstraintidx[i];
               maxabslambda = val;
            } else {
               // TODO
            }
         }
         
      }

      return minidx;
   }

   

public:  
   DevexPricing(Runtime& rt, Basis& bas, ReducedCosts& rc) : runtime(rt), basis(bas), redcosts(rc), weights(std::vector<double>(rt.instance.num_var, 1.0)) {};
   
   // B lambda = g
   // lambda = inv(B)g
   // lambda = Z'g == reduced gradient ??
   // no: lambda = Y'g !!
   // dual values updated as:
   // c_N^T  += alpha_D * a_p^T (pivotal row)
   // alpha_D = -c_q / a_pq
   int price(const Vector& x, const Vector& gradient) {
      Vector& lambda = redcosts.getReducedCosts();
		int minidx = chooseconstrainttodrop(lambda);
      return minidx;
   }

   void update_weights(const Vector& aq, const Vector& ep, int p, int q) {
      int rowindex_p = basis.getindexinfactor()[p];
      double weight_p = weights[rowindex_p];
      for (int i=0; i<runtime.instance.num_var; i++) {
         if (i == rowindex_p) {
            weights[i] = weight_p / (aq.value[rowindex_p] * aq.value[rowindex_p]);
         } else {
            weights[i] += (aq.value[i] * aq.value[i]) / (aq.value[rowindex_p] * aq.value[rowindex_p]) * weight_p * weight_p;
         }
         if (weights[i] > 10E6) {
            weights[i] = 1.0;
         }
      }
   }

};

#endif
