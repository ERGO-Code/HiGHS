#include "ratiotest.hpp"

double step(double x, double p, double l, double u, double t) {
   if (p < -t && l > -std::numeric_limits<double>::infinity()) {
      return (l - x) / p;
   } else if (p > t && u < std::numeric_limits<double>::infinity()) {
      return (u - x) / p;
   } else {
      return std::numeric_limits<double>::infinity();
   }
}

RatiotestResult ratiotest_textbook(const Vector& x, const Vector& p, const Vector& rowact, const Vector& rowmove, Instance& instance, const double alphastart, const double t) {
   RatiotestResult result;
   result.limitingconstraHighsInt = -1;
   result.alpha = alphastart;

   // check ratio towards variable bounds
   for (HighsInt j=0; j<p.num_nz; j++) {
      HighsInt i = p.index[j];
      double alpha_i = step(x.value[i], p.value[i], instance.var_lo[i], instance.var_up[i], t);
      if (alpha_i < result.alpha) {
         result.alpha = alpha_i;
         result.limitingconstraHighsInt = instance.num_con + i;
         result.nowactiveatlower = p.value[i] < 0;
      }
   }

   // check ratio towards constraHighsInt bounds
   for (HighsInt j=0; j<rowmove.num_nz; j++) {
      HighsInt i = rowmove.index[j];
      double alpha_i = step(rowact.value[i], rowmove.value[i], instance.con_lo[i], instance.con_up[i], t);
      if (alpha_i < result.alpha) {
         result.alpha = alpha_i;
         result.limitingconstraHighsInt = i;
         result.nowactiveatlower = rowmove.value[i] < 0;
      }
   }

   return result;
}

RatiotestResult ratiotest_twopass(const Vector& x, const Vector& p, const Vector& rowact, const Vector& rowmove, Instance& instance, Instance& relaxed, const double alphastart, const double t) {

   RatiotestResult res1 = ratiotest_textbook(x, p, rowact, rowmove, relaxed, alphastart, t);

   if (res1.limitingconstraHighsInt == -1) {
      return res1;
   }

   RatiotestResult result = res1;

   double max_pivot = 0;
   if (res1.limitingconstraHighsInt != -1) {
      if ((int)result.limitingconstraHighsInt < instance.num_con) {
         max_pivot = rowmove.value[result.limitingconstraint];
      } else {
         max_pivot = p.value[result.limitingconstraHighsInt - instance.num_con];
      }
   }

   for (HighsInt i=0; i<instance.num_con; i++) {
      double step_i = step(rowact.value[i], rowmove.value[i], instance.con_lo[i], instance.con_up[i], t);
      if (fabs(rowmove.value[i]) >= fabs(max_pivot) && step_i <= res1.alpha) {
         max_pivot = rowmove.value[i];
         result.limitingconstraHighsInt = i;
         result.alpha = step_i;
         result.nowactiveatlower = rowmove.value[i] < 0;
      }
   }

   for (HighsInt i=0; i<instance.num_var; i++) {
      double step_i = step(x.value[i], p.value[i], instance.var_lo[i], instance.var_up[i], t);
      if (fabs(p.value[i]) >= fabs(max_pivot) && step_i <= res1.alpha) {
         max_pivot = p.value[i];
         result.limitingconstraHighsInt = instance.num_con + i;
         result.alpha = step_i;
         result.nowactiveatlower = p.value[i] < 0;
      }
   }

   result.alpha = fmax(result.alpha, 0.0);
   return result;
}
