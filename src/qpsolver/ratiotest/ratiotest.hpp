#ifndef __SRC_LIB_RATIOTEST_HPP__
#define __SRC_LIB_RATIOTEST_HPP__

#include "../instance.hpp"
#include "../vector.hpp"

#include <limits>

struct RatiotestResult {
   double alpha;
   int limitingconstraint;
   bool nowactiveatlower;
};

RatiotestResult ratiotest_textbook(const Vector& x, const Vector& p, const Vector& rowact, const Vector& rowmove, Instance& instance, const double alphastart, const double t);

RatiotestResult ratiotest_twopass(const Vector& x, const Vector& p, const Vector& rowact, const Vector& rowmove, Instance& instance, Instance& relaxed, const double alphastart, const double t);

class Ratiotest {
public:
   virtual RatiotestResult ratiotest(const Vector& x, const Vector& p, const Vector& rowact, const Vector& rowmove, Instance& instance, const double alphastart) = 0;
};

class RatiotestTextbook : public Ratiotest {
public:
   RatiotestTextbook(double t_) : t(t_) {};

   RatiotestResult ratiotest(const Vector& x, const Vector& p, const Vector& rowact, const Vector& rowmove, Instance& instance, const double alphastart) {
      return ratiotest_textbook(x, p, rowact, rowmove, instance, alphastart, t);
   }

private:
   double t;
};

class RatiotestTwopass : public Ratiotest {
private:
   Instance relaxed_instance;
public:
   RatiotestTwopass(Instance& instance, double t_, double d) : t(t_) {
      relaxed_instance = instance;
      for (double& bound : relaxed_instance.con_lo) {
         if (bound != -std::numeric_limits<double>::infinity()) {
            bound -= d; 
         }
      }

      for (double& bound : relaxed_instance.con_up) {
         if (bound != std::numeric_limits<double>::infinity()) {
            bound += d; 
         }
      }

         for (double& bound : relaxed_instance.var_lo) {
         if (bound != -std::numeric_limits<double>::infinity()) {
            bound -= d; 
         }
      }

      for (double& bound : relaxed_instance.var_up) {
         if (bound != std::numeric_limits<double>::infinity()) {
            bound += d; 
         }
      }

   };

   RatiotestResult ratiotest(const Vector& x, const Vector& p, const Vector& rowact, const Vector& rowmove, Instance& instance, const double alphastart) {
      return ratiotest_twopass(x, p, rowact, rowmove, instance, relaxed_instance, alphastart, t);
   }

private:
   double t;
};

#endif
