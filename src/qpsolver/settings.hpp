#ifndef __SRC_LIB_SETTINGS_HPP__
#define __SRC_LIB_SETTINGS_HPP__

#include "ratiotest/ratiotest.hpp"

enum class OutputLevel {
   LIGHT,
   MEDIUM,
   HEAVY
};

struct Settings {
   bool simplexsteps = false;
   Ratiotest* ratiotest;
   double pnorm_zero_threshold = 10E-12;
   double d_zero_threshold = 10E-13;
   double lambda_zero_threshold = 10E-10;
   OutputLevel outputlevel = OutputLevel::LIGHT;
   int reportingfequency = 100;

   int reinvertfrequency = 100;
   int gradientrecomputefrequency = 1;
   int reducedgradientrecomputefrequency = std::numeric_limits<int>::infinity();
   int reducedhessianrecomputefrequency = std::numeric_limits<int>::infinity();

   int iterationlimit = std::numeric_limits<int>::infinity();
   double timelimit = std::numeric_limits<double>::infinity();

   bool rowscaling = true;
   bool varscaling = true;

   // TODO output settings
};

#endif
