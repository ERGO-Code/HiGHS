#ifndef __SRC_LIB_STATISTICS_HPP__
#define __SRC_LIB_STATISTICS_HPP__

#include <chrono>
#include <vector>

struct Statistics {
   int phase1_iterations = 0;
   int num_iterations = 0;
   std::chrono::high_resolution_clock::time_point time_start;
   std::chrono::high_resolution_clock::time_point time_end; 

   std::vector<int> iteration;
   std::vector<int> nullspacedimension;
   std::vector<double> objval;
   std::vector<double> time; 
   std::vector<double> sum_primal_infeasibilities;
   std::vector<int> num_primal_infeasibilities;
   std::vector<double> density_nullspace;
   std::vector<double> density_factor;
};

#endif
