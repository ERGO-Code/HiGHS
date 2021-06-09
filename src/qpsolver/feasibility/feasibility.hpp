#ifndef __SRC_LIB_FEASIBILITY_HPP__
#define __SRC_LIB_FEASIBILITY_HPP__

#include <cstdlib>

struct CrashSolution {
   std::vector<int> active;
   std::vector<int> inactive;
	std::vector<BasisStatus> rowstatus;
   Vector primal;
   Vector rowact;

   CrashSolution(int num_var, int num_row) : primal(Vector(num_var)), rowact(Vector(num_row)) {}
};

#endif
