#ifndef __SRC_LIB_FEASIBILITY_HPP__
#define __SRC_LIB_FEASIBILITY_HPP__

#include <cstdlib>

struct CrashSolution {
  std::vector<HighsInt> active;
  std::vector<HighsInt> inactive;
  std::vector<BasisStatus> rowstatus;
  Vector primal;
  Vector rowact;

  CrashSolution(HighsInt num_var, HighsInt num_row)
      : primal(Vector(num_var)), rowact(Vector(num_row)) {}
};

#endif
