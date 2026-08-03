#ifndef HIPO_OPTION_H
#define HIPO_OPTION_H

#include "Parameters.h"
#include "io/HighsIO.h"
#include "lp_data/HighsOptions.h"

namespace hipo {

enum ParallelType { kOff, kChoose, kOn };

struct Options {
  // Solver options
  std::string nla = kHighsChooseString;
  std::string crossover = kHighsOffString;
  std::string ordering = kHighsChooseString;
  std::string factor = kHighsChooseString;

  // Ipm parameters
  Int max_iter = kMaxIterDefault;
  double feasibility_tol = kIpmTolDefault;
  double optimality_tol = kIpmTolDefault;
  double crossover_tol = kIpmTolDefault;
  bool refine_with_ipx = true;
  double time_limit = -1.0;
  Int block_size = 0;
  Int random_seed = 0;
  std::vector<Int> parallel_type;

  // Logging
  bool display = true;
  bool timeless_log = false;
  const HighsLogOptions* log_options = nullptr;

  inline bool chooseParallel(Int bit, bool choose) const {
    if (parallel_type[bit] == kOn) return true;
    if (parallel_type[bit] == kOff) return false;
    return choose;
  }
};

enum ParallelTechnique {
  kParallelAnalyse,
  kParallelOrderNE,
  kParallelOrderAS,
  kParallelNEStruct,
  kParallelNEValues,
  kParallelTree,
  kParallelNode,
  kParallelForwardSolve,
  kParallelDiagonalSolve,
  kParallelBackwardSolve,
  kParallelCount,
};

inline bool testParallelBit(Int option, Int bit) { return option & (1 << bit); }

}  // namespace hipo

#endif