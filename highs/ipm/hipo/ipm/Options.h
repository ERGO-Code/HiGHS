#ifndef HIPO_OPTION_H
#define HIPO_OPTION_H

#include "Parameters.h"
#include "io/HighsIO.h"
#include "lp_data/HighsOptions.h"

namespace hipo {

enum ParallelType { kOff, kChoose, kOn };
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
  Int parallel[kParallelCount];

  // Logging
  bool display = true;
  bool timeless_log = false;
  const HighsLogOptions* log_options = nullptr;

  inline void chooseParallel(Int bit, bool is_on) {
    Int type_default = is_on ? kOn : kOff;
    if (parallel[bit] == kChoose) parallel[bit] = type_default;
  }
};

inline bool testParallelBit(Int option, Int bit) { return option & (1 << bit); }

}  // namespace hipo

#endif