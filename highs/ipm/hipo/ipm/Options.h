#ifndef HIPO_OPTION_H
#define HIPO_OPTION_H

#include "Parameters.h"
#include "io/HighsIO.h"

namespace hipo {

enum OptionNla {
  kOptionNlaMin = 0,
  kOptionNlaAugmented = kOptionNlaMin,
  kOptionNlaNormEq,
  kOptionNlaChoose,
  kOptionNlaMax = kOptionNlaChoose,
  kOptionNlaDefault = kOptionNlaChoose
};

enum OptionCrossover {
  kOptionCrossoverMin = 0,
  kOptionCrossoverOff = kOptionCrossoverMin,
  kOptionCrossoverOn,
  kOptionCrossoverChoose,
  kOptionCrossoverMax = kOptionCrossoverChoose,
  kOptionCrossoverDefault = kOptionCrossoverOff
};

enum OptionParallel {
  kOptionParallelMin = 0,
  kOptionParallelOff = kOptionParallelMin,  // tree off     node off
  kOptionParallelOn,                        // tree on      node on
  kOptionParallelChoose,                    // tree choose  node choose
  kOptionParallelTreeOnly,                  // tree on      node off
  kOptionParallelNodeOnly,                  // tree off     node on
  kOptionParallelMax = kOptionParallelNodeOnly,
  kOptionParallelDefault = kOptionParallelChoose
};

struct Options {
  // Solver options
  OptionNla nla = kOptionNlaDefault;
  OptionCrossover crossover = kOptionCrossoverDefault;
  OptionParallel parallel = kOptionParallelDefault;

  // Ipm parameters
  Int max_iter = kMaxIterDefault;
  double feasibility_tol = kIpmTolDefault;
  double optimality_tol = kIpmTolDefault;
  double crossover_tol = kIpmTolDefault;
  bool refine_with_ipx = true;
  double time_limit = -1.0;

  // Logging
  bool display = true;
  bool display_ipx = false;
  bool timeless_log = false;
};

}  // namespace hipo

#endif