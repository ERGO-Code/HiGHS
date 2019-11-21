
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/ICrash.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef PRESOLVE_QUADRATIC_CRASH_H_
#define PRESOLVE_QUADRATIC_CRASH_H_

#include "lp_data/HighsOptions.h"
#include "HighsTimer.h"

struct ICrashOptions {
  bool dualize;
  std::string strategy;
  double starting_weight;
  int iterations;
  int approximate_minimization_iterations;
  bool exact;
};

struct ICrashIterationDetails {
  double weight;

  double residual;
  double lp_objective;
  double quadratic_objective;
};

struct ICrashInfo {
  int num_iterations;

  double final_residual;
  double final_lp_objective;
  double final_quadratic_objective;

  double starting_weight;
  double final_weight;

  std::vector<ICrashIterationDetails> details;
  std::vector<double> x_values;
};

enum class ICrashStrategy {
  PENALTY,
  ADMM,
  ICA,
  BREAKPOINTS,
};

HighsStatus CallICrash(const HighsLp& lp, const ICrashOptions& options,
                       ICrashInfo& result, HighsTimer& timer);

#endif
