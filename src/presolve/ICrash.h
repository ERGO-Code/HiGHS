
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

#include "lp_data/HighsStatus.h"
#include "lp_data/HighsLp.h"

enum class ICrashStrategy {
  kPenalty,
  kAdmm,
  kICA
};

struct ICrashIterationDetails {
  int num;
  double weight;
  double lambda_norm_2;

  double lp_objective;
  double quadratic_objective;
  double residual_norm_2;
  
  double time;
};

struct ICrashInfo {
  int num_iterations;

  double final_lp_objective;
  double final_quadratic_objective;
  double final_residual_norm_2;

  double starting_weight;
  double final_weight;

  std::vector<ICrashIterationDetails> details;
  std::vector<double> x_values;

  double total_time;
};

struct ICrashOptions {
  bool dualize;
  ICrashStrategy strategy;
  double starting_weight;
  int iterations;
  int approximate_minimization_iterations;
  bool exact;
  bool breakpoints; // gets ignored if exact is set to true
};

HighsStatus callICrash(const HighsLp& lp, const ICrashOptions& options,
                       ICrashInfo& result);

bool parseICrashStrategy(const std::string& strategy,
                         ICrashStrategy& icrash_strategy);

#endif
