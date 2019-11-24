/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/ICrash.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "presolve/ICrash.h"
#include "lp_data/HighsLpUtils.h"

struct Quadratic {
  const HighsLp lp;
  const ICrashOptions options;
  std::vector<ICrashIterationDetails> details;

  HighsSolution xk;

  double lp_objective;
  double quadratic_objective;

  std::vector<double> residual;
  double residual_norm_1;
  double residual_norm_2;

  double mu;
  std::vector<double> lambda;
  bool ok = false;
  
  Quadratic(HighsLp lp_, ICrashOptions options_) : lp(lp_), options(options_) {}
};

bool parseICrashStrategy(const std::string& strategy,
                         ICrashStrategy& icrash_strategy) {
  std::string lower = strategy;  // todo: tolower, trim
  if (lower == "penalty")
    icrash_strategy = ICrashStrategy::kPenalty;
  else if (lower == "admm")
    icrash_strategy = ICrashStrategy::kAdmm;
  else if (lower == "ica")
    icrash_strategy = ICrashStrategy::kICA;
  else if (lower == "breakpoints")
    icrash_strategy = ICrashStrategy::kBreakpoints;
  else
    return false;
  return true;
}

bool checkOptions(const HighsLp& lp, const ICrashOptions options) {
  if (options.strategy == ICrashStrategy::kBreakpoints) {
    if (options.exact) {
      HighsPrintMessage(
          ML_ALWAYS,
          "ICrash error: exact strategy not allowed for kBreakpoints./n");
      return false;
    }
    if (options.dualize) {
      HighsPrintMessage(
          ML_ALWAYS,
          "ICrash error: kBreakpoints does not support dualize option.\n");
      return false;
    }
  }
  return true;
}

Quadratic parseOptions(const HighsLp& lp, const ICrashOptions options) {
  HighsLp ilp = lp;
  convertToMinimization(ilp);
  if (isEqualityProblem(ilp)) {
    if (options.dualize) ilp = dualizeEqualityProblem(ilp);
  } else {
    // not equality problem.
    if (options.strategy != ICrashStrategy::kBreakpoints) {
      ilp = transformIntoEqualityProblem(ilp);
      if (options.dualize) {
        // Add slacks & dualize.
        // dualizeEqualityProblem returns a minimization equality problem.
        ilp = dualizeEqualityProblem(ilp);
      }
    }
  }

  return Quadratic{ilp, options};
}

double getQuadraticObjective(const Quadratic& idata) {
  // c'x
  double quadratic = vectorProduct(idata.lp.colCost_, idata.xk.col_value);

  // lambda'x
  quadratic += vectorProduct(idata.lambda, idata.residual);

  // 1/2mu r'r
  double rtr = vectorProduct(idata.residual, idata.lambda);
  quadratic += rtr / (2 * idata.mu);

  return quadratic;
}

HighsStatus CallICrash(const HighsLp& lp, const ICrashOptions& options,
                       ICrashInfo& result) {
  if (!checkOptions(lp, options)) return HighsStatus::Error;

  Quadratic idata = parseOptions(lp, options);

  // todo: continue refactoring

  return HighsStatus::Warning;
}