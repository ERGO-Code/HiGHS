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

#include "HighsStatus.h"
#include "io/HighsIO.h"
#include "lp_data/HighsLpUtils.h"
#include "presolve/ICrashUtil.h"
#include "util/HighsUtils.h"

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

bool initialize(Quadratic& idata, const ICrashOptions& options) {
  if (!initialize(idata.lp, idata.xk, idata.lambda)) return false;

  idata.mu = options.starting_weight;

  // Maybe other values for x0.
  return true;
}

void update(Quadratic& idata) {
  idata.lp_objective = vectorProduct(idata.lp.colCost_, idata.xk.col_value);

  calculateRowValues(idata.lp, idata.xk);
  bool piecewise = (idata.options.strategy == ICrashStrategy::kBreakpoints)
                       ? true
                       : false;
  updateResidual(piecewise, idata.lp, idata.xk, idata.residual);
  idata.residual_norm_2 = getNorm2(idata.residual);
}

HighsStatus CallICrash(const HighsLp& lp, const ICrashOptions& options,
                       ICrashInfo& result) {
  if (!checkOptions(lp, options)) return HighsStatus::Error;

  Quadratic idata = parseOptions(lp, options);

  initialize(idata, options);
  update(idata);
  // todo: continue refactoring

  return HighsStatus::Warning;
}