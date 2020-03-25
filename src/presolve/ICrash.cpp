/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/ICrash.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "presolve/ICrash.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "HighsStatus.h"
#include "io/HighsIO.h"
#include "lp_data/HighsLpUtils.h"
#include "presolve/ICrashUtil.h"
#include "util/HighsUtils.h"
#include "util/stringutil.h"

constexpr double kExitTolerance = 0.00000001;

bool parseICrashStrategy(const std::string& strategy,
                         ICrashStrategy& icrash_strategy) {
  std::string lower = strategy;
  trim(lower);
  std::transform(lower.begin(), lower.end(), lower.begin(),
                 [](unsigned char c) { return std::tolower(c); });

  if (lower == "penalty")
    icrash_strategy = ICrashStrategy::kPenalty;
  else if (lower == "admm")
    icrash_strategy = ICrashStrategy::kAdmm;
  else if (lower == "ica")
    icrash_strategy = ICrashStrategy::kICA;
  else if (lower == "update_penalty")
    icrash_strategy = ICrashStrategy::kUpdatePenalty;
  else if (lower == "update_admm")
    icrash_strategy = ICrashStrategy::kUpdateAdmm;
  else
    return false;
  return true;
}

bool checkOptions(const HighsLp& lp, const ICrashOptions options) {
  if (options.exact) {
    HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                      "ICrashError: exact subproblem solution not available "
                      "at the moment.\n");
    return false;
  }

  if (options.breakpoints) {
    if (options.exact) {
      HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                        "ICrashError: exact strategy not allowed for "
                        "breakpoints minimization.\n");
      return false;
    }
    if (options.dualize) {
      HighsPrintMessage(
          options.output, options.message_level, ML_ALWAYS,
          "ICrashError: breakpoints does not support dualize option.\n");
      return false;
    }
    HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                      "ICrashError: breakpoints not implemented yet.\n");
    return false;
  }

  if (options.strategy == ICrashStrategy::kPenalty)
    HighsPrintMessage(
        options.output, options.message_level, ML_ALWAYS,
        "ICrash Warning: Using solveSubproblemICA with lambda = 0.\n");

  return true;
}

Quadratic parseOptions(const HighsLp& lp, const ICrashOptions options) {
  HighsLp ilp = lp;
  HighsLp local_lp;
  HighsStatus status;
  convertToMinimization(ilp);
  if (isEqualityProblem(ilp)) {
    if (options.dualize) {
      status = dualizeEqualityProblem(ilp, local_lp);
      if (status == HighsStatus::OK) {
        ilp = local_lp;
      } else {
        printf("Cannot dualise equality problem\n");
      }
    }
  } else {
    // not equality problem.
    assert(!options.breakpoints);  // remove when implementing breakpoints and
                                   // add if else.
    status = transformIntoEqualityProblem(ilp, local_lp);
    if (status == HighsStatus::OK) {
      ilp = local_lp;
    } else {
      printf("Cannot transform into equality problem\n");
    }
    if (options.dualize) {
      // Add slacks & dualize.
      // dualizeEqualityProblem returns a minimization equality problem.
      status = dualizeEqualityProblem(ilp, local_lp);
      if (status == HighsStatus::OK) {
        ilp = local_lp;
      } else {
        printf("Cannot dualise equality problem\n");
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
  // lp_objective
  idata.lp_objective = vectorProduct(idata.lp.colCost_, idata.xk.col_value);

  // residual & residual_norm_2
  calculateRowValues(idata.lp, idata.xk);
  updateResidual(idata.options.breakpoints, idata.lp, idata.xk, idata.residual);
  idata.residual_norm_2 = getNorm2(idata.residual);

  // quadratic_objective
  idata.quadratic_objective = idata.lp_objective;
  idata.quadratic_objective += vectorProduct(idata.lambda, idata.residual);
  idata.quadratic_objective +=
      vectorProduct(idata.residual, idata.residual) / (2 * idata.mu);
}

ICrashIterationDetails fillDetails(const int num, const Quadratic& idata) {
  double lambda_norm_2 = getNorm2(idata.lambda);
  return ICrashIterationDetails{num,
                                idata.mu,
                                lambda_norm_2,
                                idata.lp_objective,
                                idata.quadratic_objective,
                                idata.residual_norm_2,
                                0};
}

void fillICrashInfo(const int n_iterations, ICrashInfo& result) {
  assert((int)result.details.size() == n_iterations + 1);
  result.num_iterations = n_iterations;

  result.final_lp_objective = result.details[n_iterations].lp_objective;
  result.final_quadratic_objective =
      result.details[n_iterations].quadratic_objective;
  result.final_residual_norm_2 = result.details[n_iterations].residual_norm_2;

  result.starting_weight = result.details[0].weight;
  result.final_weight = result.details[n_iterations].weight;
}

void updateParameters(Quadratic& idata, const ICrashOptions& options,
                      const int iteration) {
  if (iteration == 1) return;

  // The other strategies are WIP.
  switch (options.strategy) {
    case ICrashStrategy::kPenalty: {
      idata.mu = 0.1 * idata.mu;
      break;
    }
    case ICrashStrategy::kAdmm: {
      HighsPrintMessage(
          options.output, options.message_level, ML_ALWAYS,
          "ICrash Error: ADMM parameter update not implemented\n.");
      break;
    }
    case ICrashStrategy::kICA: {
      // Update mu every third iteration, otherwise update lambda.
      if (iteration % 3 == 0) {
        idata.mu = 0.1 * idata.mu;
      } else {
        std::vector<double> residual_ica(idata.lp.numRow_, 0);
        updateResidualIca(idata.lp, idata.xk, residual_ica);
        for (int row = 0; row < idata.lp.numRow_; row++)
          idata.lambda[row] = idata.mu * residual_ica[row];
      }
      break;
    }
    case ICrashStrategy::kUpdatePenalty: {
      // Update mu every third iteration, otherwise do nothing.
      if (iteration % 3 == 0) idata.mu = 0.1 * idata.mu;
      break;
    }
    case ICrashStrategy::kUpdateAdmm: {
      // Update mu every third iteration, otherwise update lambda.
      if (iteration % 3 == 0) {
        idata.mu = 0.1 * idata.mu;
      } else {
        std::vector<double> residual_ica(idata.lp.numRow_, 0);
        updateResidualIca(idata.lp, idata.xk, residual_ica);
        for (int row = 0; row < idata.lp.numRow_; row++)
          // todo: double check clp.
          idata.lambda[row] = idata.lambda[row] + idata.mu * residual_ica[row];
      }
      break;
    }
  }
}

void solveSubproblemICA(Quadratic& idata, const ICrashOptions& options) {
  bool minor_iteration_details = false;

  std::vector<double> residual_ica(idata.lp.numRow_, 0);
  updateResidualIca(idata.lp, idata.xk, residual_ica);
  double objective_ica = 0;

  for (int k = 0; k < options.approximate_minimization_iterations; k++) {
    for (int col = 0; col < idata.lp.numCol_; col++) {
      // determine whether to minimize for col.
      // if empty skip.
      if (idata.lp.Astart_[col] == idata.lp.Astart_[col + 1]) continue;

      double delta_x = 0;
      minimizeComponentIca(col, idata.mu, idata.lambda, idata.lp, objective_ica,
                           residual_ica, idata.xk);

      if (minor_iteration_details) {
        double quadratic_objective = getQuadraticObjective(idata);
        printMinorIterationDetails(k, col, idata.xk.col_value[col] - delta_x,
                                   delta_x, objective_ica, residual_ica,
                                   quadratic_objective);
      }

      assert(std::fabs(objective_ica -
                       vectorProduct(idata.lp.colCost_, idata.xk.col_value)) <
             1e08);
    }

    // code below just for checking. Can comment out later if speed up is
    // needed.
    std::vector<double> residual_ica_check(idata.lp.numRow_, 0);
    updateResidualIca(idata.lp, idata.xk, residual_ica_check);
    double difference = getNorm2(residual_ica) - getNorm2(residual_ica_check);
    assert(std::fabs(difference) < 1e08);
    (void)difference;
  }
}

bool solveSubproblem(Quadratic& idata, const ICrashOptions& options) {
  switch (options.strategy) {
    case ICrashStrategy::kUpdatePenalty:
    case ICrashStrategy::kUpdateAdmm:
    case ICrashStrategy::kICA: {
      assert(!options.exact);
      solveSubproblemICA(idata, options);
      break;
    }
    case ICrashStrategy::kPenalty: {
      HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                        "ICrashError: Not implemented yet.\n");
      return false;
    }
    default: {
      HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                        "ICrashError: Not implemented yet.\n");
      return false;
    }
  }
  return true;
}

void reportSubproblem(const ICrashOptions options, const Quadratic& idata,
                      const int iteration) {
  std::stringstream ss;
  // Report outcome.
  if (iteration == 0) {
    ss << "Iteration " << std::setw(3) << 0 << ": objective " << std::setw(3)
       << std::fixed << std::setprecision(2) << idata.lp_objective
       << " residual " << std::setw(5) << std::scientific
       << idata.residual_norm_2 << std::endl;
  } else {
    ss << "Iter " << std::setw(3) << iteration << ", mu " << idata.mu
       << std::scientific << ", c'x " << std::setprecision(5)
       << idata.lp_objective << ", res " << idata.residual_norm_2
       << ", quad_obj " << idata.quadratic_objective << std::endl;
  }
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    ss.str().c_str());
}

std::string ICrashtrategyToString(const ICrashStrategy strategy) {
  switch (strategy) {
    case ICrashStrategy::kPenalty:
      return "Penalty";
    case ICrashStrategy::kAdmm:
      return "ADMM";
    case ICrashStrategy::kICA:
      return "ICA";
    case ICrashStrategy::kUpdatePenalty:
      return "UpdatePenalty";
    case ICrashStrategy::kUpdateAdmm:
      return "UpdateAdmm";
  }
  return "ICrashError: Unknown strategy.\n";
}

void reportOptions(const ICrashOptions& options) {
  std::stringstream ss;
  // Report outcome.
  ss << "ICrashOptions \n"
     << "dualize: " << std::boolalpha << options.dualize << "\n"
     << "strategy: " << ICrashtrategyToString(options.strategy) << "\n"
     << "starting_weight: " << std::scientific << options.starting_weight
     << "\n"
     << "iterations: " << options.iterations << "\n";
  if (!options.exact) {
    ss << "approximate_minimization_iterations: "
       << options.approximate_minimization_iterations << "\n"
       << std::boolalpha << "breakpoints: " << options.breakpoints << "\n";
  } else {
    ss << "exact: true\n";
  }
  ss << "\n";
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    ss.str().c_str());
}

HighsStatus callICrash(const HighsLp& lp, const ICrashOptions& options,
                       ICrashInfo& result) {
  if (!checkOptions(lp, options)) return HighsStatus::Error;

  // Initialize data structures and initial values.
  Quadratic idata = parseOptions(lp, options);
  reportOptions(options);
  initialize(idata, options);
  update(idata);
  reportSubproblem(options, idata, 0);
  idata.details.push_back(fillDetails(0, idata));

  // Initialize clocks.
  std::chrono::time_point<std::chrono::system_clock> start, end,
      start_iteration, end_iteration;
  std::chrono::duration<double> elapsed_seconds;
  start = std::chrono::system_clock::now();

  // Main loop.
  int iteration = 0;
  for (iteration = 1; iteration <= options.iterations; iteration++) {
    updateParameters(idata, options, iteration);

    // Solve subproblem.
    start_iteration = std::chrono::system_clock::now();
    bool success = solveSubproblem(idata, options);
    if (!success) return HighsStatus::Error;
    end_iteration = std::chrono::system_clock::now();
    elapsed_seconds = end_iteration - start_iteration;

    update(idata);

    reportSubproblem(options, idata, iteration);
    idata.details.push_back(fillDetails(iteration, idata));
    assert(iteration + 1 == (int)idata.details.size());
    idata.details[iteration].time = elapsed_seconds.count();

    // Exit if feasible.
    if (idata.residual_norm_2 < kExitTolerance) {
      HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                        "Solution feasible within exit tolerance: %g.\n",
                        kExitTolerance);
      iteration++;
      break;
    }
  }

  // Fill in return values.
  iteration--;
  result.details = std::move(idata.details);
  // reportICrashIterationDetails(result.details);
  fillICrashInfo(iteration, result);
  result.x_values = idata.xk.col_value;

  end = std::chrono::system_clock::now();
  elapsed_seconds = end - start;
  result.total_time = elapsed_seconds.count();

  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "\nICrash finished successfully after: %.2f sec.\n\n",
                    result.total_time);

  return HighsStatus::OK;
}
