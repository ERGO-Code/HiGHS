/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file ../app/RunHighs.cpp
 * @brief HiGHS main
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "HConfig.h"
#include "Highs.h"
#include "HighsIO.h"
#include "HighsMipSolver.h"
#include "HighsOptions.h"
#include "HighsRuntimeOptions.h"
#include "HighsTimer.h"
#include "presolve/HAggregator.h"

void printHighsVersionCopyright(FILE* output, const int message_level,
                                const char* message = nullptr);
void reportLpStatsOrError(FILE* output, int message_level,
                          const HighsStatus read_status, const HighsLp& lp);
void reportSolvedLpStats(FILE* output, int message_level,
                         const HighsStatus run_status, const Highs& highs);
HighsStatus callLpSolver(HighsOptions& options, const HighsLp& lp);
HighsStatus callMipSolver(HighsOptions& options, const HighsLp& lp);

int main(int argc, char** argv) {
  printHighsVersionCopyright(stdout, ML_ALWAYS);

  // Load user options.
  HighsOptions options;
  bool options_ok = loadOptions(argc, argv, options);
  if (!options_ok) return 0;

  Highs highs;
  HighsStatus read_status = highs.readModel(options.model_file);
  reportLpStatsOrError(options.output, options.message_level, read_status,
                       highs.getLp());
  if (read_status == HighsStatus::Error)
    return 1;  // todo: change to read error

  // Run LP or MIP solver.
  // todo: modify so there is no extra matrix copy.
  // requires load model to be separated from HiGHS
  const HighsLp& lp = highs.getLp();
  HighsStatus run_status = HighsStatus::Error;
  bool is_mip = false;
  for (int i = 0; i < (int) lp.integrality_.size(); i++)
    if (lp.integrality_[i] == HighsVarType::INTEGER) {
      is_mip = true;
      break;
    }

  // todo: remove options.mip
  // it was used for development of the mip solver but now
  // options.solver will control which solver is used
  // instead.
  if (options.solver == "simplex" || options.solver == "ipm" ||
      (!options.mip && !is_mip && options.presolve != "mip")) {
    run_status = callLpSolver(options, lp);
  } else {
    run_status = callMipSolver(options, lp);
  }

  return (int)run_status;
}

void printHighsVersionCopyright(FILE* output, const int message_level,
                                const char* message) {
  HighsPrintMessage(output, message_level, ML_ALWAYS,
                    "Running HiGHS %d.%d.%d [date: %s, git hash: %s]\n",
                    HIGHS_VERSION_MAJOR, HIGHS_VERSION_MINOR,
                    HIGHS_VERSION_PATCH, HIGHS_COMPILATION_DATE, HIGHS_GITHASH);
  HighsPrintMessage(output, message_level, ML_ALWAYS,
                    "Copyright (c) 2020 ERGO-Code under MIT licence terms\n\n");
#ifdef HiGHSDEV
  // Report on preprocessing macros
  if (message != nullptr) {
    HighsPrintMessage(output, message_level, ML_ALWAYS, "In %s\n", message);
  }
#ifdef OPENMP
  HighsPrintMessage(output, message_level, ML_ALWAYS,
                    "OPENMP           is     defined\n");
#else
  HighsPrintMessage(output, message_level, ML_ALWAYS,
                    "OPENMP           is not defined\n");
#endif

#ifdef SCIP_DEV
  HighsPrintMessage(output, message_level, ML_ALWAYS,
                    "SCIP_DEV         is     defined\n");
#else
  HighsPrintMessage(output, message_level, ML_ALWAYS,
                    "SCIP_DEV         is not defined\n");
#endif

#ifdef HiGHSDEV
  HighsPrintMessage(output, message_level, ML_ALWAYS,
                    "HiGHSDEV         is     defined\n");
#else
  HighsPrintMessage(output, message_level, ML_ALWAYS,
                    "HiGHSDEV         is not defined\n");
#endif
  HighsPrintMessage(output, message_level, ML_ALWAYS,
                    "Built with CMAKE_BUILD_TYPE=%s\n", CMAKE_BUILD_TYPE);
#endif
}

void reportLpStatsOrError(FILE* output, int message_level,
                          const HighsStatus read_status, const HighsLp& lp) {
  if (read_status == HighsStatus::Error) {
    HighsPrintMessage(output, message_level, ML_ALWAYS, "Error loading file\n");
  } else {
    HighsPrintMessage(output, message_level, ML_ALWAYS, "LP       : %s\n",
                      lp.model_name_.c_str());
    HighsPrintMessage(output, message_level, ML_ALWAYS, "Rows     : %d\n",
                      lp.numRow_);
    HighsPrintMessage(output, message_level, ML_ALWAYS, "Cols     : %d\n",
                      lp.numCol_);
    HighsPrintMessage(output, message_level, ML_ALWAYS, "Nonzeros : %d\n",
                      lp.Avalue_.size());
    int num_int = 0;
    for (unsigned int i = 0; i < lp.integrality_.size(); i++)
      if (lp.integrality_[i] != HighsVarType::CONTINUOUS) num_int++;
    if (num_int)
      HighsPrintMessage(output, message_level, ML_ALWAYS, "Integer  : %d\n",
                        num_int);
  }
}

void reportSolvedLpStats(FILE* output, int message_level,
                         const HighsStatus run_status, Highs& highs) {
  if (run_status == HighsStatus::Error) {
    std::string statusname = HighsStatusToString(run_status);
    HighsPrintMessage(output, message_level, ML_ALWAYS, "HiGHS status: %s\n",
                      statusname.c_str());
  } else {
    HighsPrintMessage(output, message_level, ML_ALWAYS, "\n");
    HighsModelStatus model_status = highs.getModelStatus();
    HighsModelStatus scaled_model_status = highs.getModelStatus(true);
    HighsInfo highs_info = highs.getHighsInfo();
    if (model_status != scaled_model_status) {
      if (scaled_model_status == HighsModelStatus::OPTIMAL) {
        // The scaled model has been solved to optimality, but not the
        // unscaled model, flag this up, but report the scaled model
        // status
        HighsPrintMessage(output, message_level, ML_ALWAYS,
                          "Primal infeasibility: %10.3e (%d)\n",
                          highs_info.max_primal_infeasibility,
                          highs_info.num_primal_infeasibilities);
        HighsPrintMessage(output, message_level, ML_ALWAYS,
                          "Dual   infeasibility: %10.3e (%d)\n",
                          highs_info.max_dual_infeasibility,
                          highs_info.num_dual_infeasibilities);
        model_status = scaled_model_status;
      }
    }
    HighsPrintMessage(output, message_level, ML_ALWAYS,
                      "Model   status      : %s\n",
                      highs.highsModelStatusToString(model_status).c_str());
    HighsPrintMessage(
        output, message_level, ML_ALWAYS, "Primal  status      : %s\n",
        highs.primalDualStatusToString(highs_info.primal_status).c_str());
    HighsPrintMessage(
        output, message_level, ML_ALWAYS, "Dual    status      : %s\n",
        highs.primalDualStatusToString(highs_info.dual_status).c_str());
    HighsPrintMessage(output, message_level, ML_ALWAYS,
                      "Simplex   iterations: %d\n",
                      highs_info.simplex_iteration_count);
    if (highs_info.ipm_iteration_count)
      HighsPrintMessage(output, message_level, ML_ALWAYS,
                        "IPM       iterations: %d\n",
                        highs_info.ipm_iteration_count);
    if (highs_info.crossover_iteration_count)
      HighsPrintMessage(output, message_level, ML_ALWAYS,
                        "Crossover iterations: %d\n",
                        highs_info.crossover_iteration_count);
    if (model_status == HighsModelStatus::OPTIMAL) {
      double objective_function_value;
      highs.getHighsInfoValue("objective_function_value",
                              objective_function_value);
      HighsPrintMessage(output, message_level, ML_ALWAYS,
                        "Objective value     : %17.10e\n",
                        objective_function_value);
    }
    double run_time = highs.getHighsRunTime();
    HighsPrintMessage(output, message_level, ML_ALWAYS,
                      "HiGHS run time      : %13.2f\n", run_time);
    // Possibly write the solution to a file
    const HighsOptions& options = highs.getHighsOptions();
    if (options.write_solution_to_file)
      highs.writeSolution(options.solution_file, options.write_solution_pretty);
  }
}

HighsStatus callLpSolver(HighsOptions& use_options, const HighsLp& lp) {
  FILE* output = use_options.output;
  const int message_level = use_options.message_level;

  // Solve LP case.
  Highs highs(use_options);
  const HighsOptions& options = highs.getHighsOptions();

  // // Load problem.
  highs.passModel(lp);
  // HighsStatus read_status = highs.readModel(options.model_file);
  // reportLpStatsOrError(output, message_level, read_status, highs.getLp());
  // if (read_status == HighsStatus::Error) return HighsStatus::Error;

  // Run HiGHS.
  highs.setBasis();
  HighsStatus run_status = highs.run();

  reportSolvedLpStats(output, message_level, run_status, highs);
  return run_status;
}

HighsStatus callMipSolver(HighsOptions& use_options, const HighsLp& lp) {
  FILE* output = use_options.output;
  const int message_level = use_options.message_level;

  // Highs highs(use_options);
  // const HighsOptions& options = highs.getHighsOptions();
  // HighsStatus read_status = highs.readModel(options.model_file);
  // reportLpStatsOrError(output, message_level, read_status, highs.getLp());
  // if (read_status == HighsStatus::Error) return HighsStatus::Error;
  // HighsLp lp = highs.getLp();
  // printf("running aggregator (nnz = %lu)\n", lp.Avalue_.size());
  // presolve::HAggregator aggregator(lp.rowLower_, lp.rowUpper_, lp.colCost_,
  //                                  lp.offset_, lp.integrality_, lp.colLower_,
  //                                  lp.colUpper_);

  // aggregator.fromCSC(lp.Avalue_, lp.Aindex_, lp.Astart_);
  // aggregator.run();
  // aggregator.toCSC(lp.Avalue_, lp.Aindex_, lp.Astart_);
  // printf("aggregator finished (nnz = %lu)\n", lp.Avalue_.size());

  HighsMipSolver solver(use_options, lp);
  solver.run();

  return HighsStatus::OK;
}
