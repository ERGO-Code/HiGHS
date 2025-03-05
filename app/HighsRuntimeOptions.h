/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2024 by Julian Hall, Ivet Galabova,    */
/*    Leona Gottwald and Michael Feldmeier                               */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef LP_DATA_HIGHSRUNTIMEOPTIONS_H_
#define LP_DATA_HIGHSRUNTIMEOPTIONS_H_

#include <cassert>

#include "CLI11.hpp"
#include "HConfig.h"
#include "io/HighsIO.h"
#include "io/LoadOptions.h"
#include "util/stringutil.h"

struct HighsCommandLineOptions {
  std::string model_file = "";
  std::string options_file = "";
  std::string read_solution_file = "";
  std::string cmd_presolve = "";
  std::string cmd_solver = "";
  std::string cmd_parallel = "";
  std::string cmd_crossover = "";
  std::string cmd_time_limit = "";
  std::string cmd_solution_file = "";
  std::string cmd_write_model_file = "";
  std::string cmd_random_seed = "";
  std::string cmd_ranging = "";

  bool cmd_version = false;
};

void setupCommandLineOptions(CLI::App& app,
                             HighsCommandLineOptions& cmd_options) {
  // Command line file specifications.
  app.add_option("--" + kModelFileString + "," + kModelFileString,
                 cmd_options.model_file, "File of model to solve.")
      ->required()
      ->check(CLI::ExistingFile)
      ->check([](const std::string& input) -> std::string {
        std::cout << "IG INPUT IS" << input << std::endl;
        if (input.find(' ') != std::string::npos) {
          return "Multiple files not implemented.";
        }
        return {};
      });

  app.add_option("--" + kOptionsFileString, cmd_options.options_file,
                 "File containing HiGHS options.")
      ->check(CLI::ExistingFile);

  app.add_option("--" + kReadSolutionFileString, cmd_options.read_solution_file,
                 "File of solution to read.")
      ->check(CLI::ExistingFile);

  // Command line option specifications.
  app.add_option(
      "--" + kPresolveString, cmd_options.cmd_presolve,
      "Presolve: \"choose\" by default - \"on\"/\"off\" are alternatives.");
  app.add_option("--" + kSolverString, cmd_options.cmd_solver,
                 "Solver: \"choose\" by default - \"simplex\"/\"ipm\" are "
                 "alternatives.");
  app.add_option("--" + kParallelString, cmd_options.cmd_parallel,
                 "Parallel solve: \"choose\" by default - \"on\"/\"off\" are "
                 "alternatives.");
  app.add_option("--" + kRunCrossoverString, cmd_options.cmd_crossover,
                 "Run crossover: \"on\" by default - \"choose\"/\"off\" are "
                 "alternatives.");
  app.add_option("--" + kTimeLimitString, cmd_options.cmd_time_limit,
                 "Run time limit (seconds - double).");
  app.add_option("--" + kSolutionFileString, cmd_options.cmd_solution_file,
                 "File for writing out model solution.");
  app.add_option("--" + kWriteModelFileString, cmd_options.cmd_write_model_file,
                 "File for writing out model.");
  app.add_option("--" + kRandomSeedString, cmd_options.cmd_random_seed,
                 "Seed to initialize random number generation.");
  app.add_option("--" + kRangingString, cmd_options.cmd_ranging,
                 "Compute cost, bound, RHS and basic solution ranging.");
  app.add_option("-v,--" + kVersionString, cmd_options.cmd_version,
                 "Print version.");
}

bool loadOptions(const HighsLogOptions& report_log_options,
                 const HighsCommandLineOptions& c, HighsOptions& options) {
  if (c.cmd_version) {
    std::cout << "HiGHS version " << HIGHS_VERSION_MAJOR << "."
              << HIGHS_VERSION_MINOR << "." << HIGHS_VERSION_PATCH;
    std::cout << " Githash " << HIGHS_GITHASH << ". ";
    std::cout << kHighsCopyrightStatement << std::endl;
    exit(0);
  }

  assert(c.model_file != "");

  // options file
  if (c.options_file != "") {
    switch (loadOptionsFromFile(report_log_options, options, c.options_file)) {
      case HighsLoadOptionsStatus::kError:
        return false;
      case HighsLoadOptionsStatus::kEmpty:
        writeOptionsToFile(stdout, options.log_options, options.records);
        return false;
      default:
        break;
    }
  }

  // Handle command line option specifications.

  // Presolve option.
  if (c.cmd_presolve != "") {
    if (setLocalOptionValue(report_log_options, kPresolveString,
                            options.log_options, options.records,
                            c.cmd_presolve) != OptionStatus::kOk)
      return false;
  }

  // Solver option.
  if (c.cmd_solver != "") {
    if (setLocalOptionValue(report_log_options, kSolverString,
                            options.log_options, options.records,
                            c.cmd_solver) != OptionStatus::kOk)
      return false;
  }

  // Parallel option.
  if (c.cmd_parallel != "") {
    if (setLocalOptionValue(report_log_options, kParallelString,
                            options.log_options, options.records,
                            c.cmd_parallel) != OptionStatus::kOk)
      return false;
  }

  // Run crossover option.
  if (c.cmd_crossover != "") {
    if (setLocalOptionValue(report_log_options, kRunCrossoverString,
                            options.log_options, options.records,
                            c.cmd_crossover) != OptionStatus::kOk)
      return false;
  }

  // Time limit option.
  if (c.cmd_time_limit != "") {
    double value = atof(c.cmd_time_limit.c_str());
    if (setLocalOptionValue(report_log_options, kTimeLimitString,
                            options.records, value) != OptionStatus::kOk)
      return false;
  }

  // Solution file option.
  if (c.cmd_solution_file != "") {
    if (setLocalOptionValue(report_log_options, kSolutionFileString,
                            options.log_options, options.records,
                            c.cmd_solution_file) != OptionStatus::kOk ||
        setLocalOptionValue(report_log_options, "write_solution_to_file",
                            options.records, true) != OptionStatus::kOk)
      return false;
  }

  // Write model file option.
  if (c.cmd_write_model_file != "") {
    // std::cout << "Multiple write model files not implemented.\n";
    if (setLocalOptionValue(report_log_options, kWriteModelFileString,
                            options.log_options, options.records,
                            c.cmd_write_model_file) != OptionStatus::kOk ||
        setLocalOptionValue(report_log_options, "write_model_to_file",
                            options.records, true) != OptionStatus::kOk)
      return false;
  }

  // Random seed option.
  if (c.cmd_random_seed != "") {
    HighsInt value = atof(c.cmd_time_limit.c_str());
    if (setLocalOptionValue(report_log_options, kRandomSeedString,
                            options.records, value) != OptionStatus::kOk)
      return false;
  }

  // Ranging option.
  if (c.cmd_ranging != "") {
    if (setLocalOptionValue(report_log_options, kRangingString,
                            options.log_options, options.records,
                            c.cmd_ranging) != OptionStatus::kOk)
      return false;
  }

  // If this horrible hack is still needed, please add it to RunHighs.
  // model_file is now const in this method.

  // const bool horrible_hack_for_windows_visual_studio = false;
  // if (horrible_hack_for_windows_visual_studio) {
  //   // Until I know how to debug an executable using command line
  //   // arguments in Visual Studio on Windows, this is necessary!
  //   HighsInt random_seed = -3;
  //   if (random_seed >= 0) {
  //     if (setLocalOptionValue(report_log_options, kRandomSeedString,
  //                             options.records,
  //                             random_seed) != OptionStatus::kOk)
  //       return false;
  //   }
  //   model_file = "ml.mps";
  // }

  // if (model_file.size() == 0) {
  //   std::cout << "Please specify filename in .mps|.lp|.ems format.\n";
  //   return false;
  // }

  return true;
}

#endif
