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
  bool cmd_version = false;
  double cmd_time_limit = 0;
  int cmd_random_seed = 0;

  std::string model_file = "";
  std::string cmd_read_basis_file = "";
  std::string cmd_write_basis_file = "";
  std::string options_file = "";
  std::string cmd_read_solution_file = "";
  std::string cmd_presolve = "";
  std::string cmd_solver = "";
  std::string cmd_parallel = "";
  std::string cmd_crossover = "";
  std::string cmd_write_solution_file = "";
  std::string cmd_write_model_file = "";
  std::string cmd_ranging = "";
};

void setupCommandLineOptions(CLI::App& app,
                             HighsCommandLineOptions& cmd_options) {
  // Does not work on windows with spaces in directory or file names.
  // const auto checkSingle = [](const std::string& input) -> std::string {
  //   if (input.find(' ') != std::string::npos) {
  //     return "Multiple files not implemented.";
  //   }
  //   return {};
  // };

  // Command line file specifications.
  app.add_option("--" + kModelFileString + "," + kModelFileString,
                 cmd_options.model_file, "File of model to solve.")
      // Can't use required here because it breaks version printing with -v.
      // ->required()
      // ->check(checkSingle)
      ->check(CLI::ExistingFile);

  app.add_option("--" + kOptionsFileString, cmd_options.options_file,
                 "File containing HiGHS options.")
      // ->check(checkSingle)
      ->check(CLI::ExistingFile);

  app.add_option("--" + kReadSolutionFileString,
                 cmd_options.cmd_read_solution_file,
                 "File of solution to read.")
      // ->check(checkSingle)
      ->check(CLI::ExistingFile);

  app.add_option("--" + kReadBasisFileString, cmd_options.cmd_read_basis_file,
                 "File of initial basis to read.")
      // ->check(checkSingle)
      ->check(CLI::ExistingFile);

  app.add_option("--" + kWriteModelFileString, cmd_options.cmd_write_model_file,
                 "File for writing out model.");
      // File does not need to exist
      //      ->check(CLI::ExistingFile)
      // ->check(checkSingle);

  app.add_option("--" + kWriteSolutionFileString,
                 cmd_options.cmd_write_solution_file,
                 "File for writing out solution.");
      // File does not need to exist
      //      ->check(CLI::ExistingFile)
      // ->check(checkSingle);

  app.add_option("--" + kWriteBasisFileString, cmd_options.cmd_write_basis_file,
                 "File for writing out final basis.");
      // File does not need to exist
      //      ->check(CLI::ExistingFile)
      // ->check(checkSingle);

  // Command line option specifications.
  app.add_option("--" + kPresolveString, cmd_options.cmd_presolve,
                 "Set presolve option to:\n"
                 "\"choose\" * default\n"
                 "\"on\"\n"
                 "\"off\"");

  app.add_option("--" + kSolverString, cmd_options.cmd_solver,
                 "Set solver option to:\n"
                 "\"choose\" * default\n"
                 "\"simplex\"\n"
                 "\"ipm\"");

  app.add_option("--" + kParallelString, cmd_options.cmd_parallel,
                 "Set parallel option to:\n"
                 "\"choose\" * default\n"
                 "\"on\"\n"
                 "\"off\"");

  app.add_option("--" + kRunCrossoverString, cmd_options.cmd_crossover,
                 "Set run_crossover option to:\n"
                 "\"choose\"\n"
                 "\"on\" * default\n"
                 "\"off\"");

  app.add_option("--" + kTimeLimitString, cmd_options.cmd_time_limit,
                 "Run time limit (seconds - double).");

  app.add_option("--" + kRandomSeedString, cmd_options.cmd_random_seed,
                 "Seed to initialize random number \ngeneration.");

  app.add_option("--" + kRangingString, cmd_options.cmd_ranging,
                 "Compute cost, bound, RHS and basic \nsolution ranging:\n"
                 "\"on\"\n"
                 "\"off\" * default");

  // Version.
  app.add_flag("--version,-v", cmd_options.cmd_version, "Print version.");
  app.set_help_flag("-h,--help", "Print help.");

  app.get_formatter()->column_width(33);

  app.get_formatter()->label("TEXT:FILE", "file");
  app.get_formatter()->label("TEXT", "text");
  app.get_formatter()->label("FLOAT", "float");
  app.get_formatter()->label("INT", "int");
  app.get_formatter()->label("OPTIONS", "options");
}

bool loadOptions(const CLI::App& app, const HighsLogOptions& report_log_options,
                 const HighsCommandLineOptions& c, HighsOptions& options) {
  if (c.cmd_version) {
    std::cout << "HiGHS version " << HIGHS_VERSION_MAJOR << "."
              << HIGHS_VERSION_MINOR << "." << HIGHS_VERSION_PATCH;
    std::cout << " Githash " << HIGHS_GITHASH << ". ";
    std::cout << kHighsCopyrightStatement << std::endl;
    exit(0);
  }

  // Options file.
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

  // Read solution file option.
  if (c.cmd_read_solution_file != "") {
    if (setLocalOptionValue(report_log_options, kReadSolutionFileString,
                            options.log_options, options.records,
                            c.cmd_read_solution_file) != OptionStatus::kOk)
      return false;
  }

  // Read basis file option.
  if (c.cmd_read_basis_file != "") {
    if (setLocalOptionValue(report_log_options, kReadBasisFileString,
                            options.log_options, options.records,
                            c.cmd_read_basis_file) != OptionStatus::kOk)
      return false;
  }

  // Write model file option.
  if (c.cmd_write_model_file != "") {
    if (setLocalOptionValue(report_log_options, kWriteModelFileString,
                            options.log_options, options.records,
                            c.cmd_write_model_file) != OptionStatus::kOk)
      return false;
  }

  // Write solution file option.
  if (c.cmd_write_solution_file != "") {
    if (setLocalOptionValue(report_log_options, kWriteSolutionFileString,
                            options.log_options, options.records,
                            c.cmd_write_solution_file) != OptionStatus::kOk)
      return false;
  }

  // Write basis file option.
  if (c.cmd_write_basis_file != "") {
    if (setLocalOptionValue(report_log_options, kWriteBasisFileString,
                            options.log_options, options.records,
                            c.cmd_write_basis_file) != OptionStatus::kOk)
      return false;
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
  if (app.count("--" + kTimeLimitString) > 0) {
    if (setLocalOptionValue(report_log_options, kTimeLimitString,
                            options.records,
                            c.cmd_time_limit) != OptionStatus::kOk)
      return false;
  }

  // Random seed option.
  if (app.count("--" + kRandomSeedString) > 0) {
    HighsInt value = c.cmd_random_seed;
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

  if (c.model_file.size() == 0) {
    std::cout << "Please specify filename in .mps|.lp|.ems format.\n";
    return false;
  }

  return true;
}

#endif
