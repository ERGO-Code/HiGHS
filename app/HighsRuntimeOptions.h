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

#include "cxxopts.hpp"
#include "io/HighsIO.h"
#include "io/LoadOptions.h"
#include "util/stringutil.h"

bool loadOptions(const HighsLogOptions& report_log_options, int argc,
                 char** argv, HighsOptions& options, std::string& model_file,
                 std::string& read_basis_file, std::string& output_basis_file,

                 std::string& read_solution_file) {
  try {
    cxxopts::Options cxx_options(argv[0], "HiGHS options");
    cxx_options.positional_help("[file]").show_positional_help();

    cxx_options.add_options()
        //
        // Command line file specifications
        //
        // model file
        (kModelFileString, "File of model to solve.",
         cxxopts::value<std::vector<std::string>>())
        // model file
        (kReadBasisFile, "File of initial basis to read.",
         cxxopts::value<std::vector<std::string>>())(
            kWriteBasisFile, "File of final basis to write.",
            cxxopts::value<std::vector<std::string>>())
        // read_solution file
        (kReadSolutionFileString, "File of solution to read.",
         cxxopts::value<std::vector<std::string>>())
        // options file
        (kOptionsFileString, "File containing HiGHS options.",
         cxxopts::value<std::vector<std::string>>())
        //
        // Command line option specifications
        //
        // presolve option
        (kPresolveString,
         "Presolve: \"choose\" by default - \"on\"/\"off\" are alternatives.",
         cxxopts::value<std::string>())
        // solver option
        (kSolverString,
         "Solver: \"choose\" by default - \"simplex\"/\"ipm\" are "
         "alternatives.",
         cxxopts::value<std::string>())
        // parallel option
        (kParallelString,
         "Parallel solve: \"choose\" by default - \"on\"/\"off\" are "
         "alternatives.",
         cxxopts::value<std::string>())
        // run_crossover option
        (kRunCrossoverString,
         "Run crossover: \"on\" by default - \"choose\"/\"off\" are "
         "alternatives.",
         cxxopts::value<std::string>())
        // time_limit option
        (kTimeLimitString, "Run time limit (seconds - double).",
         cxxopts::value<double>())
        // solution_file option
        (kSolutionFileString, "File for writing out model solution.",
         cxxopts::value<std::vector<std::string>>())
        // write_model_file option
        (kWriteModelFileString, "File for writing out model.",
         cxxopts::value<std::vector<std::string>>())
        // random_seed option
        (kRandomSeedString, "Seed to initialize random number generation.",
         cxxopts::value<HighsInt>())
        // ranging option
        (kRangingString, "Compute cost, bound, RHS and basic solution ranging.",
         cxxopts::value<std::string>())
        // version
        (kVersionString, "Print version.")("h, help", "Print help.");

    // Handle command line file specifications
    //
    // model file
    cxx_options.parse_positional("model_file");

    auto result = cxx_options.parse(argc, argv);

    if (result.count("help")) {
      std::cout << cxx_options.help({""}) << std::endl;
      exit(0);
    }
    if (result.count("version")) {
      std::cout << "HiGHS version " << HIGHS_VERSION_MAJOR << "."
                << HIGHS_VERSION_MINOR << "." << HIGHS_VERSION_PATCH;
      std::cout << " Githash " << HIGHS_GITHASH << ". ";
      std::cout << kHighsCopyrightStatement << std::endl;
      exit(0);
    }
    if (result.count(kModelFileString)) {
      auto& v = result[kModelFileString].as<std::vector<std::string>>();
      if (v.size() > 1) {
        HighsInt nonEmpty = 0;
        for (HighsInt i = 0; i < (HighsInt)v.size(); i++) {
          std::string arg = v[i];
          if (trim(arg).size() > 0) {
            nonEmpty++;
            model_file = arg;
          }
        }
        if (nonEmpty > 1) {
          std::cout << "Multiple files not implemented.\n";
          return false;
        }
      } else {
        model_file = v[0];
      }
    }

    read_basis_file = "";
    if (result.count(kReadBasisFile)) {
      auto& v = result[kReadBasisFile].as<std::vector<std::string>>();
      if (v.size() > 1) {
        HighsInt nonEmpty = 0;
        for (HighsInt i = 0; i < (HighsInt)v.size(); i++) {
          std::string arg = v[i];
          if (trim(arg).size() > 0) {
            nonEmpty++;
            read_basis_file = arg;
          }
        }
        if (nonEmpty > 1) {
          std::cout << "Multiple files not implemented.\n";
          return false;
        }
      } else {
        read_basis_file = v[0];
      }
    }

    output_basis_file = "";
    if (result.count(kWriteBasisFile)) {
      auto& v = result[kWriteBasisFile].as<std::vector<std::string>>();
      if (v.size() > 1) {
        HighsInt nonEmpty = 0;
        for (HighsInt i = 0; i < (HighsInt)v.size(); i++) {
          std::string arg = v[i];
          if (trim(arg).size() > 0) {
            nonEmpty++;
            output_basis_file = arg;
          }
        }
        if (nonEmpty > 1) {
          std::cout << "Multiple files not implemented.\n";
          return false;
        }
      } else {
        output_basis_file = v[0];
      }
    }

    read_solution_file = "";
    if (result.count(kReadSolutionFileString)) {
      auto& v = result[kReadSolutionFileString].as<std::vector<std::string>>();
      if (v.size() > 1) {
        HighsInt nonEmpty = 0;
        for (size_t i = 0; i < v.size(); i++) {
          std::string arg = v[i];
          if (trim(arg).size() > 0) {
            nonEmpty++;
            read_solution_file = arg;
          }
        }
        if (nonEmpty > 1) {
          std::cout << "Multiple files not implemented.\n";
          return false;
        }
      } else {
        read_solution_file = v[0];
      }
    }
    // options file
    if (result.count(kOptionsFileString)) {
      auto& v = result[kOptionsFileString].as<std::vector<std::string>>();
      if (v.size() > 1) {
        std::cout << "Multiple options files not implemented.\n";
        return false;
      }
      switch (loadOptionsFromFile(report_log_options, options, v[0])) {
        case HighsLoadOptionsStatus::kError:
          return false;
        case HighsLoadOptionsStatus::kEmpty:
          writeOptionsToFile(stdout, options.records);
          return false;
        default:
          break;
      }
    }

    // Handle command line option specifications
    //
    // presolve option
    if (result.count(kPresolveString)) {
      std::string value = result[kPresolveString].as<std::string>();
      if (setLocalOptionValue(report_log_options, kPresolveString,
                              options.log_options, options.records,
                              value) != OptionStatus::kOk)
        return false;
    }

    // solver option
    if (result.count(kSolverString)) {
      std::string value = result[kSolverString].as<std::string>();
      if (setLocalOptionValue(report_log_options, kSolverString,
                              options.log_options, options.records,
                              value) != OptionStatus::kOk)
        return false;
    }

    // parallel option
    if (result.count(kParallelString)) {
      std::string value = result[kParallelString].as<std::string>();
      if (setLocalOptionValue(report_log_options, kParallelString,
                              options.log_options, options.records,
                              value) != OptionStatus::kOk)
        return false;
    }

    // run_crossover option
    if (result.count(kRunCrossoverString)) {
      std::string value = result[kRunCrossoverString].as<std::string>();
      if (setLocalOptionValue(report_log_options, kRunCrossoverString,
                              options.log_options, options.records,
                              value) != OptionStatus::kOk)
        return false;
    }

    // time_limit option
    if (result.count(kTimeLimitString)) {
      double value = result[kTimeLimitString].as<double>();
      if (setLocalOptionValue(report_log_options, kTimeLimitString,
                              options.records, value) != OptionStatus::kOk)
        return false;
    }

    // solution_file option
    if (result.count(kSolutionFileString)) {
      auto& v = result[kSolutionFileString].as<std::vector<std::string>>();
      if (v.size() > 1) {
        std::cout << "Multiple solution files not implemented.\n";
        return false;
      }
      if (setLocalOptionValue(report_log_options, kSolutionFileString,
                              options.log_options, options.records,
                              v[0]) != OptionStatus::kOk ||
          setLocalOptionValue(report_log_options, "write_solution_to_file",
                              options.records, true) != OptionStatus::kOk)
        return false;
    }

    // write_model_file option
    if (result.count(kWriteModelFileString)) {
      auto& v = result[kWriteModelFileString].as<std::vector<std::string>>();
      if (v.size() > 1) {
        std::cout << "Multiple write model files not implemented.\n";
        return false;
      }
      if (setLocalOptionValue(report_log_options, kWriteModelFileString,
                              options.log_options, options.records,
                              v[0]) != OptionStatus::kOk ||
          setLocalOptionValue(report_log_options, "write_model_to_file",
                              options.records, true) != OptionStatus::kOk)
        return false;
    }

    // random_seed option
    if (result.count(kRandomSeedString)) {
      HighsInt value = result[kRandomSeedString].as<HighsInt>();
      if (setLocalOptionValue(report_log_options, kRandomSeedString,
                              options.records, value) != OptionStatus::kOk)
        return false;
    }

    // ranging option
    if (result.count(kRangingString)) {
      std::string value = result[kRangingString].as<std::string>();
      if (setLocalOptionValue(report_log_options, kRangingString,
                              options.log_options, options.records,
                              value) != OptionStatus::kOk)
        return false;
    }
  } catch (const cxxopts::exceptions::exception& e) {
    highsLogUser(report_log_options, HighsLogType::kError,
                 "Error parsing options: %s\n", e.what());
    return false;
  }

  const bool horrible_hack_for_windows_visual_studio = false;
  if (horrible_hack_for_windows_visual_studio) {
    // Until I know how to debug an executable using command line
    // arguments in Visual Studio on Windows, this is necessary!
    HighsInt random_seed = -3;
    if (random_seed >= 0) {
      if (setLocalOptionValue(report_log_options, kRandomSeedString,
                              options.records,
                              random_seed) != OptionStatus::kOk)
        return false;
    }
    model_file = "ml.mps";
  }

  if (model_file.size() == 0) {
    std::cout << "Please specify filename in .mps|.lp|.ems format.\n";
    return false;
  }

  return true;
}

#endif
