/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef LP_DATA_HIGHSRUNTIMEOPTIONS_H_
#define LP_DATA_HIGHSRUNTIMEOPTIONS_H_

#include "cxxopts.hpp"
#include "io/HighsIO.h"
#include "io/LoadOptions.h"
#include "util/stringutil.h"

bool loadOptions(int argc, char** argv, HighsOptions& options, std::string& model_file) {
  try {
    cxxopts::Options cxx_options(argv[0], "HiGHS options");
    cxx_options.positional_help("[file]").show_positional_help();

    std::string presolve, solver, parallel;

    cxx_options.add_options()(kModelFileString, "File of model to solve.",
                              cxxopts::value<std::vector<std::string>>())(
        kPresolveString,
        "Presolve: \"choose\" by default - \"on\"/\"off\"/\"mip\" are "
        "alternatives.",
        cxxopts::value<std::string>(presolve))(
        kSolverString,
        "Solver: \"choose\" by default - \"simplex\"/\"ipm\" are alternatives.",
        cxxopts::value<std::string>(solver))(
        kParallelString,
        "Parallel solve: \"choose\" by default - \"on\"/\"off\" are "
        "alternatives.",
        cxxopts::value<std::string>(parallel))(
        kTimeLimitString, "Run time limit (double).", cxxopts::value<double>())(
        kOptionsFileString, "File containing HiGHS options.",
        cxxopts::value<std::vector<std::string>>())("h, help", "Print help.");
    cxx_options.parse_positional("model_file");

    auto result = cxx_options.parse(argc, argv);

    if (result.count("help")) {
      std::cout << cxx_options.help({""}) << std::endl;
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

    if (result.count(kPresolveString)) {
      std::string value = result[kPresolveString].as<std::string>();
      if (setLocalOptionValue(options.log_options, kPresolveString,
                              options.records, value) != OptionStatus::kOk)
        return false;
    }

    if (result.count(kSolverString)) {
      std::string value = result[kSolverString].as<std::string>();
      if (setLocalOptionValue(options.log_options, kSolverString,
                              options.records, value) != OptionStatus::kOk)
        return false;
    }

    if (result.count(kParallelString)) {
      std::string value = result[kParallelString].as<std::string>();
      if (setLocalOptionValue(options.log_options, kParallelString,
                              options.records, value) != OptionStatus::kOk)
        return false;
    }

    if (result.count(kTimeLimitString)) {
      double value = result[kTimeLimitString].as<double>();
      if (setLocalOptionValue(options.log_options, kTimeLimitString,
                              options.records, value) != OptionStatus::kOk)
        return false;
    }

    if (result.count(kOptionsFileString)) {
      auto& v = result[kOptionsFileString].as<std::vector<std::string>>();
      if (v.size() > 1) {
        std::cout << "Multiple options files not implemented.\n";
        return false;
      }
      if (!loadOptionsFromFile(options, v[0])) return false;
    }

  } catch (const cxxopts::OptionException& e) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "Error parsing options: %s\n", e.what());
    return false;
  }

  if (model_file.size() == 0) {
    std::cout << "Please specify filename in .mps|.lp|.ems format.\n";
    return false;
  }

  return true;
}

#endif
