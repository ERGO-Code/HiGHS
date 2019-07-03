#ifndef LP_DATA_HIGHSRUNTIMEOPTIONS_H_
#define LP_DATA_HIGHSRUNTIMEOPTIONS_H_

#include "../external/cxxopts.hpp"

#include "io/HighsIO.h"
#include "io/LoadProblem.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsOptions.h"
#include "lp_data/HighsStatus.h"
#include "util/stringutil.h"

bool loadOptions(int argc, char** argv, HighsOptions& options) {
  try {
    cxxopts::Options cxx_options(argv[0], "HiGHS options");
    cxx_options.positional_help("[file]").show_positional_help();

    std::string presolve, crash, simplex, ipm, parallel, parser;

    cxx_options.add_options()(file_string, "Filename of LP to solve.",
                              cxxopts::value<std::vector<std::string>>())(
        presolve_string, "Use presolve: off by default.",
        cxxopts::value<std::string>(presolve))(
        crash_string, "Use crash to start simplex: off by default.",
        cxxopts::value<std::string>(crash))(
        parallel_string, "Use parallel solve: off by default.",
        cxxopts::value<std::string>(parallel))(
        simplex_string, "Use simplex solver: on by default.",
        cxxopts::value<std::string>(simplex))(
        ipm_string, "Use interior point method solver: off by default.",
        cxxopts::value<std::string>(ipm))(
        highs_run_time_limit_string, "Use HiGHS run time limit (double).",
        cxxopts::value<double>())(simplex_iteration_limit_string,
                                  "Use simplex iteration limit (integer).",
                                  cxxopts::value<int>())(
        options_file_string, "File containing HiGHS options.",
        cxxopts::value<std::vector<std::string>>())(
        mps_parser_type_string,
        "Mps parser type: swap back to fixed format parser.",
        cxxopts::value<std::string>(parser))("h, help", "Print help.");

    cxx_options.parse_positional("file");

    auto result = cxx_options.parse(argc, argv);

    if (result.count("help")) {
      std::cout << cxx_options.help({""}) << std::endl;
      exit(0);
    }

    if (result.count(file_string)) {
      auto& v = result[file_string].as<std::vector<std::string>>();
      if (v.size() > 1) {
        int nonEmpty = 0;
        for (int i = 0; i < (int)v.size(); i++) {
          std::string arg = v[i];
          if (trim(arg).size() > 0) {
            nonEmpty++;
            options.filename = arg;
          }
        }
        if (nonEmpty > 1) {
          std::cout << "Multiple files not implemented.\n";
          return false;
        }
      } else {
        options.filename = v[0];
      }
    }

    if (result.count(presolve_string)) {
      std::string value = result[presolve_string].as<std::string>();
      if (setPresolveValue(options, value) == OptionStatus::ILLEGAL_VALUE)
        return false;
    }

    if (result.count(crash_string)) {
      std::string value = result[crash_string].as<std::string>();
      if (setCrashValue(options, value) == OptionStatus::ILLEGAL_VALUE)
        return false;
    }

    if (result.count(parallel_string)) {
      std::string value = result[parallel_string].as<std::string>();
      if (setParallelValue(options, value) == OptionStatus::ILLEGAL_VALUE)
        return false;
    }

    if (result.count(simplex_string)) {
      std::string value = result[simplex_string].as<std::string>();
      if (setSimplexValue(options, value) == OptionStatus::ILLEGAL_VALUE)
        return false;
    }

    if (result.count(ipm_string)) {
      std::string value = result[ipm_string].as<std::string>();
      if (setIpmValue(options, value) == OptionStatus::ILLEGAL_VALUE)
        return false;
    }

    if (result.count(highs_run_time_limit_string)) {
      double value = result[highs_run_time_limit_string].as<double>();
      if (setHighsRunTimeLimitValue(options, value) ==
          OptionStatus::ILLEGAL_VALUE)
        return false;
    }

    if (result.count(simplex_iteration_limit_string)) {
      int value = result[simplex_iteration_limit_string].as<int>();
      if (setSimplexIterationLimitValue(options, value) ==
          OptionStatus::ILLEGAL_VALUE)
        return false;
    }

    if (result.count(options_file_string)) {
      auto& v = result[options_file_string].as<std::vector<std::string>>();
      if (v.size() > 1) {
        std::cout << "Multiple options files not implemented.\n";
        return false;
      }
      options.options_file = v[0];
      if (!loadOptionsFromFile(options)) return false;
    }

    // For testing of new parser
    if (result.count(mps_parser_type_string)) {
      std::string value = result[mps_parser_type_string].as<std::string>();
      if (setParserTypeValue(options, value) == OptionStatus::ILLEGAL_VALUE)
        return false;
    }
  } catch (const cxxopts::OptionException& e) {
    HighsLogMessage(HighsMessageType::ERROR, "Error parsing options: %s",
                    e.what());
    return false;
  }

  if (options.filename.size() == 0) {
    std::cout << "Please specify filename in .mps|.lp|.ems|.gz format.\n";
    return false;
  }

  return true;
}

#endif
