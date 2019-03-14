#ifndef LP_DATA_HIGHSRUNTIMEOPTIONS_H_
#define LP_DATA_HIGHSRUNTIMEOPTIONS_H_

#include "../external/cxxopts.hpp"

#include "io/LoadProblem.h"
#include "lp_data/HighsOptions.h"
#include "lp_data/HConst.h"
#include "util/stringutil.h"

bool loadOptions(int argc, char **argv, HighsOptions &options)
{
  try
  {
    cxxopts::Options cxx_options(argv[0], "HiGHS options");
    cxx_options.positional_help("[file]").show_positional_help();

    std::string presolve, crash, simplex, ipm, parallel;

    cxx_options.add_options()(
        "file",
        "Filename of LP to solve.",
        cxxopts::value<std::vector<std::string>>())(
        "presolve", "Use presolve: off by default.",
        cxxopts::value<std::string>(presolve))(
        "crash", "Use crash to start simplex: off by default.",
        cxxopts::value<std::string>(crash))(
        "simplex", "Use simplex solver: on by default.",
        cxxopts::value<std::string>(simplex))(
        "ipm", "Use interior point method solver: off by default.",
        cxxopts::value<std::string>(ipm))(
        "parallel", "Use parallel solve: off by default.",
        cxxopts::value<std::string>(parallel))(
        "time_limit", "Use time limit.",
        cxxopts::value<double>())(
        "iteration_limit", "Use iteration limit (integer).",
        cxxopts::value<int>())(
        "options_file",
        "File containing HiGHS options.",
        cxxopts::value<std::vector<std::string>>())(
        "parser_type", "Mps parser type: swap back to fixed format parser.",
        cxxopts::value<std::string>(presolve))(
        "h, help", "Print help.");

    cxx_options.parse_positional("file");

    auto result = cxx_options.parse(argc, argv);

    if (result.count("help"))
    {
      std::cout << cxx_options.help({""}) << std::endl;
      exit(0);
    }

    if (result.count("file"))
    {
      auto &v = result["file"].as<std::vector<std::string>>();
      if (v.size() > 1) {
        int nonEmpty = 0;
        for (int i=0; i<v.size(); i++) {
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

    if (result.count("presolve"))
    {
      std::string value = result["presolve"].as<std::string>();
      if (!setUserOptionValue(options, "presolve", value))
        HighsPrintMessage(ML_ALWAYS, "Unknown value for presovle option: %s. Ignored.\n", value.c_str());
    }

    if (result.count("crash"))
    {
      std::string value = result["crash"].as<std::string>();
      if (!setUserOptionValue(options, "crash", value))
        HighsPrintMessage(ML_ALWAYS, "Unknown value for crash option: %s. Ignored.\n", value.c_str());
    }

    if (result.count("parallel"))
    {
      std::string value = result["parallel"].as<std::string>();
      if (!setUserOptionValue(options, "parallel", value))
        HighsPrintMessage(ML_ALWAYS, "Unknown value for parallel option: %s. Ignored.\n", value.c_str());

      // JH will this be kept as an option?
      // else
      //  options.permutelp = true;
    }

    if (result.count("simplex"))
    {
      std::string value = result["simplex"].as<std::string>();
      if (!setUserOptionValue(options, "simplex", value))
        HighsPrintMessage(ML_ALWAYS, "Unknown value for simplex option: %s. Ignored.\n", value.c_str());
    }

    if (result.count("ipm"))
    {
      std::string value = result["ipm"].as<std::string>();
      if (!setUserOptionValue(options, "ipm", value))
        HighsPrintMessage(ML_ALWAYS, "Unknown value for ipm option: %s. Ignored.\n", value.c_str());
    }

    if (result.count("time_limit"))
    {
      double time_limit = result["time_limit"].as<double>();
      if (time_limit <= 0)
      {
        std::cout << "Time limit must be positive." << std::endl;
        std::cout << cxx_options.help({""}) << std::endl;
        exit(0);
      }
      options.highs_run_time_limit = time_limit;
    }

    if (result.count("iteration_limit"))
    {
      double iteration_limit = result["iteration-limit"].as<int>();
      if (iteration_limit <= 0)
      {
        std::cout << "Iteration limit must be positive." << std::endl;
        std::cout << cxx_options.help({""}) << std::endl;
        exit(0);
      }
      options.simplex_iteration_limit = iteration_limit;
    }

    if (result.count("options_file"))
    {
      auto &v = result["options_file"].as<std::vector<std::string>>();
      if (v.size() > 1)
      {
        std::cout << "Multiple options files not implemented.\n";
        return false;
      }
      options.options_file = v[0];
      loadOptionsFromFile(options);
    }

    // For testing of new parser
    if (result.count("parser_type"))
    {
      std::string value = result["parser_type"].as<std::string>();
      if (value == "fixed")
        options.parser_type = HighsMpsParserType::fixed;
      else if (value == "free")
        options.parser_type = HighsMpsParserType::free;
      else {
        std::cout << "Illegall value for parser." << std::endl;
        std::cout << cxx_options.help({""}) << std::endl;
        exit(0);
      }
    }
  }
  catch (const cxxopts::OptionException &e)
  {
    std::cout << "error parsing options: " << e.what() << std::endl;
    return false;
  }

  if (options.filename.size() == 0)
  {
    std::cout << "Please specify filename in .mps|.lp|.ems|.gz format.\n";
    return false;
  }

  return true;
}

#endif