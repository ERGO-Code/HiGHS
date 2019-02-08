#include "HighsOptions.h"

bool loadOptions(int argc, char **argv, HighsOptions &options)
{
  try
  {
    cxxopts::Options cxx_options(argv[0], "HiGHS options");
    cxx_options.positional_help("[filename(s)]").show_positional_help();

    bool presolve = false;
    bool crash = false;
    bool simplex = true;
    bool ipm = false;
    bool parallel = false;

    cxx_options.add_options()(
        "file",
        "Filename of LP to solve.",
        cxxopts::value<std::vector<std::string>>())(
        "presolve", "Use presolve: off by default.",
        cxxopts::value<bool>())(
        "crash", "Use crash to start simplex: off by default.",
        cxxopts::value<bool>())(
        "simplex", "Use simplex solver: on by default.",
        cxxopts::value<bool>())(
        "ipm", "Use interior point method solver: off by default.",
        cxxopts::value<bool>())(
        "time-limit", "Use time limit.", cxxopts::value<double>())(
        "h, help", "Print help.");

    cxx_options.parse_positional("file");

    auto result = cxx_options.parse(argc, argv);

    if (result.count("help"))
    {
      std::cout << cxx_options.help({""}) << std::endl;
      exit(0);
    }

    // Currently works for only one filename at a time.
    if (result.count("file"))
    {
      std::string filename = "";
      auto &v = result["file"].as<std::vector<std::string>>();
      if (v.size() > 1)
      {
        std::cout << "Multiple files not implemented.\n";
        return false;
      }
      options.filename = v[0];
    }

    if (result.count("presolve"))
    {
        options.presolve_option = PresolveOption::ON;
      }
      else
      {
        options.presolve_option = PresolveOption::OFF;
      }
      std::cout << "Presolve is set to " << data << ".\n";
    }

    if (result.count("time-limit"))
    {
      double time_limit = result["time-limit"].as<double>();
      if (time_limit <= 0)
      {
        std::cout << "Time limit must be positive." << std::endl;
        std::cout << cxx_options.help({""}) << std::endl;
        exit(0);
      }
      options.highs_run_time_limit = time_limit;
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

  // Force column permutation of the LP to be used by the solver if
  // parallel code is to be used
  //  if (options.pami || options.sip) {options.permuteLp = true;}

  return true;
}