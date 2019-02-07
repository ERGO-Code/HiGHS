#include "HighsOptions.h"

bool loadOptions(int argc, char **argv, HighsOptions &options)
{
  try
  {
    cxxopts::Options cxx_options(argv[0], "HiGHS options");
    cxx_options.positional_help("[filename(s)]").show_positional_help();

    cxx_options.add_options()(
        "f, filename",
        "Filename of LP to solve.",
        cxxopts::value<std::vector<std::string>>())(
        "p, presolve", "Presolve: on | off. On by default.",
        cxxopts::value<std::string>())(
        "c, crash",
        "Crash mode: off | ltssf | ltssf1 | ... | ltssf7 | bs | singts.",
        cxxopts::value<std::string>())(
        "e, edge-weight", "Edge weight: Dan | Dvx | DSE | DSE0 | DSE2Dvx.",
        cxxopts::value<std::string>())(
        "P, price", "Price: Row | Col | RowSw | RowSwColSw | RowUltra. ",
        cxxopts::value<std::string>())("s, sip", "Use option sip.",
                                       cxxopts::value<bool>())(
        "S, scip", "Use option SCIP (to test utilities)",
        cxxopts::value<bool>())("m, pami",
                                "Use parallel solve.",
                                cxxopts::value<bool>())(
        "t, partition", "Use pami with partition file: filename",
        cxxopts::value<std::string>())("i, ipx", "Use interior point solver.",
                                       cxxopts::value<bool>())(
        "T, time-limit", "Use time limit.", cxxopts::value<double>())(
        "h, help", "Print help.");

    cxx_options.parse_positional("filename");

    auto result = cxx_options.parse(argc, argv);

    if (result.count("help"))
    {
      std::cout << cxx_options.help({""}) << std::endl;
      exit(0);
    }

    // Currently works for only one filename at a time.
    if (result.count("filename"))
    {
      std::string filename = "";
      auto &v = result["filename"].as<std::vector<std::string>>();
      if (v.size() > 1)
      {
        std::cout << "Multiple files not implemented.\n";
        return false;
      }
      options.filename = v[0];
    }

    if (result.count("presolve"))
    {
      std::string data = result["presolve"].as<std::string>();
      std::transform(data.begin(), data.end(), data.begin(), ::tolower);
      if (data != "on" && data != "off")
      {
        std::cout << "Wrong value specified for presolve." << std::endl;
        std::cout << cxx_options.help({""}) << std::endl;
        exit(0);
      }
      if (data == "on")
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