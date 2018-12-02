/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsSetup.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef LP_DATA_HIGHS_SETUP_H_
#define LP_DATA_HIGHS_SETUP_H_

#include <algorithm>
#include <iostream>

#include "HApp.h"
#include "HighsOptions.h"
#include "cxxopts.hpp"

HModel HighsLpToHModel(const HighsLp& lp);
HighsLp HModelToHighsLp(const HModel& model);

// Class to set parameters and run HiGHS
class Highs {
 public:
  Highs() {}
  explicit Highs(const HighsOptions& opt) : options_(opt){};
  explicit Highs(const HighsStringOptions& opt) : options__(opt){};

  // The public method run(lp, solution) calls runSolver to solve problem before
  // or after presolve (or crash later?) depending on the specified options.
  HighsStatus run(const HighsLp& lp, HighsSolution& solution) const;

 private:
  // delete.
  HighsOptions options_;
  // use HighsStringOptions instead for now. Then rename to HighsOptions, once
  // previous one is gone.
  HighsStringOptions options__;
  HighsStatus runSolver(const HighsLp& lp, HighsSolution& solution) const;
};

// Checks the options calls presolve and postsolve if needed. Solvers are called
// with runSolver(..)
HighsStatus Highs::run(const HighsLp& lp, HighsSolution& solution) const {
  return runSolver(lp, solution);

  // todo
  //
  // if (!options_.presolve) {
  //  HighsSolution solution;
  //  return runSolver(lp, solution);
  //}
  // return HighsStatus::NotImplemented;

  /*
   HighsLp reduced_lp;

   // presolve(..) method below should use HPresolve now but should be
   // possible to use external presolve too. Link with ZIB presolve so clp
   // is possible later.
   status = presolve(lp, reduced_lp);
   checkStatus(status);

   switch (status) {
      case Status::ProblemReduced: {
         // Solution reduced_solution;
         // status = runSolver(reduced_lp, reduced_solution);
         // checkStatus(status);
         break;
      }
      case Status::ProblemReducedToEmpty:
         // Problem was reduced to empty so we proceed to postsolve
         break;
      default:
         checkStatus(status);
   }

   // Postsolve
   status = postsolve(lp, solution)

   // If needed set up clean up with simplex.
   */

  return HighsStatus::OK;
}

// The method below runs simplex or ipx solver on the lp.

HighsStatus Highs::runSolver(const HighsLp& lp, HighsSolution& solution) const {
  // assert(checkLp(lp) == LpError::none);

  HighsStatus status;
#ifndef IPX
  // HiGHS
  // todo: Without the presolve part, so will be
  //     = solve_simplex(options, reduced_lp, reduced_solution)
  status = solveSimplex(options_, lp, solution);
#else
  // IPX
  // todo:Check options for simplex-specific options

  status = solveIpx(options_, lp, solution);
  // If ipx crossover did not find optimality set up simplex.

#endif

  // todo:
  // assert(KktSatisfied(lp, solution));

  return status;
}

void HiGHSRun(const char* message) {
  std::cout << "Running HiGHS " << HIGHS_VERSION_MAJOR << "."
            << HIGHS_VERSION_MINOR << "." << HIGHS_VERSION_PATCH
            << " [date: " << HIGHS_COMPILATION_DATE
            << ", git hash: " << HIGHS_GITHASH << "]"
            << "\n"
            << "Copyright (c) 2018 ERGO-Code under MIT licence terms\n\n";
#ifdef HiGHSDEV
  // Report on preprocessing macros
  std::cout << "In " << message << std::endl;
  std::cout << "Built with CMAKE_BUILD_TYPE=" << CMAKE_BUILD_TYPE << std::endl;
#ifdef OLD_PARSER
  std::cout << "OLD_PARSER       is     defined" << std::endl;
#else
  std::cout << "OLD_PARSER       is not defined" << std::endl;
#endif

#ifdef OPENMP
  std::cout << "OPENMP           is     defined" << std::endl;
#else
  std::cout << "OPENMP           is not defined" << std::endl;
#endif

#ifdef SCIP_DEV
  std::cout << "SCIP_DEV         is     defined" << std::endl;
#else
  std::cout << "SCIP_DEV         is not defined" << std::endl;
#endif

#ifdef HiGHSDEV
  std::cout << "HiGHSDEV         is     defined" << std::endl;
#else
  std::cout << "HiGHSDEV         is not defined" << std::endl;
#endif

#ifdef HiGHSRELEASE
  std::cout << "HiGHSRELEASE     is     defined" << std::endl;
#else
  std::cout << "HiGHSRELEASE     is not defined" << std::endl;
#endif

#endif
};

HighsStatus loadOptions(int argc, char** argv,
                        HighsStringOptions& highs_options) {
  try {
    cxxopts::Options cxx_options(argv[0], "HiGHS options");
    cxx_options.positional_help("[filename(s)]").show_positional_help();

    cxx_options.add_options()(
        "f, filename",
        "Filename(s) of LPs to solve. The option specifier is not required.",
        cxxopts::value<std::vector<std::string>>())(
        "p, presolve", "Presolve: on | off. On by default.",
        cxxopts::value<std::string>())(
        "c, crash",
        "Crash mode: off | ltssf | ltssf1 | ... | ltssf7 | bs | singts.",
        cxxopts::value<std::string>())(
        "e, edge-weight",
        "Edge weight: Dan | Dvx | DSE | DSE0 | DSE2Dvx.",
        cxxopts::value<std::string>())(
        "P, price",
        "Price: Row | Col | RowSw | RowSwColSw | RowUltra. ",
        cxxopts::value<std::string>())(
        "s, sip",
        "Use option sip.",
        cxxopts::value<bool>())(
        "S, scip",
        "Use option SCIP (to test utilities)",
        cxxopts::value<bool>())(
        "m, pami",
        "Use pami. Cutoff optional double value.",
        cxxopts::value<double>())(
        "t, partition",
        "Use pami with partition file: filename",
        cxxopts::value<std::string>())(
        "t, time-limit",
        "Use time limit.",
        cxxopts::value<double>());
        
        
        
        ("help", "Print help.");

    cxx_options.parse_positional("filename");

    auto result = cxx_options.parse(argc, argv);

    if (result.count("help")) {
      std::cout << cxx_options.help({""}) << std::endl;
      exit(0);
    }

    if (result.count("filename")) {
      std::string filenames = "";
      std::cout << "filename = {";
      auto& v = result["filename"].as<std::vector<std::string>>();
      for (const auto& s : v) {
        std::cout << s << ", ";
        filenames = filenames + " " + s;
      }
      std::cout << "}" << std::endl;
      highs_options.setValue("filenames", filenames);
    }

    if (result.count("crash")) {
      std::string data = result["crash"].as<std::string>();
      std::transform(data.begin(), data.end(), data.begin(), ::tolower);
      if (data != "off" && 
          data != "ltssf" &&
          data != "ltssf1" &&
          data != "ltssf2" &&
          data != "ltssf3" &&
          data != "ltssf4" &&
          data != "ltssf5" &&
          data != "ltssf6" &&
          data != "ltssf7" &&
          data != "bs" &&
          data != "singts") {
        std::cout << "Wrong value specified for crash." << std::endl;
        std::cout << cxx_options.help({""}) << std::endl;
        exit(0);
      }
      highs_options.setValue("crash", data);
      std::cout << "Crash is set to " << data << ".\n";
    }

    if (result.count("edge-weight")) {
      std::string data = result["edge-weight"].as<std::string>();
      std::transform(data.begin(), data.end(), data.begin(), ::tolower);
      if (data != "dan" && 
          data != "dvx" &&
          data != "dse" &&
          data != "dse0" &&
          data != "dse2dvx") {
        std::cout << "Wrong value specified for edge-weight." << std::endl;
        std::cout << cxx_options.help({""}) << std::endl;
        exit(0);
      }
      highs_options.setValue("edge-weight", data);
      std::cout << "Edge weight is set to " << data << ".\n";
    }

    if (result.count("price")) {
      std::string data = result["price"].as<std::string>();
      std::transform(data.begin(), data.end(), data.begin(), ::tolower);
      if (data != "row" && 
          data != "col" &&
          data != "rowsw" &&
          data != "rowswcolsw" &&
          data != "rowultra") {
        std::cout << "Wrong value specified for price." << std::endl;
        std::cout << cxx_options.help({""}) << std::endl;
        exit(0);
      }
      highs_options.setValue("price", data);
      std::cout << "Price is set to " << data << ".\n";
    }

    if (result.count("presolve")) {
      std::string data = result["presolve"].as<std::string>();
      std::transform(data.begin(), data.end(), data.begin(), ::tolower);
      if (data != "on" && data != "off") {
        std::cout << "Wrong value specified for presolve." << std::endl;
        std::cout << cxx_options.help({""}) << std::endl;
        exit(0);
      }
      if (data == "off") highs_options.setValue("presolve", false);
      else highs_options.setValue("presolve", true);
      std::cout << "Presolve is set to " << data << ".\n";
    }

    if (result.count("time-limit")) {
      double time_limit = result["time-limit"].as<double>();
      if (time_limit <= 0) {
        std::cout << "Time limit must be positive." << std::endl;
        std::cout << cxx_options.help({""}) << std::endl;
        exit(0);
      }
      highs_options.setValue("time-limit", time_limit);
      std::cout << "Time limit is set to " << time_limit << ".\n";
    }

    if (result.count("partition")) {
      std::string data = result["presolve"].as<std::string>();
      std::transform(data.begin(), data.end(), data.begin(), ::tolower);
      if (data != "on" && data != "off") {
        std::cout << "Wrong value specified for presolve." << std::endl;
        std::cout << cxx_options.help({""}) << std::endl;
        exit(0);
      }
      highs_options.setValue("presolve", true);
      std::cout << "Presolve is set to " << data << ".\n";
    }

    if (result.count("presolve")) {
      std::string data = result["presolve"].as<std::string>();
      std::transform(data.begin(), data.end(), data.begin(), ::tolower);
      if (data != "on" && data != "off") {
        std::cout << "Wrong value specified for presolve." << std::endl;
        std::cout << cxx_options.help({""}) << std::endl;
        exit(0);
      }
      highs_options.setValue("presolve", true);
      std::cout << "Presolve is set to " << data << ".\n";
    }

    if (result.count("presolve")) {
      std::string data = result["presolve"].as<std::string>();
      std::transform(data.begin(), data.end(), data.begin(), ::tolower);
      if (data != "on" && data != "off") {
        std::cout << "Wrong value specified for presolve." << std::endl;
        std::cout << cxx_options.help({""}) << std::endl;
        exit(0);
      }
      highs_options.setValue("presolve", true);
      std::cout << "Presolve is set to " << data << ".\n";
    }

  } catch (const cxxopts::OptionException& e) {
    std::cout << "error parsing options: " << e.what() << std::endl;
    return HighsStatus::OptionsError;
  }

  exit(0);
}

HighsStatus loadOptions(int argc, char** argv, HighsOptions& options_) {
  // todo: replace references with options_.*
  int filename = 0;
  int presolve = 0;
  int crash = 0;
  int edgeWeight = 0;
  int price = 0;
  int pami = 0;
  int sip = 0;
  int scip = 0;
  int timeLimit = 0;

  double cut = 0;
  const char* fileName = "";
  const char* presolveMode = "";
  const char* edWtMode = "";
  const char* priceMode = "";
  const char* crashMode = "";
  const char* partitionFile = "";

  double TimeLimit_ArgV = HSOL_CONST_INF;

  if (argc == 1) {
    std::cout << "Error: No file specified. \n" << std::endl;
    printHelp(argv[0]);
    return HighsStatus::OptionsError;
  }

  if (argc == 4 && strcmp(argv[1], "-repeat") == 0) {
#ifdef HiGHSDEV
    HTester tester;
    tester.setup(argv[2]);
    tester.testUpdate(atoi(argv[3]));
#endif
    return HighsStatus::OK;
  }

  char opt;
  if (argc == 2) {
    filename = 1;
    fileName = argv[1];
  } else {
    while ((opt = getopt(argc, argv, "p:c:e:P:sSm::t:T:df:")) != EOF)
      switch (opt) {
        case 'f':
          filename = 1;
          cout << "Reading file " << optarg << endl;
          fileName = optarg;
          break;
        case 'p':
          presolveMode = optarg;
          if (presolveMode[0] == 'O' && presolveMode[1] == 'n')
            presolve = 1;
          else if (presolveMode[0] == 'E' && presolveMode[1] == 'x')
            presolve = 2;
          else
            presolve = 0;
          cout << "Presolve is set to " << optarg << endl;
          break;
        case 's':
          sip = 1;
          break;
        case 'S':
          scip = 1;
          break;
        case 'm':
          pami = 1;
          if (optarg) {
            cut = atof(optarg);
            cout << "Pami cutoff = " << cut << endl;
          }
          break;
        case 'c':
          crash = 1;
          crashMode = optarg;
          cout << "Crash is set to " << optarg << endl;
          break;
        case 'e':
          edgeWeight = 1;
          edWtMode = optarg;
          cout << "Edge weight is set to " << optarg << endl;
          break;
        case 'P':
          price = 1;
          priceMode = optarg;
          cout << "Price is set to " << optarg << endl;
          break;
        case 't':
          partitionFile = optarg;
          cout << "Partition file is set to " << optarg << endl;
          break;
        case 'T':
          timeLimit = 1;
          TimeLimit_ArgV = atof(optarg);
          cout << "Time limit is set to " << optarg << endl;
          break;
        case '?':
          if (opt == 'p')
            fprintf(stderr,
                    "Option -%c requires an argument. Current options: Off "
                    "On \n",
                    opt);
          if (opt == 'c')
            fprintf(stderr,
                    "Option -%c requires an argument. Current options: Off "
                    "LTSSF LTSSF1 LTSSF2 LTSSF3 LTSSF4 LTSSF5 LTSSF6 \n",
                    opt);
          if (opt == 'e')
            fprintf(stderr,
                    "Option -%c requires an argument. Current options: Dan Dvx "
                    "DSE DSE0 DSE2Dvx\n",
                    opt);
          if (opt == 'P')
            fprintf(stderr,
                    "Option -%c requires an argument. Current options: Row Col "
                    "RowSw RowSwColSw\n",
                    opt);
          else
            printHelp(argv[0]);
        default:
          cout << endl;
          abort();
      }
  }

  // Set defaults
  if (!filename) {
    std::cout << "No file specified. " << std::endl;
    return HighsStatus::OptionsError;
  }

  if (!presolve) {
    presolveMode = "Off";
    printf("Setting default value presolveMode = %s\n", presolveMode);
  }

  if (!crash) {
    crashMode = "Off";
    printf("Setting default value crashMode = %s\n", crashMode);
  }

  if (!edgeWeight) {
    edWtMode = "DSE2Dvx";
    printf("Setting default value edWtMode = %s\n", edWtMode);
  }

  if (!price) {
    priceMode = "RowSwColSw";
    printf("Setting default value priceMode = %s\n", priceMode);
  }
#ifdef HiGHSDEV
  printf(
      "HApp: sip = %d; scip = %d; pami = %d; presolve = %d;  crash = %d; "
      "edgeWeight = %d; price = %d; timeLimit = %d\n",
      sip, scip, pami, presolve, crash, edgeWeight, price, timeLimit);
#endif

  options_.filename = filename;
  options_.presolve = presolve;
  options_.crash = crash;
  options_.edgeWeight = edgeWeight;
  options_.price = price;
  options_.pami = pami;
  options_.sip = sip;
  options_.scip = scip;
  options_.timeLimit = TimeLimit_ArgV;

  options_.cut = cut;
  options_.fileName = fileName;
  options_.presolveMode = presolveMode;
  options_.edWtMode = edWtMode;
  options_.priceMode = priceMode;
  options_.crashMode = crashMode;
  options_.partitionFile = partitionFile;

  return HighsStatus::OK;
}

HighsStatus solveSimplex(const HighsOptions& opt, const HighsLp& lp,
                         HighsSolution& solution) {
  HModel model;
  model.load_fromArrays(lp.numCol_, lp.sense_, &lp.colCost_[0],
                        &lp.colLower_[0], &lp.colUpper_[0], lp.numRow_,
                        &lp.rowLower_[0], &lp.rowUpper_[0], lp.nnz_,
                        &lp.Astart_[0], &lp.Aindex_[0], &lp.Avalue_[0]);

  cout << "=================================================================="
          "=="
          "================"
       << endl;
  // parallel
  if (opt.sip) {
    cout << "Running solveTasks" << endl;

    solveTasks(model);
  }
  if (opt.scip) {
    cout << "Running solveSCIP" << endl;
    solveSCIP(model);
  } else if (opt.pami) {
    if (opt.partitionFile) {
      cout << "Running solveMulti" << endl;
      solveMulti(model, opt.partitionFile);
    } else if (opt.cut) {
      model.intOption[INTOPT_PRINT_FLAG] = 1;
      model.intOption[INTOPT_PERMUTE_FLAG] = 1;
      model.dblOption[DBLOPT_PAMI_CUTOFF] = opt.cut;

      model.scaleModel();

      HDual solver;
      cout << "Running solveCut" << endl;
      solver.solve(&model, HDUAL_VARIANT_MULTI, 8);

      model.util_reportSolverOutcome("Cut");
    } else {
      cout << "Running solvemulti" << endl;
      solveMulti(model);
    }
  }
  // serial
  else if (!opt.presolve && !opt.crash && !opt.edgeWeight && !opt.price &&
           opt.timeLimit == HSOL_CONST_INF) {
    cout << "Running solvePlain" << endl;
    int RtCod = solvePlain(model);
    if (RtCod != 0) {
      printf("solvePlain(API) return code is %d\n", RtCod);
    }
  }  // todo: remove case below, presolve handled elsewhere
  else if (opt.presolve && !opt.crash && !opt.edgeWeight && !opt.price &&
           opt.timeLimit == HSOL_CONST_INF) {
    if (opt.presolve == 1) {
      cout << "Running solvePlainWithPresolve" << endl;
      solvePlainWithPresolve(model);
    }
#ifdef EXT_PRESOLVE
    else if (presolve == 2) {
      cout << "Running solveExternalPresolve" << endl;
      solveExternalPresolve(fileName);
    }
#endif
  } else {
    cout << "Running solvePlainJAJH" << endl;
    solvePlainJAJH(model, opt.priceMode, opt.edWtMode, opt.crashMode,
                   opt.presolveMode, opt.timeLimit);
  }

  // todo: check what the solver outcome is and return corresponding status
  return HighsStatus::OK;
}

HighsLp HModelToHighsLp(const HModel& model) { return model.lp; }

HModel HighsLpToHModel(const HighsLp& lp) {
  HModel model;

  model.lp.numCol_ = lp.numCol_;
  model.lp.numRow_ = lp.numRow_;

  model.lp.Astart_ = lp.Astart_;
  model.lp.Aindex_ = lp.Aindex_;
  model.lp.Avalue_ = lp.Avalue_;
  model.lp.colCost_ = lp.colCost_;
  model.lp.colLower_ = lp.colLower_;
  model.lp.colUpper_ = lp.colUpper_;
  model.lp.rowLower_ = lp.rowLower_;
  model.lp.rowUpper_ = lp.rowUpper_;

  return model;
}

#endif
