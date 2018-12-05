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

#include <iostream>

#include "HApp.h"
#include "HighsLp.h"
#include "HighsOptions.h"
#include "Presolve.h"
#include "HighsModelObject.h"
#include "cxxopts.hpp"

// Class to set parameters and run HiGHS
class Highs {
 public:
  Highs() {}
  explicit Highs(const HighsOptions& opt) : options_(opt){};
  explicit Highs(const HighsStringOptions& opt) : options__(opt){};

  // Function to call just presolve. This method does not modify anything
  HighsPresolveStatus presolve(const HighsLp& lp, HighsLp& reduced_lp) const;

  HighsPresolveStatus runPresolve(PresolveInfo& presolve_info);
  HighsPostsolveStatus runPostsolve(PresolveInfo& presolve_info);
  // The public method run(lp, solution) calls runSolver to solve problem before
  // or after presolve (or crash later?) depending on the specified options.
  HighsStatus run(const HighsLp& lp, HighsSolution& solution);

 private:
  // each HighsModelObject holds a const ref to its lp_
  std::vector<HighsModelObject> lps_;

  // delete.
  HighsOptions options_;
  // use HighsStringOptions instead for now. Then rename to HighsOptions, once
  // previous one is gone.
  HighsStringOptions options__;

 public:
  // The public method run(lp, solution) calls runSolver to solve problem before
  // or after presolve (or crash later?) depending on the specified options.
  HighsStatus run(const HighsLp& lp, HighsSolution& solution) const;

  HighsStatus runSolver(HighsModelObject& model);
};

// Checks the options calls presolve and postsolve if needed. Solvers are called
// with runSolver(..)
HighsStatus Highs::run(const HighsLp& lp, HighsSolution& solution) const {
  // Make sure lp is being solved for the first time.
  if (lps_.size() > 0)
    return HighsStatus::NotImplemented;

  // todo: handle printing messages with HighsPrintMessage
/*
  // Not solved before, so create an instance of HighsModelObject.
  lps_.push_back(HighsModelObject(lps));

  // Presolve. runPresolve handles the level of presolving (0 = don't presolve).
  PresolveInfo presolve_info(options_.presolve, lp);
//  HighsPresolveStatus presolve_status = runPresolve(presolve_info);
  HighsPresolveStatus presolve_status = HighsPresolveStatus::NotReduced;

  // Run solver.
  switch (presolve_status) {
    case HighsPresolveStatus::NotReduced: {
      return runSolver(lp_[0]);
      break;
    }
    case HighsPresolveStatus::Reduced: {
      const HighsLp& reduced_lp = presolve_info.getReducedProblem();
      // Add reduced lp object to vector of HighsModelObject,
      // so the last one in lp_ is the presolved one.
      lp_.push_back(HighsModelObject(reduced_lp));
      runSolver(lp_[lp._.size() - 1], presolve_info.reduced_solution_);
      break;
    }
    case HighsPresolveStatus::ReducedToEmpty: {
      // Proceed to postsolve.
      break;
    }
    case HighsPresolveStatus::Infeasible:
    case HighsPresolveStatus::Unbounded: {
      // todo: report solver outcome.
      break;
    }
    default: {
      // case HighsPresolveStatus::Error:
      // todo: handle error.
      break;
    }
  }

  // Postsolve. Does nothing if there were no reductions during presolve.
  // HighsPostsolveStatus postsolve_status = runPostsolve(presolve_info);
  HighsPostsolveStatus postsolve_status = HighsPostsolveStatus::SolutionRecovered;
  if (postsolve_status == HighsPostsolveStatus::SolutionRecovered) {
    // todo: add finishing simplex iterations if needed.
  } else {
    // todo: handle postsolve errors.
  }
*/
  return HighsStatus::OK;
}

HighsPresolveStatus Highs::runPresolve(PresolveInfo& info) {
  if (info.presolve_ == nullptr || info.lp_ == nullptr)
    return HighsPresolveStatus::NullError;
/*
  // Initialize a new presolve class instance for the LP given in presolve info
  HPresolve* pre = new HPresolve();
  int status = pre->presolve();

  switch(status) {
    case presolve_.stat::Empty:
      info.status_ = HighsPresolveStatus::Empty;
      break;
    case presolve_.stat::Infeasible:
      info.status_ = HighsPresolveStatus::Infeasible;
      break;
    case presolve_.stat::Unbounded:
      info.status_ = HighsPresolveStatus::Unbounded;
      break;
    case presolve_.stat::Unset:
      info.status_ = HighsPresolveStatus::NotReduced;
      break;
  }
*/
  return info.presolve_status_;
}

HighsPostsolveStatus Highs::runPostsolve(PresolveInfo& info) {
/*
  if (info.presolve_ != nullptr && lp_ != nullptr) 
  {
    if (info.reduced_solution_.size() == 0)
      return HighsPostsolveStatus::ReducedSolutionEmpty; 
    if (!isSolutionConsistent(*info.lp_, info.reduced_solution_)
      return HighsPostsolveStatus::ReducedSolutionDimenionsError;
    
    // todo: error handling + see todo in run()
    pre->postsolve();

    return HighsPostsolveStatus::SolutionRecovered;
  }
  return HighsPostsolveStatus:LpOrPresolveObjectMissing;
*/
}

// The method below runs simplex or ipx solver on the lp.
HighsStatus Highs::runSolver(HighsModelObject& model) {

  // assert(checkLp(lp) == LpError::none);

  HighsStatus status;
#ifndef IPX
  // HiGHS
  // todo: Without the presolve part, so will be
  //     = solve_simplex(options, reduced_lp, reduced_solution)
  status = solveSimplex(options_, model);
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
    cxx_options.positional_help("[optional args]").show_positional_help();

    cxx_options.add_options()("p, presolve", "presolve",
                              cxxopts::value<bool>())(
        "f, filename", "Filename(s) of LPs to solve",
        cxxopts::value<std::vector<std::string>>())("help", "Print help.");

    cxx_options.parse_positional("file");

    auto result = cxx_options.parse(argc, argv);

    if (result.count("help")) {
      std::cout << cxx_options.help({""}) << std::endl;
      exit(0);
    }

    if (result.count("filename")) {
      std::cout << "filename = {";
      auto& v = result["filename"].as<std::vector<std::string>>();
      for (const auto& s : v) {
        std::cout << s << ", ";
      }
      std::cout << "}" << std::endl;
    }

    if (result.count("presolve")) {
      highs_options.setValue("presolve", true);
      std::cout << "Presolve is set to on through cxx options.";
    }

  } catch (const cxxopts::OptionException& e) {
    std::cout << "error parsing options: " << e.what() << std::endl;
    return HighsStatus::OptionsError;
  }
  return HighsStatus::OK;
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
    //printHelp(argv[0]);
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
          //else
          //  printHelp(argv[0]);
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

#endif
