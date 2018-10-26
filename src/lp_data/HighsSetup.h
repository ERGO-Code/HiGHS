#ifndef HIGHS_SETUP_
#define HIGHS_SETUP_

#include <iostream>

#include "LpData.h"
#include "HApp.h"

// Class to set parameters and run HiGHS
class Highs {
 public:
  Highs() {}
  explicit Highs(const Options& opt) : options(opt){};

  // The public method run(lp, solution) calls runSolver to solve problem before
  // or after presolve (or crash later?) depending on the specified options.
  Status run(const LpData& lp, Solution& solution) const;

  void setAllOptions(const Options& opt) { options = opt; }
  // todo: implement string based options
  void setOption(const std::string& option, int value);
  void setOption(const std::string& option, double value);
  void setOption(const std::string& option, string value);

  int getIntOption(const std::string& option);
  double getDoubleOption(const std::string& option);
  std::string getStringOption(const std::string& option);

 private:
  Options options;
  Status runSolver(const LpData& lp, Solution& solution) const;
};

// Checks the options calls presolve and postsolve if needed. Solvers are called
// with runSolver(..)
Status Highs::run(const LpData& lp, Solution& solution) const {
  if (!options.presolve) {
    Solution solution;
    return runSolver(lp, solution);
  }

  // todo
  return Status::NotImplemented;

  /*
   LpData reduced_lp;

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

  return Status::OK;
}

// The method below runs simplex or ipx solver on the lp.

Status Highs::runSolver(const LpData& lp, Solution& solution) const {
  assert(checkLp(lp) == LpError::none);

  Status status;
#ifndef IPX
  // HiGHS
  // todo: Without the presolve part, so will be
  //     = solve_simplex(options, reduced_lp, reduced_solution)
  status = solveSimplex(options, lp, solution);
#else
  // IPX
  // todo:Check options for simplex-specific options

  status = solveIpx(options, lp, solution);
  // If ipx crossover did not find optimality set up simplex.

#endif
  
  // todo: 
  // assert(KktSatisfied(lp, solution));
  
  return status;
}

// Parses the file in options.filename using the parser specified in
// options.parser
Status loadLpFromFile(const Options& options, LpData& lp) {
  // Check file exists.

  // Which parser
  // will be (options.getValue("parser") == MpsParser::new)

  // Initialize arrays

  // call MPSParser::loadProblem(arrays of LpData object)

  assert(checkLp(lp) == LpError::none);
  return Status::OK;
}

Status loadOptions(int argc, char** argv, Options& options_) {
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
  ;

  double cut = 0;
  const char* fileName = "";
  const char* presolveMode = "";
  const char* edWtMode = "";
  const char* priceMode = "";
  const char* crashMode = "";
  const char* partitionFile = "";

  double TimeLimit_ArgV = HSOL_CONST_INF;
  std::cout << "Running HiGHS " << HIGHS_VERSION_MAJOR << "."
            << HIGHS_VERSION_MINOR << "." << HIGHS_VERSION_PATCH
            << " [date: " << HIGHS_COMPILATION_DATE
            << ", git hash: " << HIGHS_GITHASH << "]"
            << "\n"
            << "Copyright (c) 2018 ERGO-Code under MIT licence terms\n\n";
#ifdef HiGHSDEV
  // Report on preprocessing macros

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

  if (argc == 1) {
    std::cout << "Error: No file specified. \n" << std::endl;
    printHelp(argv[0]);
    return Status::InputError;
  }

  if (argc == 4 && strcmp(argv[1], "-repeat") == 0) {
#ifdef HiGHSDEV
    HTester tester;
    tester.setup(argv[2]);
    tester.testUpdate(atoi(argv[3]));
#endif
    return Status::OK;
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
#ifdef HiGHSDEV
    fileName = "ml.mps";
    printf("Setting default value filenameMode = %s\n", fileName);
#else
    std::cout << "No file specified. " << std::endl;
    return Status::FileNotFound;
#endif
  }
  // Check if file exists
  else if (access(fileName, F_OK) == -1) {
    return Status::FileNotFound;
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

  return Status::OK;
}

#endif