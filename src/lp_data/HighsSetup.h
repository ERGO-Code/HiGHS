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
#include <memory>

#include "HApp.h"
#include "HighsLp.h"
#include "HighsModelObject.h"
#include "Presolve.h"
#include "cxxopts.hpp"

HModel HighsLpToHModel(const HighsLp& lp);
HighsLp HModelToHighsLp(const HModel& model);

// Class to set parameters and run HiGHS
class Highs {
 public:
  Highs() {}
  explicit Highs(const HighsOptions& opt) : options_(opt){};

  // Function to call just presolve.
  HighsPresolveStatus presolve(const HighsLp& lp, HighsLp& reduced_lp) {
    // todo: implement, from user's side.
    return HighsPresolveStatus::NullError;
  };

  // The public method run(lp, solution) calls runSolver to solve problem before
  // or after presolve (or crash later?) depending on the specified options.
  HighsStatus run(HighsLp& lp, HighsSolution& solution);
  HighsOptions options_;

 private:
  // each HighsModelObject holds a const ref to its lp_
  std::vector<HighsModelObject> lps_;

  HighsPresolveStatus runPresolve(PresolveInfo& presolve_info);
  HighsPostsolveStatus runPostsolve(PresolveInfo& presolve_info);
  HighsStatus runSolver(HighsModelObject& model);
};

// Checks the options calls presolve and postsolve if needed. Solvers are called
// with runSolver(..)
HighsStatus Highs::run(HighsLp& lp, HighsSolution& solution) {
  // todo: handle printing messages with HighsPrintMessage

  // Not solved before, so create an instance of HighsModelObject.
  lps_.push_back(HighsModelObject(lp));

  //Define clocks
  HighsTimer timer;
  int presolveClock = timer.clockDef("Presolve", "Pre");
  int scaleClock = timer.clockDef("Scale", "Scl");
  int crashClock = timer.clockDef("Crash", "Csh");
  int solveClock = timer.clockDef("Solve", "Slv");
  int postsolveClock = timer.clockDef("Postsolve", "Pst");
  //  timer.reset();

  // Presolve. runPresolve handles the level of presolving (0 = don't presolve).
  timer.start(presolveClock);

  PresolveInfo presolve_info(options_.presolveMode, lp);
  HighsPresolveStatus presolve_status = runPresolve(presolve_info);

  timer.stop(presolveClock);
 
  // Run solver.
  timer.start(solveClock);
  HighsStatus solve_status = HighsStatus::Init;
  switch (presolve_status) {
    case HighsPresolveStatus::NotReduced: {
      solve_status = runSolver(lps_[0]);
      break;
    }
    case HighsPresolveStatus::Reduced: {
      HighsLp& reduced_lp = presolve_info.getReducedProblem();
      // Add reduced lp object to vector of HighsModelObject,
      // so the last one in lp_ is the presolved one.
      lps_.push_back(HighsModelObject(reduced_lp));
      solve_status = runSolver(lps_[1]);
      break;
    }
    case HighsPresolveStatus::ReducedToEmpty: {
      // Proceed to postsolve.
      break;
    }
    case HighsPresolveStatus::Infeasible:
    case HighsPresolveStatus::Unbounded: {
      HighsStatus result = (presolve_status == HighsPresolveStatus::Infeasible) ?
               HighsStatus::Infeasible : HighsStatus::Unbounded;

      std::cout << "Problem status detected on presolve: "
                << HighsStatusToString(result);
      return result;
    }
    default: {
      // case HighsPresolveStatus::Error:a
      std::cout << "Error during presolve.";
      return HighsStatus::PresolveError;
    }
  }
  timer.stop(solveClock);

  timer.start(postsolveClock);
  // Postsolve. Does nothing if there were no reductions during presolve.
  if (solve_status == HighsStatus::Optimal) {
    if (presolve_status == HighsPresolveStatus::Reduced) {
      presolve_info.reduced_solution_ = lps_[1].solution_;
      presolve_info.presolve_[0].setBasisInfo(
          lps_[1].basis_info_.basis_index, lps_[1].basis_info_.nonbasic_flag,
          lps_[1].basis_info_.nonbasic_move);
    }

    HighsPostsolveStatus postsolve_status = runPostsolve(presolve_info);
    if (postsolve_status == HighsPostsolveStatus::SolutionRecovered) {
      std::cout << "Postsolve finished.\n";

      // Set solution and basis info for simplex clean up.
      // Original LP is in lp_[0] so we set the basis information there.
      lps_[0].basis_info_.basis_index =
          presolve_info.presolve_[0].getBasisIndex();
      lps_[0].basis_info_.nonbasic_flag =
          presolve_info.presolve_[0].getNonbasicFlag();
      lps_[0].basis_info_.nonbasic_move =
          presolve_info.presolve_[0].getNonbasicMove();

      options_.clean_up = true;

      solve_status = runSolver(lps_[0]);
    }
  }
  timer.stop(postsolveClock);

  if (solve_status != HighsStatus::Optimal) {
    if (solve_status == HighsStatus::Infeasible ||
        solve_status == HighsStatus::Unbounded) {
      if (options_.presolveMode == "on") {
        std::cout << "Reduced problem status: "
                  << HighsStatusToString(solve_status);
        // todo: handle case. Try to solve again with no presolve?
        return HighsStatus::NotImplemented;
      } else {
        std::cout << "Solver terminated with a non-optimal status: "
                  << HighsStatusToString(solve_status) << std::endl;
        lps_[0].hmodel_[0].intOption[INTOPT_PRINT_FLAG] = 1;
        lps_[0].hmodel_[0].util_reportSolverOutcome("Run");
      }
    }
  } else {
    // Report in old way so tests pass.
    lps_[0].hmodel_[0].intOption[INTOPT_PRINT_FLAG] = 1;
    lps_[0].hmodel_[0].util_reportSolverOutcome("Run");
  }
  // Report times
  std::vector<int> clockList{presolveClock, scaleClock, crashClock, solveClock, postsolveClock};
  timer.report(clockList);

  return HighsStatus::OK;
}

HighsPresolveStatus Highs::runPresolve(PresolveInfo& info) {
  if (options_.presolveMode != "on") return HighsPresolveStatus::NotReduced;

  if (info.lp_ == nullptr) return HighsPresolveStatus::NullError;

  if (info.presolve_.size() == 0) return HighsPresolveStatus::NotReduced;

  info.presolve_[0].load(*(info.lp_));

  // Initialize a new presolve class instance for the LP given in presolve info
  return info.presolve_[0].presolve();
}

HighsPostsolveStatus Highs::runPostsolve(PresolveInfo& info) {
  if (info.presolve_.size() != 0) {
    bool solution_ok =
        isSolutionConsistent(info.getReducedProblem(), info.reduced_solution_);
    if (!solution_ok)
      return HighsPostsolveStatus::ReducedSolutionDimenionsError;

    // todo: error handling + see todo in run()
    info.presolve_[0].postsolve(info.reduced_solution_,
                                info.recovered_solution_);

    return HighsPostsolveStatus::SolutionRecovered;
  } else {
    return HighsPostsolveStatus::NoPostsolve;
  }
}

// The method below runs simplex or ipx solver on the lp.
HighsStatus Highs::runSolver(HighsModelObject& model) {
  assert(checkLp(model.lp_) == HighsInputStatus::OK);

  HighsStatus status = HighsStatus::Init;
#ifndef IPX
  // HiGHS
  // todo: Without the presolve part, so will be
  //     = solve_simplex(options, reduced_lp, reduced_solution)
  status = runSimplexSolver(options_, model);
#else
  // IPX
  // todo:Check options for simplex-specific options
  // use model.lp_, model.solution_ and model.hmodel_ remains empty.
  status = runIpxSolver(options_, lp, solution);
  // If ipx crossover did not find optimality set up simplex.

#endif

  if (status != HighsStatus::Optimal) return status;

  // Check.
  if (!isSolutionConsistent(model.lp_, model.solution_)) {
    std::cout << "Error: Inconsistent solution returned from solver.\n";
  }

  // todo:
  // assert(KktSatisfied(lp, solution));

  return status;
}

void HiGHSRun(const char* message = nullptr) {
  std::cout << "Running HiGHS " << HIGHS_VERSION_MAJOR << "."
            << HIGHS_VERSION_MINOR << "." << HIGHS_VERSION_PATCH
            << " [date: " << HIGHS_COMPILATION_DATE
            << ", git hash: " << HIGHS_GITHASH << "]"
            << "\n"
            << "Copyright (c) 2018 ERGO-Code under MIT licence terms.\n\n";
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

#endif
};

HighsStatus loadOptions(int argc, char** argv, HighsOptions& options) {
  try {
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

    if (result.count("help")) {
      std::cout << cxx_options.help({""}) << std::endl;
      exit(0);
    }

    // Currently works for only one filename at a time.
    if (result.count("filename")) {
      std::string filenames = "";
      auto& v = result["filename"].as<std::vector<std::string>>();
      if (v.size() > 1) {
        std::cout << "Multiple files not implemented yet.\n";
        return HighsStatus::LpError;
      }
      for (const auto& s : v) {
        filenames = filenames + s;
      }
      options.filenames = filenames;
    }

    if (result.count("crash")) {
      std::string data = result["crash"].as<std::string>();
      std::transform(data.begin(), data.end(), data.begin(), ::tolower);
      if (data != "off" && data != "ltssf" && data != "ltssf1" &&
          data != "ltssf2" && data != "ltssf3" && data != "ltssf4" &&
          data != "ltssf5" && data != "ltssf6" && data != "ltssf7" &&
          data != "bs" && data != "singts") {
        std::cout << "Wrong value specified for crash." << std::endl;
        std::cout << cxx_options.help({""}) << std::endl;
        exit(0);
      }
      options.crashMode = data;
      std::cout << "Crash is set to " << data << ".\n";
    }

    if (result.count("edge-weight")) {
      std::string data = result["edge-weight"].as<std::string>();
      std::transform(data.begin(), data.end(), data.begin(), ::tolower);
      if (data != "dan" && data != "dvx" && data != "dse" && data != "dse0" &&
          data != "dse2dvx") {
        std::cout << "Wrong value specified for edge-weight." << std::endl;
        std::cout << cxx_options.help({""}) << std::endl;
        exit(0);
      }
      options.edWtMode = data;
      std::cout << "Edge weight is set to " << data << ".\n";
    }

    if (result.count("price")) {
      std::string data = result["price"].as<std::string>();
      std::transform(data.begin(), data.end(), data.begin(), ::tolower);
      if (data != "row" && data != "col" && data != "rowsw" &&
          data != "rowswcolsw" && data != "rowultra") {
        std::cout << "Wrong value specified for price." << std::endl;
        std::cout << cxx_options.help({""}) << std::endl;
        exit(0);
      }
      options.priceMode = data;
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
      options.presolveMode = data;
      std::cout << "Presolve is set to " << data << ".\n";
    }

    if (result.count("time-limit")) {
      double time_limit = result["time-limit"].as<double>();
      if (time_limit <= 0) {
        std::cout << "Time limit must be positive." << std::endl;
        std::cout << cxx_options.help({""}) << std::endl;
        exit(0);
      }
      options.timeLimit = time_limit;
    }

    if (result.count("partition")) {
      std::string data = result["partition"].as<std::string>();
      std::transform(data.begin(), data.end(), data.begin(), ::tolower);
      //      highs_options.setValue("partition", data);
      std::cout << "Partition is set to " << data << ".\n";
    }

    if (result.count("sip")) {
      options.sip = true;
      std::cout << "Option sip enabled."
                << ".\n";
    }

    if (result.count("scip")) {
      options.scip = true;
      std::cout << "Option scip enabled."
                << ".\n";
    }

    if (result.count("pami")) {
      options.pami = true;
      std::cout << "Option pami enabled (parallel solve).\n";
    }

  } catch (const cxxopts::OptionException& e) {
    std::cout << "error parsing options: " << e.what() << std::endl;
    return HighsStatus::OptionsError;
  }

  if (options.filenames.size() == 0) {
    std::cout << "Please specify filename in .mps or .gz format.\n";
    return HighsStatus::LpError;
  }

  return HighsStatus::OK;
}

#endif
