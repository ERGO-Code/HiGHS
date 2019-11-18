/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/Highs.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "Highs.h"

#include <algorithm>
#include <iostream>
#include <memory>
#include <sstream>

#include "HConfig.h"
#include "io/Filereader.h"
#include "io/HighsIO.h"
#include "io/LoadOptions.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelUtils.h"
#include "lp_data/HighsStatus.h"
#include "presolve/FindFeasibility.h"
#include "presolve/Presolve.h"
#include "simplex/HApp.h"
#include "simplex/HighsSimplexInterface.h"

// until add_row_.. functions are moved to HighsLpUtils.h
#include "simplex/HSimplex.h"

#ifdef IPX_ON
#include "interior_point/IpxWrapper.h"
#else
#include "interior_point/IpxWrapperEmpty.h"
#endif

Highs::Highs() {
  hmos_.clear();
  HighsModelObject* hmo = new HighsModelObject(lp_, options_, timer_);
  hmos_.push_back(*hmo);
  //   hmos_.push_back(HighsModelObject(lp_, options_, timer_));
  //  allow_presolve_ = true;
}

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const bool value) {
  if (setOptionValue(option, options_.records, value) == OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const int value) {
  if (setOptionValue(option, options_.records, value) == OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const double value) {
  if (setOptionValue(option, options_.records, value) == OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const std::string value) {
  if (setOptionValue(option, options_.records, value) == OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const char* value) {
  if (setOptionValue(option, options_.records, value) == OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::readHighsOptions(const std::string filename) {
  if (filename.size() <= 0) {
    HighsLogMessage(HighsMessageType::WARNING,
                    "Empty file name so not reading options");
    return HighsStatus::Error;
  }
  options_.options_file = filename;
  if (!loadOptionsFromFile(options_)) return HighsStatus::Error;
  return HighsStatus::OK;    
}

HighsStatus Highs::passHighsOptions(const HighsOptions& options) {

  if (passOptions(options, options_) == OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::getHighsOptionValue(const std::string& option, bool& value) {
  if (getOptionValue(option, options_.records, value) == OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::getHighsOptionValue(const std::string& option, int& value) {
  if (getOptionValue(option, options_.records, value) == OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::getHighsOptionValue(const std::string& option,
                                       double& value) {
  if (getOptionValue(option, options_.records, value) == OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::getHighsOptionValue(const std::string& option,
                                       std::string& value) {
  if (getOptionValue(option, options_.records, value) == OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::writeHighsOptions(const std::string filename) {
  return reportOptionsToFile(filename, options_.records);
}

const HighsOptions& Highs::getHighsOptions() const { return options_; }

const HighsInfo& Highs::getHighsInfo() const { return info_; }

HighsStatus Highs::getHighsInfoValue(const std::string& info, int& value) {
  if (getInfoValue(info, info_.records, value) == InfoStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::getHighsInfoValue(const std::string& info, double& value) {
  if (getInfoValue(info, info_.records, value) == InfoStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::writeHighsInfo(const std::string filename) {
  return reportInfoToFile(filename, info_.records);
}

HighsStatus Highs::passModel(const HighsLp& lp) {
  // Copy the LP to the internal LP
  lp_ = lp;
  // Check validity of the LP, normalising its values (by default).
  HighsStatus return_status = assessLp(lp_, options_);
  if (return_status != HighsStatus::OK) return return_status;
  // hmos_[0] is the HighsModelObject corresponding to the original LP
  hmos_.clear();
  hmos_.push_back(HighsModelObject(lp_, options_, timer_));
  //  allow_presolve_ = true;
  return HighsStatus::OK;
}

HighsStatus Highs::readModel(const std::string filename) {
  Filereader* reader = Filereader::getFilereader(filename.c_str());
  HighsLp model;
  this->options_.model_file = filename;

  FilereaderRetcode retcode = reader->readModelFromFile(this->options_, model);
  if (retcode != FilereaderRetcode::OK) {
    return HighsStatus::Error;
  }

  return this->passModel(model);
}

HighsStatus Highs::writeModel(const std::string filename) {
  HighsLp model = this->lp_;

  if (filename == "") {
    // Empty file name: report model on stdout
    HighsLogMessage(HighsMessageType::WARNING,
                    "Empty file name so reporting model on stdout");
    reportLp(model, 2);
    return HighsStatus::Warning;
  } else {
    Filereader* writer = Filereader::getFilereader(filename.c_str());
    return writer->writeModelToFile(filename.c_str(), model);
  }
}

// Checks the options calls presolve and postsolve if needed. Solvers are called
// with runSolver(..)
HighsStatus Highs::run() {
    /*
  if (options_.message_level >= 0) {
    printf("\n!! Actually solving an LP with %d cols, %d rows", lp_.numCol_, lp_.numRow_);
    if (lp_.numCol_) printf(" and %d nonzeros", lp_.Astart_[lp_.numCol_]);
    printf(":basis.valid_ = %d: basis_.valid_ = %d: simplex_lp_status_.has_basis = %d!!\n\n",
	   basis_.valid_,
	   hmos_[0].basis_.valid_,
	   hmos_[0].simplex_lp_status_.has_basis);
    if (basis_.valid_ != hmos_[0].basis_.valid_) {
      printf("NB %d = basis_.valid_ != hmos_[0].basis_.valid_ = %d\n", basis_.valid_, hmos_[0].basis_.valid_);
    }
  }
    */
  // If running as hsol, reset any changed options
  if (options_.run_as_hsol) setHsolOptions(options_);
  // Initialise the HiGHS model status values
  model_status_ = HighsModelStatus::NOTSET;
  scaled_model_status_ = HighsModelStatus::NOTSET;
  
#ifdef HIGHSDEV
  // Shouldn't have to check validity of the LP since this is done when it is
  // loaded or modified
  //  bool normalise = true;
  HighsStatus return_status = assessLp(lp_, options_);  //, normalise);
  assert(return_status == HighsStatus::OK);
  if (return_status != HighsStatus::OK) return HighsStatus::Error;
#endif

  // For the moment runFeasibility as standalone.
  if (options_.find_feasibility) {
    // use when you do something with solution depending on whether we have
    // dualized or not.
    // HighsSolution& solution = solution_;

    // options_.message_level = HighsPrintMessageLevel::ML_DETAILED;
    // HighsSetIO(options_);

    // Add slacks and make sure a minimization problem is passed to
    // runFeasibility.
    HighsLp primal = transformIntoEqualityProblem(lp_);
    if (options_.feasibility_strategy_dualize) {
      // Add slacks & dualize.
      HighsLp dual = dualizeEqualityProblem(primal);
      // dualizeEqualityProblem returns a minimization problem.
      passModel(dual);
    } else {
      // If maximization, minimize before calling runFeasibility.
      if (primal.sense_ != OBJSENSE_MINIMIZE) {
        for (int col = 0; col < primal.numCol_; col++)
          primal.colCost_[col] = -primal.colCost_[col];
      }
      passModel(primal);
    }

    if (options_.feasibility_strategy ==
        FEASIBILITY_STRATEGY_kApproxComponentWise)
      return runFeasibility(lp_, solution_, MinimizationType::kComponentWise);
    else if (options_.feasibility_strategy == FEASIBILITY_STRATEGY_kApproxExact)
      return runFeasibility(lp_, solution_, MinimizationType::kExact);
    else if (options_.feasibility_strategy ==
             FEASIBILITY_STRATEGY_kDirectSolve) {
      // Proceed to normal exection of run().
      // If dualize has been called replace LP is replaced with dual in code
      // above.
    }
  }

  // Return immediately if the LP has no columns
  if (!lp_.numCol_) {
    hmos_[0].model_status_ = HighsModelStatus::MODEL_EMPTY;
    model_status_ = hmos_[0].model_status_;
    return highsStatusFromHighsModelStatus(hmos_[0].model_status_);
  }

  HighsSetIO(options_);
#ifdef HiGHSDEV
  if (checkOptions(options_.records) != OptionStatus::OK) return HighsStatus::Error;
#endif
  // Report all the options to an options file
  //  reportOptionsToFile("Highs.set", options_.records);
  // Report all the options as HTML
  //  reportOptionsToFile("Highs.html", options_.records);
  // Possibly report options settings
  reportOptions(stdout, options_.records);  //, true);
  HighsPrintMessage(ML_VERBOSE, "Solving %s", lp_.model_name_.c_str());

  // IPX with no presolve yet.
  if (options_.solver == "ipm") {
    int ipx_iteration_count = 0;
    return callRunSolver(hmos_[0], ipx_iteration_count, "IPX");
  }
  if (options_.mip) return runBnb();

  // Running as LP solver: start the HiGHS clock unless it's already running
  bool run_highs_clock_already_running = timer_.runningRunHighsClock();
  if (!run_highs_clock_already_running) timer_.startRunHighsClock();
  // Record the initial time and zero the overall iteration count
  double initial_time = timer_.readRunHighsClock();
  int solve_iteration_count = 0;
  int postsolve_iteration_count = 0;
  // Define identifiers to refer to the HMO of the original LP
  // (0) and the HMO created when using presolve (1)
  const int original_hmo = 0;
  const int presolve_hmo = 1;
  // Keep track of the hmo that is the most recently solved. By default it's the
  // original LP
  int solved_hmo = original_hmo;
  // Initial solve. Presolve, choose solver (simplex, ipx), postsolve.
  int iteration_count;
  //  printf("\nHighs::run() 1: basis_.valid_ = %d\n",
  //  basis_.valid_);fflush(stdout);
  if (!basis_.valid_) {
    // No HiGHS basis so consider presolve
    // Presolve. runPresolve handles the level of presolving (0 = don't
    // presolve).
    timer_.start(timer_.presolve_clock);
    PresolveInfo presolve_info(options_.presolve, lp_);
    HighsPresolveStatus presolve_status = runPresolve(presolve_info);
    timer_.stop(timer_.presolve_clock);
    //    printf("\nHighs::run() 2: presolve status = %d\n",
    //    (int)presolve_status);fflush(stdout);

    // Run solver.
    switch (presolve_status) {
      case HighsPresolveStatus::NotPresolved: {
        hmos_[solved_hmo].lp_.lp_name_ = "Original LP";
        HighsStatus return_status =
            callRunSolver(hmos_[solved_hmo], iteration_count,
                          "Not presolved: solving the LP");
        solve_iteration_count += iteration_count;
        if (return_status != HighsStatus::OK) return return_status;
        break;
      }
      case HighsPresolveStatus::NotReduced: {
        hmos_[solved_hmo].lp_.lp_name_ = "Unreduced LP";
        HighsStatus return_status =
            callRunSolver(hmos_[solved_hmo], iteration_count,
                          "Problem not reduced by presolve: solving the LP");
        solve_iteration_count += iteration_count;
        if (return_status != HighsStatus::OK) return return_status;
        break;
      }
      case HighsPresolveStatus::Reduced: {
        HighsLp& reduced_lp = presolve_info.getReducedProblem();
        // Validate the reduced LP
        assert(assessLp(reduced_lp, options_) == HighsStatus::OK);
        // Add reduced lp object to vector of HighsModelObject,
        // so the last one in lp_ is the presolved one.

        hmos_.push_back(HighsModelObject(reduced_lp, options_, timer_));
        // Log the presolve reductions
        logPresolveReductions(hmos_[original_hmo].lp_, hmos_[presolve_hmo].lp_);
        // Record the HMO to be solved
        solved_hmo = presolve_hmo;
        hmos_[solved_hmo].lp_.lp_name_ = "Presolved LP";
        HighsStatus return_status = callRunSolver(
            hmos_[solved_hmo], iteration_count, "Solving the presolved LP");
        solve_iteration_count += iteration_count;
        if (return_status != HighsStatus::OK) return return_status;
        break;
      }
      case HighsPresolveStatus::ReducedToEmpty: {
        hmos_[0].model_status_ = HighsModelStatus::OPTIMAL;
        // Proceed to postsolve.
        break;
      }
        //	printf("\nHighs::run() 3: presolve status = %d\n",
        //(int)presolve_status);fflush(stdout);
      case HighsPresolveStatus::Infeasible:
      case HighsPresolveStatus::Unbounded: {
        if (presolve_status == HighsPresolveStatus::Infeasible) {
          hmos_[original_hmo].model_status_ =
              HighsModelStatus::PRIMAL_INFEASIBLE;
        } else {
          hmos_[original_hmo].model_status_ =
              HighsModelStatus::PRIMAL_UNBOUNDED;
        }
        HighsLogMessage(
            HighsMessageType::INFO, "Problem status detected on presolve: %s",
            highsModelStatusToString(hmos_[original_hmo].model_status_)
                .c_str());

        // Report this way for the moment. May modify after merge with
        // OSIinterface branch which has new way of setting up a
        // HighsModelObject and can support multiple calls to run(). Stop and
        // read the HiGHS clock, then work out time for this call
        if (!run_highs_clock_already_running) timer_.stopRunHighsClock();

	/*
        double lp_solve_final_time = timer_.readRunHighsClock();
        std::stringstream message_not_opt;
        message_not_opt << std::endl;
        message_not_opt << "Run status : "
                        << highsModelStatusToString(
                               hmos_[original_hmo].model_status_)
                        << std::endl;
        message_not_opt << "Time       : " << std::fixed << std::setprecision(3)
                        << lp_solve_final_time - initial_time << std::endl;

        message_not_opt << std::endl;

        HighsPrintMessage(ML_MINIMAL, message_not_opt.str().c_str());
	*/
	model_status_ = hmos_[original_hmo].model_status_;
        return HighsStatus::OK;
      }
      default: {
        // case HighsPresolveStatus::Error
        HighsPrintMessage(ML_ALWAYS, "Presolve failed.");
        if (!run_highs_clock_already_running) timer_.stopRunHighsClock();
	hmos_[original_hmo].model_status_ = HighsModelStatus::PRESOLVE_ERROR;
	model_status_ = hmos_[original_hmo].model_status_;
        return HighsStatus::Error;
      }
    }
    // Postsolve. Does nothing if there were no reductions during presolve.
    if (hmos_[solved_hmo].model_status_ == HighsModelStatus::OPTIMAL) {
      if (presolve_status == HighsPresolveStatus::Reduced ||
          presolve_status == HighsPresolveStatus::ReducedToEmpty) {
        // If presolve is nontrivial, extract the optimal solution
        // and basis for the presolved problem in order to generate
        // the solution and basis for postsolve to use to generate a
        // solution(?) and basis that is, hopefully, optimal. This is
        // confirmed or corrected by hot-starting the simplex solver
        presolve_info.reduced_solution_ = hmos_[solved_hmo].solution_;
        presolve_info.presolve_[0].setBasisInfo(
            hmos_[solved_hmo].basis_.col_status,
            hmos_[solved_hmo].basis_.row_status);
        // Run postsolve
        timer_.start(timer_.postsolve_clock);
        HighsPostsolveStatus postsolve_status = runPostsolve(presolve_info);
        timer_.stop(timer_.postsolve_clock);
        if (postsolve_status == HighsPostsolveStatus::SolutionRecovered) {
          HighsPrintMessage(ML_VERBOSE, "Postsolve finished.");
          // Set solution(?) and basis to hot-start the simplex solver
          // for the original_hmo
          hmos_[original_hmo].solution_ = presolve_info.recovered_solution_;

          hmos_[original_hmo].basis_.col_status =
              presolve_info.presolve_[0].getColStatus();
          hmos_[original_hmo].basis_.row_status =
              presolve_info.presolve_[0].getRowStatus();
          hmos_[original_hmo].basis_.valid_ = true;
          // Analyse the Highs basic solution returned from postsolve
          int report_level = -1;
#ifdef HiGHSDEV
          report_level = 1;
#endif
	  HighsSolutionParams solution_params;
	  copyToSolutionParams(solution_params,
			       hmos_[original_hmo].options_,
			       hmos_[original_hmo].simplex_info_);
	  hmos_[original_hmo].model_status_ = 
	    analyseHighsBasicSolution(hmos_[original_hmo].lp_,
				 hmos_[original_hmo].basis_,
				 hmos_[original_hmo].solution_,
				 solution_params, report_level, "after returning from postsolve");
	  copyFromSolutionParams(hmos_[original_hmo].simplex_info_, solution_params);

          // Now hot-start the simplex solver for the original_hmo
          solved_hmo = original_hmo;
          // Save the options to allow the best simplex strategy to
          // be used
          HighsOptions save_options = hmos_[solved_hmo].options_;
          HighsOptions& options = hmos_[solved_hmo].options_;
          options.simplex_strategy = SIMPLEX_STRATEGY_CHOOSE;
          // Set the message level to ML_ALWAYS so that data for
          // individual iterations are reported
          bool full_iteration_logging = false;
          if (full_iteration_logging) HighsSetMessagelevel(ML_ALWAYS);
          hmos_[solved_hmo].lp_.lp_name_ = "Postsolve LP";
          HighsStatus return_status = callRunSolver(
              hmos_[solved_hmo], iteration_count,
              "Solving the original LP from the solution after postsolve");
          postsolve_iteration_count = iteration_count;
          solve_iteration_count += iteration_count;
          // Recover the options
          options = save_options;
          // Reset the message level
          if (full_iteration_logging)
            HighsSetMessagelevel(options_.message_level);
          if (return_status != HighsStatus::OK) return return_status;
        }
      }
    } else {
      // Optimal solution of presolved problem has not been found
      // The original model inherits the solved model's status and iteration count
      hmos_[original_hmo].model_status_ = hmos_[solved_hmo].model_status_;
      hmos_[original_hmo].simplex_info_.iteration_count = hmos_[solved_hmo].simplex_info_.iteration_count;
    }
  } else {
    // The problem has been solved before so we ignore presolve/postsolve/ipx.
    solved_hmo = original_hmo;
    hmos_[solved_hmo].lp_.lp_name_ = "Re-solved LP";
    HighsStatus return_status =
        callRunSolver(hmos_[solved_hmo], iteration_count, "Re-solving the LP");
    solve_iteration_count += iteration_count;
    if (return_status != HighsStatus::OK) return return_status;
  }
  // else if (reduced problem failed to solve) {
  //   todo: handle case when presolved problem failed to solve. Try to solve
  //   again with no presolve.
  // }

  assert(solved_hmo == original_hmo);

  int hmos_size = hmos_.size();
  assert(hmos_size > 0);
  // Copy HMO solution/basis to HiGHS solution/basis: this resizes solution_ and basis_
  // ToDo: make sure the model_status values are corrected
  model_status_ = hmos_[original_hmo].model_status_;
  scaled_model_status_ = hmos_[original_hmo].scaled_model_status_;
  
  info_.objective_function_value = hmos_[original_hmo].simplex_info_.dual_objective_value;
  info_.simplex_iteration_count = 0;
  for (int k = 0; k < hmos_size; k++) {
    info_.simplex_iteration_count += hmos_[k].simplex_info_.iteration_count;
  }
  info_.primal_status = hmos_[original_hmo].simplex_info_.primal_status;
  info_.dual_status = hmos_[original_hmo].simplex_info_.dual_status;
  solution_ = hmos_[original_hmo].solution_;
  basis_ = hmos_[original_hmo].basis_;
  // Report times
  if (hmos_[original_hmo].report_model_operations_clock) {
    std::vector<int> clockList{timer_.presolve_clock, timer_.solve_clock,
                               timer_.postsolve_clock};
    timer_.report("ModelOperations", clockList);
  }
  // Stop and read the HiGHS clock, then work out time for this call
  if (!run_highs_clock_already_running) timer_.stopRunHighsClock();

  double lp_solve_final_time = timer_.readRunHighsClock();
  /*
  std::stringstream message;
  message << std::endl;
  message << "Run status : "
          << highsModelStatusToString(hmos_[solved_hmo].model_status_)
          << std::endl;
  message
      << "Iterations : "
      << solve_iteration_count  // hmos_[solved_hmo].simplex_info_.iteration_count
      << std::endl;

  if (hmos_[solved_hmo].model_status_ == HighsModelStatus::OPTIMAL)
    message << "Objective  : " << std::scientific
            << hmos_[original_hmo].simplex_info_.dual_objective_value
            << std::endl;

  message << "Time       : " << std::fixed << std::setprecision(3)
          << lp_solve_final_time - initial_time << std::endl;

  message << "Postsolve  : " << postsolve_iteration_count << std::endl;

  message << std::endl;

  HighsPrintMessage(ML_MINIMAL, message.str().c_str());
  */
  HighsPrintMessage(ML_MINIMAL, "Postsolve  : %d\n", postsolve_iteration_count);
  HighsPrintMessage(ML_MINIMAL, "Time       : %0.3g\n", lp_solve_final_time - initial_time);
  model_status_ = hmos_[original_hmo].model_status_;
  return highsStatusFromHighsModelStatus(hmos_[original_hmo].model_status_);
}

const HighsLp& Highs::getLp() const { return lp_; }

const HighsSolution& Highs::getSolution() const { return solution_; }

const HighsBasis& Highs::getBasis() const { return basis_; }

const HighsModelStatus& Highs::getModelStatus(const bool scaled_model) const {
  if (scaled_model) {
    return scaled_model_status_;
  } else {
    return model_status_;
  }
}

HighsStatus Highs::getBasicVariables(int* basic_variables) {
  if (hmos_.size() == 0) return HighsStatus::Error;
  if (!hmos_[0].simplex_lp_status_.has_basis) {
    HighsLogMessage(HighsMessageType::ERROR,
                    "No basis available in getBasicVariables");
    return HighsStatus::Error;
  }
  int numRow = hmos_[0].lp_.numRow_;
  int numCol = hmos_[0].lp_.numCol_;
  if (numRow != hmos_[0].simplex_lp_.numRow_) {
    HighsLogMessage(
        HighsMessageType::ERROR,
        "Model LP and simplex LP row dimension difference (%d-%d=%d", numRow,
        hmos_[0].simplex_lp_.numRow_, numRow - hmos_[0].simplex_lp_.numRow_);
    return HighsStatus::Error;
  }
  for (int row = 0; row < numRow; row++) {
    int var = hmos_[0].simplex_basis_.basicIndex_[row];
    if (var < numCol) {
      basic_variables[row] = var;
    } else {
      basic_variables[row] = -(1 + var - numCol);
    }
  }
  return HighsStatus::OK;
}

HighsStatus Highs::getBasisInverseRow(const int row, double* row_vector,
                                      int* row_num_nz, int* row_indices) {
  if (hmos_.size() == 0) return HighsStatus::Error;
  int numRow = hmos_[0].lp_.numRow_;
  if (row < 0 || row >= numRow) {
    HighsLogMessage(HighsMessageType::ERROR,
                    "Row index %d out of range [0, %d] in getBasisInverseRow",
                    row, numRow - 1);
    return HighsStatus::Error;
  }
  if (!hmos_[0].simplex_lp_status_.has_invert) {
    HighsLogMessage(HighsMessageType::ERROR,
                    "No invertible representation for getBasisInverseRow");
    return HighsStatus::Error;
  }
  // Compute a row i of the inverse of the basis matrix by solving B^Tx=e_i
  vector<double> rhs;
  rhs.assign(numRow, 0);
  rhs[row] = 1;
  HighsSimplexInterface simplex_interface(hmos_[0]);
  simplex_interface.basisSolve(rhs, row_vector, row_num_nz, row_indices, true);
  return HighsStatus::OK;
}

HighsStatus Highs::getBasisInverseCol(const int col, double* col_vector,
                                      int* col_num_nz, int* col_indices) {
  if (hmos_.size() == 0) return HighsStatus::Error;
  int numRow = hmos_[0].lp_.numRow_;
  if (col < 0 || col >= numRow) {
    HighsLogMessage(
        HighsMessageType::ERROR,
        "Column index %d out of range [0, %d] in getBasisInverseCol", col,
        numRow - 1);
    return HighsStatus::Error;
  }
  if (!hmos_[0].simplex_lp_status_.has_invert) {
    HighsLogMessage(HighsMessageType::ERROR,
                    "No invertible representation for getBasisInverseCol");
    return HighsStatus::Error;
  }
  // Compute a col i of the inverse of the basis matrix by solving Bx=e_i
  vector<double> rhs;
  rhs.assign(numRow, 0);
  rhs[col] = 1;
  HighsSimplexInterface simplex_interface(hmos_[0]);
  simplex_interface.basisSolve(rhs, col_vector, col_num_nz, col_indices, false);
  return HighsStatus::OK;
}

HighsStatus Highs::getBasisSolve(const double* Xrhs, double* solution_vector,
                                 int* solution_num_nz, int* solution_indices) {
  if (hmos_.size() == 0) return HighsStatus::Error;
  if (!hmos_[0].simplex_lp_status_.has_invert) {
    HighsLogMessage(HighsMessageType::ERROR,
                    "No invertible representation for getBasisSolve");
    return HighsStatus::Error;
  }
  int numRow = hmos_[0].lp_.numRow_;
  vector<double> rhs;
  rhs.assign(numRow, 0);
  for (int row = 0; row < numRow; row++) rhs[row] = Xrhs[row];
  HighsSimplexInterface simplex_interface(hmos_[0]);
  simplex_interface.basisSolve(rhs, solution_vector, solution_num_nz,
                               solution_indices, false);
  return HighsStatus::OK;
}

HighsStatus Highs::getBasisTransposeSolve(const double* Xrhs,
                                          double* solution_vector,
                                          int* solution_num_nz,
                                          int* solution_indices) {
  if (hmos_.size() == 0) return HighsStatus::Error;
  if (!hmos_[0].simplex_lp_status_.has_invert) {
    HighsLogMessage(HighsMessageType::ERROR,
                    "No invertible representation for getBasisTransposeSolve");
    return HighsStatus::Error;
  }
  int numRow = hmos_[0].lp_.numRow_;
  vector<double> rhs;
  rhs.assign(numRow, 0);
  for (int row = 0; row < numRow; row++) rhs[row] = Xrhs[row];
  HighsSimplexInterface simplex_interface(hmos_[0]);
  simplex_interface.basisSolve(rhs, solution_vector, solution_num_nz,
                               solution_indices, true);
  return HighsStatus::OK;
}

HighsStatus Highs::getReducedRow(const int row, double* row_vector,
                                 int* row_num_nz, int* row_indices) {
  if (hmos_.size() == 0) return HighsStatus::Error;
  if (row < 0 || row >= hmos_[0].lp_.numRow_) {
    HighsLogMessage(HighsMessageType::ERROR,
                    "Row index %d out of range [0, %d] in getReducedRow", row,
                    hmos_[0].lp_.numRow_ - 1);
    return HighsStatus::Error;
  }
  if (!hmos_[0].simplex_lp_status_.has_invert) {
    HighsLogMessage(HighsMessageType::ERROR,
                    "No invertible representation for getReducedRow");
    return HighsStatus::Error;
  }
  HighsLp& lp = hmos_[0].lp_;
  int numRow = lp.numRow_;
  vector<double> rhs;
  vector<double> col_vector;
  vector<int> col_indices;
  int col_num_nz;
  rhs.assign(numRow, 0);
  rhs[row] = 1;
  col_vector.resize(numRow, 0);
  col_indices.resize(numRow, 0);
  HighsSimplexInterface simplex_interface(hmos_[0]);
  // Form B^{-T}e_{row}
  simplex_interface.basisSolve(rhs, &col_vector[0], &col_num_nz,
                               &col_indices[0], true);
  bool return_indices = row_num_nz != NULL;
  if (return_indices) *row_num_nz = 0;
  for (int col = 0; col < lp.numCol_; col++) {
    double value = 0;
    for (int el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
      int row = lp.Aindex_[el];
      value += lp.Avalue_[el] * col_vector[row];
    }
    row_vector[col] = 0;
    if (fabs(value) > HIGHS_CONST_TINY) {
      if (return_indices) row_indices[(*row_num_nz)++] = col;
      row_vector[col] = value;
    }
  }
  return HighsStatus::OK;
}

HighsStatus Highs::getReducedColumn(const int col, double* col_vector,
                                    int* col_num_nz, int* col_indices) {
  if (hmos_.size() == 0) return HighsStatus::Error;
  if (col < 0 || col >= hmos_[0].lp_.numCol_) {
    HighsLogMessage(HighsMessageType::ERROR,
                    "Column index %d out of range [0, %d] in getReducedColumn",
                    col, hmos_[0].lp_.numCol_ - 1);
    return HighsStatus::Error;
  }
  if (!hmos_[0].simplex_lp_status_.has_invert) {
    HighsLogMessage(HighsMessageType::ERROR,
                    "No invertible representation for getReducedColumn");
    return HighsStatus::Error;
  }
  HighsLp& lp = hmos_[0].lp_;
  int numRow = lp.numRow_;
  vector<double> rhs;
  rhs.assign(numRow, 0);
  for (int el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++)
    rhs[lp.Aindex_[el]] = lp.Avalue_[el];
  HighsSimplexInterface simplex_interface(hmos_[0]);
  simplex_interface.basisSolve(rhs, col_vector, col_num_nz, col_indices, false);
  return HighsStatus::OK;
}

HighsStatus Highs::setSolution(const HighsSolution& solution) {
  underDevelopmentLogMessage("setSolution");
  // Check if solution is valid.
  assert((int)solution_.col_value.size() != 0 ||
         (int)solution_.col_value.size() != lp_.numCol_);
  assert((int)solution.col_dual.size() == 0 ||
         (int)solution.col_dual.size() == lp_.numCol_);
  assert((int)solution.row_dual.size() == 0 ||
         (int)solution.row_dual.size() == lp_.numRow_);

  if (solution.col_value.size()) solution_.col_value = solution.col_value;
  if (solution.col_dual.size()) solution_.col_dual = solution.col_dual;
  if (solution.row_dual.size()) solution_.row_dual = solution.row_dual;

  if (solution.col_value.size() > 0) {
    HighsStatus return_status = calculateRowValues(lp_, solution_);
    if (return_status != HighsStatus::OK) return return_status;
  }
  if (solution.row_dual.size() > 0) {
    HighsStatus return_status = calculateColDuals(lp_, solution_);
    if (return_status != HighsStatus::OK) return return_status;
  }

  return HighsStatus::OK;
}

HighsStatus Highs::setBasis(const HighsBasis& basis) {
  underDevelopmentLogMessage("setBasis");
  if (!basisOk(lp_, basis)) {
    HighsLogMessage(HighsMessageType::ERROR, "setBasis: invalid basis");
    return HighsStatus::Error;
  }
  basis_ = basis;
  basis_.valid_ = true;
  return HighsStatus::OK;
}

bool Highs::addRow(const double lower_bound, const double upper_bound,
                   const int num_new_nz, const int* indices,
                   const double* values) {
  int starts = 0;
  return addRows(1, &lower_bound, &upper_bound, num_new_nz, &starts, indices,
                 values);
}

bool Highs::addRows(const int num_new_row, const double* lower_bounds,
                    const double* upper_bounds, const int num_new_nz,
                    const int* starts, const int* indices,
                    const double* values) {
  underDevelopmentLogMessage("addRows");
  HighsStatus return_status = HighsStatus::Error;
  assert(hmos_.size() > 0);
  HighsSimplexInterface interface(hmos_[0]);
  return_status = interface.addRows(num_new_row, lower_bounds, upper_bounds,
				    num_new_nz, starts, indices, values);
  updateHighsSolutionBasis();
  return return_status == HighsStatus::OK;
}

bool Highs::addCol(const double cost, const double lower_bound,
                   const double upper_bound, const int num_new_nz,
                   const int* indices, const double* values) {
  int starts = 0;
  return addCols(1, &cost, &lower_bound, &upper_bound, num_new_nz, &starts,
                 indices, values);
}

bool Highs::addCols(const int num_new_col, const double* costs,
                    const double* lower_bounds, const double* upper_bounds,
                    const int num_new_nz, const int* starts, const int* indices,
                    const double* values) {
  underDevelopmentLogMessage("addCols");
  HighsStatus return_status = HighsStatus::Error;
  assert(hmos_.size() > 0);
  HighsSimplexInterface interface(hmos_[0]);
  return_status =
      interface.addCols(num_new_col, costs, lower_bounds, upper_bounds,
                        num_new_nz, starts, indices, values);
  updateHighsSolutionBasis();
  return return_status == HighsStatus::OK;
}

bool Highs::changeObjectiveSense(const int sense) {
  underDevelopmentLogMessage("changeObjectiveSense");
  HighsStatus return_status = HighsStatus::Error;
  assert(hmos_.size() > 0);
  HighsSimplexInterface interface(hmos_[0]);
  return_status = interface.changeObjectiveSense(sense);
  return return_status == HighsStatus::OK;
}

bool Highs::changeColCost(const int col, const double cost) {
  return changeColsCost(1, &col, &cost);
}

bool Highs::changeColsCost(const int num_set_entries, const int* set,
                           const double* cost) {
  underDevelopmentLogMessage("changeColsCost");
  HighsStatus return_status = HighsStatus::Error;
  assert(hmos_.size() > 0);
  HighsSimplexInterface interface(hmos_[0]);
  return_status = interface.changeCosts(num_set_entries, set, cost);
  return return_status == HighsStatus::OK;
}

bool Highs::changeColsCost(const int* mask, const double* cost) {
  underDevelopmentLogMessage("changeColsCost");
  HighsStatus return_status = HighsStatus::Error;
  assert(hmos_.size() > 0);
  HighsSimplexInterface interface(hmos_[0]);
  return_status = interface.changeCosts(mask, cost);
  return return_status == HighsStatus::OK;
}

bool Highs::changeColBounds(const int col, const double lower,
                            const double upper) {
  return changeColsBounds(1, &col, &lower, &upper);
}

bool Highs::changeColsBounds(const int num_set_entries, const int* set,
                             const double* lower, const double* upper) {
  underDevelopmentLogMessage("changeColsBounds");
  HighsStatus return_status = HighsStatus::Error;
  assert(hmos_.size() > 0);
  HighsSimplexInterface interface(hmos_[0]);
  return_status = interface.changeColBounds(num_set_entries, set, lower, upper);
  return return_status == HighsStatus::OK;
}

bool Highs::changeColsBounds(const int from_col, const int to_col,
                             const double* lower, const double* upper) {
  underDevelopmentLogMessage("changeColsBounds");
  HighsStatus return_status = HighsStatus::Error;
  assert(hmos_.size() > 0);
  HighsSimplexInterface interface(hmos_[0]);
  return_status = interface.changeColBounds(from_col, to_col, lower, upper);
  return return_status == HighsStatus::OK;
}

bool Highs::changeColsBounds(const int* mask, const double* lower,
                             const double* upper) {
  underDevelopmentLogMessage("changeColsBounds");
  HighsStatus return_status = HighsStatus::Error;
  assert(hmos_.size() > 0);
  HighsSimplexInterface interface(hmos_[0]);
  return_status = interface.changeColBounds(mask, lower, upper);
  return return_status == HighsStatus::OK;
}

bool Highs::changeRowBounds(const int row, const double lower,
                            const double upper) {
  return changeRowsBounds(1, &row, &lower, &upper);
}

bool Highs::changeRowsBounds(const int num_set_entries, const int* set,
                             const double* lower, const double* upper) {
  underDevelopmentLogMessage("changeRowsBounds");
  HighsStatus return_status = HighsStatus::Error;
  assert(hmos_.size() > 0);
  HighsSimplexInterface interface(hmos_[0]);
  return_status = interface.changeRowBounds(num_set_entries, set, lower, upper);
  return return_status == HighsStatus::OK;
}

bool Highs::changeRowsBounds(const int* mask, const double* lower,
                             const double* upper) {
  underDevelopmentLogMessage("changeRowsBounds");
  HighsStatus return_status = HighsStatus::Error;
  assert(hmos_.size() > 0);
  HighsSimplexInterface interface(hmos_[0]);
  return_status = interface.changeRowBounds(mask, lower, upper);
  return return_status == HighsStatus::OK;
}

bool Highs::changeCoeff(const int row, const int col, const double value) {
  underDevelopmentLogMessage("changeCoeff");
  HighsStatus return_status = HighsStatus::Error;
  assert(hmos_.size() > 0);
  HighsSimplexInterface interface(hmos_[0]);
  return_status = interface.changeCoefficient(row, col, value);
  return return_status == HighsStatus::OK;
}

bool Highs::getCols(const int from_col, const int to_col, int& num_col,
                    double* costs, double* lower, double* upper, int& num_nz,
                    int* start, int* index, double* value) {
  underDevelopmentLogMessage("getCols");
  HighsStatus return_status = HighsStatus::Error;
  assert(hmos_.size() > 0);
  HighsSimplexInterface interface(hmos_[0]);
  return_status = interface.getCols(from_col, to_col, num_col, costs, lower,
                                    upper, num_nz, start, index, value);
  return return_status == HighsStatus::OK;
}

bool Highs::getCols(const int n, const int* set, int& num_col, double* costs,
                    double* lower, double* upper, int& num_nz, int* start,
                    int* index, double* value) {
  underDevelopmentLogMessage("getCols");
  HighsStatus return_status = HighsStatus::Error;
  assert(hmos_.size() > 0);
  HighsSimplexInterface interface(hmos_[0]);
  return_status = interface.getCols(n, set, num_col, costs, lower, upper,
                                    num_nz, start, index, value);
  return return_status == HighsStatus::OK;
}

bool Highs::getCols(const int* col_mask, int& num_col, double* costs,
                    double* lower, double* upper, int& num_nz, int* start,
                    int* index, double* value) {
  underDevelopmentLogMessage("getCols");
  HighsStatus return_status = HighsStatus::Error;
  assert(hmos_.size() > 0);
  HighsSimplexInterface interface(hmos_[0]);
  return_status = interface.getCols(col_mask, num_col, costs, lower, upper,
                                    num_nz, start, index, value);
  return return_status == HighsStatus::OK;
}

bool Highs::getRows(const int from_row, const int to_row, int& num_row,
                    double* lower, double* upper, int& num_nz, int* start,
                    int* index, double* value) {
  underDevelopmentLogMessage("getRows");
  HighsStatus return_status = HighsStatus::Error;
  assert(hmos_.size() > 0);
  HighsSimplexInterface interface(hmos_[0]);
  return_status = interface.getRows(from_row, to_row, num_row, lower, upper,
                                    num_nz, start, index, value);
  return return_status == HighsStatus::OK;
}

bool Highs::getRows(const int num_set_entries, const int* set, int& num_row,
                    double* lower, double* upper, int& num_nz, int* start,
                    int* index, double* value) {
  underDevelopmentLogMessage("getRows");
  HighsStatus return_status = HighsStatus::Error;
  assert(hmos_.size() > 0);
  HighsSimplexInterface interface(hmos_[0]);
  return_status = interface.getRows(num_set_entries, set, num_row, lower, upper,
                                    num_nz, start, index, value);
  return return_status == HighsStatus::OK;
}

bool Highs::getRows(const int* mask, int& num_row, double* lower, double* upper,
                    int& num_nz, int* start, int* index, double* value) {
  underDevelopmentLogMessage("getRows");
  HighsStatus return_status = HighsStatus::Error;
  assert(hmos_.size() > 0);
  HighsSimplexInterface interface(hmos_[0]);
  return_status = interface.getRows(mask, num_row, lower, upper, num_nz, start,
                                    index, value);
  return return_status == HighsStatus::OK;
}

bool Highs::getCoeff(const int row, const int col, double& value) {
  underDevelopmentLogMessage("getCoeff");
  HighsStatus return_status = HighsStatus::Error;
  assert(hmos_.size() > 0);
  HighsSimplexInterface interface(hmos_[0]);

  return_status = interface.getCoefficient(row, col, value);
  return return_status == HighsStatus::OK;
}

bool Highs::deleteCols(const int from_col, const int to_col) {
  underDevelopmentLogMessage("deleteCols");
  HighsStatus return_status = HighsStatus::Error;
  assert(hmos_.size() > 0);
  HighsSimplexInterface interface(hmos_[0]);
  return_status = interface.deleteCols(from_col, to_col);
  updateHighsSolutionBasis();
  return return_status == HighsStatus::OK;
}

bool Highs::deleteCols(const int num_set_entries, const int* set) {
  underDevelopmentLogMessage("deleteCols");
  HighsStatus return_status = HighsStatus::Error;
  assert(hmos_.size() > 0);
  HighsSimplexInterface interface(hmos_[0]);
  return_status = interface.deleteCols(num_set_entries, set);
  updateHighsSolutionBasis();
  return return_status == HighsStatus::OK;
}

bool Highs::deleteCols(int* mask) {
  underDevelopmentLogMessage("deleteCols");
  HighsStatus return_status = HighsStatus::Error;
  assert(hmos_.size() > 0);
  HighsSimplexInterface interface(hmos_[0]);
  return_status = interface.deleteCols(mask);
  updateHighsSolutionBasis();
  return return_status == HighsStatus::OK;
}

bool Highs::deleteRows(const int from_row, const int to_row) {
  underDevelopmentLogMessage("deleteRows");
  HighsStatus return_status = HighsStatus::Error;
  assert(hmos_.size() > 0);
  HighsSimplexInterface interface(hmos_[0]);
  return_status = interface.deleteRows(from_row, to_row);
  updateHighsSolutionBasis();
  return return_status == HighsStatus::OK;
}

bool Highs::deleteRows(const int num_set_entries, const int* set) {
  underDevelopmentLogMessage("deleteRows");
  HighsStatus return_status = HighsStatus::Error;
  assert(hmos_.size() > 0);
  HighsSimplexInterface interface(hmos_[0]);
  return_status = interface.deleteRows(num_set_entries, set);
  updateHighsSolutionBasis();
  return return_status == HighsStatus::OK;
}

bool Highs::deleteRows(int* mask) {
  underDevelopmentLogMessage("deleteRows");
  HighsStatus return_status = HighsStatus::Error;
  assert(hmos_.size() > 0);
  HighsSimplexInterface interface(hmos_[0]);
  return_status = interface.deleteRows(mask);
  updateHighsSolutionBasis();
  return return_status == HighsStatus::OK;
}

HighsStatus Highs::clearSolver() {
  underDevelopmentLogMessage("clearSolver");
  basis_.valid_ = false;
  return HighsStatus::OK;
}

void Highs::reportModelStatusSolutionBasis(const std::string message, const HighsModelStatus model_status, const HighsLp &lp, const HighsSolution &solution, const HighsBasis &basis) {
  printf("\n%s\nModelStatus = %s; LP(%d, %d); solution (%d, %d; %d, %d); basis %d (%d, %d)\n\n",
	 message.c_str(), utilHighsModelStatusToString(model_status).c_str(), lp.numCol_, lp.numRow_,
	 (int)solution.col_value.size(), (int)solution.row_value.size(), (int)solution.col_dual.size(), (int)solution.row_dual.size(),
	 basis.valid_, (int)basis.col_status.size(), (int)basis.row_status.size());
}

std::string Highs::highsModelStatusToString(const HighsModelStatus model_status) { return utilHighsModelStatusToString(model_status); }

std::string Highs::highsPrimalDualStatusToString(const int primal_dual_status) { return utilPrimalDualStatusToString(primal_dual_status); }

// Private methods
HighsPresolveStatus Highs::runPresolve(PresolveInfo& info) {
  if (options_.presolve == off_string) return HighsPresolveStatus::NotPresolved;

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
HighsStatus Highs::callRunSolver(HighsModelObject& model, int& iteration_count,
                                 const string message) {
  HighsLogMessage(HighsMessageType::INFO, message.c_str());

  if (options_.solver == "ipm") {
    HighsPrintMessage(ML_ALWAYS, "Starting IPX...\n");
    IpxStatus ipx_return = solveModelWithIpx(lp_, options_, model_status_,
					     info_, solution_, basis_);
    if (ipx_return != IpxStatus::OK) {
      // todo:
      return HighsStatus::Error;
    }
    return HighsStatus::OK;
  }

  HighsStatus solver_return_status;
  if (!model.lp_.numRow_) {
    // Handle the case of unconstrained LPs here
    HighsSimplexInterface simplex_interface(model);
    solver_return_status = solveUnconstrainedLp(model);
    iteration_count = 0;
  } else {
    int initial_iteration_count = model.simplex_info_.iteration_count;
    solver_return_status = runSolver(model);
    int final_iteration_count = model.simplex_info_.iteration_count;
    iteration_count = final_iteration_count - initial_iteration_count;
  }
  return solver_return_status;
}

// The method below runs simplex or ipx solver on the lp.
HighsStatus Highs::runSolver(HighsModelObject& model) {
#ifdef HIGHSDEV
  // Shouldn't have to check validity of the LP since this is done when it is
  // loaded or modified
  //  bool normalise = true;
  HighsStatus assess_lp_status = assessLp(lp_, options_);  //, normalise);
  assert(assess_lp_status == HighsStatus::OK);
  if (assess_lp_status != HighsStatus::OK) return HighsStatus::Error;
#endif
#ifndef IPX
  // HiGHS
  // todo: Without the presolve part, so will be
  //     = solve_simplex(options, reduced_lp, reduced_solution)
  HighsStatus return_status = solveModelSimplex(model);
  if (return_status == HighsStatus::Error) return HighsStatus::Error;
    //  allow_presolve_ = false;
#else
  // IPX
  // todo:Check options for simplex-specific options
  // use model.lp_, model.solution_
  // status = runIpxSolver(options_, lp_, solution_);
  // If ipx crossover did not find optimality set up simplex.

#endif

  if (model.model_status_ != HighsModelStatus::OPTIMAL)
    return highsStatusFromHighsModelStatus(model.model_status_);

  // Check.
  if (!isSolutionConsistent(model.lp_, model.solution_)) {
    std::cout << "Error: Inconsistent solution returned from solver.\n";
  }

  // todo:
  // assert(KktSatisfied(lp, solution));

  return highsStatusFromHighsModelStatus(model.model_status_);
}

// Branch-and-bound code below here:
// Solve a mixed integer problem using branch and bound.
HighsStatus Highs::runBnb() {
  HighsPrintMessage(ML_ALWAYS, "Using branch and bound solver\n");

  // Need to start the HiGHS clock unless it's already running
  bool run_highs_clock_already_running = timer_.runningRunHighsClock();
  if (!run_highs_clock_already_running) timer_.startRunHighsClock();
  double mip_solve_initial_time = timer_.readRunHighsClock();

  // Start tree by making root node.
  std::unique_ptr<Node> root = std::unique_ptr<Node>(new Node(-1, 0, 0));

  root->integer_variables = lp_.integrality_;
  root->col_lower_bound = lp_.colLower_;
  root->col_upper_bound = lp_.colUpper_;

  HighsStatus return_status = solveRootNode(*(root.get()));
  if (return_status != HighsStatus::OK) return return_status;
  if (hmos_[0].model_status_ != HighsModelStatus::OPTIMAL) {
    HighsPrintMessage(ML_ALWAYS,
                      "Root note not solved to optimality. Status: %s\n",
                      utilHighsModelStatusToString(hmos_[0].model_status_).c_str());
    return highsStatusFromHighsModelStatus(hmos_[0].model_status_);
  }

  // The method branch(...) below calls chooseBranchingVariable(..) which
  // currently returns the first violated one. If a branching variable is found
  // children are added to the stack. If there are no more violated integrality
  // constraints we have a feasible solution, if it is best than current best,
  // the current best is updated.
  int message_level = ML_DETAILED | ML_VERBOSE;
  options_.message_level = message_level;
  Tree tree(*(root.get()));
  tree.branch(*(root.get()));

  // While stack not empty.
  //   Solve node.
  //   Branch.
  while (!tree.empty()) {
    Node& node = tree.next();
    HighsStatus return_status = solveNode(node);
    if (return_status != HighsStatus::OK) return return_status;
    tree.pop();

    if (hmos_[0].model_status_ == HighsModelStatus::PRIMAL_INFEASIBLE) continue;

    options_.message_level = message_level;
    tree.branch(node);
  }

  // Stop and read the HiGHS clock, then work out time for this call
  if (!run_highs_clock_already_running) timer_.stopRunHighsClock();
  double mip_solve_final_time = timer_.readRunHighsClock();

  if (tree.getBestSolution().size() > 0) {
    hmos_[0].model_status_ = HighsModelStatus::OPTIMAL;
    std::stringstream message;
    message << std::endl;
    message << "Optimal solution found.";
    message << std::endl;
    message << "Run status : "
            << highsModelStatusToString(hmos_[0].model_status_) << std::endl;
    message << "Objective  : " << std::scientific << tree.getBestObjective()
            << std::endl;
    message << "Time       : " << std::fixed << std::setprecision(3)
            << mip_solve_final_time - mip_solve_initial_time << std::endl;
    message << std::endl;

    HighsPrintMessage(ML_MINIMAL, message.str().c_str());
  } else {
    hmos_[0].model_status_ = HighsModelStatus::PRIMAL_INFEASIBLE;
    HighsPrintMessage(ML_ALWAYS, "No feasible solution found.\n");
  }

  return HighsStatus::OK;
}

HighsStatus Highs::solveNode(Node& node) {
  // Apply column bounds from node to LP.
  const bool check_call = false;
  const bool call_changeColsBounds = true;
  if (call_changeColsBounds) {
    changeColsBounds(0, lp_.numCol_, &node.col_lower_bound[0],
                     &node.col_upper_bound[0]);
  } else {
    // Change the LP directly and invalidate the simplex information
    lp_.colLower_ = node.col_lower_bound;
    lp_.colUpper_ = node.col_upper_bound;
    hmos_[0].simplex_lp_status_.valid = false;
  }

  // Call warm start.
  //  HighsStatus status = run();
  // call works but simply calling run() should be enough and will call hot
  // start in the same way as a user would call it from the outside

  int iteration_count0;
  int iteration_count1;
  int solve0_iteration_count;
  int solve1_iteration_count;
  double solve0_objective_value;
  double solve1_objective_value;
  int solve0_status;
  int solve1_status;

  iteration_count0 = hmos_[0].simplex_info_.iteration_count;

  HighsStatus return_status = solveModelSimplex(hmos_[0]);
  if (return_status == HighsStatus::Error) return HighsStatus::Error;
  //  allow_presolve_ = false;

  iteration_count1 = hmos_[0].simplex_info_.iteration_count;
  solve0_iteration_count = iteration_count1 - iteration_count0;
  solve0_objective_value = hmos_[0].simplex_info_.dual_objective_value;
  solve0_status = (int)hmos_[0].model_status_;
  printf("Solve0: Obj = %12g; Iter =%6d; Status =%2d\n", solve0_objective_value,
         solve0_iteration_count, solve0_status);

  if (check_call) {
    // Generate a fresh model object for the LP at this node
    hmos_[0].simplex_lp_status_.has_basis = false;
    hmos_[0].basis_.valid_ = false;
    iteration_count0 = hmos_[0].simplex_info_.iteration_count;
    HighsStatus return_status = solveModelSimplex(hmos_[0]);
    if (return_status != HighsStatus::OK) return return_status;
    iteration_count1 = hmos_[0].simplex_info_.iteration_count;
    solve1_iteration_count = iteration_count1 - iteration_count0;
    solve1_objective_value = hmos_[0].simplex_info_.dual_objective_value;
    solve1_status = (int)hmos_[0].model_status_;
    printf("Solve1: Obj = %12g; Iter =%6d; Status =%2d\n",
           solve1_objective_value, solve1_iteration_count, solve1_status);
    double rlv_objective_value_difference =
        fabs(solve1_objective_value - solve0_objective_value) /
        max(1.0, fabs(solve1_objective_value));
    if (solve0_status != solve1_status) {
      // Look for unequal status
      printf(
          "!! NodeSolveInequality: Status difference: Status0=%2d; Status1=%2d "
          "!!\n",
          solve0_status, solve1_status);
    } else if (solve0_status != (int)HighsModelStatus::PRIMAL_INFEASIBLE) {
      // Unless infeasible, look for unequal objective
      if (rlv_objective_value_difference > 1e-12)
        printf(
            "!! NodeSolveInequality: Relative objective difference = %12g !!\n",
            rlv_objective_value_difference);
    }
  }

  // Set solution.
  if (hmos_[0].model_status_ == HighsModelStatus::OPTIMAL) {
    node.primal_solution = hmos_[0].solution_.col_value;
    node.objective_value = hmos_[0].simplex_info_.dual_objective_value;
  }

  // Solve with a new hmo (replace with code above)
  // passModel(lp_);
  // lp_.colLower_ = node.col_lower_bound;
  // lp_.colUpper_ = node.col_upper_bound;

  // HighsStatus status = solveModelSimplex(hmos_[0]);

  // // Set solution.
  // if (status == HighsStatus::Optimal) {
  //   node.primal_solution = hmos_[0].solution_.col_value;
  //   node.objective_value = hmos_[0].simplex_info_.dual_objective_value;
  // }

  return highsStatusFromHighsModelStatus(hmos_[0].model_status_);
}

HighsStatus Highs::solveRootNode(Node& root) {
  // No presolve for the moment.
  options_.message_level = ML_NONE;
  // HighsStatus status = run();
  // call works but simply calling run() should be enough.
  HighsStatus return_status = solveModelSimplex(hmos_[0]);
  if (return_status == HighsStatus::Error) return HighsStatus::Error;
  //  allow_presolve_ = false;

  if (hmos_[0].model_status_ == HighsModelStatus::OPTIMAL) {
    root.primal_solution = hmos_[0].solution_.col_value;
    root.objective_value = hmos_[0].simplex_info_.dual_objective_value;
  }

  return highsStatusFromHighsModelStatus(hmos_[0].model_status_);
}

HighsStatus Highs::writeSolution(const std::string filename, const bool pretty) {
  HighsLp lp = this->lp_;
  FILE* file;
  if (filename == "") {
    // Empty file name: report model on stdout
    HighsLogMessage(HighsMessageType::WARNING, "Empty file name so reporting solution on stdout");
    file = stdout;
  } else {
    file = fopen(filename.c_str(), "w");
    if (file == 0) {
      HighsLogMessage(HighsMessageType::ERROR, "writeSolution: cannot open file");
      return HighsStatus::Error;
    }
  }
  if (pretty) {
    reportModelBoundSol(file,
			true, lp_.numCol_, lp_.colLower_, lp_.colUpper_,
			lp_.col_names_, solution_.col_value, solution_.col_dual,
			basis_.col_status);
    reportModelBoundSol(file,
			false, lp_.numRow_, lp_.rowLower_, lp_.rowUpper_,
			lp_.row_names_, solution_.row_value, solution_.row_dual,
			basis_.row_status);
  } else {
    fprintf(file, "%d %d : Number of columns and rows for primal and dual solution and basis\n", lp.numCol_, lp.numRow_);
    const bool with_basis = basis_.valid_;
    if (with_basis) {
      fprintf(file, "T\n");
    } else {
      fprintf(file, "F\n");
    }
    for (int iCol = 0; iCol < lp.numCol_; iCol++) {
      fprintf(file, "%g %g", solution_.col_value[iCol], solution_.col_dual[iCol]);
      if (with_basis) fprintf(file, " %d", (int)basis_.col_status[iCol]);
      fprintf(file, " \n");
    }
    for (int iRow = 0; iRow < lp.numRow_; iRow++) {
      fprintf(file, "%g %g", solution_.row_value[iRow], solution_.row_dual[iRow]);
      if (with_basis) fprintf(file, " %d", (int)basis_.row_status[iRow]);
      fprintf(file, " \n");
    }
  }
  if (file == stdout) return HighsStatus::Warning;
  return HighsStatus::OK;
}

void Highs::updateHighsSolutionBasis() {
  assert(hmos_.size() > 0);
  solution_.col_value.resize(lp_.numCol_);
  solution_.row_value.resize(lp_.numRow_);
  solution_.col_dual.resize(lp_.numCol_);
  solution_.row_dual.resize(lp_.numRow_);
  hmos_[0].solution_.col_value.resize(lp_.numCol_);
  hmos_[0].solution_.row_value.resize(lp_.numRow_);
  hmos_[0].solution_.col_dual.resize(lp_.numCol_);
  hmos_[0].solution_.row_dual.resize(lp_.numRow_);

  if (hmos_[0].basis_.valid_) {
    basis_ = hmos_[0].basis_;
  } else {
    basis_.valid_ = false;
    basis_.col_status.resize(lp_.numCol_);
    basis_.row_status.resize(lp_.numRow_);
  }
  //  if (hmos[0].simplex_lp_status_.has_simplex_lp_.
}  

void Highs::underDevelopmentLogMessage(const string method_name) {
  HighsLogMessage(
      HighsMessageType::WARNING,
      "Method %s is still under development and behaviour may be unpredictable",
      method_name.c_str());
}
