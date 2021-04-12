/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
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
//#include <string>
#include <iostream>
#include <memory>
#include <sstream>

#include "HConfig.h"
#include "io/Filereader.h"
#include "io/HighsIO.h"
#include "io/LoadOptions.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelUtils.h"
#include "lp_data/HighsSolution.h"
#include "lp_data/HighsSolve.h"
#include "mip/HighsMipSolver.h"
#include "simplex/HSimplexDebug.h"
#include "util/HighsMatrixPic.h"

#ifdef OPENMP
#include "omp.h"
#endif

Highs::Highs() {
  hmos_.clear();
  hmos_.push_back(HighsModelObject(lp_, options_, timer_));
}

HighsStatus Highs::setOptionValue(const std::string& option,
                                       const bool value) {
  if (setLocalOptionValue(options_.log_options, option, options_.records, value) ==
      OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::setOptionValue(const std::string& option,
                                       const HighsInt value) {
  if (setLocalOptionValue(options_.log_options, option, options_.records, value) ==
      OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::setOptionValue(const std::string& option,
                                       const double value) {
  if (setLocalOptionValue(options_.log_options, option, options_.records, value) ==
      OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::setOptionValue(const std::string& option,
                                       const std::string value) {
  if (setLocalOptionValue(options_.log_options, option, options_.records, value) ==
      OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::setOptionValue(const std::string& option,
                                       const char* value) {
  if (setLocalOptionValue(options_.log_options, option, options_.records, value) ==
      OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::readOptions(const std::string filename) {
  if (filename.size() <= 0) {
    highsLogUser(options_.log_options, HighsLogType::WARNING,
                 "Empty file name so not reading options\n");
    return HighsStatus::Warning;
  }
  options_.options_file = filename;
  if (!loadOptionsFromFile(options_)) return HighsStatus::Error;
  return HighsStatus::OK;
}

HighsStatus Highs::passOptions(const HighsOptions& options) {
  if (passLocalOptions(options_.log_options, options, options_) == OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

const HighsOptions& Highs::getOptions() { return options_; }

HighsStatus Highs::getOptionValue(const std::string& option, bool& value) {
  if (getLocalOptionValue(options_.log_options, option, options_.records, value) ==
      OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::getOptionValue(const std::string& option,
                                       HighsInt& value) {
  if (getLocalOptionValue(options_.log_options, option, options_.records, value) ==
      OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::getOptionValue(const std::string& option,
                                       double& value) {
  if (getLocalOptionValue(options_.log_options, option, options_.records, value) ==
      OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::getOptionValue(const std::string& option,
                                       std::string& value) {
  if (getLocalOptionValue(options_.log_options, option, options_.records, value) ==
      OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::getHighsOptionType(const std::string& option,
                                      HighsOptionType& type) {
  if (getOptionType(options_.log_options, option, options_.records, type) ==
      OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::resetHighsOptions() {
  resetOptions(options_.records);
  return HighsStatus::OK;
}

HighsStatus Highs::writeHighsOptions(
    const std::string filename, const bool report_only_non_default_values) {
  HighsStatus return_status = HighsStatus::OK;
  FILE* file;
  bool html;
  return_status = interpretCallStatus(
      openWriteFile(filename, "writeHighsOptions", file, html), return_status,
      "openWriteFile");
  if (return_status == HighsStatus::Error) return return_status;

  return_status = interpretCallStatus(
      writeOptionsToFile(file, options_.records, report_only_non_default_values,
                         html),
      return_status, "writeOptionsToFile");
  return return_status;
}

const HighsOptions& Highs::getHighsOptions() const { return options_; }

const HighsInfo& Highs::getHighsInfo() const { return info_; }

HighsStatus Highs::getHighsInfoValue(const std::string& info, HighsInt& value) {
  if (getInfoValue(options_, info, info_.records, value) == InfoStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::getHighsInfoValue(const std::string& info,
                                     double& value) const {
  if (getInfoValue(options_, info, info_.records, value) == InfoStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::writeHighsInfo(const std::string filename) {
  HighsStatus return_status = HighsStatus::OK;
  FILE* file;
  bool html;
  return_status =
      interpretCallStatus(openWriteFile(filename, "writeHighsInfo", file, html),
                          return_status, "openWriteFile");
  if (return_status == HighsStatus::Error) return return_status;

  return_status =
      interpretCallStatus(writeInfoToFile(file, info_.records, html),
                          return_status, "writeInfoToFile");
  return return_status;
}

// Methods below change the incumbent model or solver infomation
// associated with it. Hence returnFromHighs is called at the end of
// each
HighsStatus Highs::reset() {
  HighsStatus return_status = HighsStatus::OK;
  // Clear the status, solution, basis and info associated with any previous
  // model
  return_status =
      interpretCallStatus(clearSolver(), return_status, "clearSolver");
  if (return_status == HighsStatus::Error) return return_status;
  // Clear any HiGHS model object
  hmos_.clear();
  // Create a HiGHS model object for this LP
  hmos_.push_back(HighsModelObject(lp_, options_, timer_));

  presolve_.clear();

  return returnFromHighs(return_status);
}

HighsStatus Highs::passModel(const HighsLp lp) {
  HighsStatus return_status = HighsStatus::OK;
  // move the copy of the LP to the internal LP
  lp_ = std::move(lp);
  // Ensure that the LP is column-wise
  setOrientation(lp_);
  // Check validity of the LP, normalising its values
  return_status =
      interpretCallStatus(assessLp(lp_, options_), return_status, "assessLp");
  if (return_status == HighsStatus::Error) return return_status;
  // Clear solver status, solution, basis and info associated with any
  // previous model; clear any HiGHS model object; create a HiGHS
  // model object for this LP
  return_status = interpretCallStatus(reset(), return_status, "reset");
  return returnFromHighs(return_status);
}

HighsStatus Highs::passModel(const HighsInt num_col, const HighsInt num_row,
                             const HighsInt num_nz, const double* costs,
                             const double* col_lower, const double* col_upper,
                             const double* row_lower, const double* row_upper,
                             const HighsInt* astart, const HighsInt* aindex,
                             const double* avalue,
                             const HighsInt* integrality) {
  HighsLp lp;
  lp.numCol_ = num_col;
  lp.numRow_ = num_row;
  if (num_col > 0) {
    assert(costs != NULL);
    assert(col_lower != NULL);
    assert(col_upper != NULL);
    lp.colCost_.assign(costs, costs + num_col);
    lp.colLower_.assign(col_lower, col_lower + num_col);
    lp.colUpper_.assign(col_upper, col_upper + num_col);
  }
  if (num_row > 0) {
    assert(row_lower != NULL);
    assert(row_upper != NULL);
    lp.rowLower_.assign(row_lower, row_lower + num_row);
    lp.rowUpper_.assign(row_upper, row_upper + num_row);
  }
  if (num_nz > 0) {
    assert(num_col > 0);
    assert(num_row > 0);
    assert(astart != NULL);
    assert(aindex != NULL);
    assert(avalue != NULL);
    lp.Astart_.assign(astart, astart + num_col);
    lp.Aindex_.assign(aindex, aindex + num_nz);
    lp.Avalue_.assign(avalue, avalue + num_nz);
  }
  lp.Astart_.resize(num_col + 1);
  lp.Astart_[num_col] = num_nz;
  lp.orientation_ = MatrixOrientation::COLWISE;
  if (num_col > 0 && integrality != NULL) {
    lp.integrality_.resize(num_col);
    for (HighsInt iCol = 0; iCol < num_col; iCol++) {
      HighsInt integrality_status = integrality[iCol];
      assert(integrality_status == (HighsInt)HighsVarType::CONTINUOUS ||
             integrality_status == (HighsInt)HighsVarType::INTEGER);
      lp.integrality_[iCol] = (HighsVarType)integrality_status;
    }
  }
  return passModel(std::move(lp));
}

HighsStatus Highs::readModel(const std::string filename) {
  HighsStatus return_status = HighsStatus::OK;
  Filereader* reader = Filereader::getFilereader(filename);
  if (reader == NULL) {
    highsLogUser(options_.log_options, HighsLogType::ERROR,
                 "Model file %s not supported\n", filename.c_str());
    return HighsStatus::Error;
  }

  HighsLp model;
  options_.model_file = filename;

  FilereaderRetcode call_code = reader->readModelFromFile(options_, model);
  delete reader;
  if (call_code != FilereaderRetcode::OK) {
    interpretFilereaderRetcode(options_.log_options, filename.c_str(),
                               call_code);
    return_status = interpretCallStatus(HighsStatus::Error, return_status,
                                        "readModelFromFile");
    if (return_status == HighsStatus::Error) return return_status;
  }
  model.model_name_ = extractModelName(filename);
  return_status =
      interpretCallStatus(passModel(model), return_status, "passModel");
  return returnFromHighs(return_status);
}

HighsStatus Highs::clearModel() {
  HighsStatus return_status = HighsStatus::OK;
  // Remove all HighsModelObject entries
  hmos_.clear();
  lp_.clear();
  hmos_.push_back(HighsModelObject(lp_, options_, timer_));
  return_status =
      interpretCallStatus(clearSolver(), return_status, "clearSolver");
  if (return_status == HighsStatus::Error) return return_status;
  return returnFromHighs(return_status);
}

HighsStatus Highs::readBasis(const std::string filename) {
  HighsStatus return_status = HighsStatus::OK;
  // Try to read basis file into read_basis
  HighsBasis read_basis = basis_;
  return_status = interpretCallStatus(
      readBasisFile(options_.log_options, read_basis, filename), return_status,
      "readBasis");
  if (return_status != HighsStatus::OK) return return_status;
  // Basis read OK: check whether it's consistent with the LP
  if (!isBasisConsistent(lp_, read_basis)) {
    highsLogUser(options_.log_options, HighsLogType::ERROR,
                 "readBasis: invalid basis\n");
    return HighsStatus::Error;
  }
  // Update the HiGHS basis and invalidate any simplex basis for the model
  basis_ = read_basis;
  basis_.valid_ = true;
  if (hmos_.size() > 0) {
    clearBasisInterface();
  }
  // Can't use returnFromHighs since...
  return HighsStatus::OK;
}

HighsStatus Highs::writeModel(const std::string filename) {
  HighsStatus return_status = HighsStatus::OK;
  HighsLp model = lp_;

  // Ensure that the LP is column-wise
  setOrientation(model);
  if (filename == "") {
    // Empty file name: report model on stdout
    reportLp(options_.log_options, model, HighsLogType::VERBOSE);
    return_status = HighsStatus::OK;
  } else {
    Filereader* writer = Filereader::getFilereader(filename);
    if (writer == NULL) {
      highsLogUser(options_.log_options, HighsLogType::ERROR,
                   "Model file %s not supported\n", filename.c_str());
      return HighsStatus::Error;
    }
    return_status =
        interpretCallStatus(writer->writeModelToFile(options_, filename, model),
                            return_status, "writeModelToFile");
    delete writer;
  }
  return returnFromHighs(return_status);
}

HighsStatus Highs::writeBasis(const std::string filename) {
  HighsStatus return_status = HighsStatus::OK;
  return_status = interpretCallStatus(
      writeBasisFile(options_.log_options, basis_, filename), return_status,
      "writeBasis");
  return returnFromHighs(return_status);
}

// Checks the options calls presolve and postsolve if needed. Solvers are called
// with callSolveLp(..)
HighsStatus Highs::run() {
  if (!haveHmo("run")) return HighsStatus::Error;
  // Ensure that there is exactly one Highs model object
  assert((HighsInt)hmos_.size() == 1);
  HighsInt min_highs_debug_level = HIGHS_DEBUG_LEVEL_MIN;
  //  HIGHS_DEBUG_LEVEL_CHEAP;
  //  HIGHS_DEBUG_LEVEL_COSTLY;
  //  HIGHS_DEBUG_LEVEL_EXPENSIVE;
  //  HIGHS_DEBUG_LEVEL_MAX;
#ifdef HiGHSDEV
  min_highs_debug_level =  // HIGHS_DEBUG_LEVEL_MIN;
                           //  HIGHS_DEBUG_LEVEL_CHEAP;
      HIGHS_DEBUG_LEVEL_COSTLY;
  //  HIGHS_DEBUG_LEVEL_EXPENSIVE;
  //  HIGHS_DEBUG_LEVEL_MAX;
  if (options_.highs_debug_level < min_highs_debug_level)
    printf(
        "Highs::run() HiGHSDEV define so switching options_.highs_debug_level "
        "from %" HIGHSINT_FORMAT " to %" HIGHSINT_FORMAT "\n",
        options_.highs_debug_level, min_highs_debug_level);
    //  writeModel("HighsRunModel.mps");
    //  if (lp_.numRow_>0 && lp_.numCol_>0) writeLpMatrixPicToFile(options_,
    //  "LpMatrix", lp_);
#endif
  if (options_.highs_debug_level < min_highs_debug_level)
    options_.highs_debug_level = min_highs_debug_level;

#ifdef OPENMP
  omp_max_threads = omp_get_max_threads();
  assert(omp_max_threads > 0);
#ifdef HiGHSDEV
  if (omp_max_threads <= 0)
    printf("WARNING: omp_get_max_threads() returns %" HIGHSINT_FORMAT "\n",
           omp_max_threads);
  printf("Running with %" HIGHSINT_FORMAT " OMP thread(s)\n", omp_max_threads);
#endif
#endif
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Zero the HiGHS iteration counts
  zeroHighsIterationCounts(info_);

  // Initialise the HiGHS model status values
  hmos_[0].scaled_model_status_ = HighsModelStatus::NOTSET;
  hmos_[0].unscaled_model_status_ = HighsModelStatus::NOTSET;
  model_status_ = hmos_[0].scaled_model_status_;
  scaled_model_status_ = hmos_[0].unscaled_model_status_;
  // Return immediately if the LP has no columns
  if (!lp_.numCol_) {
    model_status_ = HighsModelStatus::MODEL_EMPTY;
    scaled_model_status_ = model_status_;
    hmos_[0].unscaled_model_status_ = model_status_;
    hmos_[0].scaled_model_status_ = model_status_;
    return_status = highsStatusFromHighsModelStatus(model_status_);
    return returnFromRun(return_status);
  }

  // Ensure that the LP (and any simplex LP) has the matrix column-wise
  setOrientation(lp_);
  if (hmos_[0].ekk_instance_.simplex_lp_status_.valid)
    setOrientation(hmos_[0].ekk_instance_.simplex_lp_);
#ifdef HIGHSDEV
  // Shouldn't have to check validity of the LP since this is done when it is
  // loaded or modified
  call_status = assessLp(lp_, options_);
  // If any errors have been found or normalisation carried out,
  // call_status will be ERROR or WARNING, so only valid return is OK.
  assert(call_status == HighsStatus::OK);
  return_status = interpretCallStatus(call_status, return_status, "assessLp");
  if (return_status == HighsStatus::Error) return returnFromRun(return_status);
#endif

  highsSetLogCallback(options_);
#ifdef HiGHSDEV
  if (checkOptions(options_.log_options, options_.records) !=
      OptionStatus::OK) {
    return_status = HighsStatus::Error;
    return returnFromRun(return_status);
  }
#endif
  if (lp_.model_name_.compare(""))
    highsLogDev(options_.log_options, HighsLogType::VERBOSE,
                "Solving model: %s\n", lp_.model_name_.c_str());

  // Start the HiGHS clock unless it's already running
  bool run_highs_clock_already_running = timer_.runningRunHighsClock();
  if (!run_highs_clock_already_running) timer_.startRunHighsClock();

  if (!options_.solver.compare(choose_string) && isMip(lp_)) {
    // Solve the model as a MIP
    call_status = callSolveMip();
    return_status =
        interpretCallStatus(call_status, return_status, "callSolveMip");
    if (!run_highs_clock_already_running) timer_.stopRunHighsClock();
    return returnFromRun(return_status);
  }

  // Solve the model as an LP
  //
  // Record the initial time and set the component times and postsolve
  // iteration count to -1 to identify whether they are not required
  double initial_time = timer_.readRunHighsClock();
  double this_presolve_time = -1;
  double this_solve_presolved_lp_time = -1;
  double this_postsolve_time = -1;
  double this_solve_original_lp_time = -1;
  HighsInt postsolve_iteration_count = -1;

  // Define identifiers to refer to the HMO of the original LP (0) and
  // the HMO created when using presolve. The index of this HMO is 1
  // when solving a one-off LP, but greater than one if presolve has
  // been called multiple times. It's equal to the size of HMO
  const HighsInt original_hmo = 0;
  const HighsInt presolve_hmo = hmos_.size();
  // Keep track of the hmo that is the most recently solved. By default it's the
  // original LP
  HighsInt solved_hmo = original_hmo;

  if (!basis_.valid_ && options_.presolve != off_string) {
    // No HiGHS basis so consider presolve
    //
    // If using IPX to solve the reduced LP, crossover must be run
    // since a basic solution is required by postsolve
    if (options_.solver == ipm_string && !options_.run_crossover) {
      highsLogUser(options_.log_options, HighsLogType::WARNING,
                   "Forcing IPX to use crossover after presolve\n");
      options_.run_crossover = true;
    }

    hmos_[original_hmo].scaled_model_status_ = HighsModelStatus::NOTSET;
    // Possibly presolve - according to option_.presolve
    const double from_presolve_time = timer_.read(timer_.presolve_clock);
    this_presolve_time = -from_presolve_time;
    timer_.start(timer_.presolve_clock);
    HighsPresolveStatus presolve_status = runPresolve();
    timer_.stop(timer_.presolve_clock);
    const double to_presolve_time = timer_.read(timer_.presolve_clock);
    this_presolve_time += to_presolve_time;
    presolve_.info_.presolve_time = this_presolve_time;

    // Set an illegal local pivot threshold value that's updated after
    // solving the presolved LP - if simplex is used
    double factor_pivot_threshold = -1;

    // Run solver.
    switch (presolve_status) {
      case HighsPresolveStatus::NotPresolved: {
        hmos_[solved_hmo].lp_.lp_name_ = "Original LP";
        this_solve_original_lp_time = -timer_.read(timer_.solve_clock);
        timer_.start(timer_.solve_clock);
        call_status = callSolveLp(solved_hmo, "Not presolved: solving the LP");
        timer_.stop(timer_.solve_clock);
        this_solve_original_lp_time += timer_.read(timer_.solve_clock);
        return_status =
            interpretCallStatus(call_status, return_status, "callSolveLp");
        if (return_status == HighsStatus::Error)
          return returnFromRun(return_status);
        break;
      }
      case HighsPresolveStatus::NotReduced: {
        hmos_[solved_hmo].lp_.lp_name_ = "Unreduced LP";
        // Log the presolve reductions
        reportPresolveReductions(hmos_[original_hmo].options_.log_options,
                                 hmos_[original_hmo].lp_, false);
        this_solve_original_lp_time = -timer_.read(timer_.solve_clock);
        timer_.start(timer_.solve_clock);
        call_status = callSolveLp(
            solved_hmo, "Problem not reduced by presolve: solving the LP");
        timer_.stop(timer_.solve_clock);
        this_solve_original_lp_time += timer_.read(timer_.solve_clock);
        return_status =
            interpretCallStatus(call_status, return_status, "callSolveLp");
        if (return_status == HighsStatus::Error)
          return returnFromRun(return_status);
        break;
      }
      case HighsPresolveStatus::Reduced: {
        HighsLp& reduced_lp = presolve_.getReducedProblem();
        // Validate the reduced LP
        assert(assessLp(reduced_lp, options_) == HighsStatus::OK);
        call_status = cleanBounds(options_, reduced_lp);
        // Ignore any warning from clean bounds since the original LP
        // is still solved after presolve
        if (interpretCallStatus(call_status, return_status, "cleanBounds") ==
            HighsStatus::Error)
          return HighsStatus::Error;
        // Add reduced lp object to vector of HighsModelObject,
        // so the last one in lp_ is the presolved one.

        hmos_.push_back(HighsModelObject(reduced_lp, options_, timer_));
        // Log the presolve reductions
        reportPresolveReductions(hmos_[original_hmo].options_.log_options,
                                 hmos_[original_hmo].lp_,
                                 hmos_[presolve_hmo].lp_);
        // Record the HMO to be solved
        solved_hmo = presolve_hmo;
        hmos_[solved_hmo].lp_.lp_name_ = "Presolved LP";
        // Don't try dual cut-off when solving the presolved LP, as the
        // objective values aren't correct
        //	HighsOptions& options = hmos_[solved_hmo].options_;
        //	HighsOptions save_options = options;
        const double save_dual_objective_value_upper_bound =
            options_.dual_objective_value_upper_bound;
        options_.dual_objective_value_upper_bound = HIGHS_CONST_INF;
        this_solve_presolved_lp_time = -timer_.read(timer_.solve_clock);
        timer_.start(timer_.solve_clock);
        call_status = callSolveLp(solved_hmo, "Solving the presolved LP");
        timer_.stop(timer_.solve_clock);
        this_solve_presolved_lp_time += timer_.read(timer_.solve_clock);
        if (hmos_[solved_hmo].ekk_instance_.simplex_lp_status_.valid) {
          // Record the pivot threshold resulting from solving the presolved LP
          // with simplex
          factor_pivot_threshold =
              hmos_[solved_hmo]
                  .ekk_instance_.simplex_info_.factor_pivot_threshold;
        }
        // Restore the dual objective cut-off
        options_.dual_objective_value_upper_bound =
            save_dual_objective_value_upper_bound;
        return_status =
            interpretCallStatus(call_status, return_status, "callSolveLp");
        if (return_status == HighsStatus::Error)
          return returnFromRun(return_status);

        break;
      }
      case HighsPresolveStatus::ReducedToEmpty: {
        reportPresolveReductions(hmos_[original_hmo].options_.log_options,
                                 hmos_[original_hmo].lp_, true);
        hmos_[original_hmo].unscaled_model_status_ = HighsModelStatus::OPTIMAL;
        hmos_[original_hmo].scaled_model_status_ =
            hmos_[original_hmo].unscaled_model_status_;
        // Proceed to postsolve.
        break;
      }
        //	printf("\nHighs::run() 3: presolve status = %" HIGHSINT_FORMAT
        //"\n", (HighsInt)presolve_status);fflush(stdout);
      case HighsPresolveStatus::Infeasible:
      case HighsPresolveStatus::Unbounded: {
        if (presolve_status == HighsPresolveStatus::Infeasible) {
          model_status_ = HighsModelStatus::PRIMAL_INFEASIBLE;
        } else {
          // If presolve returns (primal) unbounded, the problem may
          // not be feasible, in which case
          // HighsModelStatus::PRIMAL_INFEASIBLE_OR_UNBOUNDED rather
          // than HighsModelStatus::PRIMAL_UNBOUNDED should be
          // returned
          model_status_ = HighsModelStatus::PRIMAL_INFEASIBLE_OR_UNBOUNDED;
        }
        highsLogUser(options_.log_options, HighsLogType::INFO,
                     "Problem status detected on presolve: %s\n",
                     modelStatusToString(model_status_).c_str());

        // Report this way for the moment. May modify after merge with
        // OSIinterface branch which has new way of setting up a
        // HighsModelObject and can support multiple calls to run(). Stop and
        // read the HiGHS clock, then work out time for this call
        if (!run_highs_clock_already_running) timer_.stopRunHighsClock();

        // Transfer the model status to the scaled model status and orriginal
        // HMO statuses;
        scaled_model_status_ = model_status_;
        hmos_[original_hmo].unscaled_model_status_ = model_status_;
        hmos_[original_hmo].scaled_model_status_ = model_status_;
        return_status = HighsStatus::OK;
        return returnFromRun(return_status);
      }
      case HighsPresolveStatus::Timeout: {
        model_status_ = HighsModelStatus::PRESOLVE_ERROR;
        highsLogDev(options_.log_options, HighsLogType::ERROR,
                    "Presolve reached timeout\n");
        if (run_highs_clock_already_running) timer_.stopRunHighsClock();
        return HighsStatus::Warning;
      }
      case HighsPresolveStatus::OptionsError: {
        model_status_ = HighsModelStatus::PRESOLVE_ERROR;
        highsLogDev(options_.log_options, HighsLogType::ERROR,
                    "Presolve options error.\n");
        if (run_highs_clock_already_running) timer_.stopRunHighsClock();
        return HighsStatus::Warning;
      }
      default: {
        // case HighsPresolveStatus::Error
        model_status_ = HighsModelStatus::PRESOLVE_ERROR;
        highsLogDev(options_.log_options, HighsLogType::ERROR,
                    "Presolve failed.\n");
        if (run_highs_clock_already_running) timer_.stopRunHighsClock();
        // Transfer the model status to the scaled model status and orriginal
        // HMO statuses;
        scaled_model_status_ = model_status_;
        hmos_[original_hmo].unscaled_model_status_ = model_status_;
        hmos_[original_hmo].scaled_model_status_ = model_status_;
        return_status = HighsStatus::Error;
        return returnFromRun(return_status);
      }
    }
    // Postsolve. Does nothing if there were no reductions during presolve.
    if (hmos_[solved_hmo].scaled_model_status_ == HighsModelStatus::OPTIMAL) {
      if (presolve_status == HighsPresolveStatus::Reduced ||
          presolve_status == HighsPresolveStatus::ReducedToEmpty) {
        // If presolve is nontrivial, extract the optimal solution
        // and basis for the presolved problem in order to generate
        // the solution and basis for postsolve to use to generate a
        // solution(?) and basis that is, hopefully, optimal. This is
        // confirmed or corrected by hot-starting the simplex solver
        if (presolve_status == HighsPresolveStatus::ReducedToEmpty) {
          clearSolutionUtil(hmos_[solved_hmo].solution_);
          clearBasisUtil(hmos_[solved_hmo].basis_);
        }

        presolve_.data_.recovered_solution_ = hmos_[solved_hmo].solution_;
        presolve_.data_.recovered_basis_.col_status =
            hmos_[solved_hmo].basis_.col_status;
        presolve_.data_.recovered_basis_.row_status =
            hmos_[solved_hmo].basis_.row_status;

        this_postsolve_time = -timer_.read(timer_.postsolve_clock);
        timer_.start(timer_.postsolve_clock);
        HighsPostsolveStatus postsolve_status = runPostsolve();
        timer_.stop(timer_.postsolve_clock);
        this_postsolve_time += -timer_.read(timer_.postsolve_clock);
        presolve_.info_.postsolve_time = this_postsolve_time;

        if (postsolve_status == HighsPostsolveStatus::SolutionRecovered) {
          highsLogDev(options_.log_options, HighsLogType::VERBOSE,
                      "Postsolve finished\n");
          //
          // Now hot-start the simplex solver for the original_hmo:
          //
          // The original model hasn't been solved, so set up its solution
          // parameters
          resetModelStatusAndSolutionParams(hmos_[original_hmo]);
          // Set solution and its status
          hmos_[original_hmo].solution_ = presolve_.data_.recovered_solution_;

          // Set basis and its status
          hmos_[original_hmo].basis_.valid_ = true;
          hmos_[original_hmo].basis_.col_status =
              presolve_.data_.recovered_basis_.col_status;
          hmos_[original_hmo].basis_.row_status =
              presolve_.data_.recovered_basis_.row_status;

          // Possibly force debug to perform KKT check on what's
          // returned from postsolve
          const bool force_debug = false;
          HighsInt save_highs_debug_level = options_.highs_debug_level;
          if (force_debug)
            options_.highs_debug_level = HIGHS_DEBUG_LEVEL_COSTLY;
          debugHighsBasicSolution("After returning from postsolve", options_,
                                  lp_, hmos_[original_hmo].basis_,
                                  hmos_[original_hmo].solution_);
          options_.highs_debug_level = save_highs_debug_level;

          // Now hot-start the simplex solver for the original_hmo
          solved_hmo = original_hmo;
          // Save the options to allow the best simplex strategy to
          // be used
          HighsOptions& options = hmos_[solved_hmo].options_;
          HighsOptions save_options = options;
          const bool full_logging = false;
          if (full_logging) options.log_dev_level = LOG_DEV_LEVEL_VERBOSE;
          // Force the use of simplex to clean up if IPM has been used
          // to solve the presolved problem
          if (options.solver == ipm_string) options.solver = simplex_string;
          options.simplex_strategy = SIMPLEX_STRATEGY_CHOOSE;
          // Ensure that the parallel solver isn't used
          options.highs_min_threads = 1;
          options.highs_max_threads = 1;
          // Use any pivot threshold resulting from solving the presolved LP
          if (factor_pivot_threshold > 0)
            options.factor_pivot_threshold = factor_pivot_threshold;

          // The basis returned from postsolve is just basic/nonbasic
          // and EKK expects a refined basis, so set it up now
          refineBasis(lp_, hmos_[original_hmo].solution_,
                      hmos_[original_hmo].basis_);

          hmos_[solved_hmo].lp_.lp_name_ = "Postsolve LP";
          HighsInt iteration_count0 = info_.simplex_iteration_count;
          this_solve_original_lp_time = -timer_.read(timer_.solve_clock);
          timer_.start(timer_.solve_clock);
          call_status = callSolveLp(
              solved_hmo,
              "Solving the original LP from the solution after postsolve");
          timer_.stop(timer_.solve_clock);
          this_solve_original_lp_time += timer_.read(timer_.solve_clock);
          return_status =
              interpretCallStatus(call_status, return_status, "callSolveLp");
          // Recover the options
          options = save_options;
          if (return_status == HighsStatus::Error)
            return returnFromRun(return_status);
          postsolve_iteration_count =
              info_.simplex_iteration_count - iteration_count0;
        } else {
          highsLogUser(options_.log_options, HighsLogType::ERROR,
                       "Postsolve return status is %" HIGHSINT_FORMAT "\n",
                       (HighsInt)postsolve_status);
          model_status_ = HighsModelStatus::POSTSOLVE_ERROR;
          scaled_model_status_ = model_status_;
          hmos_[0].unscaled_model_status_ = model_status_;
          hmos_[0].scaled_model_status_ = model_status_;
          return_status = HighsStatus::Error;
          return returnFromRun(return_status);
        }
      }
    } else {
      // Optimal solution of presolved problem has not been found
      // The original model inherits the solved model's status
      hmos_[original_hmo].unscaled_model_status_ =
          hmos_[solved_hmo].unscaled_model_status_;
      hmos_[original_hmo].scaled_model_status_ =
          hmos_[solved_hmo].scaled_model_status_;
    }
  } else {
    // There is a valid basis for the problem or presolve is off
    solved_hmo = original_hmo;
    hmos_[solved_hmo].lp_.lp_name_ = "LP without presolve or with basis";
    if (basis_.valid_) {
      // There is a valid HiGHS basis, so use it to initialise the basis
      // in the HMO to be solved after refining any status values that
      // are simply HighsBasisStatus::NONBASIC
      refineBasis(hmos_[solved_hmo].lp_, solution_, basis_);
      hmos_[solved_hmo].basis_ = basis_;
    }
    this_solve_original_lp_time = -timer_.read(timer_.solve_clock);
    timer_.start(timer_.solve_clock);
    call_status =
        callSolveLp(solved_hmo, "Solving LP without presolve or with basis");
    timer_.stop(timer_.solve_clock);
    this_solve_original_lp_time += timer_.read(timer_.solve_clock);
    return_status =
        interpretCallStatus(call_status, return_status, "callSolveLp");
    if (return_status == HighsStatus::Error)
      return returnFromRun(return_status);
  }
  // else if (reduced problem failed to solve) {
  //   todo: handle case when presolved problem failed to solve. Try to solve
  //   again with no presolve.
  // }

  //   assert(solved_hmo == original_hmo);
  // solved_hmo will be original_hmo unless the presolved LP is found to be
  // infeasible or unbounded

  if (!getHighsModelStatusAndInfo(solved_hmo)) {
    return_status = HighsStatus::Error;
    return returnFromRun(return_status);
  }

  // Copy HMO solution/basis to HiGHS solution/basis: this resizes solution_ and
  // basis_ The HiGHS solution and basis have to come from the original_hmo for
  // them to have the right dimension.
  solution_ = hmos_[original_hmo].solution_;
  basis_ = hmos_[original_hmo].basis_;
  // Stop and read the HiGHS clock, then work out time for this call
  if (!run_highs_clock_already_running) timer_.stopRunHighsClock();

  double lp_solve_final_time = timer_.readRunHighsClock();
  double this_solve_time = lp_solve_final_time - initial_time;
  if (postsolve_iteration_count < 0) {
    highsLogDev(options_.log_options, HighsLogType::INFO, "Postsolve  : \n");
  } else {
    highsLogDev(options_.log_options, HighsLogType::INFO,
                "Postsolve  : %" HIGHSINT_FORMAT "\n",
                postsolve_iteration_count);
  }
  highsLogDev(options_.log_options, HighsLogType::INFO, "Time       : %8.2f\n",
              this_solve_time);
  highsLogDev(options_.log_options, HighsLogType::INFO, "Time Pre   : %8.2f\n",
              this_presolve_time);
  highsLogDev(options_.log_options, HighsLogType::INFO, "Time PreLP : %8.2f\n",
              this_solve_presolved_lp_time);
  highsLogDev(options_.log_options, HighsLogType::INFO, "Time PostLP: %8.2f\n",
              this_solve_original_lp_time);
  if (this_solve_time > 0) {
    highsLogDev(options_.log_options, HighsLogType::INFO, "For LP %16s",
                hmos_[original_hmo].lp_.model_name_.c_str());
    double sum_time = 0;
    if (this_presolve_time > 0) {
      sum_time += this_presolve_time;
      HighsInt pct = (100 * this_presolve_time) / this_solve_time;
      highsLogDev(options_.log_options, HighsLogType::INFO,
                  ": Presolve %8.2f (%3" HIGHSINT_FORMAT "%%)",
                  this_presolve_time, pct);
    }
    if (this_solve_presolved_lp_time > 0) {
      sum_time += this_solve_presolved_lp_time;
      HighsInt pct = (100 * this_solve_presolved_lp_time) / this_solve_time;
      highsLogDev(options_.log_options, HighsLogType::INFO,
                  ": Solve presolved LP %8.2f (%3" HIGHSINT_FORMAT "%%)",
                  this_solve_presolved_lp_time, pct);
    }
    if (this_postsolve_time > 0) {
      sum_time += this_postsolve_time;
      HighsInt pct = (100 * this_postsolve_time) / this_solve_time;
      highsLogDev(options_.log_options, HighsLogType::INFO,
                  ": Postsolve %8.2f (%3" HIGHSINT_FORMAT "%%)",
                  this_postsolve_time, pct);
    }
    if (this_solve_original_lp_time > 0) {
      sum_time += this_solve_original_lp_time;
      HighsInt pct = (100 * this_solve_original_lp_time) / this_solve_time;
      highsLogDev(options_.log_options, HighsLogType::INFO,
                  ": Solve original LP %8.2f (%3" HIGHSINT_FORMAT "%%)",
                  this_solve_original_lp_time, pct);
    }
    highsLogDev(options_.log_options, HighsLogType::INFO, "\n");
    double rlv_time_difference =
        fabs(sum_time - this_solve_time) / this_solve_time;
    if (rlv_time_difference > 0.1)
      highsLogDev(options_.log_options, HighsLogType::INFO,
                  "Strange: Solve time = %g; Sum times = %g: relative "
                  "difference = %g\n",
                  this_solve_time, sum_time, rlv_time_difference);
  }
  // Assess success according to the scaled model status, unless
  // something worse has happened earlier
  call_status = highsStatusFromHighsModelStatus(scaled_model_status_);
  return_status = interpretCallStatus(call_status, return_status);
  return returnFromRun(return_status);
}

const HighsModelStatus& Highs::getModelStatus(const bool scaled_model) const {
  if (scaled_model) {
    return scaled_model_status_;
  } else {
    return model_status_;
  }
}

HighsStatus Highs::getDualRay(bool& has_dual_ray, double* dual_ray_value) {
  if (!haveHmo("getDualRay")) return HighsStatus::Error;
  return getDualRayInterface(has_dual_ray, dual_ray_value);
}

HighsStatus Highs::getPrimalRay(bool& has_primal_ray,
                                double* primal_ray_value) {
  underDevelopmentLogMessage("getPrimalRay");
  if (!haveHmo("getPrimalRay")) return HighsStatus::Error;
  return getPrimalRayInterface(has_primal_ray, primal_ray_value);
}

HighsStatus Highs::getRanging(HighsRanging& ranging) {
  underDevelopmentLogMessage("getRanging");
  if (!haveHmo("getRanging")) return HighsStatus::Error;
  return getRangingData(ranging, hmos_[0]);
}

HighsStatus Highs::getBasicVariables(HighsInt* basic_variables) {
  if (!haveHmo("getBasicVariables")) return HighsStatus::Error;
  if (basic_variables == NULL) {
    highsLogUser(options_.log_options, HighsLogType::ERROR,
                 "getBasicVariables: basic_variables is NULL\n");
    return HighsStatus::Error;
  }
  return getBasicVariablesInterface(basic_variables);
}

HighsStatus Highs::getBasisInverseRow(const HighsInt row, double* row_vector,
                                      HighsInt* row_num_nz,
                                      HighsInt* row_indices) {
  if (!haveHmo("getBasisInverseRow")) return HighsStatus::Error;
  if (row_vector == NULL) {
    highsLogUser(options_.log_options, HighsLogType::ERROR,
                 "getBasisInverseRow: row_vector is NULL\n");
    return HighsStatus::Error;
  }
  // row_indices can be NULL - it's the trigger that determines
  // whether they are identified or not
  HighsInt numRow = lp_.numRow_;
  if (row < 0 || row >= numRow) {
    highsLogUser(options_.log_options, HighsLogType::ERROR,
                 "Row index %" HIGHSINT_FORMAT
                 " out of range [0, %" HIGHSINT_FORMAT
                 "] in getBasisInverseRow\n",
                 row, numRow - 1);
    return HighsStatus::Error;
  }
  bool has_invert = hmos_[0].ekk_instance_.simplex_lp_status_.has_invert;
  if (!has_invert) {
    highsLogUser(options_.log_options, HighsLogType::ERROR,
                 "No invertible representation for getBasisInverseRow\n");
    return HighsStatus::Error;
  }
  // Compute a row i of the inverse of the basis matrix by solving B^Tx=e_i
  vector<double> rhs;
  rhs.assign(numRow, 0);
  rhs[row] = 1;
  basisSolveInterface(rhs, row_vector, row_num_nz, row_indices, true);
  return HighsStatus::OK;
}

HighsStatus Highs::getBasisInverseCol(const HighsInt col, double* col_vector,
                                      HighsInt* col_num_nz,
                                      HighsInt* col_indices) {
  if (!haveHmo("getBasisInverseCol")) return HighsStatus::Error;
  if (col_vector == NULL) {
    highsLogUser(options_.log_options, HighsLogType::ERROR,
                 "getBasisInverseCol: col_vector is NULL\n");
    return HighsStatus::Error;
  }
  // col_indices can be NULL - it's the trigger that determines
  // whether they are identified or not
  HighsInt numRow = lp_.numRow_;
  if (col < 0 || col >= numRow) {
    highsLogUser(options_.log_options, HighsLogType::ERROR,
                 "Column index %" HIGHSINT_FORMAT
                 " out of range [0, %" HIGHSINT_FORMAT
                 "] in getBasisInverseCol\n",
                 col, numRow - 1);
    return HighsStatus::Error;
  }
  bool has_invert = hmos_[0].ekk_instance_.simplex_lp_status_.has_invert;
  if (!has_invert) {
    highsLogUser(options_.log_options, HighsLogType::ERROR,
                 "No invertible representation for getBasisInverseCol\n");
    return HighsStatus::Error;
  }
  // Compute a col i of the inverse of the basis matrix by solving Bx=e_i
  vector<double> rhs;
  rhs.assign(numRow, 0);
  rhs[col] = 1;
  basisSolveInterface(rhs, col_vector, col_num_nz, col_indices, false);
  return HighsStatus::OK;
}

HighsStatus Highs::getBasisSolve(const double* Xrhs, double* solution_vector,
                                 HighsInt* solution_num_nz,
                                 HighsInt* solution_indices) {
  if (!haveHmo("getBasisSolve")) return HighsStatus::Error;
  if (Xrhs == NULL) {
    highsLogUser(options_.log_options, HighsLogType::ERROR,
                 "getBasisSolve: Xrhs is NULL\n");
    return HighsStatus::Error;
  }
  if (solution_vector == NULL) {
    highsLogUser(options_.log_options, HighsLogType::ERROR,
                 "getBasisSolve: solution_vector is NULL\n");
    return HighsStatus::Error;
  }
  // solution_indices can be NULL - it's the trigger that determines
  // whether they are identified or not
  bool has_invert = hmos_[0].ekk_instance_.simplex_lp_status_.has_invert;
  if (!has_invert) {
    highsLogUser(options_.log_options, HighsLogType::ERROR,
                 "No invertible representation for getBasisSolve\n");
    return HighsStatus::Error;
  }
  HighsInt numRow = lp_.numRow_;
  vector<double> rhs;
  rhs.assign(numRow, 0);
  for (HighsInt row = 0; row < numRow; row++) rhs[row] = Xrhs[row];
  basisSolveInterface(rhs, solution_vector, solution_num_nz, solution_indices,
                      false);
  return HighsStatus::OK;
}

HighsStatus Highs::getBasisTransposeSolve(const double* Xrhs,
                                          double* solution_vector,
                                          HighsInt* solution_num_nz,
                                          HighsInt* solution_indices) {
  if (!haveHmo("getBasisTransposeSolve")) return HighsStatus::Error;
  if (Xrhs == NULL) {
    highsLogUser(options_.log_options, HighsLogType::ERROR,
                 "getBasisTransposeSolve: Xrhs is NULL\n");
    return HighsStatus::Error;
  }
  if (solution_vector == NULL) {
    highsLogUser(options_.log_options, HighsLogType::ERROR,
                 "getBasisTransposeSolve: solution_vector is NULL\n");
    return HighsStatus::Error;
  }
  // solution_indices can be NULL - it's the trigger that determines
  // whether they are identified or not
  bool has_invert = hmos_[0].ekk_instance_.simplex_lp_status_.has_invert;
  if (!has_invert) {
    highsLogUser(options_.log_options, HighsLogType::ERROR,
                 "No invertible representation for getBasisTransposeSolve\n");
    return HighsStatus::Error;
  }
  HighsInt numRow = lp_.numRow_;
  vector<double> rhs;
  rhs.assign(numRow, 0);
  for (HighsInt row = 0; row < numRow; row++) rhs[row] = Xrhs[row];
  basisSolveInterface(rhs, solution_vector, solution_num_nz, solution_indices,
                      true);
  return HighsStatus::OK;
}

HighsStatus Highs::getReducedRow(const HighsInt row, double* row_vector,
                                 HighsInt* row_num_nz, HighsInt* row_indices,
                                 const double* pass_basis_inverse_row_vector) {
  if (!haveHmo("getReducedRow")) return HighsStatus::Error;
  // Ensure that the LP is column-wise
  setOrientation(lp_);
  HighsLp& lp = lp_;
  if (row_vector == NULL) {
    highsLogUser(options_.log_options, HighsLogType::ERROR,
                 "getReducedRow: row_vector is NULL\n");
    return HighsStatus::Error;
  }
  // row_indices can be NULL - it's the trigger that determines
  // whether they are identified or not pass_basis_inverse_row_vector
  // NULL - it's the trigger to determine whether it's computed or not
  if (row < 0 || row >= lp.numRow_) {
    highsLogUser(options_.log_options, HighsLogType::ERROR,
                 "Row index %" HIGHSINT_FORMAT
                 " out of range [0, %" HIGHSINT_FORMAT "] in getReducedRow\n",
                 row, lp.numRow_ - 1);
    return HighsStatus::Error;
  }
  bool has_invert = hmos_[0].ekk_instance_.simplex_lp_status_.has_invert;
  if (!has_invert) {
    highsLogUser(options_.log_options, HighsLogType::ERROR,
                 "No invertible representation for getReducedRow\n");
    return HighsStatus::Error;
  }
  HighsInt numRow = lp.numRow_;
  vector<double> basis_inverse_row;
  double* basis_inverse_row_vector = (double*)pass_basis_inverse_row_vector;
  if (basis_inverse_row_vector == NULL) {
    vector<double> rhs;
    vector<HighsInt> col_indices;
    rhs.assign(numRow, 0);
    rhs[row] = 1;
    basis_inverse_row.resize(numRow, 0);
    // Form B^{-T}e_{row}
    basisSolveInterface(rhs, &basis_inverse_row[0], NULL, NULL, true);
    basis_inverse_row_vector = &basis_inverse_row[0];
  }
  bool return_indices = row_num_nz != NULL;
  if (return_indices) *row_num_nz = 0;
  for (HighsInt col = 0; col < lp.numCol_; col++) {
    double value = 0;
    for (HighsInt el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
      HighsInt row = lp.Aindex_[el];
      value += lp.Avalue_[el] * basis_inverse_row_vector[row];
    }
    row_vector[col] = 0;
    if (fabs(value) > HIGHS_CONST_TINY) {
      if (return_indices) row_indices[(*row_num_nz)++] = col;
      row_vector[col] = value;
    }
  }
  return HighsStatus::OK;
}

HighsStatus Highs::getReducedColumn(const HighsInt col, double* col_vector,
                                    HighsInt* col_num_nz,
                                    HighsInt* col_indices) {
  if (!haveHmo("getReducedColumn")) return HighsStatus::Error;
  // Ensure that the LP is column-wise
  setOrientation(lp_);
  HighsLp& lp = lp_;
  if (col_vector == NULL) {
    highsLogUser(options_.log_options, HighsLogType::ERROR,
                 "getReducedColumn: col_vector is NULL\n");
    return HighsStatus::Error;
  }
  // col_indices can be NULL - it's the trigger that determines
  // whether they are identified or not
  if (col < 0 || col >= lp.numCol_) {
    highsLogUser(options_.log_options, HighsLogType::ERROR,
                 "Column index %" HIGHSINT_FORMAT
                 " out of range [0, %" HIGHSINT_FORMAT
                 "] in getReducedColumn\n",
                 col, lp.numCol_ - 1);
    return HighsStatus::Error;
  }
  bool has_invert = hmos_[0].ekk_instance_.simplex_lp_status_.has_invert;
  if (!has_invert) {
    highsLogUser(options_.log_options, HighsLogType::ERROR,
                 "No invertible representation for getReducedColumn\n");
    return HighsStatus::Error;
  }
  HighsInt numRow = lp.numRow_;
  vector<double> rhs;
  rhs.assign(numRow, 0);
  for (HighsInt el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++)
    rhs[lp.Aindex_[el]] = lp.Avalue_[el];
  basisSolveInterface(rhs, col_vector, col_num_nz, col_indices, false);
  return HighsStatus::OK;
}

HighsStatus Highs::setSolution(const HighsSolution& solution) {
  HighsStatus return_status = HighsStatus::OK;
  // Check if solution is valid.
  assert((HighsInt)solution_.col_value.size() != 0 ||
         (HighsInt)solution_.col_value.size() != lp_.numCol_);
  assert((HighsInt)solution.col_dual.size() == 0 ||
         (HighsInt)solution.col_dual.size() == lp_.numCol_);
  assert((HighsInt)solution.row_dual.size() == 0 ||
         (HighsInt)solution.row_dual.size() == lp_.numRow_);

  if (solution.col_value.size()) solution_.col_value = solution.col_value;
  if (solution.col_dual.size()) solution_.col_dual = solution.col_dual;
  if (solution.row_dual.size()) solution_.row_dual = solution.row_dual;

  if (solution.col_value.size() > 0) {
    return_status = interpretCallStatus(calculateRowValues(lp_, solution_),
                                        return_status, "calculateRowValues");
    if (return_status == HighsStatus::Error) return return_status;
  }
  if (solution.row_dual.size() > 0) {
    return_status = interpretCallStatus(calculateColDuals(lp_, solution_),
                                        return_status, "calculateColDuals");
    if (return_status == HighsStatus::Error) return return_status;
  }
  return returnFromHighs(return_status);
}

HighsStatus Highs::setBasis(const HighsBasis& basis) {
  // Check the user-supplied basis
  if (!isBasisConsistent(lp_, basis)) {
    highsLogUser(options_.log_options, HighsLogType::ERROR,
                 "setBasis: invalid basis\n");
    return HighsStatus::Error;
  }
  // Update the HiGHS basis
  basis_ = basis;
  basis_.valid_ = true;
  // Follow implications of a new HiGHS basis
  newHighsBasis();
  // Can't use returnFromHighs since...
  return HighsStatus::OK;
}

HighsStatus Highs::setBasis() {
  // Invalidate the basis for HiGHS Don't set to logical basis since
  // that causes presolve to be skipped
  basis_.valid_ = false;
  // Follow implications of a new HiGHS basis
  newHighsBasis();
  // Can't use returnFromHighs since...
  return HighsStatus::OK;
}

bool Highs::addRow(const double lower_bound, const double upper_bound,
                   const HighsInt num_new_nz, const HighsInt* indices,
                   const double* values) {
  HighsInt starts = 0;
  return addRows(1, &lower_bound, &upper_bound, num_new_nz, &starts, indices,
                 values);
}

bool Highs::addRows(const HighsInt num_new_row, const double* lower_bounds,
                    const double* upper_bounds, const HighsInt num_new_nz,
                    const HighsInt* starts, const HighsInt* indices,
                    const double* values) {
  HighsStatus return_status = HighsStatus::OK;
  // Check that there is a HighsModelObject
  if (!haveHmo("addRows")) return false;
  return_status = interpretCallStatus(
      addRowsInterface(num_new_row, lower_bounds, upper_bounds, num_new_nz,
                       starts, indices, values),
      return_status, "addRows");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::addCol(const double cost, const double lower_bound,
                   const double upper_bound, const HighsInt num_new_nz,
                   const HighsInt* indices, const double* values) {
  HighsInt starts = 0;
  return addCols(1, &cost, &lower_bound, &upper_bound, num_new_nz, &starts,
                 indices, values);
}

bool Highs::addCols(const HighsInt num_new_col, const double* costs,
                    const double* lower_bounds, const double* upper_bounds,
                    const HighsInt num_new_nz, const HighsInt* starts,
                    const HighsInt* indices, const double* values) {
  HighsStatus return_status = HighsStatus::OK;
  if (!haveHmo("addCols")) return false;
  return_status = interpretCallStatus(
      addColsInterface(num_new_col, costs, lower_bounds, upper_bounds,
                       num_new_nz, starts, indices, values),
      return_status, "addCols");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::changeObjectiveSense(const ObjSense sense) {
  HighsStatus return_status = HighsStatus::OK;
  if (!haveHmo("changeObjectiveSense")) return false;
  HighsStatus call_status;
  call_status = changeObjectiveSenseInterface(sense);
  return_status =
      interpretCallStatus(call_status, return_status, "changeObjectiveSense");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::changeColCost(const HighsInt col, const double cost) {
  return changeColsCost(1, &col, &cost);
}

bool Highs::changeColsCost(const HighsInt from_col, const HighsInt to_col,
                           const double* cost) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numCol_;
  index_collection.is_interval_ = true;
  index_collection.from_ = from_col;
  index_collection.to_ = to_col;
  if (!haveHmo("changeColsCost")) return false;
  call_status = changeCostsInterface(index_collection, cost);
  return_status =
      interpretCallStatus(call_status, return_status, "changeCosts");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::changeColsCost(const HighsInt num_set_entries, const HighsInt* set,
                           const double* cost) {
  if (num_set_entries <= 0) return true;
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Create a local set that is not const since index_collection.set_
  // cannot be const as it may change if the set is not ordered
  vector<HighsInt> local_set{set, set + num_set_entries};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numCol_;
  index_collection.is_set_ = true;
  index_collection.set_ = &local_set[0];
  index_collection.set_num_entries_ = num_set_entries;
  if (!haveHmo("changeColsCost")) return false;
  call_status = changeCostsInterface(index_collection, cost);
  return_status =
      interpretCallStatus(call_status, return_status, "changeCosts");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::changeColsCost(const HighsInt* mask, const double* cost) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Create a local mask that is not const since
  // index_collection.mask_ cannot be const as it changes when
  // deleting rows/columns
  vector<HighsInt> local_mask{mask, mask + lp_.numCol_};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numCol_;
  index_collection.is_mask_ = true;
  index_collection.mask_ = &local_mask[0];
  if (!haveHmo("changeColsCost")) return false;
  call_status = changeCostsInterface(index_collection, cost);
  return_status =
      interpretCallStatus(call_status, return_status, "changeCosts");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::changeColBounds(const HighsInt col, const double lower,
                            const double upper) {
  return changeColsBounds(1, &col, &lower, &upper);
}

bool Highs::changeColsBounds(const HighsInt from_col, const HighsInt to_col,
                             const double* lower, const double* upper) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numCol_;
  index_collection.is_interval_ = true;
  index_collection.from_ = from_col;
  index_collection.to_ = to_col;
  if (!haveHmo("changeColsBounds")) return false;
  call_status = changeColBoundsInterface(index_collection, lower, upper);
  return_status =
      interpretCallStatus(call_status, return_status, "changeColBounds");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::changeColsBounds(const HighsInt num_set_entries,
                             const HighsInt* set, const double* lower,
                             const double* upper) {
  if (num_set_entries <= 0) return true;
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Create a local set that is not const since index_collection.set_
  // cannot be const as it may change if the set is not ordered
  vector<HighsInt> local_set{set, set + num_set_entries};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numCol_;
  index_collection.is_set_ = true;
  index_collection.set_ = &local_set[0];
  index_collection.set_num_entries_ = num_set_entries;
  if (!haveHmo("changeColsBounds")) return false;
  call_status = changeColBoundsInterface(index_collection, lower, upper);
  return_status =
      interpretCallStatus(call_status, return_status, "changeColBounds");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::changeColsBounds(const HighsInt* mask, const double* lower,
                             const double* upper) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Create a local mask that is not const since
  // index_collection.mask_ cannot be const as it changes when
  // deleting rows/columns
  vector<HighsInt> local_mask{mask, mask + lp_.numCol_};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numCol_;
  index_collection.is_mask_ = true;
  index_collection.mask_ = &local_mask[0];
  if (!haveHmo("changeColsBounds")) return false;
  call_status = changeColBoundsInterface(index_collection, lower, upper);
  return_status =
      interpretCallStatus(call_status, return_status, "changeColBounds");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::changeRowBounds(const HighsInt row, const double lower,
                            const double upper) {
  return changeRowsBounds(1, &row, &lower, &upper);
}

bool Highs::changeRowsBounds(const HighsInt from_row, const HighsInt to_row,
                             const double* lower, const double* upper) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numRow_;
  index_collection.is_interval_ = true;
  index_collection.from_ = from_row;
  index_collection.to_ = to_row;
  if (!haveHmo("changeRowsBounds")) return false;
  call_status = changeRowBoundsInterface(index_collection, lower, upper);
  return_status =
      interpretCallStatus(call_status, return_status, "changeRowBounds");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::changeRowsBounds(const HighsInt num_set_entries,
                             const HighsInt* set, const double* lower,
                             const double* upper) {
  if (num_set_entries <= 0) return true;
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Create a local set that is not const since index_collection.set_
  // cannot be const as it may change if the set is not ordered
  vector<HighsInt> local_set{set, set + num_set_entries};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numRow_;
  index_collection.is_set_ = true;
  index_collection.set_ = &local_set[0];
  index_collection.set_num_entries_ = num_set_entries;
  if (!haveHmo("changeRowsBounds")) return false;
  call_status = changeRowBoundsInterface(index_collection, lower, upper);
  return_status =
      interpretCallStatus(call_status, return_status, "changeRowBounds");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::changeRowsBounds(const HighsInt* mask, const double* lower,
                             const double* upper) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Create a local mask that is not const since
  // index_collection.mask_ cannot be const as it changes when
  // deleting rows/columns
  vector<HighsInt> local_mask{mask, mask + lp_.numRow_};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numRow_;
  index_collection.is_mask_ = true;
  index_collection.mask_ = &local_mask[0];
  if (!haveHmo("changeRowsBounds")) return false;
  call_status = changeRowBoundsInterface(index_collection, lower, upper);
  return_status =
      interpretCallStatus(call_status, return_status, "changeRowBounds");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::changeCoeff(const HighsInt row, const HighsInt col,
                        const double value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  if (!haveHmo("changeCoeff")) return false;
  call_status = changeCoefficientInterface(row, col, value);
  return_status =
      interpretCallStatus(call_status, return_status, "changeCoefficient");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::getObjectiveSense(ObjSense& sense) {
  if (!haveHmo("getObjectiveSense")) return false;
  sense = lp_.sense_;
  return true;
}

bool Highs::getCols(const HighsInt from_col, const HighsInt to_col,
                    HighsInt& num_col, double* costs, double* lower,
                    double* upper, HighsInt& num_nz, HighsInt* start,
                    HighsInt* index, double* value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numCol_;
  index_collection.is_interval_ = true;
  index_collection.from_ = from_col;
  index_collection.to_ = to_col;
  if (!haveHmo("getCols")) return false;
  call_status = getColsInterface(index_collection, num_col, costs, lower, upper,
                                 num_nz, start, index, value);
  return_status = interpretCallStatus(call_status, return_status, "getCols");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::getCols(const HighsInt num_set_entries, const HighsInt* set,
                    HighsInt& num_col, double* costs, double* lower,
                    double* upper, HighsInt& num_nz, HighsInt* start,
                    HighsInt* index, double* value) {
  if (num_set_entries <= 0) return true;
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Create a local set that is not const since index_collection.set_
  // cannot be const as it may change if the set is not ordered
  vector<HighsInt> local_set{set, set + num_set_entries};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numCol_;
  index_collection.is_set_ = true;
  index_collection.set_ = &local_set[0];
  index_collection.set_num_entries_ = num_set_entries;
  if (!haveHmo("getCols")) return false;
  call_status = getColsInterface(index_collection, num_col, costs, lower, upper,
                                 num_nz, start, index, value);
  return_status = interpretCallStatus(call_status, return_status, "getCols");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::getCols(const HighsInt* mask, HighsInt& num_col, double* costs,
                    double* lower, double* upper, HighsInt& num_nz,
                    HighsInt* start, HighsInt* index, double* value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Create a local mask that is not const since
  // index_collection.mask_ cannot be const as it changes when
  // deleting rows/columns
  vector<HighsInt> local_mask{mask, mask + lp_.numCol_};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numCol_;
  index_collection.is_mask_ = true;
  index_collection.mask_ = &local_mask[0];
  if (!haveHmo("getCols")) return false;
  call_status = getColsInterface(index_collection, num_col, costs, lower, upper,
                                 num_nz, start, index, value);
  return_status = interpretCallStatus(call_status, return_status, "getCols");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::getRows(const HighsInt from_row, const HighsInt to_row,
                    HighsInt& num_row, double* lower, double* upper,
                    HighsInt& num_nz, HighsInt* start, HighsInt* index,
                    double* value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numRow_;
  index_collection.is_interval_ = true;
  index_collection.from_ = from_row;
  index_collection.to_ = to_row;
  if (!haveHmo("getRows")) return false;
  call_status = getRowsInterface(index_collection, num_row, lower, upper,
                                 num_nz, start, index, value);
  return_status = interpretCallStatus(call_status, return_status, "getRows");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::getRows(const HighsInt num_set_entries, const HighsInt* set,
                    HighsInt& num_row, double* lower, double* upper,
                    HighsInt& num_nz, HighsInt* start, HighsInt* index,
                    double* value) {
  if (num_set_entries <= 0) return true;
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Create a local set that is not const since index_collection.set_
  // cannot be const as it may change if the set is not ordered
  vector<HighsInt> local_set{set, set + num_set_entries};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numRow_;
  index_collection.is_set_ = true;
  index_collection.set_ = &local_set[0];
  index_collection.set_num_entries_ = num_set_entries;
  if (!haveHmo("getRows")) return false;
  call_status = getRowsInterface(index_collection, num_row, lower, upper,
                                 num_nz, start, index, value);
  return_status = interpretCallStatus(call_status, return_status, "getRows");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::getRows(const HighsInt* mask, HighsInt& num_row, double* lower,
                    double* upper, HighsInt& num_nz, HighsInt* start,
                    HighsInt* index, double* value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Create a local mask that is not const since
  // index_collection.mask_ cannot be const as it changes when
  // deleting rows/columns
  vector<HighsInt> local_mask{mask, mask + lp_.numRow_};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numRow_;
  index_collection.is_mask_ = true;
  index_collection.mask_ = &local_mask[0];
  if (!haveHmo("getRows")) return false;
  call_status = getRowsInterface(index_collection, num_row, lower, upper,
                                 num_nz, start, index, value);
  return_status = interpretCallStatus(call_status, return_status, "getRows");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::getCoeff(const HighsInt row, const HighsInt col, double& value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  if (!haveHmo("getCoeff")) return false;
  call_status = getCoefficientInterface(row, col, value);
  return_status =
      interpretCallStatus(call_status, return_status, "getCoefficient");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::deleteCols(const HighsInt from_col, const HighsInt to_col) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numCol_;
  index_collection.is_interval_ = true;
  index_collection.from_ = from_col;
  index_collection.to_ = to_col;
  if (!haveHmo("deleteCols")) return false;
  call_status = deleteColsInterface(index_collection);
  return_status = interpretCallStatus(call_status, return_status, "deleteCols");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::deleteCols(const HighsInt num_set_entries, const HighsInt* set) {
  if (num_set_entries <= 0) return true;
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Create a local set that is not const since index_collection.set_
  // cannot be const as it may change if the set is not ordered
  vector<HighsInt> local_set{set, set + num_set_entries};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numCol_;
  index_collection.is_set_ = true;
  index_collection.set_ = &local_set[0];
  index_collection.set_num_entries_ = num_set_entries;
  if (!haveHmo("deleteCols")) return false;
  call_status = deleteColsInterface(index_collection);
  return_status = interpretCallStatus(call_status, return_status, "deleteCols");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::deleteCols(HighsInt* mask) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numCol_;
  index_collection.is_mask_ = true;
  index_collection.mask_ = &mask[0];
  if (!haveHmo("deleteCols")) return false;
  call_status = deleteColsInterface(index_collection);
  return_status = interpretCallStatus(call_status, return_status, "deleteCols");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::deleteRows(const HighsInt from_row, const HighsInt to_row) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numRow_;
  index_collection.is_interval_ = true;
  index_collection.from_ = from_row;
  index_collection.to_ = to_row;
  if (!haveHmo("deleteRows")) return false;
  call_status = deleteRowsInterface(index_collection);
  return_status = interpretCallStatus(call_status, return_status, "deleteRows");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::deleteRows(const HighsInt num_set_entries, const HighsInt* set) {
  if (num_set_entries <= 0) return true;
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Create a local set that is not const since index_collection.set_
  // cannot be const as it may change if the set is not ordered
  vector<HighsInt> local_set{set, set + num_set_entries};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numRow_;
  index_collection.is_set_ = true;
  index_collection.set_ = &local_set[0];
  index_collection.set_num_entries_ = num_set_entries;
  if (!haveHmo("deleteRows")) return false;
  call_status = deleteRowsInterface(index_collection);
  return_status = interpretCallStatus(call_status, return_status, "deleteRows");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::deleteRows(HighsInt* mask) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numRow_;
  index_collection.is_mask_ = true;
  index_collection.mask_ = &mask[0];
  if (!haveHmo("deleteRows")) return false;
  call_status = deleteRowsInterface(index_collection);
  return_status = interpretCallStatus(call_status, return_status, "deleteRows");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::scaleCol(const HighsInt col, const double scaleval) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  if (!haveHmo("scaleCol")) return false;
  call_status = scaleColInterface(col, scaleval);
  return_status = interpretCallStatus(call_status, return_status, "scaleCol");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::scaleRow(const HighsInt row, const double scaleval) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  if (!haveHmo("scaleRow")) return false;
  call_status = scaleRowInterface(row, scaleval);
  return_status = interpretCallStatus(call_status, return_status, "scaleRow");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

double Highs::getHighsInfinity() { return HIGHS_CONST_INF; }

double Highs::getHighsRunTime() { return timer_.readRunHighsClock(); }

HighsStatus Highs::clearSolver() {
  clearModelStatus();
  clearSolution();
  clearBasis();
  clearInfo();
  return HighsStatus::OK;
}

#ifdef HiGHSDEV
void Highs::reportModelStatusSolutionBasis(const std::string message,
                                           const HighsInt hmo_ix) {
  HighsModelStatus& model_status = model_status_;
  HighsModelStatus& scaled_model_status = scaled_model_status_;
  HighsSolution& solution = solution_;
  HighsBasis& basis = basis_;
  HighsInt unscaled_primal_status = info_.primal_status;
  HighsInt unscaled_dual_status = info_.dual_status;
  HighsLp& lp = lp_;
  if (hmo_ix >= 0) {
    assert(hmo_ix < (HighsInt)hmos_.size());
    model_status = hmos_[hmo_ix].unscaled_model_status_;
    scaled_model_status = hmos_[hmo_ix].scaled_model_status_;
    solution = hmos_[hmo_ix].solution_;
    basis = hmos_[hmo_ix].basis_;
    unscaled_primal_status = hmos_[hmo_ix].solution_params_.primal_status;
    unscaled_dual_status = hmos_[hmo_ix].solution_params_.dual_status;
    lp = hmos_[hmo_ix].lp_;
  }
  printf(
      "\n%s\nModel status = %s; Scaled model status = %s; LP(%" HIGHSINT_FORMAT
      ", %" HIGHSINT_FORMAT
      "); solution "
      "([%" HIGHSINT_FORMAT "] %" HIGHSINT_FORMAT ", %" HIGHSINT_FORMAT
      "; [%" HIGHSINT_FORMAT "] %" HIGHSINT_FORMAT ", %" HIGHSINT_FORMAT
      "); basis %" HIGHSINT_FORMAT
      " "
      "(%" HIGHSINT_FORMAT ", %" HIGHSINT_FORMAT ")\n\n",
      message.c_str(), modelStatusToString(model_status).c_str(),
      modelStatusToString(scaled_model_status).c_str(), lp.numCol_, lp.numRow_,
      unscaled_primal_status, (HighsInt)solution.col_value.size(),
      (HighsInt)solution.row_value.size(), unscaled_dual_status,
      (HighsInt)solution.col_dual.size(), (HighsInt)solution.row_dual.size(),
      basis.valid_, (HighsInt)basis.col_status.size(),
      (HighsInt)basis.row_status.size());
}
#endif

std::string Highs::modelStatusToString(
    const HighsModelStatus model_status) const {
  return utilModelStatusToString(model_status);
}

std::string Highs::primalDualStatusToString(const HighsInt primal_dual_status) {
  return utilPrimalDualStatusToString(primal_dual_status);
}

void Highs::setMatrixOrientation(const MatrixOrientation& desired_orientation) {
  setOrientation(lp_, desired_orientation);
}

// Private methods
HighsPresolveStatus Highs::runPresolve() {
  presolve_.clear();
  // Exit if the problem is empty or if presolve is set to off.
  if (options_.presolve == off_string) return HighsPresolveStatus::NotPresolved;

  // Ensure that the LP is column-wise
  // setOrientation(lp_);

  if (lp_.numCol_ == 0 && lp_.numRow_ == 0)
    return HighsPresolveStatus::NullError;

  // Clear info from previous runs if lp_ has been modified.
  double start_presolve = timer_.readRunHighsClock();

  // Set time limit.
  if (options_.time_limit > 0 && options_.time_limit < HIGHS_CONST_INF) {
    double left = options_.time_limit - start_presolve;
    if (left <= 0) {
      highsLogDev(options_.log_options, HighsLogType::ERROR,
                  "Time limit reached while reading in matrix\n");
      return HighsPresolveStatus::Timeout;
    }

    highsLogDev(options_.log_options, HighsLogType::VERBOSE,
                "Time limit set: reading matrix took %.2g, presolve "
                "time left: %.2g\n",
                start_presolve, left);
  }

  // Presolve.
  presolve_.init(lp_, timer_);
  presolve_.options_ = &options_;
  if (options_.time_limit > 0 && options_.time_limit < HIGHS_CONST_INF) {
    double current = timer_.readRunHighsClock();
    double time_init = current - start_presolve;
    double left = presolve_.options_->time_limit - time_init;
    if (left <= 0) {
      highsLogDev(options_.log_options, HighsLogType::ERROR,
                  "Time limit reached while copying matrix into presolve.\n");
      return HighsPresolveStatus::Timeout;
    }

    highsLogDev(options_.log_options, HighsLogType::VERBOSE,
                "Time limit set: copying matrix took %.2g, presolve "
                "time left: %.2g\n",
                time_init, left);
  }

  HighsPresolveStatus presolve_return_status = presolve_.run();

  // Update reduction counts.
  switch (presolve_.presolve_status_) {
    case HighsPresolveStatus::Reduced: {
      HighsLp& reduced_lp = presolve_.getReducedProblem();
      presolve_.info_.n_cols_removed = lp_.numCol_ - reduced_lp.numCol_;
      presolve_.info_.n_rows_removed = lp_.numRow_ - reduced_lp.numRow_;
      presolve_.info_.n_nnz_removed =
          (HighsInt)lp_.Avalue_.size() - (HighsInt)reduced_lp.Avalue_.size();
      break;
    }
    case HighsPresolveStatus::ReducedToEmpty: {
      presolve_.info_.n_cols_removed = lp_.numCol_;
      presolve_.info_.n_rows_removed = lp_.numRow_;
      presolve_.info_.n_nnz_removed = (HighsInt)lp_.Avalue_.size();
      break;
    }
    default:
      break;
  }
  return presolve_return_status;
}

HighsPostsolveStatus Highs::runPostsolve() {
  // assert(presolve_.has_run_);
  bool solution_ok = isSolutionRightSize(presolve_.getReducedProblem(),
                                         presolve_.data_.recovered_solution_);
  if (!solution_ok) return HighsPostsolveStatus::ReducedSolutionDimenionsError;

  presolve_.data_.postSolveStack.undo(options_,
                                      presolve_.data_.recovered_solution_,
                                      presolve_.data_.recovered_basis_);

  if (lp_.sense_ == ObjSense::MAXIMIZE) presolve_.negateReducedLpColDuals(true);

  return HighsPostsolveStatus::SolutionRecovered;
}

// The method below runs calls solveLp to solve the LP associated with
// a particular model, integrating the iteration counts into the
// overall values in HighsInfo
HighsStatus Highs::callSolveLp(const HighsInt model_index,
                               const string message) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;

  // Check that the model index is OK
  bool model_index_ok =
      model_index >= 0 && model_index < (HighsInt)hmos_.size();
  assert(model_index_ok);
  if (!model_index_ok) return HighsStatus::Error;

  HighsModelObject& model = hmos_[model_index];
  // Check that the model isn't row-wise
  assert(model.lp_.orientation_ != MatrixOrientation::ROWWISE);

  // Transfer the LP solver iteration counts to this model
  HighsIterationCounts& iteration_counts = hmos_[model_index].iteration_counts_;
  copyHighsIterationCounts(info_, iteration_counts);

  // Solve the LP
  call_status = solveLp(model, message);
  return_status = interpretCallStatus(call_status, return_status, "solveLp");
  if (return_status == HighsStatus::Error) return return_status;

  // Transfer this model's LP solver iteration counts to HiGHS
  copyHighsIterationCounts(iteration_counts, info_);

  return return_status;
}

HighsStatus Highs::callSolveMip() {
  HighsStatus return_status = HighsStatus::OK;
  // Run the MIP solver
  options_.log_dev_level = LOG_DEV_LEVEL_INFO;
  // Check that the model isn't row-wise
  assert(lp_.orientation_ != MatrixOrientation::ROWWISE);
  HighsMipSolver solver(options_, lp_);
  solver.run();
  // Cheating now, but need to set this honestly!
  HighsStatus call_status = HighsStatus::OK;
  return_status =
      interpretCallStatus(call_status, return_status, "HighsMipSolver::solver");
  if (return_status == HighsStatus::Error) return return_status;
  // Cheating now, but need to set this honestly!
  scaled_model_status_ = HighsModelStatus::OPTIMAL;
  scaled_model_status_ = solver.modelstatus_;
  model_status_ = scaled_model_status_;
  // Set the values in HighsInfo instance info_
  info_.mip_node_count = solver.node_count_;
  info_.simplex_iteration_count = -1;    // Not known
  info_.ipm_iteration_count = -1;        // Not known
  info_.crossover_iteration_count = -1;  // Not known
  info_.primal_status = PrimalDualStatus::STATUS_NO_SOLUTION;
  info_.dual_status = PrimalDualStatus::STATUS_NO_SOLUTION;
  info_.objective_function_value = solver.solution_objective_;
  info_.mip_dual_bound = solver.dual_bound_;
  info_.mip_gap =
      100 * std::abs(info_.objective_function_value - info_.mip_dual_bound) /
      std::max(1.0, std::abs(info_.objective_function_value));
  info_.num_primal_infeasibilities = -1;  // Not known
  // Are the violations max or sum?
  info_.max_primal_infeasibility =
      std::max({solver.bound_violation_, solver.row_violation_,
                solver.integrality_violation_});
  info_.sum_primal_infeasibilities = -1;  // Not known
  info_.num_dual_infeasibilities = -1;    // Not known
  info_.max_dual_infeasibility = -1;      // Not known
  info_.sum_dual_infeasibilities = -1;    // Not known
  // The solution needs to be here, but just resize it for now
  if (solver.solution_objective_ != HIGHS_CONST_INF) {
    info_.primal_status = PrimalDualStatus::STATUS_FEASIBLE_POINT;
    HighsInt solver_solution_size = solver.solution_.size();
    assert(solver_solution_size >= lp_.numCol_);
    solution_.col_value.resize(lp_.numCol_);
    for (HighsInt iCol = 0; iCol < lp_.numCol_; iCol++)
      solution_.col_value[iCol] = solver.solution_[iCol];
  }

  //  assert((HighsInt)solution_.col_value.size() == lp_.numCol_);
  //  assert((HighsInt)solution_.row_value.size() == lp_.numRow_);
  //  solution_.row_value.resize(lp_.numRow_);

  return return_status;
}

HighsStatus Highs::writeSolution(const std::string filename,
                                 const bool pretty) const {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsLp lp = lp_;
  HighsBasis basis = basis_;
  HighsSolution solution = solution_;
  FILE* file;
  bool html;
  call_status = openWriteFile(filename, "writeSolution", file, html);
  return_status =
      interpretCallStatus(call_status, return_status, "openWriteFile");
  if (return_status == HighsStatus::Error) return return_status;

  //  std::cout << "warning: Feature under development" << std::endl;

  writeSolutionToFile(file, lp, basis, solution, pretty);

  //  return HighsStatus::Warning;
  return HighsStatus::OK;
}

// Actions to take if there is a new Highs basis
void Highs::newHighsBasis() {
  if (hmos_.size() > 0) {
    // Copy this basis to the HMO basis
    hmos_[0].basis_ = basis_;
    // Clear any simplex basis
    clearBasisInterface();
  }
}

// Ensure that the HiGHS solution and basis have the same size as the
// model, and that the HiGHS basis is kept up-to-date with any solved
// basis
void Highs::forceHighsSolutionBasisSize() {
  // Ensure that the HiGHS solution vectors are the right size
  solution_.col_value.resize(lp_.numCol_);
  solution_.row_value.resize(lp_.numRow_);
  solution_.col_dual.resize(lp_.numCol_);
  solution_.row_dual.resize(lp_.numRow_);
  // Ensure that the HiGHS basis vectors are the right size,
  // invalidating the basis if they aren't
  if ((HighsInt)basis_.col_status.size() != lp_.numCol_) {
    basis_.col_status.resize(lp_.numCol_);
    basis_.valid_ = false;
  }
  if ((HighsInt)basis_.row_status.size() != lp_.numRow_) {
    basis_.row_status.resize(lp_.numRow_);
    basis_.valid_ = false;
  }
}

bool Highs::getHighsModelStatusAndInfo(const HighsInt solved_hmo) {
  if (!haveHmo("getHighsModelStatusAndInfo")) return false;

  model_status_ = hmos_[solved_hmo].unscaled_model_status_;
  scaled_model_status_ = hmos_[solved_hmo].scaled_model_status_;

  HighsSolutionParams& solution_params = hmos_[solved_hmo].solution_params_;

  info_.primal_status = solution_params.primal_status;
  info_.dual_status = solution_params.dual_status;
  info_.objective_function_value = solution_params.objective_function_value;
  info_.num_primal_infeasibilities = solution_params.num_primal_infeasibility;
  info_.max_primal_infeasibility = solution_params.max_primal_infeasibility;
  info_.sum_primal_infeasibilities = solution_params.sum_primal_infeasibility;
  info_.num_dual_infeasibilities = solution_params.num_dual_infeasibility;
  info_.max_dual_infeasibility = solution_params.max_dual_infeasibility;
  info_.sum_dual_infeasibilities = solution_params.sum_dual_infeasibility;
  return true;
}

HighsStatus Highs::openWriteFile(const string filename,
                                 const string method_name, FILE*& file,
                                 bool& html) const {
  html = false;
  if (filename == "") {
    // Empty file name: use stdout
    file = stdout;
  } else {
    file = fopen(filename.c_str(), "w");
    if (file == 0) {
      highsLogUser(options_.log_options, HighsLogType::ERROR,
                   "Cannot open writeable file \"%s\" in %s\n",
                   filename.c_str(), method_name.c_str());
      return HighsStatus::Error;
    }
    const char* dot = strrchr(filename.c_str(), '.');
    if (dot && dot != filename) html = strcmp(dot + 1, "html") == 0;
  }
  return HighsStatus::OK;
}

HighsStatus Highs::getUseModelStatus(
    HighsModelStatus& use_model_status,
    const double unscaled_primal_feasibility_tolerance,
    const double unscaled_dual_feasibility_tolerance,
    const bool rerun_from_logical_basis) {
  if (model_status_ != HighsModelStatus::NOTSET) {
    use_model_status = model_status_;
  } else {
    // Handle the case where the status of the unscaled model is not set
    HighsStatus return_status = HighsStatus::OK;
    HighsStatus call_status;
    const double report = false;  // true;//
    if (unscaledOptimal(unscaled_primal_feasibility_tolerance,
                        unscaled_dual_feasibility_tolerance, report)) {
      use_model_status = HighsModelStatus::OPTIMAL;
    } else if (rerun_from_logical_basis) {
      std::string save_presolve = options_.presolve;
      basis_.valid_ = false;
      options_.presolve = on_string;
      call_status = run();
      return_status = interpretCallStatus(call_status, return_status, "run()");
      options_.presolve = save_presolve;
      if (return_status == HighsStatus::Error) return return_status;

      if (report)
        printf(
            "Unscaled model status was NOTSET: after running from logical "
            "basis it is %s\n",
            modelStatusToString(model_status_).c_str());

      if (model_status_ != HighsModelStatus::NOTSET) {
        use_model_status = model_status_;
      } else if (unscaledOptimal(unscaled_primal_feasibility_tolerance,
                                 unscaled_dual_feasibility_tolerance, report)) {
        use_model_status = HighsModelStatus::OPTIMAL;
      }
    } else {
      // Nothing to be done: use original unscaled model status
      use_model_status = model_status_;
    }
  }
  return HighsStatus::OK;
}

bool Highs::unscaledOptimal(const double unscaled_primal_feasibility_tolerance,
                            const double unscaled_dual_feasibility_tolerance,
                            const bool report) {
  if (scaled_model_status_ == HighsModelStatus::OPTIMAL) {
    const double max_primal_infeasibility = info_.max_primal_infeasibility;
    const double max_dual_infeasibility = info_.max_dual_infeasibility;
    if (report)
      printf(
          "Scaled model status is OPTIMAL: max unscaled (primal / dual) "
          "infeasibilities are (%g / %g)\n",
          max_primal_infeasibility, max_dual_infeasibility);
    if ((max_primal_infeasibility > unscaled_primal_feasibility_tolerance) ||
        (max_dual_infeasibility > unscaled_dual_feasibility_tolerance)) {
      printf(
          "Use model status of NOTSET since max unscaled (primal / dual) "
          "infeasibilities are (%g / %g)\n",
          max_primal_infeasibility, max_dual_infeasibility);
    } else {
      if (report)
        printf(
            "Set unscaled model status to OPTIMAL since unscaled "
            "infeasibilities are tolerable\n");
      return true;
    }
  }
  return false;
}

bool Highs::haveHmo(const string method_name) const {
  bool have_hmo = hmos_.size() > 0;
#ifdef HiGHSDEV
  if (!have_hmo)
    highsLogUser(options_.log_options, HighsLogType::ERROR,
                 "Method %s called without any HighsModelObject\n",
                 method_name.c_str());
#endif
  assert(have_hmo);
  return have_hmo;
}

void Highs::clearModelStatus() {
  model_status_ = HighsModelStatus::NOTSET;
  scaled_model_status_ = HighsModelStatus::NOTSET;
}

void Highs::clearSolution() {
  info_.primal_status = (HighsInt)PrimalDualStatus::STATUS_NOTSET;
  info_.dual_status = (HighsInt)PrimalDualStatus::STATUS_NOTSET;
  clearSolutionUtil(solution_);
}

void Highs::clearBasis() { clearBasisUtil(basis_); }

void Highs::clearInfo() { info_.clear(); }

// Applies checks before returning from run()
HighsStatus Highs::returnFromRun(const HighsStatus run_return_status) {
  bool have_solution = false;
  HighsStatus return_status = run_return_status;
  if (hmos_.size() == 0) {
    // No model has been loaded: ensure that the status, solution,
    // basis and info associated with any previous model are cleared
    clearSolver();
    return returnFromHighs(return_status);
  } else {
    // A model has been loaded: remove any additional HMO created when solving
    if (hmos_.size() > 1) hmos_.pop_back();
    // There should be only one entry in hmos_
    assert((HighsInt)hmos_.size() == 1);
    // Make sure that the unscaled status, solution, basis and info
    // are consistent with the scaled status
#ifdef HiGHSDEV
    reportModelStatusSolutionBasis("returnFromRun(HiGHS)");
    reportModelStatusSolutionBasis("returnFromRun(HMO_0)", 0);
#endif
    switch (scaled_model_status_) {
      // First consider the error returns
      case HighsModelStatus::NOTSET:
      case HighsModelStatus::LOAD_ERROR:
      case HighsModelStatus::MODEL_ERROR:
      case HighsModelStatus::PRESOLVE_ERROR:
      case HighsModelStatus::SOLVE_ERROR:
      case HighsModelStatus::POSTSOLVE_ERROR:
        clearSolver();
        assert(return_status == HighsStatus::Error);
        break;

      // Then consider the OK returns
      case HighsModelStatus::MODEL_EMPTY:
        clearSolution();
        clearBasis();
        clearInfo();
        assert(model_status_ == scaled_model_status_);
        assert(return_status == HighsStatus::OK);
        break;

      case HighsModelStatus::OPTIMAL:
        have_solution = true;
        // The following is an aspiration
        //        assert(info_.primal_status ==
        //                   (HighsInt)PrimalDualStatus::STATUS_FEASIBLE_POINT);
        //        assert(info_.dual_status ==
        //                   (HighsInt)PrimalDualStatus::STATUS_FEASIBLE_POINT);
        assert(model_status_ == HighsModelStatus::NOTSET ||
               model_status_ == HighsModelStatus::OPTIMAL);
        assert(return_status == HighsStatus::OK);
        break;

      case HighsModelStatus::PRIMAL_INFEASIBLE:
        clearSolution();
        // May have a basis, according to whether infeasibility was
        // detected in presolve or solve
        assert(model_status_ == scaled_model_status_);
        assert(return_status == HighsStatus::OK);
        break;

      case HighsModelStatus::PRIMAL_INFEASIBLE_OR_UNBOUNDED:
        clearSolution();
        // May have a basis, according to whether infeasibility was
        // detected in presolve or solve
        clearInfo();
        assert(model_status_ == scaled_model_status_);
        assert(return_status == HighsStatus::OK);
        break;

      case HighsModelStatus::PRIMAL_UNBOUNDED:
        clearSolution();
        // May have a basis, according to whether infeasibility was
        // detected in presolve or solve
        clearInfo();
        assert(model_status_ == scaled_model_status_);
        assert(return_status == HighsStatus::OK);
        break;

      case HighsModelStatus::PRIMAL_DUAL_INFEASIBLE:
        clearSolution();
        // May have a basis, according to whether infeasibility was
        // detected in presolve or solve
        clearInfo();
        assert(model_status_ == scaled_model_status_);
        assert(return_status == HighsStatus::OK);
        break;

      case HighsModelStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND:
        clearSolution();
        clearBasis();
        clearInfo();
        assert(model_status_ == scaled_model_status_);
        assert(return_status == HighsStatus::OK);
        break;

      // Finally consider the warning returns
      case HighsModelStatus::REACHED_TIME_LIMIT:
      case HighsModelStatus::REACHED_ITERATION_LIMIT:
        clearSolution();
        clearBasis();
        clearInfo();
        assert(model_status_ == scaled_model_status_);
        assert(return_status == HighsStatus::Warning);
        break;
      case HighsModelStatus::DUAL_INFEASIBLE:
        clearSolution();
        // May have a basis, according to whether infeasibility was
        // detected in presolve or solve
        clearInfo();
        assert(model_status_ == scaled_model_status_);
        assert(return_status == HighsStatus::Warning);
        break;
    }
  }
  if (have_solution) debugSolutionRightSize(options_, lp_, solution_);
  bool have_basis = basis_.valid_;
  if (have_basis) {
    if (debugBasisRightSize(options_, lp_, basis_) ==
        HighsDebugStatus::LOGICAL_ERROR)
      return_status = HighsStatus::Error;
  }

  if (have_solution && have_basis) {
    if (debugHighsBasicSolution("Return from run()", options_, lp_, basis_,
                                solution_, info_, model_status_) ==
        HighsDebugStatus::LOGICAL_ERROR)
      return_status = HighsStatus::Error;
  }
  return returnFromHighs(return_status);
}

// Applies checks before returning from HiGHS
HighsStatus Highs::returnFromHighs(HighsStatus highs_return_status) {
  HighsStatus return_status = highs_return_status;

  forceHighsSolutionBasisSize();

  const bool consistent = debugBasisConsistent(options_, lp_, basis_) !=
                          HighsDebugStatus::LOGICAL_ERROR;
  if (!consistent) {
    highsLogUser(
        options_.log_options, HighsLogType::ERROR,
        "returnFromHighs: Supposed to be a HiGHS basis, but not consistent\n");
    assert(consistent);
    return_status = HighsStatus::Error;
  }

  if (hmos_.size()) {
    bool simplex_lp_ok =
        ekkDebugSimplexLp(hmos_[0]) != HighsDebugStatus::LOGICAL_ERROR;
    if (!simplex_lp_ok) {
      highsLogUser(options_.log_options, HighsLogType::ERROR,
                   "returnFromHighs: Simplex LP not OK\n");
      assert(simplex_lp_ok);
      return_status = HighsStatus::Error;
    }
  }

  return return_status;
}
void Highs::underDevelopmentLogMessage(const string method_name) {
  highsLogUser(options_.log_options, HighsLogType::WARNING,
               "Method %s is still under development and behaviour may be "
               "unpredictable\n",
               method_name.c_str());
}
