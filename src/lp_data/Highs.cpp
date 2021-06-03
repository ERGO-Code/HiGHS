/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/Highs.cpp
 * @brief
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
#include "lp_data/HighsInfoDebug.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelUtils.h"
#include "lp_data/HighsSolution.h"
#include "lp_data/HighsSolve.h"
#include "mip/HighsMipSolver.h"
#include "model/HighsHessianUtils.h"
#include "simplex/HSimplexDebug.h"
#include "util/HighsMatrixPic.h"

#ifdef OPENMP
#include "omp.h"
#endif

Highs::Highs() {
  hmos_.clear();
  hmos_.push_back(HighsModelObject(model_.lp_, options_, timer_));
}

HighsStatus Highs::clear() {
  resetOptions();
  return clearModel();
}

HighsStatus Highs::clearModel() {
  model_.clear();
  return clearSolver();
}

HighsStatus Highs::clearSolver() {
  HighsStatus return_status = HighsStatus::kOk;
  clearPresolve();
  clearUserSolverData();
  hmos_.clear();
  // Clear any HighsModelObject instances and create a fresh one for
  // the incumbent model
  hmos_.push_back(HighsModelObject(model_.lp_, options_, timer_));
  return returnFromHighs(return_status);
}

HighsStatus Highs::setOptionValue(const std::string& option, const bool value) {
  if (setLocalOptionValue(options_.log_options, option, options_.records,
                          value) == OptionStatus::kOk)
    return HighsStatus::kOk;
  return HighsStatus::kError;
}

HighsStatus Highs::setOptionValue(const std::string& option,
                                  const HighsInt value) {
  if (setLocalOptionValue(options_.log_options, option, options_.records,
                          value) == OptionStatus::kOk)
    return HighsStatus::kOk;
  return HighsStatus::kError;
}

HighsStatus Highs::setOptionValue(const std::string& option,
                                  const double value) {
  if (setLocalOptionValue(options_.log_options, option, options_.records,
                          value) == OptionStatus::kOk)
    return HighsStatus::kOk;
  return HighsStatus::kError;
}

HighsStatus Highs::setOptionValue(const std::string& option,
                                  const std::string value) {
  if (setLocalOptionValue(options_.log_options, option, options_.records,
                          value) == OptionStatus::kOk)
    return HighsStatus::kOk;
  return HighsStatus::kError;
}

HighsStatus Highs::setOptionValue(const std::string& option,
                                  const char* value) {
  if (setLocalOptionValue(options_.log_options, option, options_.records,
                          value) == OptionStatus::kOk)
    return HighsStatus::kOk;
  return HighsStatus::kError;
}

HighsStatus Highs::readOptions(const std::string filename) {
  if (filename.size() <= 0) {
    highsLogUser(options_.log_options, HighsLogType::kWarning,
                 "Empty file name so not reading options\n");
    return HighsStatus::kWarning;
  }
  if (!loadOptionsFromFile(options_, filename)) return HighsStatus::kError;
  return HighsStatus::kOk;
}

HighsStatus Highs::passOptions(const HighsOptions& options) {
  if (passLocalOptions(options_.log_options, options, options_) ==
      OptionStatus::kOk)
    return HighsStatus::kOk;
  return HighsStatus::kError;
}

const HighsOptions& Highs::getOptions() { return options_; }

HighsStatus Highs::getOptionValue(const std::string& option, bool& value) {
  if (getLocalOptionValue(options_.log_options, option, options_.records,
                          value) == OptionStatus::kOk)
    return HighsStatus::kOk;
  return HighsStatus::kError;
}

HighsStatus Highs::getOptionValue(const std::string& option, HighsInt& value) {
  if (getLocalOptionValue(options_.log_options, option, options_.records,
                          value) == OptionStatus::kOk)
    return HighsStatus::kOk;
  return HighsStatus::kError;
}

HighsStatus Highs::getOptionValue(const std::string& option, double& value) {
  if (getLocalOptionValue(options_.log_options, option, options_.records,
                          value) == OptionStatus::kOk)
    return HighsStatus::kOk;
  return HighsStatus::kError;
}

HighsStatus Highs::getOptionValue(const std::string& option,
                                  std::string& value) {
  if (getLocalOptionValue(options_.log_options, option, options_.records,
                          value) == OptionStatus::kOk)
    return HighsStatus::kOk;
  return HighsStatus::kError;
}

HighsStatus Highs::getOptionType(const std::string& option,
                                 HighsOptionType& type) {
  if (getLocalOptionType(options_.log_options, option, options_.records,
                         type) == OptionStatus::kOk)
    return HighsStatus::kOk;
  return HighsStatus::kError;
}

HighsStatus Highs::resetOptions() {
  resetLocalOptions(options_.records);
  return HighsStatus::kOk;
}

HighsStatus Highs::writeOptions(const std::string filename,
                                const bool report_only_non_default_values) {
  HighsStatus return_status = HighsStatus::kOk;
  FILE* file;
  bool html;
  return_status =
      interpretCallStatus(openWriteFile(filename, "writeOptions", file, html),
                          return_status, "openWriteFile");
  if (return_status == HighsStatus::kError) return return_status;

  return_status = interpretCallStatus(
      writeOptionsToFile(file, options_.records, report_only_non_default_values,
                         html),
      return_status, "writeOptionsToFile");
  if (file != stdout) fclose(file);
  return return_status;
}

const HighsOptions& Highs::getOptions() const { return options_; }

const HighsInfo& Highs::getInfo() const { return info_; }

HighsStatus Highs::getInfoValue(const std::string& info, HighsInt& value) {
  InfoStatus status =
      getLocalInfoValue(options_, info, info_.valid, info_.records, value);
  if (status == InfoStatus::kOk) {
    return HighsStatus::kOk;
  } else if (status == InfoStatus::kUnavailable) {
    return HighsStatus::kWarning;
  } else {
    return HighsStatus::kError;
  }
}

HighsStatus Highs::getInfoValue(const std::string& info, double& value) const {
  InfoStatus status =
      getLocalInfoValue(options_, info, info_.valid, info_.records, value);
  if (status == InfoStatus::kOk) {
    return HighsStatus::kOk;
  } else if (status == InfoStatus::kUnavailable) {
    return HighsStatus::kWarning;
  } else {
    return HighsStatus::kError;
  }
}

HighsStatus Highs::writeInfo(const std::string filename) {
  HighsStatus return_status = HighsStatus::kOk;
  FILE* file;
  bool html;
  return_status =
      interpretCallStatus(openWriteFile(filename, "writeInfo", file, html),
                          return_status, "openWriteFile");
  if (return_status == HighsStatus::kError) return return_status;

  return_status = interpretCallStatus(
      writeInfoToFile(file, info_.valid, info_.records, html), return_status,
      "writeInfoToFile");
  if (file != stdout) fclose(file);
  return return_status;
}

// Methods below change the incumbent model or solver infomation
// associated with it. Hence returnFromHighs is called at the end of
// each
HighsStatus Highs::passModel(const HighsModel model) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsLp& lp = model_.lp_;
  HighsHessian& hessian = model_.hessian_;
  // Move the model's LP and Hessian to the internal LP and Hessian
  lp = std::move(model.lp_);
  hessian = std::move(model.hessian_);
  // Ensure that the LP is column-wise
  setOrientation(lp);
  // Check validity of the LP, normalising its values
  return_status =
      interpretCallStatus(assessLp(lp, options_), return_status, "assessLp");
  if (return_status == HighsStatus::kError) return return_status;
  // Check validity of any Hessian, normalising its entries
  return_status = interpretCallStatus(assessHessian(hessian, options_),
                                      return_status, "assessHessian");
  if (return_status == HighsStatus::kError) return return_status;
  if (hessian.dim_) {
    // Clear any zero Hessian
    if (hessian.q_start_[hessian.dim_] == 0) {
      highsLogUser(options_.log_options, HighsLogType::kInfo,
                   "Hessian has dimension %" HIGHSINT_FORMAT
                   " but no nonzeros, so is ignored\n",
                   hessian.dim_);
      hessian.clear();
    }
  }

  // Clear solver status, solution, basis and info associated with any
  // previous model; clear any HiGHS model object; create a HiGHS
  // model object for this LP
  return_status =
      interpretCallStatus(clearSolver(), return_status, "clearSolver");
  return returnFromHighs(return_status);
}

HighsStatus Highs::passModel(const HighsLp lp) {
  HighsModel model;
  model.lp_ = std::move(lp);
  return passModel(model);
}

HighsStatus Highs::passModel(const HighsInt num_col, const HighsInt num_row,
                             const HighsInt num_nz, const bool rowwise,
                             const HighsInt hessian_num_nz, const double* costs,
                             const double* col_lower, const double* col_upper,
                             const double* row_lower, const double* row_upper,
                             const HighsInt* astart, const HighsInt* aindex,
                             const double* avalue, const HighsInt* q_start,
                             const HighsInt* q_index, const double* q_value,
                             const HighsInt* integrality) {
  HighsModel model;
  HighsLp& lp = model.lp_;
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
    if (rowwise) {
      lp.Astart_.assign(astart, astart + num_row);
    } else {
      lp.Astart_.assign(astart, astart + num_col);
    }
    lp.Aindex_.assign(aindex, aindex + num_nz);
    lp.Avalue_.assign(avalue, avalue + num_nz);
  }
  if (rowwise) {
    lp.Astart_.resize(num_row + 1);
    lp.Astart_[num_row] = num_nz;
    lp.orientation_ = MatrixOrientation::kRowwise;
  } else {
    lp.Astart_.resize(num_col + 1);
    lp.Astart_[num_col] = num_nz;
    lp.orientation_ = MatrixOrientation::kColwise;
  }
  if (num_col > 0 && integrality != NULL) {
    lp.integrality_.resize(num_col);
    for (HighsInt iCol = 0; iCol < num_col; iCol++) {
      HighsInt integrality_status = integrality[iCol];
      assert(integrality_status == (HighsInt)HighsVarType::kContinuous ||
             integrality_status == (HighsInt)HighsVarType::kInteger);
      lp.integrality_[iCol] = (HighsVarType)integrality_status;
    }
  }
  if (hessian_num_nz > 0) {
    assert(num_col > 0);
    assert(q_start != NULL);
    assert(q_index != NULL);
    assert(q_value != NULL);
    HighsHessian& hessian = model.hessian_;
    hessian.dim_ = num_col;
    hessian.q_start_.assign(q_start, q_start + num_col);
    hessian.q_start_.resize(num_col + 1);
    hessian.q_start_[num_col] = hessian_num_nz;
    hessian.q_index_.assign(q_index, q_index + hessian_num_nz);
    hessian.q_value_.assign(q_value, q_value + hessian_num_nz);
  }
  return passModel(std::move(model));
}

HighsStatus Highs::passModel(const HighsInt num_col, const HighsInt num_row,
                             const HighsInt num_nz, const bool rowwise,
                             const double* costs, const double* col_lower,
                             const double* col_upper, const double* row_lower,
                             const double* row_upper, const HighsInt* astart,
                             const HighsInt* aindex, const double* avalue,
                             const HighsInt* integrality) {
  return passModel(num_col, num_row, num_nz, rowwise, 0, costs, col_lower,
                   col_upper, row_lower, row_upper, astart, aindex, avalue,
                   NULL, NULL, NULL, integrality);
}

HighsStatus Highs::readModel(const std::string filename) {
  HighsStatus return_status = HighsStatus::kOk;
  Filereader* reader = Filereader::getFilereader(filename);
  if (reader == NULL) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Model file %s not supported\n", filename.c_str());
    return HighsStatus::kError;
  }

  HighsModel model;
  FilereaderRetcode call_code =
      reader->readModelFromFile(options_, filename, model);
  delete reader;
  if (call_code != FilereaderRetcode::kOk) {
    interpretFilereaderRetcode(options_.log_options, filename.c_str(),
                               call_code);
    return_status = interpretCallStatus(HighsStatus::kError, return_status,
                                        "readModelFromFile");
    if (return_status == HighsStatus::kError) return return_status;
  }
  model.lp_.model_name_ = extractModelName(filename);
  return_status =
      interpretCallStatus(passModel(model), return_status, "passModel");
  return returnFromHighs(return_status);
}

HighsStatus Highs::readBasis(const std::string filename) {
  HighsStatus return_status = HighsStatus::kOk;
  // Try to read basis file into read_basis
  HighsBasis read_basis = basis_;
  return_status = interpretCallStatus(
      readBasisFile(options_.log_options, read_basis, filename), return_status,
      "readBasis");
  if (return_status != HighsStatus::kOk) return return_status;
  // Basis read OK: check whether it's consistent with the LP
  if (!isBasisConsistent(model_.lp_, read_basis)) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "readBasis: invalid basis\n");
    return HighsStatus::kError;
  }
  // Update the HiGHS basis and invalidate any simplex basis for the model
  basis_ = read_basis;
  basis_.valid = true;
  if (hmos_.size() > 0) {
    clearBasisInterface();
  }
  // Can't use returnFromHighs since...
  return HighsStatus::kOk;
}

HighsStatus Highs::writeModel(const std::string filename) {
  HighsStatus return_status = HighsStatus::kOk;

  // Ensure that the LP is column-wise
  setOrientation(model_.lp_);
  if (filename == "") {
    // Empty file name: report model on logging stream
    reportModel();
    return_status = HighsStatus::kOk;
  } else {
    Filereader* writer = Filereader::getFilereader(filename);
    if (writer == NULL) {
      highsLogUser(options_.log_options, HighsLogType::kError,
                   "Model file %s not supported\n", filename.c_str());
      return HighsStatus::kError;
    }
    return_status = interpretCallStatus(
        writer->writeModelToFile(options_, filename, model_), return_status,
        "writeModelToFile");
    delete writer;
  }
  return returnFromHighs(return_status);
}

HighsStatus Highs::writeBasis(const std::string filename) {
  HighsStatus return_status = HighsStatus::kOk;
  return_status = interpretCallStatus(
      writeBasisFile(options_.log_options, basis_, filename), return_status,
      "writeBasis");
  return returnFromHighs(return_status);
}

// Checks the options calls presolve and postsolve if needed. Solvers are called
// with callSolveLp(..)
HighsStatus Highs::run() {
  if (!haveHmo("run")) return HighsStatus::kError;
  // Ensure that there is exactly one Highs model object
  assert((HighsInt)hmos_.size() == 1);
  HighsInt min_highs_debug_level = kHighsDebugLevelMin;
  //    kHighsDebugLevelCostly;
#ifdef HiGHSDEV
  min_highs_debug_level =  // kHighsDebugLevelMin;
                           //  kHighsDebugLevelCheap;
      kHighsDebugLevelCostly;
  //  kHighsDebugLevelExpensive;
  //  kHighsDebugLevelMax;
  if (options_.highs_debug_level < min_highs_debug_level)
    highsLogDev(options_.log_options, HighsLogType::kWarning,
                "Highs::run() HiGHSDEV defined, so switching "
                "options_.highs_debug_level "
                "from %" HIGHSINT_FORMAT " to %" HIGHSINT_FORMAT "\n",
                options_.highs_debug_level, min_highs_debug_level);
    //  writeModel("HighsRunModel.mps");
    //  if (model_.lp_.numRow_>0 && model_.lp_.numCol_>0)
    //  writeLpMatrixPicToFile(options_, "LpMatrix", model_.lp_);
#endif
  if (options_.highs_debug_level < min_highs_debug_level)
    options_.highs_debug_level = min_highs_debug_level;

#ifdef OPENMP
  omp_max_threads = omp_get_max_threads();
  assert(omp_max_threads > 0);
#ifdef HiGHSDEV
  if (omp_max_threads <= 0)
    highsLogDev(options_.log_options, HighsLogType::kWarning,
                "WARNING: omp_get_max_threads() returns %" HIGHSINT_FORMAT "\n",
                omp_max_threads);

  highsLogDev(options_.log_options, HighsLogType::kDetailed,
              "Running with %" HIGHSINT_FORMAT " OMP thread(s)\n",
              omp_max_threads);
#endif
#endif
  assert(called_return_from_run);
  if (!called_return_from_run) {
    highsLogDev(options_.log_options, HighsLogType::kError,
                "Highs::run() called with called_return_from_run false\n");
    return HighsStatus::kError;
  }
  // Set this so that calls to returnFromRun() can be checked
  called_return_from_run = false;
  // From here all return statements execute returnFromRun()
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  // Initialise the HiGHS model status values
  hmos_[0].scaled_model_status_ = HighsModelStatus::kNotset;
  hmos_[0].unscaled_model_status_ = HighsModelStatus::kNotset;
  model_status_ = hmos_[0].scaled_model_status_;
  scaled_model_status_ = hmos_[0].unscaled_model_status_;
  // Clear the run info
  clearInfo();
  // Zero the HiGHS iteration counts
  zeroHighsIterationCounts(iteration_counts_);
  // Start the HiGHS run clock
  timer_.startRunHighsClock();
  // Return immediately if the model has no columns
  if (!model_.lp_.numCol_) {
    setHighsModelStatusAndInfo(HighsModelStatus::kModelEmpty);
    return returnFromRun(HighsStatus::kOk);
  }
  // Return immediately if the model is infeasible due to inconsistent bounds
  if (isBoundInfeasible(options_.log_options, model_.lp_)) {
    setHighsModelStatusAndInfo(HighsModelStatus::kInfeasible);
    return returnFromRun(return_status);
  }
  if (!options_.solver.compare(kHighsChooseString) && model_.isQp()) {
    // Solve the model as a QP
    call_status = callSolveQp();
    return_status =
        interpretCallStatus(call_status, return_status, "callSolveQp");
    return returnFromRun(return_status);
  }

  // Ensure that the LP (and any simplex LP) has the matrix column-wise
  setOrientation(model_.lp_);
  if (hmos_[0].ekk_instance_.status_.valid)
    setOrientation(hmos_[0].ekk_instance_.lp_);
#ifdef HIGHSDEV
  // Shouldn't have to check validity of the LP since this is done when it is
  // loaded or modified
  call_status = assessLp(model_.lp_, options_);
  // If any errors have been found or normalisation carried out,
  // call_status will be kError or kWarning, so only valid return is OK.
  assert(call_status == HighsStatus::kOk);
  return_status = interpretCallStatus(call_status, return_status, "assessLp");
  if (return_status == HighsStatus::kError) return returnFromRun(return_status);
#endif

  highsSetLogCallback(options_);
#ifdef HiGHSDEV
  if (checkOptions(options_.log_options, options_.records) !=
      OptionStatus::kOk) {
    return_status = HighsStatus::kError;
    return returnFromRun(return_status);
  }
#endif
  if (model_.lp_.model_name_.compare(""))
    highsLogDev(options_.log_options, HighsLogType::kVerbose,
                "Solving model: %s\n", model_.lp_.model_name_.c_str());

  if (!options_.solver.compare(kHighsChooseString) && model_.isMip()) {
    // Solve the model as a MIP
    call_status = callSolveMip();
    return_status =
        interpretCallStatus(call_status, return_status, "callSolveMip");
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

  if (basis_.valid || options_.presolve == kHighsOffString) {
    // There is a valid basis for the problem or presolve is off
    solved_hmo = original_hmo;
    hmos_[solved_hmo].ekk_instance_.lp_name_ =
        "LP without presolve or with basis";
    if (basis_.valid) {
      // There is a valid HiGHS basis, so use it to initialise the basis
      // in the HMO to be solved after refining any status values that
      // are simply HighsBasisStatus::kNonbasic
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
    if (return_status == HighsStatus::kError)
      return returnFromRun(return_status);
  } else {
    // No HiGHS basis so consider presolve
    //
    // If using IPX to solve the reduced LP, crossover must be run
    // since a basic solution is required by postsolve
    if (options_.solver == kIpmString && !options_.run_crossover) {
      highsLogUser(options_.log_options, HighsLogType::kWarning,
                   "Forcing IPX to use crossover after presolve\n");
      options_.run_crossover = true;
    }
    // Possibly presolve - according to option_.presolve
    const double from_presolve_time = timer_.read(timer_.presolve_clock);
    this_presolve_time = -from_presolve_time;
    timer_.start(timer_.presolve_clock);
    model_presolve_status_ = runPresolve();
    timer_.stop(timer_.presolve_clock);
    const double to_presolve_time = timer_.read(timer_.presolve_clock);
    this_presolve_time += to_presolve_time;
    presolve_.info_.presolve_time = this_presolve_time;

    // Set an illegal local pivot threshold value that's updated after
    // solving the presolved LP - if simplex is used
    double factor_pivot_threshold = -1;

    // Run solver.
    bool have_optimal_solution = false;
    switch (model_presolve_status_) {
      case HighsPresolveStatus::kNotPresolved: {
        hmos_[solved_hmo].ekk_instance_.lp_name_ = "Original LP";
        this_solve_original_lp_time = -timer_.read(timer_.solve_clock);
        timer_.start(timer_.solve_clock);
        call_status = callSolveLp(solved_hmo, "Not presolved: solving the LP");
        timer_.stop(timer_.solve_clock);
        this_solve_original_lp_time += timer_.read(timer_.solve_clock);
        return_status =
            interpretCallStatus(call_status, return_status, "callSolveLp");
        if (return_status == HighsStatus::kError)
          return returnFromRun(return_status);
        break;
      }
      case HighsPresolveStatus::kNotReduced: {
        hmos_[solved_hmo].ekk_instance_.lp_name_ = "Unreduced LP";
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
        if (return_status == HighsStatus::kError)
          return returnFromRun(return_status);
        break;
      }
      case HighsPresolveStatus::kReduced: {
        HighsLp& reduced_lp = presolve_.getReducedProblem();
        // Validate the reduced LP
        assert(assessLp(reduced_lp, options_) == HighsStatus::kOk);
        call_status = cleanBounds(options_, reduced_lp);
        // Ignore any warning from clean bounds since the original LP
        // is still solved after presolve
        if (interpretCallStatus(call_status, return_status, "cleanBounds") ==
            HighsStatus::kError)
          return HighsStatus::kError;
        // Add reduced lp object to vector of HighsModelObject,
        // so the last one in lp_ is the presolved one.

        hmos_.push_back(HighsModelObject(reduced_lp, options_, timer_));
        // Log the presolve reductions
        reportPresolveReductions(hmos_[original_hmo].options_.log_options,
                                 hmos_[original_hmo].lp_,
                                 hmos_[presolve_hmo].lp_);
        // Record the HMO to be solved
        solved_hmo = presolve_hmo;
        hmos_[solved_hmo].ekk_instance_.lp_name_ = "Presolved LP";
        // Don't try dual cut-off when solving the presolved LP, as the
        // objective values aren't correct
        //	HighsOptions& options = hmos_[solved_hmo].options_;
        //	HighsOptions save_options = options;
        const double save_objective_bound = options_.objective_bound;
        options_.objective_bound = kHighsInf;
        this_solve_presolved_lp_time = -timer_.read(timer_.solve_clock);
        timer_.start(timer_.solve_clock);
        call_status = callSolveLp(solved_hmo, "Solving the presolved LP");
        timer_.stop(timer_.solve_clock);
        this_solve_presolved_lp_time += timer_.read(timer_.solve_clock);
        if (hmos_[solved_hmo].ekk_instance_.status_.valid) {
          // Record the pivot threshold resulting from solving the presolved LP
          // with simplex
          factor_pivot_threshold =
              hmos_[solved_hmo].ekk_instance_.info_.factor_pivot_threshold;
        }
        // Restore the dual objective cut-off
        options_.objective_bound = save_objective_bound;
        return_status =
            interpretCallStatus(call_status, return_status, "callSolveLp");
        if (return_status == HighsStatus::kError)
          return returnFromRun(return_status);
        have_optimal_solution = hmos_[solved_hmo].scaled_model_status_ ==
                                HighsModelStatus::kOptimal;
        break;
      }
      case HighsPresolveStatus::kReducedToEmpty: {
        reportPresolveReductions(hmos_[original_hmo].options_.log_options,
                                 hmos_[original_hmo].lp_, true);
        // Create a trivial optimal solution for postsolve to use
        clearSolutionUtil(hmos_[original_hmo].solution_);
        clearBasisUtil(hmos_[original_hmo].basis_);
        have_optimal_solution = true;
        break;
      }
      case HighsPresolveStatus::kInfeasible: {
        setHighsModelStatusAndInfo(HighsModelStatus::kInfeasible);
        highsLogUser(options_.log_options, HighsLogType::kInfo,
                     "Problem status detected on presolve: %s\n",
                     modelStatusToString(model_status_).c_str());
        return returnFromRun(return_status);
      }
      case HighsPresolveStatus::kUnboundedOrInfeasible: {
        if (options_.allow_unbounded_or_infeasible) {
          setHighsModelStatusAndInfo(HighsModelStatus::kUnboundedOrInfeasible);
          highsLogUser(options_.log_options, HighsLogType::kInfo,
                       "Problem status detected on presolve: %s\n",
                       modelStatusToString(model_status_).c_str());
          return returnFromRun(return_status);
        }
        // Presolve has returned kUnboundedOrInfeasible, but HiGHS
        // can't reurn this. Use primal simplex solver on the original
        // LP
        std::string solver = options_.solver;
        HighsInt simplex_strategy = options_.simplex_strategy;
        options_.solver = "simplex";
        options_.simplex_strategy = kSimplexStrategyPrimal;
        this_solve_original_lp_time = -timer_.read(timer_.solve_clock);
        timer_.start(timer_.solve_clock);
        call_status = callSolveLp(original_hmo,
                                  "Solving the original LP with primal simplex "
                                  "to determine infeasible or unbounded");
        timer_.stop(timer_.solve_clock);
        this_solve_original_lp_time += timer_.read(timer_.solve_clock);
        if (return_status == HighsStatus::kError)
          return returnFromRun(return_status);
        setHighsModelStatusBasisSolutionAndInfo();
        assert(model_status_ == HighsModelStatus::kInfeasible ||
               model_status_ == HighsModelStatus::kUnbounded);
        return returnFromRun(return_status);
      }
      case HighsPresolveStatus::kTimeout: {
        setHighsModelStatusAndInfo(HighsModelStatus::kTimeLimit);
        highsLogDev(options_.log_options, HighsLogType::kError,
                    "Presolve reached timeout\n");
        return returnFromRun(HighsStatus::kWarning);
      }
      case HighsPresolveStatus::kOptionsError: {
        setHighsModelStatusAndInfo(HighsModelStatus::kPresolveError);
        highsLogDev(options_.log_options, HighsLogType::kError,
                    "Presolve options error\n");
        return returnFromRun(HighsStatus::kError);
      }
      default: {
        // case HighsPresolveStatus::kError
        setHighsModelStatusAndInfo(HighsModelStatus::kPresolveError);
        highsLogDev(options_.log_options, HighsLogType::kError,
                    "Presolve returned status %d\n",
                    (int)model_presolve_status_);
        return returnFromRun(HighsStatus::kError);
      }
    }
    // End of presolve
    assert(model_presolve_status_ == HighsPresolveStatus::kNotPresolved ||
           model_presolve_status_ == HighsPresolveStatus::kNotReduced ||
           model_presolve_status_ == HighsPresolveStatus::kReduced ||
           model_presolve_status_ == HighsPresolveStatus::kReducedToEmpty);

    // Postsolve. Does nothing if there were no reductions during presolve.

    if (have_optimal_solution) {
      assert(hmos_[solved_hmo].scaled_model_status_ ==
                 HighsModelStatus::kOptimal ||
             model_presolve_status_ == HighsPresolveStatus::kReducedToEmpty);
      if (model_presolve_status_ == HighsPresolveStatus::kReduced ||
          model_presolve_status_ == HighsPresolveStatus::kReducedToEmpty) {
        // If presolve is nontrivial, extract the optimal solution
        // and basis for the presolved problem in order to generate
        // the solution and basis for postsolve to use to generate a
        // solution(?) and basis that is, hopefully, optimal. This is
        // confirmed or corrected by hot-starting the simplex solver
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

        if (postsolve_status == HighsPostsolveStatus::kSolutionRecovered) {
          highsLogDev(options_.log_options, HighsLogType::kVerbose,
                      "Postsolve finished\n");
          //
          // Now hot-start the simplex solver for the original_hmo:
          //
          // The original model hasn't been solved, so set up its solution
          // parameters
          resetModelStatusAndSolutionParams(hmos_[original_hmo]);
          // Set solution and its status
          hmos_[original_hmo].solution_ = presolve_.data_.recovered_solution_;
          hmos_[original_hmo].solution_.value_valid = true;
          hmos_[original_hmo].solution_.dual_valid = true;

          // Set basis and its status
          hmos_[original_hmo].basis_.valid = true;
          hmos_[original_hmo].basis_.col_status =
              presolve_.data_.recovered_basis_.col_status;
          hmos_[original_hmo].basis_.row_status =
              presolve_.data_.recovered_basis_.row_status;

          // Possibly force debug to perform KKT check on what's
          // returned from postsolve
          const bool force_debug = false;
          HighsInt save_highs_debug_level = options_.highs_debug_level;
          if (force_debug) options_.highs_debug_level = kHighsDebugLevelCostly;
          if (debugHighsSolution("After returning from postsolve", options_,
                                 model_.lp_, hmos_[original_hmo].solution_,
                                 hmos_[original_hmo].basis_) ==
              HighsDebugStatus::kLogicalError)
            return returnFromRun(HighsStatus::kError);
          options_.highs_debug_level = save_highs_debug_level;

          // Now hot-start the simplex solver for the original_hmo
          solved_hmo = original_hmo;
          // Save the options to allow the best simplex strategy to
          // be used
          HighsOptions& options = hmos_[solved_hmo].options_;
          HighsOptions save_options = options;
          const bool full_logging = false;
          if (full_logging) options.log_dev_level = kHighsLogDevLevelVerbose;
          // Force the use of simplex to clean up if IPM has been used
          // to solve the presolved problem
          if (options.solver == kIpmString) options.solver = kSimplexString;
          options.simplex_strategy = kSimplexStrategyChoose;
          // Ensure that the parallel solver isn't used
          options.highs_min_threads = 1;
          options.highs_max_threads = 1;
          // Use any pivot threshold resulting from solving the presolved LP
          if (factor_pivot_threshold > 0)
            options.factor_pivot_threshold = factor_pivot_threshold;

          // The basis returned from postsolve is just basic/nonbasic
          // and EKK expects a refined basis, so set it up now
          refineBasis(model_.lp_, hmos_[original_hmo].solution_,
                      hmos_[original_hmo].basis_);

          hmos_[solved_hmo].ekk_instance_.lp_name_ = "Postsolve LP";
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
          if (return_status == HighsStatus::kError)
            return returnFromRun(return_status);
          postsolve_iteration_count =
              info_.simplex_iteration_count - iteration_count0;
        } else {
          highsLogUser(options_.log_options, HighsLogType::kError,
                       "Postsolve return status is %d\n",
                       (int)postsolve_status);
          setHighsModelStatusAndInfo(HighsModelStatus::kPostsolveError);
          return returnFromRun(HighsStatus::kError);
        }
      } else {
        // LP was not reduced by presolve, so have simply solved the original LP
        assert(model_presolve_status_ == HighsPresolveStatus::kNotReduced);
        assert(solved_hmo == original_hmo);
      }
    } else {
      // Optimal solution of presolved problem has not been found
      // The original model inherits the solved model's status
      //      hmos_[original_hmo].unscaled_model_status_ =
      //          hmos_[solved_hmo].unscaled_model_status_;
      //      hmos_[original_hmo].scaled_model_status_ =
      //          hmos_[solved_hmo].scaled_model_status_;
    }
  }
  // solved_hmo will be original_hmo unless the presolved LP is found to be
  // infeasible or unbounded, or if the time/iteration limit is reached
  if (solved_hmo != original_hmo) {
    HighsModelStatus solved_model_status =
        hmos_[solved_hmo].unscaled_model_status_;
    assert(solved_model_status == HighsModelStatus::kInfeasible ||
           solved_model_status == HighsModelStatus::kUnbounded ||
           solved_model_status == HighsModelStatus::kUnboundedOrInfeasible ||
           solved_model_status == HighsModelStatus::kTimeLimit ||
           solved_model_status == HighsModelStatus::kIterationLimit);
    setHighsModelStatusAndInfo(solved_model_status);
  } else {
    setHighsModelStatusBasisSolutionAndInfo();
  }
  double lp_solve_final_time = timer_.readRunHighsClock();
  double this_solve_time = lp_solve_final_time - initial_time;
  if (postsolve_iteration_count < 0) {
    highsLogDev(options_.log_options, HighsLogType::kInfo, "Postsolve  : \n");
  } else {
    highsLogDev(options_.log_options, HighsLogType::kInfo,
                "Postsolve  : %" HIGHSINT_FORMAT "\n",
                postsolve_iteration_count);
  }
  highsLogDev(options_.log_options, HighsLogType::kInfo, "Time       : %8.2f\n",
              this_solve_time);
  highsLogDev(options_.log_options, HighsLogType::kInfo, "Time Pre   : %8.2f\n",
              this_presolve_time);
  highsLogDev(options_.log_options, HighsLogType::kInfo, "Time PreLP : %8.2f\n",
              this_solve_presolved_lp_time);
  highsLogDev(options_.log_options, HighsLogType::kInfo, "Time PostLP: %8.2f\n",
              this_solve_original_lp_time);
  if (this_solve_time > 0) {
    highsLogDev(options_.log_options, HighsLogType::kInfo, "For LP %16s",
                hmos_[original_hmo].lp_.model_name_.c_str());
    double sum_time = 0;
    if (this_presolve_time > 0) {
      sum_time += this_presolve_time;
      HighsInt pct = (100 * this_presolve_time) / this_solve_time;
      highsLogDev(options_.log_options, HighsLogType::kInfo,
                  ": Presolve %8.2f (%3" HIGHSINT_FORMAT "%%)",
                  this_presolve_time, pct);
    }
    if (this_solve_presolved_lp_time > 0) {
      sum_time += this_solve_presolved_lp_time;
      HighsInt pct = (100 * this_solve_presolved_lp_time) / this_solve_time;
      highsLogDev(options_.log_options, HighsLogType::kInfo,
                  ": Solve presolved LP %8.2f (%3" HIGHSINT_FORMAT "%%)",
                  this_solve_presolved_lp_time, pct);
    }
    if (this_postsolve_time > 0) {
      sum_time += this_postsolve_time;
      HighsInt pct = (100 * this_postsolve_time) / this_solve_time;
      highsLogDev(options_.log_options, HighsLogType::kInfo,
                  ": Postsolve %8.2f (%3" HIGHSINT_FORMAT "%%)",
                  this_postsolve_time, pct);
    }
    if (this_solve_original_lp_time > 0) {
      sum_time += this_solve_original_lp_time;
      HighsInt pct = (100 * this_solve_original_lp_time) / this_solve_time;
      highsLogDev(options_.log_options, HighsLogType::kInfo,
                  ": Solve original LP %8.2f (%3" HIGHSINT_FORMAT "%%)",
                  this_solve_original_lp_time, pct);
    }
    highsLogDev(options_.log_options, HighsLogType::kInfo, "\n");
    double rlv_time_difference =
        fabs(sum_time - this_solve_time) / this_solve_time;
    if (rlv_time_difference > 0.1)
      highsLogDev(options_.log_options, HighsLogType::kInfo,
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
  if (!haveHmo("getDualRay")) return HighsStatus::kError;
  return getDualRayInterface(has_dual_ray, dual_ray_value);
}

HighsStatus Highs::getPrimalRay(bool& has_primal_ray,
                                double* primal_ray_value) {
  underDevelopmentLogMessage("getPrimalRay");
  if (!haveHmo("getPrimalRay")) return HighsStatus::kError;
  return getPrimalRayInterface(has_primal_ray, primal_ray_value);
}

HighsStatus Highs::getRanging(HighsRanging& ranging) {
  underDevelopmentLogMessage("getRanging");
  if (!haveHmo("getRanging")) return HighsStatus::kError;
  return getRangingData(ranging, hmos_[0]);
}

HighsStatus Highs::getBasicVariables(HighsInt* basic_variables) {
  if (!haveHmo("getBasicVariables")) return HighsStatus::kError;
  if (basic_variables == NULL) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "getBasicVariables: basic_variables is NULL\n");
    return HighsStatus::kError;
  }
  return getBasicVariablesInterface(basic_variables);
}

HighsStatus Highs::getBasisInverseRow(const HighsInt row, double* row_vector,
                                      HighsInt* row_num_nz,
                                      HighsInt* row_indices) {
  if (!haveHmo("getBasisInverseRow")) return HighsStatus::kError;
  if (row_vector == NULL) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "getBasisInverseRow: row_vector is NULL\n");
    return HighsStatus::kError;
  }
  // row_indices can be NULL - it's the trigger that determines
  // whether they are identified or not
  HighsInt numRow = model_.lp_.numRow_;
  if (row < 0 || row >= numRow) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Row index %" HIGHSINT_FORMAT
                 " out of range [0, %" HIGHSINT_FORMAT
                 "] in getBasisInverseRow\n",
                 row, numRow - 1);
    return HighsStatus::kError;
  }
  bool has_invert = hmos_[0].ekk_instance_.status_.has_invert;
  if (!has_invert) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "No invertible representation for getBasisInverseRow\n");
    return HighsStatus::kError;
  }
  // Compute a row i of the inverse of the basis matrix by solving B^Tx=e_i
  vector<double> rhs;
  rhs.assign(numRow, 0);
  rhs[row] = 1;
  basisSolveInterface(rhs, row_vector, row_num_nz, row_indices, true);
  return HighsStatus::kOk;
}

HighsStatus Highs::getBasisInverseCol(const HighsInt col, double* col_vector,
                                      HighsInt* col_num_nz,
                                      HighsInt* col_indices) {
  if (!haveHmo("getBasisInverseCol")) return HighsStatus::kError;
  if (col_vector == NULL) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "getBasisInverseCol: col_vector is NULL\n");
    return HighsStatus::kError;
  }
  // col_indices can be NULL - it's the trigger that determines
  // whether they are identified or not
  HighsInt numRow = model_.lp_.numRow_;
  if (col < 0 || col >= numRow) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Column index %" HIGHSINT_FORMAT
                 " out of range [0, %" HIGHSINT_FORMAT
                 "] in getBasisInverseCol\n",
                 col, numRow - 1);
    return HighsStatus::kError;
  }
  bool has_invert = hmos_[0].ekk_instance_.status_.has_invert;
  if (!has_invert) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "No invertible representation for getBasisInverseCol\n");
    return HighsStatus::kError;
  }
  // Compute a col i of the inverse of the basis matrix by solving Bx=e_i
  vector<double> rhs;
  rhs.assign(numRow, 0);
  rhs[col] = 1;
  basisSolveInterface(rhs, col_vector, col_num_nz, col_indices, false);
  return HighsStatus::kOk;
}

HighsStatus Highs::getBasisSolve(const double* Xrhs, double* solution_vector,
                                 HighsInt* solution_num_nz,
                                 HighsInt* solution_indices) {
  if (!haveHmo("getBasisSolve")) return HighsStatus::kError;
  if (Xrhs == NULL) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "getBasisSolve: Xrhs is NULL\n");
    return HighsStatus::kError;
  }
  if (solution_vector == NULL) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "getBasisSolve: solution_vector is NULL\n");
    return HighsStatus::kError;
  }
  // solution_indices can be NULL - it's the trigger that determines
  // whether they are identified or not
  bool has_invert = hmos_[0].ekk_instance_.status_.has_invert;
  if (!has_invert) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "No invertible representation for getBasisSolve\n");
    return HighsStatus::kError;
  }
  HighsInt numRow = model_.lp_.numRow_;
  vector<double> rhs;
  rhs.assign(numRow, 0);
  for (HighsInt row = 0; row < numRow; row++) rhs[row] = Xrhs[row];
  basisSolveInterface(rhs, solution_vector, solution_num_nz, solution_indices,
                      false);
  return HighsStatus::kOk;
}

HighsStatus Highs::getBasisTransposeSolve(const double* Xrhs,
                                          double* solution_vector,
                                          HighsInt* solution_num_nz,
                                          HighsInt* solution_indices) {
  if (!haveHmo("getBasisTransposeSolve")) return HighsStatus::kError;
  if (Xrhs == NULL) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "getBasisTransposeSolve: Xrhs is NULL\n");
    return HighsStatus::kError;
  }
  if (solution_vector == NULL) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "getBasisTransposeSolve: solution_vector is NULL\n");
    return HighsStatus::kError;
  }
  // solution_indices can be NULL - it's the trigger that determines
  // whether they are identified or not
  bool has_invert = hmos_[0].ekk_instance_.status_.has_invert;
  if (!has_invert) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "No invertible representation for getBasisTransposeSolve\n");
    return HighsStatus::kError;
  }
  HighsInt numRow = model_.lp_.numRow_;
  vector<double> rhs;
  rhs.assign(numRow, 0);
  for (HighsInt row = 0; row < numRow; row++) rhs[row] = Xrhs[row];
  basisSolveInterface(rhs, solution_vector, solution_num_nz, solution_indices,
                      true);
  return HighsStatus::kOk;
}

HighsStatus Highs::getReducedRow(const HighsInt row, double* row_vector,
                                 HighsInt* row_num_nz, HighsInt* row_indices,
                                 const double* pass_basis_inverse_row_vector) {
  if (!haveHmo("getReducedRow")) return HighsStatus::kError;
  // Ensure that the LP is column-wise
  setOrientation(model_.lp_);
  HighsLp& lp = model_.lp_;
  if (row_vector == NULL) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "getReducedRow: row_vector is NULL\n");
    return HighsStatus::kError;
  }
  // row_indices can be NULL - it's the trigger that determines
  // whether they are identified or not pass_basis_inverse_row_vector
  // NULL - it's the trigger to determine whether it's computed or not
  if (row < 0 || row >= lp.numRow_) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Row index %" HIGHSINT_FORMAT
                 " out of range [0, %" HIGHSINT_FORMAT "] in getReducedRow\n",
                 row, lp.numRow_ - 1);
    return HighsStatus::kError;
  }
  bool has_invert = hmos_[0].ekk_instance_.status_.has_invert;
  if (!has_invert) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "No invertible representation for getReducedRow\n");
    return HighsStatus::kError;
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
    if (fabs(value) > kHighsTiny) {
      if (return_indices) row_indices[(*row_num_nz)++] = col;
      row_vector[col] = value;
    }
  }
  return HighsStatus::kOk;
}

HighsStatus Highs::getReducedColumn(const HighsInt col, double* col_vector,
                                    HighsInt* col_num_nz,
                                    HighsInt* col_indices) {
  if (!haveHmo("getReducedColumn")) return HighsStatus::kError;
  // Ensure that the LP is column-wise
  setOrientation(model_.lp_);
  HighsLp& lp = model_.lp_;
  if (col_vector == NULL) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "getReducedColumn: col_vector is NULL\n");
    return HighsStatus::kError;
  }
  // col_indices can be NULL - it's the trigger that determines
  // whether they are identified or not
  if (col < 0 || col >= lp.numCol_) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Column index %" HIGHSINT_FORMAT
                 " out of range [0, %" HIGHSINT_FORMAT
                 "] in getReducedColumn\n",
                 col, lp.numCol_ - 1);
    return HighsStatus::kError;
  }
  bool has_invert = hmos_[0].ekk_instance_.status_.has_invert;
  if (!has_invert) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "No invertible representation for getReducedColumn\n");
    return HighsStatus::kError;
  }
  HighsInt numRow = lp.numRow_;
  vector<double> rhs;
  rhs.assign(numRow, 0);
  for (HighsInt el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++)
    rhs[lp.Aindex_[el]] = lp.Avalue_[el];
  basisSolveInterface(rhs, col_vector, col_num_nz, col_indices, false);
  return HighsStatus::kOk;
}

HighsStatus Highs::setSolution(const HighsSolution& solution) {
  HighsStatus return_status = HighsStatus::kOk;
  // Check if primal solution is valid.
  if (model_.lp_.numCol_ > 0 &&
      solution.col_value.size() >= model_.lp_.numCol_) {
    // Worth considering the column values
    solution_.col_value = solution.col_value;
    if (model_.lp_.numRow_ > 0) {
      // Worth computing the row values
      solution_.row_value.resize(model_.lp_.numRow_);
      return_status =
          interpretCallStatus(calculateRowValues(model_.lp_, solution_),
                              return_status, "calculateRowValues");
      if (return_status == HighsStatus::kError) return return_status;
    }
    solution_.value_valid = true;
  } else {
    // Primal solution not valid
    solution_.value_valid = false;
  }
  // Check if dual solution is valid.
  if (model_.lp_.numRow_ > 0 &&
      solution.row_dual.size() >= model_.lp_.numRow_) {
    // Worth considering the row duals
    solution_.row_dual = solution.row_dual;
    if (model_.lp_.numCol_ > 0) {
      // Worth computing the column duals
      solution_.col_dual.resize(model_.lp_.numCol_);
      return_status =
          interpretCallStatus(calculateColDuals(model_.lp_, solution_),
                              return_status, "calculateColDuals");
      if (return_status == HighsStatus::kError) return return_status;
    }
    solution_.dual_valid = true;
  } else {
    // Dual solution not valid
    solution_.dual_valid = false;
  }
  return returnFromHighs(return_status);
}

HighsStatus Highs::setBasis(const HighsBasis& basis) {
  // Check the user-supplied basis
  if (!isBasisConsistent(model_.lp_, basis)) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "setBasis: invalid basis\n");
    return HighsStatus::kError;
  }
  // Update the HiGHS basis
  basis_ = basis;
  basis_.valid = true;
  // Follow implications of a new HiGHS basis
  newHighsBasis();
  // Can't use returnFromHighs since...
  return HighsStatus::kOk;
}

HighsStatus Highs::setBasis() {
  // Invalidate the basis for HiGHS Don't set to logical basis since
  // that causes presolve to be skipped
  basis_.valid = false;
  // Follow implications of a new HiGHS basis
  newHighsBasis();
  // Can't use returnFromHighs since...
  return HighsStatus::kOk;
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
  HighsStatus return_status = HighsStatus::kOk;
  clearPresolve();
  // Check that there is a HighsModelObject
  if (!haveHmo("addRows")) return false;
  return_status = interpretCallStatus(
      addRowsInterface(num_new_row, lower_bounds, upper_bounds, num_new_nz,
                       starts, indices, values),
      return_status, "addRows");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
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
  HighsStatus return_status = HighsStatus::kOk;
  clearPresolve();
  if (!haveHmo("addCols")) return false;
  return_status = interpretCallStatus(
      addColsInterface(num_new_col, costs, lower_bounds, upper_bounds,
                       num_new_nz, starts, indices, values),
      return_status, "addCols");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

bool Highs::changeObjectiveSense(const ObjSense sense) {
  HighsStatus return_status = HighsStatus::kOk;
  clearPresolve();
  if (!haveHmo("changeObjectiveSense")) return false;
  HighsStatus call_status;
  call_status = changeObjectiveSenseInterface(sense);
  return_status =
      interpretCallStatus(call_status, return_status, "changeObjectiveSense");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

bool Highs::changeColIntegrality(const HighsInt col,
                                 const HighsVarType integrality) {
  return changeColsIntegrality(1, &col, &integrality);
}

bool Highs::changeColsIntegrality(const HighsInt from_col,
                                  const HighsInt to_col,
                                  const HighsVarType* integrality) {
  HighsStatus return_status = HighsStatus::kOk;
  clearPresolve();
  HighsStatus call_status;
  HighsIndexCollection index_collection;
  index_collection.dimension_ = model_.lp_.numCol_;
  index_collection.is_interval_ = true;
  index_collection.from_ = from_col;
  index_collection.to_ = to_col;
  if (!haveHmo("changeColsIntegrality")) return false;
  call_status = changeIntegralityInterface(index_collection, integrality);
  return_status =
      interpretCallStatus(call_status, return_status, "changeIntegrality");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

bool Highs::changeColsIntegrality(const HighsInt num_set_entries,
                                  const HighsInt* set,
                                  const HighsVarType* integrality) {
  if (num_set_entries <= 0) return true;
  HighsStatus return_status = HighsStatus::kOk;
  clearPresolve();
  HighsStatus call_status;
  // Create a local set that is not const since index_collection.set_
  // cannot be const as it may change if the set is not ordered
  vector<HighsInt> local_set{set, set + num_set_entries};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = model_.lp_.numCol_;
  index_collection.is_set_ = true;
  index_collection.set_ = &local_set[0];
  index_collection.set_num_entries_ = num_set_entries;
  if (!haveHmo("changeColsIntegrality")) return false;
  call_status = changeIntegralityInterface(index_collection, integrality);
  return_status =
      interpretCallStatus(call_status, return_status, "changeIntegrality");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

bool Highs::changeColsIntegrality(const HighsInt* mask,
                                  const HighsVarType* integrality) {
  HighsStatus return_status = HighsStatus::kOk;
  clearPresolve();
  HighsStatus call_status;
  // Create a local mask that is not const since
  // index_collection.mask_ cannot be const as it changes when
  // deleting rows/columns
  vector<HighsInt> local_mask{mask, mask + model_.lp_.numCol_};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = model_.lp_.numCol_;
  index_collection.is_mask_ = true;
  index_collection.mask_ = &local_mask[0];
  if (!haveHmo("changeColsIntegrality")) return false;
  call_status = changeIntegralityInterface(index_collection, integrality);
  return_status =
      interpretCallStatus(call_status, return_status, "changeIntegrality");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

bool Highs::changeColCost(const HighsInt col, const double cost) {
  return changeColsCost(1, &col, &cost);
}

bool Highs::changeColsCost(const HighsInt from_col, const HighsInt to_col,
                           const double* cost) {
  HighsStatus return_status = HighsStatus::kOk;
  clearPresolve();
  HighsStatus call_status;
  HighsIndexCollection index_collection;
  index_collection.dimension_ = model_.lp_.numCol_;
  index_collection.is_interval_ = true;
  index_collection.from_ = from_col;
  index_collection.to_ = to_col;
  if (!haveHmo("changeColsCost")) return false;
  call_status = changeCostsInterface(index_collection, cost);
  return_status =
      interpretCallStatus(call_status, return_status, "changeCosts");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

bool Highs::changeColsCost(const HighsInt num_set_entries, const HighsInt* set,
                           const double* cost) {
  if (num_set_entries <= 0) return true;
  HighsStatus return_status = HighsStatus::kOk;
  clearPresolve();
  HighsStatus call_status;
  // Create a local set that is not const since index_collection.set_
  // cannot be const as it may change if the set is not ordered
  vector<HighsInt> local_set{set, set + num_set_entries};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = model_.lp_.numCol_;
  index_collection.is_set_ = true;
  index_collection.set_ = &local_set[0];
  index_collection.set_num_entries_ = num_set_entries;
  if (!haveHmo("changeColsCost")) return false;
  call_status = changeCostsInterface(index_collection, cost);
  return_status =
      interpretCallStatus(call_status, return_status, "changeCosts");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

bool Highs::changeColsCost(const HighsInt* mask, const double* cost) {
  HighsStatus return_status = HighsStatus::kOk;
  clearPresolve();
  HighsStatus call_status;
  // Create a local mask that is not const since
  // index_collection.mask_ cannot be const as it changes when
  // deleting rows/columns
  vector<HighsInt> local_mask{mask, mask + model_.lp_.numCol_};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = model_.lp_.numCol_;
  index_collection.is_mask_ = true;
  index_collection.mask_ = &local_mask[0];
  if (!haveHmo("changeColsCost")) return false;
  call_status = changeCostsInterface(index_collection, cost);
  return_status =
      interpretCallStatus(call_status, return_status, "changeCosts");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

bool Highs::changeColBounds(const HighsInt col, const double lower,
                            const double upper) {
  return changeColsBounds(1, &col, &lower, &upper);
}

bool Highs::changeColsBounds(const HighsInt from_col, const HighsInt to_col,
                             const double* lower, const double* upper) {
  HighsStatus return_status = HighsStatus::kOk;
  clearPresolve();
  HighsStatus call_status;
  HighsIndexCollection index_collection;
  index_collection.dimension_ = model_.lp_.numCol_;
  index_collection.is_interval_ = true;
  index_collection.from_ = from_col;
  index_collection.to_ = to_col;
  if (!haveHmo("changeColsBounds")) return false;
  call_status = changeColBoundsInterface(index_collection, lower, upper);
  return_status =
      interpretCallStatus(call_status, return_status, "changeColBounds");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

bool Highs::changeColsBounds(const HighsInt num_set_entries,
                             const HighsInt* set, const double* lower,
                             const double* upper) {
  if (num_set_entries <= 0) return true;
  HighsStatus return_status = HighsStatus::kOk;
  clearPresolve();
  HighsStatus call_status;
  // Create a local set that is not const since index_collection.set_
  // cannot be const as it may change if the set is not ordered
  vector<HighsInt> local_set{set, set + num_set_entries};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = model_.lp_.numCol_;
  index_collection.is_set_ = true;
  index_collection.set_ = &local_set[0];
  index_collection.set_num_entries_ = num_set_entries;
  if (!haveHmo("changeColsBounds")) return false;
  call_status = changeColBoundsInterface(index_collection, lower, upper);
  return_status =
      interpretCallStatus(call_status, return_status, "changeColBounds");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

bool Highs::changeColsBounds(const HighsInt* mask, const double* lower,
                             const double* upper) {
  HighsStatus return_status = HighsStatus::kOk;
  clearPresolve();
  HighsStatus call_status;
  // Create a local mask that is not const since
  // index_collection.mask_ cannot be const as it changes when
  // deleting rows/columns
  vector<HighsInt> local_mask{mask, mask + model_.lp_.numCol_};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = model_.lp_.numCol_;
  index_collection.is_mask_ = true;
  index_collection.mask_ = &local_mask[0];
  if (!haveHmo("changeColsBounds")) return false;
  call_status = changeColBoundsInterface(index_collection, lower, upper);
  return_status =
      interpretCallStatus(call_status, return_status, "changeColBounds");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

bool Highs::changeRowBounds(const HighsInt row, const double lower,
                            const double upper) {
  return changeRowsBounds(1, &row, &lower, &upper);
}

bool Highs::changeRowsBounds(const HighsInt from_row, const HighsInt to_row,
                             const double* lower, const double* upper) {
  HighsStatus return_status = HighsStatus::kOk;
  clearPresolve();
  HighsStatus call_status;
  HighsIndexCollection index_collection;
  index_collection.dimension_ = model_.lp_.numRow_;
  index_collection.is_interval_ = true;
  index_collection.from_ = from_row;
  index_collection.to_ = to_row;
  if (!haveHmo("changeRowsBounds")) return false;
  call_status = changeRowBoundsInterface(index_collection, lower, upper);
  return_status =
      interpretCallStatus(call_status, return_status, "changeRowBounds");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

bool Highs::changeRowsBounds(const HighsInt num_set_entries,
                             const HighsInt* set, const double* lower,
                             const double* upper) {
  if (num_set_entries <= 0) return true;
  HighsStatus return_status = HighsStatus::kOk;
  clearPresolve();
  HighsStatus call_status;
  // Create a local set that is not const since index_collection.set_
  // cannot be const as it may change if the set is not ordered
  vector<HighsInt> local_set{set, set + num_set_entries};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = model_.lp_.numRow_;
  index_collection.is_set_ = true;
  index_collection.set_ = &local_set[0];
  index_collection.set_num_entries_ = num_set_entries;
  if (!haveHmo("changeRowsBounds")) return false;
  call_status = changeRowBoundsInterface(index_collection, lower, upper);
  return_status =
      interpretCallStatus(call_status, return_status, "changeRowBounds");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

bool Highs::changeRowsBounds(const HighsInt* mask, const double* lower,
                             const double* upper) {
  HighsStatus return_status = HighsStatus::kOk;
  clearPresolve();
  HighsStatus call_status;
  // Create a local mask that is not const since
  // index_collection.mask_ cannot be const as it changes when
  // deleting rows/columns
  vector<HighsInt> local_mask{mask, mask + model_.lp_.numRow_};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = model_.lp_.numRow_;
  index_collection.is_mask_ = true;
  index_collection.mask_ = &local_mask[0];
  if (!haveHmo("changeRowsBounds")) return false;
  call_status = changeRowBoundsInterface(index_collection, lower, upper);
  return_status =
      interpretCallStatus(call_status, return_status, "changeRowBounds");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

bool Highs::changeCoeff(const HighsInt row, const HighsInt col,
                        const double value) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  if (!haveHmo("changeCoeff")) return false;
  call_status = changeCoefficientInterface(row, col, value);
  return_status =
      interpretCallStatus(call_status, return_status, "changeCoefficient");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

bool Highs::getObjectiveSense(ObjSense& sense) {
  if (!haveHmo("getObjectiveSense")) return false;
  sense = model_.lp_.sense_;
  return true;
}

bool Highs::getCols(const HighsInt from_col, const HighsInt to_col,
                    HighsInt& num_col, double* costs, double* lower,
                    double* upper, HighsInt& num_nz, HighsInt* start,
                    HighsInt* index, double* value) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  HighsIndexCollection index_collection;
  index_collection.dimension_ = model_.lp_.numCol_;
  index_collection.is_interval_ = true;
  index_collection.from_ = from_col;
  index_collection.to_ = to_col;
  if (!haveHmo("getCols")) return false;
  call_status = getColsInterface(index_collection, num_col, costs, lower, upper,
                                 num_nz, start, index, value);
  return_status = interpretCallStatus(call_status, return_status, "getCols");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

bool Highs::getCols(const HighsInt num_set_entries, const HighsInt* set,
                    HighsInt& num_col, double* costs, double* lower,
                    double* upper, HighsInt& num_nz, HighsInt* start,
                    HighsInt* index, double* value) {
  if (num_set_entries <= 0) return true;
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  // Create a local set that is not const since index_collection.set_
  // cannot be const as it may change if the set is not ordered
  vector<HighsInt> local_set{set, set + num_set_entries};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = model_.lp_.numCol_;
  index_collection.is_set_ = true;
  index_collection.set_ = &local_set[0];
  index_collection.set_num_entries_ = num_set_entries;
  if (!haveHmo("getCols")) return false;
  call_status = getColsInterface(index_collection, num_col, costs, lower, upper,
                                 num_nz, start, index, value);
  return_status = interpretCallStatus(call_status, return_status, "getCols");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

bool Highs::getCols(const HighsInt* mask, HighsInt& num_col, double* costs,
                    double* lower, double* upper, HighsInt& num_nz,
                    HighsInt* start, HighsInt* index, double* value) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  // Create a local mask that is not const since
  // index_collection.mask_ cannot be const as it changes when
  // deleting rows/columns
  vector<HighsInt> local_mask{mask, mask + model_.lp_.numCol_};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = model_.lp_.numCol_;
  index_collection.is_mask_ = true;
  index_collection.mask_ = &local_mask[0];
  if (!haveHmo("getCols")) return false;
  call_status = getColsInterface(index_collection, num_col, costs, lower, upper,
                                 num_nz, start, index, value);
  return_status = interpretCallStatus(call_status, return_status, "getCols");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

bool Highs::getRows(const HighsInt from_row, const HighsInt to_row,
                    HighsInt& num_row, double* lower, double* upper,
                    HighsInt& num_nz, HighsInt* start, HighsInt* index,
                    double* value) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  HighsIndexCollection index_collection;
  index_collection.dimension_ = model_.lp_.numRow_;
  index_collection.is_interval_ = true;
  index_collection.from_ = from_row;
  index_collection.to_ = to_row;
  if (!haveHmo("getRows")) return false;
  call_status = getRowsInterface(index_collection, num_row, lower, upper,
                                 num_nz, start, index, value);
  return_status = interpretCallStatus(call_status, return_status, "getRows");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

bool Highs::getRows(const HighsInt num_set_entries, const HighsInt* set,
                    HighsInt& num_row, double* lower, double* upper,
                    HighsInt& num_nz, HighsInt* start, HighsInt* index,
                    double* value) {
  if (num_set_entries <= 0) return true;
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  // Create a local set that is not const since index_collection.set_
  // cannot be const as it may change if the set is not ordered
  vector<HighsInt> local_set{set, set + num_set_entries};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = model_.lp_.numRow_;
  index_collection.is_set_ = true;
  index_collection.set_ = &local_set[0];
  index_collection.set_num_entries_ = num_set_entries;
  if (!haveHmo("getRows")) return false;
  call_status = getRowsInterface(index_collection, num_row, lower, upper,
                                 num_nz, start, index, value);
  return_status = interpretCallStatus(call_status, return_status, "getRows");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

bool Highs::getRows(const HighsInt* mask, HighsInt& num_row, double* lower,
                    double* upper, HighsInt& num_nz, HighsInt* start,
                    HighsInt* index, double* value) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  // Create a local mask that is not const since
  // index_collection.mask_ cannot be const as it changes when
  // deleting rows/columns
  vector<HighsInt> local_mask{mask, mask + model_.lp_.numRow_};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = model_.lp_.numRow_;
  index_collection.is_mask_ = true;
  index_collection.mask_ = &local_mask[0];
  if (!haveHmo("getRows")) return false;
  call_status = getRowsInterface(index_collection, num_row, lower, upper,
                                 num_nz, start, index, value);
  return_status = interpretCallStatus(call_status, return_status, "getRows");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

bool Highs::getCoeff(const HighsInt row, const HighsInt col, double& value) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  if (!haveHmo("getCoeff")) return false;
  call_status = getCoefficientInterface(row, col, value);
  return_status =
      interpretCallStatus(call_status, return_status, "getCoefficient");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

bool Highs::deleteCols(const HighsInt from_col, const HighsInt to_col) {
  HighsStatus return_status = HighsStatus::kOk;
  clearPresolve();
  HighsStatus call_status;
  HighsIndexCollection index_collection;
  index_collection.dimension_ = model_.lp_.numCol_;
  index_collection.is_interval_ = true;
  index_collection.from_ = from_col;
  index_collection.to_ = to_col;
  if (!haveHmo("deleteCols")) return false;
  call_status = deleteColsInterface(index_collection);
  return_status = interpretCallStatus(call_status, return_status, "deleteCols");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

bool Highs::deleteCols(const HighsInt num_set_entries, const HighsInt* set) {
  if (num_set_entries <= 0) return true;
  HighsStatus return_status = HighsStatus::kOk;
  clearPresolve();
  HighsStatus call_status;
  // Create a local set that is not const since index_collection.set_
  // cannot be const as it may change if the set is not ordered
  vector<HighsInt> local_set{set, set + num_set_entries};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = model_.lp_.numCol_;
  index_collection.is_set_ = true;
  index_collection.set_ = &local_set[0];
  index_collection.set_num_entries_ = num_set_entries;
  if (!haveHmo("deleteCols")) return false;
  call_status = deleteColsInterface(index_collection);
  return_status = interpretCallStatus(call_status, return_status, "deleteCols");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

bool Highs::deleteCols(HighsInt* mask) {
  HighsStatus return_status = HighsStatus::kOk;
  clearPresolve();
  HighsStatus call_status;
  HighsIndexCollection index_collection;
  index_collection.dimension_ = model_.lp_.numCol_;
  index_collection.is_mask_ = true;
  index_collection.mask_ = &mask[0];
  if (!haveHmo("deleteCols")) return false;
  call_status = deleteColsInterface(index_collection);
  return_status = interpretCallStatus(call_status, return_status, "deleteCols");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

bool Highs::deleteRows(const HighsInt from_row, const HighsInt to_row) {
  HighsStatus return_status = HighsStatus::kOk;
  clearPresolve();
  HighsStatus call_status;
  HighsIndexCollection index_collection;
  index_collection.dimension_ = model_.lp_.numRow_;
  index_collection.is_interval_ = true;
  index_collection.from_ = from_row;
  index_collection.to_ = to_row;
  if (!haveHmo("deleteRows")) return false;
  call_status = deleteRowsInterface(index_collection);
  return_status = interpretCallStatus(call_status, return_status, "deleteRows");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

bool Highs::deleteRows(const HighsInt num_set_entries, const HighsInt* set) {
  if (num_set_entries <= 0) return true;
  HighsStatus return_status = HighsStatus::kOk;
  clearPresolve();
  HighsStatus call_status;
  // Create a local set that is not const since index_collection.set_
  // cannot be const as it may change if the set is not ordered
  vector<HighsInt> local_set{set, set + num_set_entries};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = model_.lp_.numRow_;
  index_collection.is_set_ = true;
  index_collection.set_ = &local_set[0];
  index_collection.set_num_entries_ = num_set_entries;
  if (!haveHmo("deleteRows")) return false;
  call_status = deleteRowsInterface(index_collection);
  return_status = interpretCallStatus(call_status, return_status, "deleteRows");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

bool Highs::deleteRows(HighsInt* mask) {
  HighsStatus return_status = HighsStatus::kOk;
  clearPresolve();
  HighsStatus call_status;
  HighsIndexCollection index_collection;
  index_collection.dimension_ = model_.lp_.numRow_;
  index_collection.is_mask_ = true;
  index_collection.mask_ = &mask[0];
  if (!haveHmo("deleteRows")) return false;
  call_status = deleteRowsInterface(index_collection);
  return_status = interpretCallStatus(call_status, return_status, "deleteRows");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

bool Highs::scaleCol(const HighsInt col, const double scaleval) {
  HighsStatus return_status = HighsStatus::kOk;
  clearPresolve();
  HighsStatus call_status;
  if (!haveHmo("scaleCol")) return false;
  call_status = scaleColInterface(col, scaleval);
  return_status = interpretCallStatus(call_status, return_status, "scaleCol");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

bool Highs::scaleRow(const HighsInt row, const double scaleval) {
  HighsStatus return_status = HighsStatus::kOk;
  clearPresolve();
  HighsStatus call_status;
  if (!haveHmo("scaleRow")) return false;
  call_status = scaleRowInterface(row, scaleval);
  return_status = interpretCallStatus(call_status, return_status, "scaleRow");
  if (return_status == HighsStatus::kError) return false;
  return returnFromHighs(return_status) != HighsStatus::kError;
}

double Highs::getInfinity() { return kHighsInf; }

double Highs::getRunTime() { return timer_.readRunHighsClock(); }

void Highs::deprecationMessage(const std::string method_name,
                               const std::string alt_method_name) const {
  if (alt_method_name.compare("None") == 0) {
    highsLogUser(options_.log_options, HighsLogType::kWarning,
                 "Method %s is deprecated: no alternative method\n",
                 method_name.c_str());
  } else {
    highsLogUser(options_.log_options, HighsLogType::kWarning,
                 "Method %s is deprecated: alternative method is %s\n",
                 method_name.c_str(), alt_method_name.c_str());
  }
}

#ifdef HiGHSDEV
void Highs::reportModelStatusSolutionBasis(const std::string message,
                                           const HighsInt hmo_ix) {
  HighsModelStatus& model_status = model_status_;
  HighsModelStatus& scaled_model_status = scaled_model_status_;
  HighsSolution& solution = solution_;
  HighsBasis& basis = basis_;
  HighsInt unscaled_primal_solution_status = info_.primal_solution_status;
  HighsInt unscaled_dual_solution_status = info_.dual_solution_status;
  HighsLp& lp = model_.lp_;
  if (hmo_ix >= 0) {
    assert(hmo_ix < (HighsInt)hmos_.size());
    model_status = hmos_[hmo_ix].unscaled_model_status_;
    scaled_model_status = hmos_[hmo_ix].scaled_model_status_;
    solution = hmos_[hmo_ix].solution_;
    basis = hmos_[hmo_ix].basis_;
    unscaled_primal_solution_status =
        hmos_[hmo_ix].solution_params_.primal_solution_status;
    unscaled_dual_solution_status =
        hmos_[hmo_ix].solution_params_.dual_solution_status;
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
      unscaled_primal_solution_status, (HighsInt)solution.col_value.size(),
      (HighsInt)solution.row_value.size(), unscaled_dual_solution_status,
      (HighsInt)solution.col_dual.size(), (HighsInt)solution.row_dual.size(),
      basis.valid, (HighsInt)basis.col_status.size(),
      (HighsInt)basis.row_status.size());
}
#endif

std::string Highs::modelStatusToString(
    const HighsModelStatus model_status) const {
  return utilModelStatusToString(model_status);
}

std::string Highs::solutionStatusToString(const HighsInt solution_status) {
  return utilSolutionStatusToString(solution_status);
}

void Highs::setMatrixOrientation(const MatrixOrientation& desired_orientation) {
  setOrientation(model_.lp_, desired_orientation);
}

// Private methods
HighsPresolveStatus Highs::runPresolve() {
  presolve_.clear();
  // Exit if the problem is empty or if presolve is set to off.
  if (options_.presolve == kHighsOffString)
    return HighsPresolveStatus::kNotPresolved;

  // @FlipRowDual Side-stpe presolve until @leona has fixed it wrt row dual flip
  const bool force_no_presolve = false;
  if (force_no_presolve) {
    printf("Forcing no presolve!!\n");
    return HighsPresolveStatus::kNotPresolved;
  }

  // Ensure that the LP is column-wise
  // setOrientation(model_.lp_);

  if (model_.lp_.numCol_ == 0 && model_.lp_.numRow_ == 0)
    return HighsPresolveStatus::kNullError;

  // Clear info from previous runs if model_.lp_ has been modified.
  double start_presolve = timer_.readRunHighsClock();

  // Set time limit.
  if (options_.time_limit > 0 && options_.time_limit < kHighsInf) {
    double left = options_.time_limit - start_presolve;
    if (left <= 0) {
      highsLogDev(options_.log_options, HighsLogType::kError,
                  "Time limit reached while reading in matrix\n");
      return HighsPresolveStatus::kTimeout;
    }

    highsLogDev(options_.log_options, HighsLogType::kVerbose,
                "Time limit set: reading matrix took %.2g, presolve "
                "time left: %.2g\n",
                start_presolve, left);
  }

  // Presolve.
  presolve_.init(model_.lp_, timer_);
  presolve_.options_ = &options_;
  if (options_.time_limit > 0 && options_.time_limit < kHighsInf) {
    double current = timer_.readRunHighsClock();
    double time_init = current - start_presolve;
    double left = presolve_.options_->time_limit - time_init;
    if (left <= 0) {
      highsLogDev(options_.log_options, HighsLogType::kError,
                  "Time limit reached while copying matrix into presolve.\n");
      return HighsPresolveStatus::kTimeout;
    }
    highsLogDev(options_.log_options, HighsLogType::kVerbose,
                "Time limit set: copying matrix took %.2g, presolve "
                "time left: %.2g\n",
                time_init, left);
  }

  HighsPresolveStatus presolve_return_status = presolve_.run();

  highsLogDev(options_.log_options, HighsLogType::kVerbose,
              "presolve_.run() returns status: %s\n",
              presolve_.presolveStatusToString(presolve_return_status).c_str());

  // Update reduction counts.
  switch (presolve_.presolve_status_) {
    case HighsPresolveStatus::kReduced: {
      HighsLp& reduced_lp = presolve_.getReducedProblem();
      presolve_.info_.n_cols_removed = model_.lp_.numCol_ - reduced_lp.numCol_;
      presolve_.info_.n_rows_removed = model_.lp_.numRow_ - reduced_lp.numRow_;
      presolve_.info_.n_nnz_removed = (HighsInt)model_.lp_.Avalue_.size() -
                                      (HighsInt)reduced_lp.Avalue_.size();
      break;
    }
    case HighsPresolveStatus::kReducedToEmpty: {
      presolve_.info_.n_cols_removed = model_.lp_.numCol_;
      presolve_.info_.n_rows_removed = model_.lp_.numRow_;
      presolve_.info_.n_nnz_removed = (HighsInt)model_.lp_.Avalue_.size();
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
  if (!solution_ok) return HighsPostsolveStatus::kReducedSolutionDimenionsError;

  presolve_.data_.postSolveStack.undo(options_,
                                      presolve_.data_.recovered_solution_,
                                      presolve_.data_.recovered_basis_);

  if (model_.lp_.sense_ == ObjSense::kMaximize)
    presolve_.negateReducedLpColDuals(true);

  return HighsPostsolveStatus::kSolutionRecovered;
}

void Highs::clearPresolve() {
  model_presolve_status_ = HighsPresolveStatus::kNotPresolved;
  presolve_.clear();
}

void Highs::clearUserSolverData() {
  clearModelStatus();
  clearSolution();
  clearBasis();
  clearInfo();
}

void Highs::clearModelStatus() {
  model_status_ = HighsModelStatus::kNotset;
  scaled_model_status_ = HighsModelStatus::kNotset;
}

void Highs::clearSolution() {
  info_.primal_solution_status = kSolutionStatusNone;
  info_.dual_solution_status = kSolutionStatusNone;
  clearSolutionUtil(solution_);
}

void Highs::clearBasis() { clearBasisUtil(basis_); }

void Highs::clearInfo() { info_.clear(); }

// The method below runs calls solveLp to solve the LP associated with
// a particular model, integrating the iteration counts into the
// overall values in HighsInfo
HighsStatus Highs::callSolveLp(const HighsInt model_index,
                               const string message) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;

  // Check that the model index is OK
  bool model_index_ok =
      model_index >= 0 && model_index < (HighsInt)hmos_.size();
  assert(model_index_ok);
  if (!model_index_ok) return HighsStatus::kError;

  HighsModelObject& model = hmos_[model_index];
  // Check that the model isn't row-wise
  assert(model_.lp_.orientation_ != MatrixOrientation::kRowwise);

  // Copy the LP solver iteration counts to this model so that they
  // are updated
  hmos_[model_index].iteration_counts_ = iteration_counts_;

  // Solve the LP
  call_status = solveLp(model, message);
  return_status = interpretCallStatus(call_status, return_status, "solveLp");
  if (return_status == HighsStatus::kError) return return_status;

  // Copy this model's iteration counts to the LP solver iteration counts so
  // that they are updated
  iteration_counts_ = hmos_[model_index].iteration_counts_;
  return return_status;
}

HighsStatus Highs::callSolveQp() {
  HighsStatus return_status = HighsStatus::kOk;
  // Check that the model isn't row-wise - not yet in master
  HighsLp& lp = model_.lp_;
  HighsHessian& hessian = model_.hessian_;
  assert(lp.orientation_ != MatrixOrientation::kRowwise);
  //
  // Run the QP solver
  return_status = HighsStatus::kError;
  /*
  Instance instance(lp.numCol_, lp.numRow_);


  instance.num_con = lp.numRow_;
  instance.num_var = lp.numCol_;

  instance.A.mat.num_col = lp.numCol_;
  instance.A.mat.num_row = lp.numRow_;
  instance.A.mat.start = lp.Astart_;
  instance.A.mat.index = lp.Aindex_;
  instance.A.mat.value = lp.Avalue_;
  instance.c.value = lp.colCost_;
  instance.con_lo = lp.rowLower_;
  instance.con_up = lp.rowUpper_;
  instance.var_lo = lp.colLower_;
  instance.var_up = lp.colUpper_;
  instance.Q.mat.num_col = lp.numCol_;
  instance.Q.mat.num_row = lp.numCol_;
  instance.Q.mat.start = hessian.q_start_;
  instance.Q.mat.index = hessian.q_index_;
  instance.Q.mat.value = hessian.q_value_;

  if (lp.sense_ != ObjSense::kMinimize) {
    for (double& i : instance.c.value) {
      i *= -1.0;
    }
  }

  Runtime runtime(instance);

  runtime.settings.reportingfequency = 1000;
  runtime.endofiterationevent.subscribe(reportIteration);
  runtime.settings.iterationlimit = std::numeric_limits<int>::max();
  runtime.settings.ratiotest = new RatiotestTwopass(instance, 0.000000001,
  0.000001); Solver solver(runtime); solver.solve();





  //
  // Cheating now, but need to set this honestly!
  HighsStatus call_status = HighsStatus::kOk;
  return_status =
      interpretCallStatus(call_status, return_status, "QpSolver");
  if (return_status == HighsStatus::kError) return return_status;
  // Cheating now, but need to set this honestly!
  scaled_model_status_ = runtime.status == ProblemStatus::OPTIMAL ?
  HighsModelStatus::kOptimal : runtime.status == ProblemStatus::UNBOUNDED ?
  HighsModelStatus::kUnbounded : HighsModelStatus::kInfeasible; model_status_ =
  scaled_model_status_;
  // Set the values in HighsInfo instance info_
  info_.qp_iteration_count = runtime.statistics.num_iterations;
  info_.simplex_iteration_count = runtime.statistics.phase1_iterations;
  info_.ipm_iteration_count = -1;
  info_.crossover_iteration_count = -1;
  info_.primal_status = runtime.status == ProblemStatus::OPTIMAL ?
  SolutionStatus::kSolutionStatusFeasible :
  SolutionStatus::kSolutionStatusInfeasible; info_.dual_status = runtime.status
  == ProblemStatus::OPTIMAL ? SolutionStatus::kSolutionStatusFeasible :
  SolutionStatus::kSolutionStatusInfeasible; info_.objective_function_value =
  runtime.instance.objval(runtime.primal); info_.num_primal_infeasibilities =
  -1;  // Not known
  // Are the violations max or sum?
  info_.max_primal_infeasibility =0.0; //
  info_.sum_primal_infeasibilities = -1;  // Not known
  info_.num_dual_infeasibilities = -1;    // Not known
  info_.max_dual_infeasibility = -1;      // Not known
  info_.sum_dual_infeasibilities = -1;    // Not known
  // The solution needs to be here, but just resize it for now

  info_.primal_status = SolutionStatus::kSolutionStatusFeasible;
  solution_.col_value.resize(lp.numCol_);
  solution_.col_dual.resize(lp.numCol_);
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    solution_.col_value[iCol] = runtime.primal.value[iCol]; //
    solution_.col_dual[iCol] = runtime.dualvar.value[iCol];
  }

  solution_.row_value.resize(lp.numRow_);
  solution_.row_dual.resize(lp.numRow_);
  for (int iRow = 0; iRow < lp.numRow_; iRow++) {
    solution_.row_value[iRow] = runtime.rowactivity.value[iRow];
    solution_.row_dual[iRow] = runtime.dualcon.value[iRow];
  }
  */
  return return_status;
}

HighsStatus Highs::callSolveMip() {
  HighsStatus return_status = HighsStatus::kOk;
  // Ensure that any solver data for users in Highs class members are
  // cleared
  clearUserSolverData();
  // Run the MIP solver
  HighsInt log_dev_level = options_.log_dev_level;
  //  options_.log_dev_level = kHighsLogDevLevelInfo;
  // Check that the model isn't row-wise
  assert(model_.lp_.orientation_ != MatrixOrientation::kRowwise);
  HighsMipSolver solver(options_, model_.lp_, solution_);
  solver.run();
  options_.log_dev_level = log_dev_level;
  HighsStatus call_status = HighsStatus::kOk;
  return_status =
      interpretCallStatus(call_status, return_status, "HighsMipSolver::solver");
  if (return_status == HighsStatus::kError) return return_status;
  scaled_model_status_ = solver.modelstatus_;
  model_status_ = scaled_model_status_;
  // Use generic method to set data required for info
  HighsSolutionParams solution_params;
  solution_params.primal_feasibility_tolerance =
      options_.primal_feasibility_tolerance;
  solution_params.dual_feasibility_tolerance =
      options_.dual_feasibility_tolerance;
  if (solver.solution_objective_ != kHighsInf) {
    // There is a primal solution
    HighsInt solver_solution_size = solver.solution_.size();
    assert(solver_solution_size >= model_.lp_.numCol_);
    solution_.col_value.resize(model_.lp_.numCol_);
    solution_.row_value.assign(model_.lp_.numRow_, 0);
    for (HighsInt iCol = 0; iCol < model_.lp_.numCol_; iCol++) {
      double value = solver.solution_[iCol];
      for (HighsInt iEl = model_.lp_.Astart_[iCol];
           iEl < model_.lp_.Astart_[iCol + 1]; iEl++) {
        HighsInt iRow = model_.lp_.Aindex_[iEl];
        solution_.row_value[iRow] += value * model_.lp_.Avalue_[iEl];
      }
      solution_.col_value[iCol] = value;
    }
    solution_.value_valid = true;
  } else {
    // There is no primal solution: should be so by default
    assert(!solution_.value_valid);
  }
  // There is no dual solution: should be so by default
  assert(!solution_.dual_valid);
  // There is no basis: should be so by default
  assert(!basis_.valid);
  getKktFailures(model_.lp_, solution_, basis_, solution_params);
  // Set the values in HighsInfo instance info_.
  solution_params.objective_function_value = solver.solution_objective_;
  //  Most come from solution_params...
  copyFromSolutionParams(info_, solution_params);
  // ... but others are MIP-specific.
  info_.mip_node_count = solver.node_count_;
  info_.mip_dual_bound = solver.dual_bound_;
  info_.mip_gap =
      100 * std::abs(info_.objective_function_value - info_.mip_dual_bound) /
      std::max(1.0, std::abs(info_.objective_function_value));
  info_.valid = true;
  return return_status;
}

HighsStatus Highs::writeSolution(const std::string filename,
                                 const bool pretty) const {
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  FILE* file;
  bool html;
  call_status = openWriteFile(filename, "writeSolution", file, html);
  return_status =
      interpretCallStatus(call_status, return_status, "openWriteFile");
  if (return_status == HighsStatus::kError) return return_status;
  writeSolutionToFile(file, model_.lp_, basis_, solution_, pretty);
  if (file != stdout) fclose(file);
  return HighsStatus::kOk;
}

void Highs::reportModel() {
  reportLp(options_.log_options, model_.lp_, HighsLogType::kVerbose);
  //  reportHessian();
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
  solution_.col_value.resize(model_.lp_.numCol_);
  solution_.row_value.resize(model_.lp_.numRow_);
  solution_.col_dual.resize(model_.lp_.numCol_);
  solution_.row_dual.resize(model_.lp_.numRow_);
  // Ensure that the HiGHS basis vectors are the right size,
  // invalidating the basis if they aren't
  if ((HighsInt)basis_.col_status.size() != model_.lp_.numCol_) {
    basis_.col_status.resize(model_.lp_.numCol_);
    basis_.valid = false;
  }
  if ((HighsInt)basis_.row_status.size() != model_.lp_.numRow_) {
    basis_.row_status.resize(model_.lp_.numRow_);
    basis_.valid = false;
  }
}

void Highs::setHighsModelStatusAndInfo(const HighsModelStatus model_status) {
  clearUserSolverData();
  model_status_ = model_status;
  scaled_model_status_ = model_status_;
  info_.simplex_iteration_count = iteration_counts_.simplex;
  info_.ipm_iteration_count = iteration_counts_.ipm;
  info_.crossover_iteration_count = iteration_counts_.crossover;
  info_.valid = true;
}

void Highs::setHighsModelStatusBasisSolutionAndInfo() {
  assert(haveHmo("setHighsModelStatusBasisSolutionAndInfo"));
  clearUserSolverData();

  model_status_ = hmos_[0].unscaled_model_status_;
  scaled_model_status_ = hmos_[0].scaled_model_status_;

  basis_ = hmos_[0].basis_;
  solution_ = hmos_[0].solution_;

  info_.simplex_iteration_count = iteration_counts_.simplex;
  info_.ipm_iteration_count = iteration_counts_.ipm;
  info_.crossover_iteration_count = iteration_counts_.crossover;

  HighsSolutionParams& solution_params = hmos_[0].solution_params_;
  info_.primal_solution_status = solution_params.primal_solution_status;
  info_.dual_solution_status = solution_params.dual_solution_status;
  info_.objective_function_value = solution_params.objective_function_value;
  info_.num_primal_infeasibilities = solution_params.num_primal_infeasibility;
  info_.max_primal_infeasibility = solution_params.max_primal_infeasibility;
  info_.sum_primal_infeasibilities = solution_params.sum_primal_infeasibility;
  info_.num_dual_infeasibilities = solution_params.num_dual_infeasibility;
  info_.max_dual_infeasibility = solution_params.max_dual_infeasibility;
  info_.sum_dual_infeasibilities = solution_params.sum_dual_infeasibility;
  info_.valid = true;
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
      highsLogUser(options_.log_options, HighsLogType::kError,
                   "Cannot open writeable file \"%s\" in %s\n",
                   filename.c_str(), method_name.c_str());
      return HighsStatus::kError;
    }
    const char* dot = strrchr(filename.c_str(), '.');
    if (dot && dot != filename) html = strcmp(dot + 1, "html") == 0;
  }
  return HighsStatus::kOk;
}

HighsStatus Highs::getUseModelStatus(
    HighsModelStatus& use_model_status,
    const double unscaled_primal_feasibility_tolerance,
    const double unscaled_dual_feasibility_tolerance,
    const bool rerun_from_logical_basis) {
  if (model_status_ != HighsModelStatus::kNotset) {
    use_model_status = model_status_;
  } else {
    // Handle the case where the status of the unscaled model is not set
    HighsStatus return_status = HighsStatus::kOk;
    HighsStatus call_status;
    const double report = false;  // true;//
    if (unscaledOptimal(unscaled_primal_feasibility_tolerance,
                        unscaled_dual_feasibility_tolerance, report)) {
      use_model_status = HighsModelStatus::kOptimal;
    } else if (rerun_from_logical_basis) {
      std::string save_presolve = options_.presolve;
      basis_.valid = false;
      options_.presolve = kHighsOnString;
      call_status = run();
      return_status = interpretCallStatus(call_status, return_status, "run()");
      options_.presolve = save_presolve;
      if (return_status == HighsStatus::kError) return return_status;

      if (report)
        printf(
            "Unscaled model status was NOTSET: after running from logical "
            "basis it is %s\n",
            modelStatusToString(model_status_).c_str());

      if (model_status_ != HighsModelStatus::kNotset) {
        use_model_status = model_status_;
      } else if (unscaledOptimal(unscaled_primal_feasibility_tolerance,
                                 unscaled_dual_feasibility_tolerance, report)) {
        use_model_status = HighsModelStatus::kOptimal;
      }
    } else {
      // Nothing to be done: use original unscaled model status
      use_model_status = model_status_;
    }
  }
  return HighsStatus::kOk;
}

bool Highs::unscaledOptimal(const double unscaled_primal_feasibility_tolerance,
                            const double unscaled_dual_feasibility_tolerance,
                            const bool report) {
  if (scaled_model_status_ == HighsModelStatus::kOptimal) {
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
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Method %s called without any HighsModelObject\n",
                 method_name.c_str());
#endif
  assert(have_hmo);
  return have_hmo;
}

// Applies checks before returning from run()
HighsStatus Highs::returnFromRun(const HighsStatus run_return_status) {
  assert(!called_return_from_run);
  HighsStatus return_status =
      highsStatusFromHighsModelStatus(scaled_model_status_);
  assert(return_status == run_return_status);
  //  return_status = run_return_status;
  if (hmos_.size() == 0) {
    // No model has been loaded: ensure that any solver data for users
    // in Highs class members are cleared
    clearUserSolverData();
    // Record that returnFromRun() has been called, and stop the Highs
    // run clock
    called_return_from_run = true;
    return returnFromHighs(return_status);
  }
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
  // ToDo: Outcome of Run() should be driven by model_status_, not
  // scaled_model_status_. This is currently done because latter may
  // be optimal but tolerances not satisfied for unscaled model.
  switch (scaled_model_status_) {
      // First consider the error returns
    case HighsModelStatus::kNotset:
    case HighsModelStatus::kLoadError:
    case HighsModelStatus::kModelError:
    case HighsModelStatus::kPresolveError:
    case HighsModelStatus::kSolveError:
    case HighsModelStatus::kPostsolveError:
      clearUserSolverData();
      assert(return_status == HighsStatus::kError);
      break;

      // Then consider the OK returns
    case HighsModelStatus::kModelEmpty:
      clearInfo();
      clearSolution();
      clearBasis();
      assert(model_status_ == scaled_model_status_);
      assert(return_status == HighsStatus::kOk);
      break;

    case HighsModelStatus::kOptimal:
      // The following is an aspiration
      //
      // assert(info_.primal_solution_status == kSolutionStatusFeasible);
      //
      // assert(info_.dual_solution_status == kSolutionStatusFeasible);
      assert(model_status_ == HighsModelStatus::kNotset ||
             model_status_ == HighsModelStatus::kOptimal);
      assert(return_status == HighsStatus::kOk);
      break;

    case HighsModelStatus::kInfeasible:
    case HighsModelStatus::kUnbounded:
    case HighsModelStatus::kObjectiveBound:
    case HighsModelStatus::kObjectiveTarget:
      // For kInfeasible, will not have a basis, if infeasibility was
      // detected in presolve or by IPX without crossover
      assert(model_status_ == scaled_model_status_);
      assert(return_status == HighsStatus::kOk);
      break;

    case HighsModelStatus::kUnboundedOrInfeasible:
      if (options_.allow_unbounded_or_infeasible ||
          (options_.solver == kIpmString && options_.run_crossover)) {
        assert(model_status_ == scaled_model_status_);
        assert(return_status == HighsStatus::kOk);
      } else {
        // This model status is not permitted unless IPM is run without
        // crossover
        highsLogUser(
            options_.log_options, HighsLogType::kError,
            "returnFromHighs: HighsModelStatus::kUnboundedOrInfeasible is not "
            "permitted\n");
        assert(options_.allow_unbounded_or_infeasible);
        return_status = HighsStatus::kError;
      }
      break;

      // Finally consider the warning returns
    case HighsModelStatus::kTimeLimit:
    case HighsModelStatus::kIterationLimit:
    case HighsModelStatus::kUnknown:
      assert(model_status_ == scaled_model_status_);
      assert(return_status == HighsStatus::kWarning);
      break;
    default:
      // All cases should have been considered so assert on reaching here
      assert(1 == 0);
  }
  // Now to check what's available with each model status
  //
  const bool have_info = info_.valid;
  const bool have_primal_solution = solution_.value_valid;
  const bool have_dual_solution = solution_.dual_valid;
  // Can't have a dual solution without a primal solution
  assert(have_primal_solution || !have_dual_solution);
  //  const bool have_solution = have_primal_solution && have_dual_solution;
  const bool have_basis = basis_.valid;
  switch (scaled_model_status_) {
    case HighsModelStatus::kNotset:
    case HighsModelStatus::kLoadError:
    case HighsModelStatus::kModelError:
    case HighsModelStatus::kPresolveError:
    case HighsModelStatus::kSolveError:
    case HighsModelStatus::kPostsolveError:
    case HighsModelStatus::kModelEmpty:
      // No info, primal solution or basis
      assert(have_info == false);
      assert(have_primal_solution == false);
      assert(have_basis == false);
      break;
    case HighsModelStatus::kOptimal:
    case HighsModelStatus::kInfeasible:
    case HighsModelStatus::kUnbounded:
    case HighsModelStatus::kObjectiveBound:
    case HighsModelStatus::kObjectiveTarget:
    case HighsModelStatus::kUnboundedOrInfeasible:
    case HighsModelStatus::kTimeLimit:
    case HighsModelStatus::kIterationLimit:
    case HighsModelStatus::kUnknown:
      // Have info and primal solution (unless infeasible). No primal solution
      // in some other case, too!
      assert(have_info == true);
      //      if (have_primal_solution == true || scaled_model_status_ ==
      //      HighsModelStatus::kInfeasible);
      break;
    default:
      // All cases should have been considered so assert on reaching here
      assert(1 == 0);
  }
  if (have_primal_solution) {
    if (debugPrimalSolutionRightSize(options_, model_.lp_, solution_) ==
        HighsDebugStatus::kLogicalError)
      return_status = HighsStatus::kError;
  }
  if (have_dual_solution) {
    if (debugDualSolutionRightSize(options_, model_.lp_, solution_) ==
        HighsDebugStatus::kLogicalError)
      return_status = HighsStatus::kError;
  }
  if (have_basis) {
    if (debugBasisRightSize(options_, model_.lp_, basis_) ==
        HighsDebugStatus::kLogicalError)
      return_status = HighsStatus::kError;
  }
  if (debugHighsSolution("Return from run()", options_, model_.lp_, solution_,
                         basis_, model_status_,
                         info_) == HighsDebugStatus::kLogicalError)
    return_status = HighsStatus::kError;
  //  getReportKktFailures(options_, model_.lp_, solution_, basis_);
  if (debugInfo(options_, model_.lp_, basis_, solution_, info_,
                scaled_model_status_) == HighsDebugStatus::kLogicalError)
    return_status = HighsStatus::kError;
  // Record that returnFromRun() has been called, and stop the Highs
  // run clock
  called_return_from_run = true;

  return returnFromHighs(return_status);
}

// Applies checks before returning from HiGHS
HighsStatus Highs::returnFromHighs(HighsStatus highs_return_status) {
  HighsStatus return_status = highs_return_status;

  forceHighsSolutionBasisSize();

  const bool consistent = debugBasisConsistent(options_, model_.lp_, basis_) !=
                          HighsDebugStatus::kLogicalError;
  if (!consistent) {
    highsLogUser(
        options_.log_options, HighsLogType::kError,
        "returnFromHighs: Supposed to be a HiGHS basis, but not consistent\n");
    assert(consistent);
    return_status = HighsStatus::kError;
  }

  if (hmos_.size()) {
    bool simplex_lp_ok =
        ekkDebugSimplexLp(hmos_[0]) != HighsDebugStatus::kLogicalError;
    if (!simplex_lp_ok) {
      highsLogUser(options_.log_options, HighsLogType::kError,
                   "returnFromHighs: Simplex LP not OK\n");
      assert(simplex_lp_ok);
      return_status = HighsStatus::kError;
    }
  }
  // Check that returnFromRun() has been called
  if (!called_return_from_run) {
    highsLogDev(
        options_.log_options, HighsLogType::kError,
        "Highs::returnFromHighs() called with called_return_from_run false\n");
    assert(called_return_from_run);
  }
  // Stop the HiGHS run clock if it is running
  if (timer_.runningRunHighsClock()) timer_.stopRunHighsClock();
  return return_status;
}
void Highs::underDevelopmentLogMessage(const std::string method_name) {
  highsLogUser(options_.log_options, HighsLogType::kWarning,
               "Method %s is still under development and behaviour may be "
               "unpredictable\n",
               method_name.c_str());
}
