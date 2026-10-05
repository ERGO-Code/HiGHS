/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "presolve/HPresolve.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>

#include "../extern/pdqsort/pdqsort.h"
#include "Highs.h"
#include "io/HighsIO.h"
#include "lp_data/HConst.h"
#include "lp_data/HStruct.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelUtils.h"
#include "lp_data/HighsSolution.h"
#include "mip/HighsCliqueTable.h"
#include "mip/HighsImplications.h"
#include "mip/HighsMipSolverData.h"
#include "mip/HighsObjectiveFunction.h"
#include "mip/MipTimer.h"
#include "presolve/HPresolveInitialSweep.h"
#include "presolve/HPresolveUtils.h"
#include "presolve/HighsPostsolveStack.h"
#include "presolve/PresolveTimer.h"
#include "util/HFactor.h"
#include "util/HighsCDouble.h"
#include "util/HighsIntegers.h"
#include "util/HighsLinearSumBounds.h"
#include "util/HighsMemoryAllocation.h"
#include "util/HighsSplay.h"
#include "util/HighsUtils.h"

namespace presolve {

#ifndef NDEBUG
void HPresolve::debugPrintRow(HighsPostsolveStack& postsolve_stack,
                              HighsInt row) {
  printf("(row %" HIGHSINT_FORMAT ") %.15g (impl: %.15g) <= ",
         postsolve_stack.getOrigRowIndex(row), model->row_lower_[row],
         impliedRowBounds.getSumLower(row));

  for (const HighsSliceNonzero& nonzero : getSortedRowVector(row)) {
    // for (HighsInt rowiter = rowhead[row]; rowiter != -1; rowiter =
    // ARnext[rowiter]) {
    char colchar =
        model->integrality_[nonzero.index()] == HighsVarType::kInteger ? 'y'
                                                                       : 'x';
    char signchar = nonzero.value() < 0 ? '-' : '+';
    printf("%c%g %c%" HIGHSINT_FORMAT " ", signchar, std::abs(nonzero.value()),
           colchar, postsolve_stack.getOrigColIndex(nonzero.index()));
  }

  printf("<= %.15g (impl: %.15g)\n", model->row_upper_[row],
         impliedRowBounds.getSumUpper(row));
}
#endif

void HPresolve::setInput(HighsLp& model_, const HighsOptions& options_,
                         const HighsInt presolve_reduction_limit,
                         HighsTimer* timer) {
  this->model = &model_;
  this->options = &options_;
  this->timer = timer;
  // Set up the logic to allow presolve rules
  this->chooseRules();
  // Set up profiling for presolve rules and logging for their effectiveness
  analysis_.setup(this->model, this->options, this->numDeletedRows,
                  this->numDeletedCols, this->timer);
  analysis_.presolveTimerStart(kPresolveClockPresolve);

  if (mipsolver == nullptr) {
    this->primal_feastol = options->primal_feasibility_tolerance;
    model->integrality_.assign(model->num_col_, HighsVarType::kContinuous);
  } else
    this->primal_feastol = options->mip_feasibility_tolerance;

  // Take value passed in as reduction limit, allowing different
  // values to be used for initial presolve, and after restart
  this->reductionLimit =
      presolve_reduction_limit < 0 ? kHighsSize_tInf : presolve_reduction_limit;
  if (options->presolve != kHighsOffString &&
      reductionLimit < kHighsSize_tInf) {
    highsLogDev(options->log_options, HighsLogType::kInfo,
                "HPresolve::setInput reductionLimit = %d\n",
                static_cast<int>(this->reductionLimit));
  }
  // last_reduction_ is used to identify when HPresolve::checkLimits
  // is called for the first time following a reduction
  this->last_reduction_ = 0;
}

// for MIP presolve
void HPresolve::setInput(HighsMipSolver& mipsolver,
                         const HighsInt presolve_reduction_limit) {
  this->mipsolver = &mipsolver;
  presolveCliqueTable.rebuild(mipsolver.mipdata_->cliquetable);

  probingContingent = 1000;
  probingNumDelCol = 0;
  numProbed = 0;
  numProbes.assign(mipsolver.numCol(), 0);

  if (mipsolver.model_ != &mipsolver.mipdata_->presolvedModel) {
    mipsolver.mipdata_->presolvedModel = *mipsolver.model_;
    mipsolver.model_ = &mipsolver.mipdata_->presolvedModel;
  } else {
    mipsolver.mipdata_->presolvedModel.col_lower_ =
        mipsolver.mipdata_->getDomain().col_lower_;
    mipsolver.mipdata_->presolvedModel.col_upper_ =
        mipsolver.mipdata_->getDomain().col_upper_;
  }

  mipsolver.mipdata_->implications.resize(mipsolver.numCol(),
                                          mipsolver.numRow());

  setInput(mipsolver.mipdata_->presolvedModel, *mipsolver.options_mip_,
           presolve_reduction_limit, &mipsolver.timer_);
}

bool HPresolve::okSetupPresolveDataStructures() {
  analysis_.presolveTimerStart(kPresolveClockSetupResize);
  if (!okResize(colLowerSource, model->num_col_, HighsInt{-1})) return false;
  if (!okResize(colUpperSource, model->num_col_, HighsInt{-1})) return false;
  if (!okResize(implColLower, model->num_col_, -kHighsInf)) return false;
  if (!okResize(implColUpper, model->num_col_, kHighsInf)) return false;
  if (!okResize(colImplSourceByRow, model->num_row_)) return false;
  if (!okResize(implRowDualSourceByCol, model->num_col_)) return false;
  if (!okResize(rowDualLower, model->num_row_, -kHighsInf)) return false;
  if (!okResize(rowDualUpper, model->num_row_, kHighsInf)) return false;
  if (!okResize(implRowDualLower, model->num_row_, -kHighsInf)) return false;
  if (!okResize(implRowDualUpper, model->num_row_, kHighsInf)) return false;
  if (!okResize(rowDualUpperSource, model->num_row_, HighsInt{-1}))
    return false;
  if (!okResize(rowDualLowerSource, model->num_row_, HighsInt{-1}))
    return false;

  for (HighsInt i = 0; i != model->num_row_; ++i) {
    if (model->row_lower_[i] == -kHighsInf) rowDualUpper[i] = 0;
    if (model->row_upper_[i] == kHighsInf) rowDualLower[i] = 0;
  }

  analysis_.presolveTimerStop(kPresolveClockSetupResize);

  analysis_.presolveTimerStart(kPresolveClockSetupToCsc);
  if (model->a_matrix_.isRowwise()) {
    // Does this even happen?
    assert(model->a_matrix_.isColwise());
    if (!okFromCSR(model->a_matrix_.value_, model->a_matrix_.index_,
                   model->a_matrix_.start_))
      return false;
  } else {
    if (!okFromCSC(model->a_matrix_.value_, model->a_matrix_.index_,
                   model->a_matrix_.start_))
      return false;
  }
  analysis_.presolveTimerStop(kPresolveClockSetupToCsc);

  // numDeletedCols and numDeletedRows are not cumulative through the
  // whole of presolve, but "since the last time the model had no
  // deleted columns or rows" - ie from initialisation here, or from a
  // call to shrinkProblem
  numDeletedCols = 0;
  numDeletedRows = 0;
  // number of equations checked by the single-equation extension of dual
  // fixing
  numSingleEquationChecked = 0;
  // Need to reset current number of deleted rows and columns in logging
  analysis_.resetNumDeleted();

  analysis_.presolveTimerStart(kPresolveClockSetupResize);
  // initialize everything as changed, but do not add all indices
  // since the first thing presolve will do is a scan for easy reductions
  // of each row and column and set the flag of processed columns to false
  // from then on they are added to the vector whenever there are changes
  if (!okResize(changedRowFlag, model->num_row_, uint8_t{1})) return false;
  if (!okResize(rowDeleted, model->num_row_)) return false;
  if (!okReserve(changedRowIndices, model->num_row_)) return false;
  if (!okResize(changedColFlag, model->num_col_, uint8_t{1})) return false;
  if (!okResize(colDeleted, model->num_col_)) return false;
  if (!okReserve(changedColIndices, model->num_col_)) return false;
  if (!okReserve(liftingOpportunities, model->num_row_)) return false;
  if (!okResize(singleEquationChecked, model->num_row_)) return false;
  analysis_.presolveTimerStop(kPresolveClockSetupResize);
  return true;
}

void HPresolve::setupSubstitutionOpportunities() {
  for (HighsInt row = 0; row != model->num_row_; ++row) {
    if (!isDualImpliedFree(row)) continue;
    for (const HighsSliceNonzero& nonzero : getRowVector(row)) {
      if (isImpliedFree(nonzero.index()))
        substitutionOpportunities.emplace_back(row, nonzero.index());
    }
  }
}

bool HPresolve::rowCoefficientsIntegral(HighsInt row, double scale) const {
  for (const HighsSliceNonzero& nz : getRowVector(row)) {
    if (fractionality(nz.value() * scale) > options->small_matrix_value)
      return false;
  }

  return true;
}

bool HPresolve::isLowerImplied(HighsInt col) const {
  return (model->col_lower_[col] == -kHighsInf ||
          implColLower[col] >= model->col_lower_[col] - primal_feastol);
}

bool HPresolve::isLowerStrictlyImplied(HighsInt col, double* tolerance) const {
  return (model->col_lower_[col] == -kHighsInf ||
          implColLower[col] >
              model->col_lower_[col] +
                  (tolerance != nullptr ? *tolerance : primal_feastol));
}

bool HPresolve::isUpperImplied(HighsInt col) const {
  return (model->col_upper_[col] == kHighsInf ||
          implColUpper[col] <= model->col_upper_[col] + primal_feastol);
}

bool HPresolve::isUpperStrictlyImplied(HighsInt col, double* tolerance) const {
  return (model->col_upper_[col] == kHighsInf ||
          implColUpper[col] <
              model->col_upper_[col] -
                  (tolerance != nullptr ? *tolerance : primal_feastol));
}

bool HPresolve::isBinary(HighsInt col) const {
  return model->integrality_[col] == HighsVarType::kInteger &&
         model->col_lower_[col] == 0.0 && model->col_upper_[col] == 1.0;
}

bool HPresolve::isImpliedFree(HighsInt col) const {
  return isLowerImplied(col) && isUpperImplied(col);
}

bool HPresolve::isDualImpliedFree(HighsInt row) const {
  return isEquation(row) ||
         (model->row_upper_[row] != kHighsInf &&
          implRowDualUpper[row] <= options->dual_feasibility_tolerance) ||
         (model->row_lower_[row] != -kHighsInf &&
          implRowDualLower[row] >= -options->dual_feasibility_tolerance);
}

void HPresolve::dualImpliedFreeGetRhsAndRowType(
    HighsInt row, double& rhs, HighsPostsolveStack::RowType& rowType,
    bool relaxRowDualBounds) {
  assert(isDualImpliedFree(row));
  if (isEquation(row)) {
    rowType = HighsPostsolveStack::RowType::kEq;
    rhs = model->row_upper_[row];
  } else if (model->row_upper_[row] != kHighsInf &&
             implRowDualUpper[row] <= options->dual_feasibility_tolerance) {
    rowType = HighsPostsolveStack::RowType::kLeq;
    rhs = model->row_upper_[row];
    if (relaxRowDualBounds) changeRowDualUpper(row, kHighsInf);
  } else {
    rowType = HighsPostsolveStack::RowType::kGeq;
    rhs = model->row_lower_[row];
    if (relaxRowDualBounds) changeRowDualLower(row, -kHighsInf);
  }
}

bool HPresolve::isEquation(HighsInt row) const {
  return (model->row_lower_[row] == model->row_upper_[row]);
}

bool HPresolve::isRanged(HighsInt row) const {
  return (model->row_lower_[row] != -kHighsInf &&
          model->row_upper_[row] != kHighsInf);
}

bool HPresolve::isRedundant(HighsInt row, double sumLower,
                            double sumUpper) const {
  return sumLower >= model->row_lower_[row] - primal_feastol &&
         sumUpper <= model->row_upper_[row] + primal_feastol;
}

bool HPresolve::isRedundant(HighsInt row) const {
  return isRedundant(row, impliedRowBounds.getSumLower(row),
                     impliedRowBounds.getSumUpper(row));
}

bool HPresolve::yieldsImpliedLowerBound(HighsInt row, double val) const {
  return ((val < 0 && model->row_upper_[row] != kHighsInf) ||
          (val > 0 && model->row_lower_[row] != -kHighsInf));
}

bool HPresolve::yieldsImpliedUpperBound(HighsInt row, double val) const {
  return yieldsImpliedLowerBound(row, -val);
}

bool HPresolve::isImpliedEquationAtLower(HighsInt row) const {
  // if the implied lower bound on a row dual is strictly positive then the row
  // is an implied equation (using its lower bound) due to complementary
  // slackness
  bool isLbndPositive =
      implRowDualLower[row] > options->dual_feasibility_tolerance;
  assert(!isLbndPositive || model->row_lower_[row] != -kHighsInf);
  return isLbndPositive;
}

bool HPresolve::isImpliedEquationAtUpper(HighsInt row) const {
  // if the implied upper bound on a row dual is strictly negative then the row
  // is an implied equation (using its upper bound) due to complementary
  // slackness
  bool isUbndNegative =
      implRowDualUpper[row] < -options->dual_feasibility_tolerance;
  assert(!isUbndNegative || model->row_upper_[row] != kHighsInf);
  return isUbndNegative;
}

HPresolve::StatusResult HPresolve::isImpliedIntegral(HighsInt col) {
  // check if the integer constraint on a variable is implied by the model
  assert(model->integrality_[col] == HighsVarType::kInteger);

  bool runDualDetection = true;

  for (const HighsSliceNonzero& nz : getColumnVector(col)) {
    // if not all other columns are integer, skip row and also do not try the
    // dual detection in the second loop as it must hold for all rows
    if (rowsize[nz.index()] < 2 ||
        rowsizeInteger[nz.index()] < rowsize[nz.index()]) {
      runDualDetection = false;
      continue;
    }

    double rowLower = isImpliedEquationAtUpper(nz.index())
                          ? model->row_upper_[nz.index()]
                          : model->row_lower_[nz.index()];

    double rowUpper = isImpliedEquationAtLower(nz.index())
                          ? model->row_lower_[nz.index()]
                          : model->row_upper_[nz.index()];

    if (rowUpper == rowLower) {
      // if there is an equation the dual detection does not need to be tried
      runDualDetection = false;
      double scale = 1.0 / nz.value();

      if (!rowCoefficientsIntegral(nz.index(), scale)) continue;

      if (fractionality(rowLower * scale) > primal_feastol)
        return StatusResult(Result::kPrimalInfeasible);

      return StatusResult(true);
    }
  }

  if (!runDualDetection) return StatusResult(false);

  for (const HighsSliceNonzero& nz : getColumnVector(col)) {
    double scale = 1.0 / nz.value();
    // if row coefficients are not integral, variable is not (implied) integral
    if (!rowCoefficientsIntegral(nz.index(), scale)) return StatusResult(false);
    if (model->row_upper_[nz.index()] != kHighsInf) {
      // right-hand side: scale, round down and unscale again
      double rUpper =
          std::abs(nz.value()) *
          std::floor(model->row_upper_[nz.index()] * std::abs(scale) +
                     primal_feastol);
      // check if modification is large enough
      if (std::abs(model->row_upper_[nz.index()] - rUpper) >
          options->small_matrix_value)
        changeRowUpper(nz.index(), rUpper);
    }
    if (model->row_lower_[nz.index()] != -kHighsInf) {
      // left-hand side: scale, round up and unscale again
      double rLower =
          std::abs(nz.value()) *
          std::ceil(model->row_lower_[nz.index()] * std::abs(scale) -
                    primal_feastol);
      // check if modification is large enough
      if (std::abs(model->row_lower_[nz.index()] - rLower) >
          options->small_matrix_value)
        changeRowLower(nz.index(), rLower);
    }
  }

  return StatusResult(true);
}

HPresolve::StatusResult HPresolve::isImpliedInteger(HighsInt col) const {
  // check if a continuous variable is implied integer
  assert(model->integrality_[col] == HighsVarType::kContinuous);

  bool runDualDetection = true;

  for (const HighsSliceNonzero& nz : getColumnVector(col)) {
    // if not all other columns are integer, skip row and also do not try the
    // dual detection in the second loop as it must hold for all rows
    if (rowsize[nz.index()] < 2 ||
        rowsizeInteger[nz.index()] + rowsizeImplInt[nz.index()] <
            rowsize[nz.index()] - 1) {
      runDualDetection = false;
      continue;
    }

    double rowLower = isImpliedEquationAtUpper(nz.index())
                          ? model->row_upper_[nz.index()]
                          : model->row_lower_[nz.index()];

    double rowUpper = isImpliedEquationAtLower(nz.index())
                          ? model->row_lower_[nz.index()]
                          : model->row_upper_[nz.index()];

    if (rowUpper == rowLower) {
      // if there is an equation the dual detection does not need to be tried
      runDualDetection = false;
      double scale = 1.0 / nz.value();

      if (fractionality(rowLower * scale) > primal_feastol) continue;

      if (!rowCoefficientsIntegral(nz.index(), scale)) continue;

      return StatusResult(true);
    }
  }

  if (!runDualDetection) return StatusResult(false);

  if ((model->col_lower_[col] != -kHighsInf &&
       fractionality(model->col_lower_[col]) > options->small_matrix_value) ||
      (model->col_upper_[col] != kHighsInf &&
       fractionality(model->col_upper_[col]) > options->small_matrix_value))
    return StatusResult(false);

  for (const HighsSliceNonzero& nz : getColumnVector(col)) {
    double scale = 1.0 / nz.value();
    if (model->row_upper_[nz.index()] != kHighsInf &&
        fractionality(model->row_upper_[nz.index()] * scale) > primal_feastol)
      return StatusResult(false);

    if (model->row_lower_[nz.index()] != -kHighsInf &&
        fractionality(model->row_lower_[nz.index()] * scale) > primal_feastol)
      return StatusResult(false);

    if (!rowCoefficientsIntegral(nz.index(), scale)) return StatusResult(false);
  }

  return StatusResult(true);
}

HPresolve::StatusResult HPresolve::convertImpliedInteger(HighsInt col,
                                                         HighsInt row,
                                                         bool skipInputChecks) {
  // return if column was deleted
  if (colDeleted[col]) return StatusResult(false);

  // return if column is not continuous or cannot be converted to an implied
  // integer
  if (!skipInputChecks) {
    if (model->integrality_[col] != HighsVarType::kContinuous)
      return StatusResult(false);
    StatusResult impliedInteger = isImpliedInteger(col);
    if (!impliedInteger) return impliedInteger;
  }

  // convert to implied integer
  model->integrality_[col] = HighsVarType::kImplicitInteger;

  if (row != -1) {
    // use row index supplied by caller (e.g. singleton)
    ++rowsizeImplInt[row];
  } else {
    // iterate over rows
    for (const HighsSliceNonzero& nonzero : getColumnVector(col))
      ++rowsizeImplInt[nonzero.index()];
  }

  // Potentially strengthen existing bound using implied bounds
  // If not done then there may be a stronger fractional
  // implied bound for a non-continuous column, which causes
  // errors in rules downstream that assume integrality.
  // changeColBounds will perform rounding
  double newLower = model->col_lower_[col];
  double newUpper = model->col_upper_[col];
  if (implColLower[col] > newLower + primal_feastol)
    newLower = implColLower[col];
  if (implColUpper[col] < newUpper - primal_feastol)
    newUpper = implColUpper[col];
  // round and update bounds
  return StatusResult(changeColBounds(col, newLower, newUpper));
}

void HPresolve::chooseRules() {
  const bool silent = silentLog();
  // By default all presolve rules are allowed
  this->allow_rule_.assign(kPresolveRuleCount, true);

  // All rules except kPresolveRuleDominatedCol can be switched off,
  // although the most fundamental can only be switched off in initial
  // sweep
  auto allowedOffInInitialSweep = [&](const HighsInt rule_type) {
    if (rule_type == kPresolveRuleEmptyRow) return true;
    if (rule_type == kPresolveRuleSingletonRow) return true;
    if (rule_type == kPresolveRuleRedundantRow) return true;
    if (rule_type == kPresolveRuleEmptyCol) return true;
    if (rule_type == kPresolveRuleFixedCol) return true;
    return false;
  };

  auto logRule = [&](const HighsInt rule_type) {
    highsLogUser(options->log_options, HighsLogType::kInfo,
                 " %1s Rule %2d (set bit %2d = %7d): %s\n",
                 allowedOffInInitialSweep(rule_type) ? "*" : " ",
                 int(rule_type), int(rule_type), int(1 << rule_type),
                 utilPresolveRuleTypeToString(rule_type).c_str());
  };

  if (!silent && options->log_dev_level) {
    // State which rules can be off, and what bit to set
    highsLogUser(options->log_options, HighsLogType::kInfo,
                 "Permitted suppression of presolve rules via "
                 "presolve_rule_off option:\n");
    for (HighsInt rule_type = kPresolveRuleMin; rule_type < kPresolveRuleCount;
         rule_type++) {
      // Only kPresolveRuleDominatedCol cannot be switched off
      if (rule_type == kPresolveRuleDominatedCol) continue;
      logRule(rule_type);
    }
    highsLogUser(options->log_options, HighsLogType::kInfo,
                 " * Only in initial sweep\n");
  }

  if (options->presolve_light == kHighsOnString) {
    // Switch off the rules not used in presolve_light mode
    allow_rule_[kPresolveRuleDependentEquations] = false;
    allow_rule_[kPresolveRuleDependentFreeCols] = false;
    allow_rule_[kPresolveRuleAggregator] = false;
    allow_rule_[kPresolveRuleParallelRowsAndCols] = false;
    allow_rule_[kPresolveRuleSparsify] = false;
    allow_rule_[kPresolveRuleProbing] = false;
    allow_rule_[kPresolveRuleEnumeration] = false;
    allow_rule_[kPresolveRuleDualFixing] = false;
    allow_rule_[kPresolveRuleColStuffing] = false;
    allow_rule_[kPresolveRuleFourierMotzkin] = false;
  }

  if (options->presolve_rule_off) {
    // Some presolve rules are off
    //
    // Transform options->presolve_rule_off into logical settings in
    // allow_rule_[*], commenting on the rules switched off
    if (!silent)
      highsLogUser(options->log_options, HighsLogType::kInfo,
                   "Presolve rules switched off:\n");
    bool off_in_initial_sweep = false;
    for (HighsInt rule_type = kPresolveRuleMin; rule_type < kPresolveRuleCount;
         rule_type++) {
      HighsInt bit = 1 << rule_type;
      // Identify whether this rule is allowed
      if (!(options->presolve_rule_off & bit)) continue;
      if (rule_type != kPresolveRuleDominatedCol) {
        off_in_initial_sweep =
            allowedOffInInitialSweep(rule_type) || off_in_initial_sweep;
        allow_rule_[rule_type] = false;
        if (!silent) logRule(rule_type);
      } else {
        // This is a rule that cannot be switched off so, if an
        // attempt is made, don't allow it to be off and possibly
        // comment negatively
        if (!silent)
          highsLogUser(options->log_options, HighsLogType::kWarning,
                       "Cannot disallow rule %2d (bit %2d = %5d): %s\n",
                       int(rule_type), int(rule_type), int(bit),
                       utilPresolveRuleTypeToString(rule_type).c_str());
      }
    }
    if (!silent && off_in_initial_sweep)
      highsLogUser(options->log_options, HighsLogType::kInfo,
                   " * Only in initial sweep\n");
  }

  may_require_basis_postsolve_ = false;
  may_require_primal_dual_postsolve_ = false;
  if (mipsolver == nullptr) {
    // Record whether one or both of basis and primal-dual postsolve
    // may be required, and switch off rules that are correspondingly
    // incompatible
    if (mayRequireBasisPostsolve(options)) {
      may_require_basis_postsolve_ = true;
      // Cannot use sparsify or any parallel row rule using
      // EqualityRowAddition(s). Latter are avoided logically in
      // HPresolve::detectParallelRowsAndCols
      allow_rule_[kPresolveRuleSparsify] = false;
    }
    if (mayRequirePrimalDualPostsolve(options)) {
      may_require_primal_dual_postsolve_ = true;
      // Cannot use weakly dominated column rule (see
      // test-weakly-dominated-column-primal-dual-postsolve)
      allow_rule_[kPresolveRuleWeaklyDominatedCol] = false;
      // Cannot use kPresolveRuleDualFixing rule (see
      // test-dual-fixing-primal-dual-postsolve)
      allow_rule_[kPresolveRuleDualFixing] = false;
    }
  }
}

void HPresolve::link(HighsInt pos) {
  Anext[pos] = colhead[Acol[pos]];
  Aprev[pos] = -1;
  colhead[Acol[pos]] = pos;
  if (Anext[pos] != -1) Aprev[Anext[pos]] = pos;

  ++colsize[Acol[pos]];

  ARleft[pos] = -1;
  ARright[pos] = -1;
  auto get_row_left = [&](HighsInt pos) -> HighsInt& { return ARleft[pos]; };
  auto get_row_right = [&](HighsInt pos) -> HighsInt& { return ARright[pos]; };
  auto get_row_key = [&](HighsInt pos) { return Acol[pos]; };
  highs_splay_link(pos, rowroot[Arow[pos]], get_row_left, get_row_right,
                   get_row_key);

  ++rowsize[Arow[pos]];
  if (model->integrality_[Acol[pos]] == HighsVarType::kInteger)
    ++rowsizeInteger[Arow[pos]];
  else if (model->integrality_[Acol[pos]] == HighsVarType::kImplicitInteger)
    ++rowsizeImplInt[Arow[pos]];
}

void HPresolve::unlink(HighsInt pos) {
  HighsInt next = Anext[pos];
  HighsInt prev = Aprev[pos];

  if (next != -1) Aprev[next] = prev;

  if (prev != -1)
    Anext[prev] = next;
  else
    colhead[Acol[pos]] = next;
  --colsize[Acol[pos]];

  if (!colDeleted[Acol[pos]]) {
    if (colsize[Acol[pos]] == 1)
      singletonColumns.push_back(Acol[pos]);
    else
      markChangedCol(Acol[pos]);
  }

  auto get_row_left = [&](HighsInt pos) -> HighsInt& { return ARleft[pos]; };
  auto get_row_right = [&](HighsInt pos) -> HighsInt& { return ARright[pos]; };
  auto get_row_key = [&](HighsInt pos) { return Acol[pos]; };
  highs_splay_unlink(pos, rowroot[Arow[pos]], get_row_left, get_row_right,
                     get_row_key);
  --rowsize[Arow[pos]];
  if (model->integrality_[Acol[pos]] == HighsVarType::kInteger)
    --rowsizeInteger[Arow[pos]];
  else if (model->integrality_[Acol[pos]] == HighsVarType::kImplicitInteger)
    --rowsizeImplInt[Arow[pos]];

  if (!rowDeleted[Arow[pos]]) {
    if (rowsize[Arow[pos]] == 1)
      singletonRows.push_back(Arow[pos]);
    else
      markChangedRow(Arow[pos]);
  }

  matrixNonZeroChanged(Arow[pos], Acol[pos], Avalue[pos], 0.0,
                       rowDeleted[Arow[pos]], colDeleted[Acol[pos]]);

  // remove non-zero
  Avalue[pos] = 0;

  freeslots.push_back(pos);
}

void HPresolve::markChangedRow(HighsInt row) {
  if (!changedRowFlag[row]) {
    changedRowIndices.push_back(row);
    changedRowFlag[row] = true;
  }
  resetSingleEquationChecked(row);
}

void HPresolve::markChangedCol(HighsInt col) {
  if (!changedColFlag[col]) {
    changedColIndices.push_back(col);
    changedColFlag[col] = true;
  }
  for (const auto& nz : getColumnVector(col))
    if (!resetSingleEquationChecked(nz.index())) break;
}

bool HPresolve::resetSingleEquationChecked(HighsInt row) {
  // clears the flag of the row and returns whether any row is still flagged
  if (singleEquationChecked[row]) {
    singleEquationChecked[row] = false;
    --numSingleEquationChecked;
  }
  return numSingleEquationChecked > 0;
}

double HPresolve::getMaxAbsColVal(HighsInt col) const {
  double maxVal = 0.0;
  for (const auto& nz : getColumnVector(col))
    maxVal = std::max(std::abs(nz.value()), maxVal);
  return maxVal;
}

double HPresolve::getMaxAbsRowVal(HighsInt row) const {
  double maxVal = 0.0;

  for (const auto& nz : getRowVector(row))
    maxVal = std::max(std::abs(nz.value()), maxVal);

  return maxVal;
}

bool HPresolve::checkUpdateRowDualImpliedBounds(HighsInt col,
                                                double* dualRowLower,
                                                double* dualRowUpper) const {
  // check if implied bounds of row duals in given column can be updated (i.e.
  // dual row has finite bounds and number of infinite contributions to
  // corresponding activity bounds is at most one)

  // if the column has an infinite lower bound the reduced cost cannot be
  // positive, i.e. the column corresponds to a <= constraint in the dual with
  // right hand side -cost which becomes a >= constraint with side +cost.
  // Furthermore, we can ignore strictly redundant primal
  // column bounds and treat them as if they are infinite
  double impliedMargin = colsize[col] != 1 ? primal_feastol : -primal_feastol;

  double myDualRowLower = isLowerStrictlyImplied(col, &impliedMargin)
                              ? model->col_cost_[col]
                              : -kHighsInf;

  double myDualRowUpper = isUpperStrictlyImplied(col, &impliedMargin)
                              ? model->col_cost_[col]
                              : kHighsInf;

  if (dualRowLower != nullptr) *dualRowLower = myDualRowLower;
  if (dualRowUpper != nullptr) *dualRowUpper = myDualRowUpper;

  return (myDualRowLower != -kHighsInf &&
          impliedDualRowBounds.getNumInfSumUpperOrig(col) <= 1) ||
         (myDualRowUpper != kHighsInf &&
          impliedDualRowBounds.getNumInfSumLowerOrig(col) <= 1);
}

void HPresolve::updateRowDualImpliedBounds(HighsInt row, HighsInt col,
                                           double val) {
  // propagate implied row dual bound
  double dualRowLower, dualRowUpper;
  if (!checkUpdateRowDualImpliedBounds(col, &dualRowLower, &dualRowUpper))
    return;

  const double threshold = 1000 * options->dual_feasibility_tolerance;

  auto checkImpliedBound = [&](HighsInt row, HighsInt col, double val,
                               double dualRowBnd, double residualAct,
                               HighsInt direction) {
    if (direction * residualAct <= -kHighsInf) return;
    double impliedBound = static_cast<double>(
        (static_cast<HighsCDouble>(dualRowBnd) - residualAct) / val);

    if (std::abs(impliedBound) * kHighsTiny >
        options->dual_feasibility_tolerance)
      return;

    if (direction * val > 0) {
      // only tighten bound if it is tighter by a wide enough margin
      if (impliedBound < implRowDualUpper[row] - threshold)
        changeImplRowDualUpper(row, impliedBound, col);
    } else {
      if (impliedBound > implRowDualLower[row] + threshold)
        changeImplRowDualLower(row, impliedBound, col);
    }
  };

  if (dualRowUpper != kHighsInf)
    checkImpliedBound(
        row, col, val, dualRowUpper,
        impliedDualRowBounds.getResidualSumLowerOrig(col, row, val),
        HighsInt{1});

  if (dualRowLower != -kHighsInf)
    checkImpliedBound(
        row, col, val, dualRowLower,
        impliedDualRowBounds.getResidualSumUpperOrig(col, row, val),
        HighsInt{-1});
}

void HPresolve::updateRowDualImpliedBounds(HighsInt col) {
  // update dual implied bounds of all rows in given column
  assert(col >= 0 && col < model->num_col_);
  if (!checkUpdateRowDualImpliedBounds(col)) return;
  for (const HighsSliceNonzero& nonzero : getColumnVector(col))
    updateRowDualImpliedBounds(nonzero.index(), col, nonzero.value());
}

bool HPresolve::checkUpdateColImpliedBounds(HighsInt row, double* rowLower,
                                            double* rowUpper) const {
  // check if implied bounds of columns in given row can be updated (i.e. if
  // row's left-hand or right-hand side is finite and number of infinite
  // contributions to corresponding activity bounds is at most one)
  double myRowLower = isImpliedEquationAtUpper(row) ? model->row_upper_[row]
                                                    : model->row_lower_[row];
  double myRowUpper = isImpliedEquationAtLower(row) ? model->row_lower_[row]
                                                    : model->row_upper_[row];
  assert(myRowLower != kHighsInf);
  assert(myRowUpper != -kHighsInf);

  if (rowLower != nullptr) *rowLower = myRowLower;
  if (rowUpper != nullptr) *rowUpper = myRowUpper;

  return (myRowLower != -kHighsInf &&
          impliedRowBounds.getNumInfSumUpperOrig(row) <= 1) ||
         (myRowUpper != kHighsInf &&
          impliedRowBounds.getNumInfSumLowerOrig(row) <= 1);
}

HPresolve::Result HPresolve::updateColImpliedBounds(HighsInt row, HighsInt col,
                                                    double val) {
  // propagate implied column bound upper bound if row has an upper bound
  double rowLower, rowUpper;
  if (!checkUpdateColImpliedBounds(row, &rowLower, &rowUpper))
    return Result::kOk;

  const double threshold = 1000 * primal_feastol;

  auto checkImpliedBound = [&](HighsInt row, HighsInt col, double val,
                               double rowBnd, double residualAct,
                               HighsInt direction) {
    if (direction * residualAct <= -kHighsInf) return Result::kOk;
    double impliedBound = static_cast<double>(
        (static_cast<HighsCDouble>(rowBnd) - residualAct) / val);

    if (std::abs(impliedBound) * kHighsTiny > primal_feastol)
      return Result::kOk;

    // do not use the implied bound if this a not a model row, since the
    // row can be removed and should not be used, e.g., to identify a
    // column as implied free
    bool useImplBound = mipsolver == nullptr ||
                        !mipsolver->mipdata_->postSolveStack.isCutRow(row);

    if (direction * val > 0) {
      // upper bound
      // when solving a MIP, keep (1) tighter bounds on integer variables or (2)
      // implied bound (as column bound)
      bool updateColBound =
          mipsolver != nullptr &&
          ((model->integrality_[col] != HighsVarType::kContinuous &&
            impliedBound < model->col_upper_[col] - primal_feastol) ||
           (!useImplBound &&
            impliedBound < model->col_upper_[col] - threshold));

      // only tighten bound if it is tighter by a wide enough margin
      if (useImplBound && impliedBound < implColUpper[col] - threshold)
        changeImplColUpper(col, impliedBound, row);
      if (updateColBound)
        HPRESOLVE_CHECKED_CALL(changeColUpper(col, impliedBound));
    } else {
      // lower bound
      // when solving a MIP, keep (1) tighter bounds on integer variables or (2)
      // implied bound (as column bound)
      bool updateColBound =
          mipsolver != nullptr &&
          ((model->integrality_[col] != HighsVarType::kContinuous &&
            impliedBound > model->col_lower_[col] + primal_feastol) ||
           (!useImplBound &&
            impliedBound > model->col_lower_[col] + threshold));

      // only tighten bound if it is tighter by a wide enough margin
      if (useImplBound && impliedBound > implColLower[col] + threshold)
        changeImplColLower(col, impliedBound, row);
      if (updateColBound)
        HPRESOLVE_CHECKED_CALL(changeColLower(col, impliedBound));
    }
    return Result::kOk;
  };

  if (rowUpper != kHighsInf)
    HPRESOLVE_CHECKED_CALL(checkImpliedBound(
        row, col, val, rowUpper,
        impliedRowBounds.getResidualSumLowerOrig(row, col, val), HighsInt{1}));

  if (rowLower != -kHighsInf)
    HPRESOLVE_CHECKED_CALL(checkImpliedBound(
        row, col, val, rowLower,
        impliedRowBounds.getResidualSumUpperOrig(row, col, val), HighsInt{-1}));

  return Result::kOk;
}

HPresolve::Result HPresolve::updateColImpliedBounds(HighsInt row) {
  // update implied bounds of all columns in given row
  assert(row >= 0 && row < model->num_row_);
  if (!checkUpdateColImpliedBounds(row)) return Result::kOk;
  for (const HighsSliceNonzero& nonzero : getRowVector(row))
    HPRESOLVE_CHECKED_CALL(
        updateColImpliedBounds(row, nonzero.index(), nonzero.value()));
  return Result::kOk;
}

void HPresolve::resetColImpliedBounds(HighsInt col, HighsInt row) {
  assert(row == -1 || colLowerSource[col] == row || colUpperSource[col] == row);
  if (!colDeleted[col]) {
    // set implied bounds to infinite values if (1) they were deduced from the
    // given row or (2) no row was given
    if (colLowerSource[col] != -1 && (row == -1 || colLowerSource[col] == row))
      changeImplColLower(col, -kHighsInf, -1);
    if (colUpperSource[col] != -1 && (row == -1 || colUpperSource[col] == row))
      changeImplColUpper(col, kHighsInf, -1);
  } else if (row != -1 && !rowDeleted[row]) {
    // remove column from row-wise implied bound storage
    colImplSourceByRow[row].erase(col);
  }
}

void HPresolve::resetRowDualImpliedBounds(HighsInt row, HighsInt col) {
  assert(col == -1 || rowDualLowerSource[row] == col ||
         rowDualUpperSource[row] == col);
  if (!rowDeleted[row]) {
    // set implied bounds to infinite values if (1) they were deduced from the
    // given column or (2) no column was given
    if (rowDualLowerSource[row] != -1 &&
        (col == -1 || rowDualLowerSource[row] == col))
      changeImplRowDualLower(row, -kHighsInf, -1);
    if (rowDualUpperSource[row] != -1 &&
        (col == -1 || rowDualUpperSource[row] == col))
      changeImplRowDualUpper(row, kHighsInf, -1);
  } else if (col != -1 && !colDeleted[col]) {
    // remove row from column-wise implied bound storage
    implRowDualSourceByCol[col].erase(row);
  }
}

void HPresolve::resetColImpliedBoundsDerivedFromRow(HighsInt row) {
  // reset implied column bounds affected by a modification in a row
  // (removed / added non-zeros, etc.)
  if (colImplSourceByRow[row].empty()) return;
  std::set<HighsInt> affectedCols;
  affectedCols.swap(colImplSourceByRow[row]);
  for (const HighsInt& col : affectedCols) {
    // set implied bounds to infinite values if they were deduced from the
    // given row
    resetColImpliedBounds(col, row);
  }
}

void HPresolve::resetRowDualImpliedBoundsDerivedFromCol(HighsInt col) {
  // reset implied row dual bounds affected by a modification in a column
  // (removed / added non-zeros, etc.)
  if (implRowDualSourceByCol[col].empty()) return;
  std::set<HighsInt> affectedRows;
  affectedRows.swap(implRowDualSourceByCol[col]);
  for (const HighsInt& row : affectedRows) {
    // set implied bounds to infinite values if they were deduced from the
    // given column
    resetRowDualImpliedBounds(row, col);
  }
}

void HPresolve::matrixNonZeroChanged(HighsInt row, HighsInt col, double oldCoef,
                                     double newCoef, bool rowIsDeleted,
                                     bool colIsDeleted) {
  if (!colIsDeleted) {
    // update dual activity; relaxed activity bound invalidates derived implied
    // bounds
    if (!impliedDualRowBounds.implBndsValidAfterCoefChange(col, row, oldCoef,
                                                           newCoef))
      resetRowDualImpliedBoundsDerivedFromCol(col);
    // not relaxed; only reset if sourced from this column
    else if (rowDualLowerSource[row] == col || rowDualUpperSource[row] == col)
      resetRowDualImpliedBounds(row, col);
  } else {
    // col is deleted; unconditionally reset derived implied bounds
    resetRowDualImpliedBoundsDerivedFromCol(col);
  }

  if (!rowIsDeleted) {
    // update primal activity; relaxed activity bound invalidates derived
    // implied bounds
    if (!impliedRowBounds.implBndsValidAfterCoefChange(row, col, oldCoef,
                                                       newCoef))
      resetColImpliedBoundsDerivedFromRow(row);
    // not relaxed; only reset if sourced from this row
    else if (colLowerSource[col] == row || colUpperSource[col] == row)
      resetColImpliedBounds(col, row);
  } else {
    // row is deleted; unconditionally reset derived implied bounds
    resetColImpliedBoundsDerivedFromRow(row);
  }

  // modifications to row invalidate lifting opportunities
  clearLiftingOpportunities(row);

  // inform implication table about row modification
  if (mipsolver != nullptr) mipsolver->mipdata_->implications.rowModified(row);
}

void HPresolve::changeRowLower(HighsInt row, double newLower,
                               bool skipRowDualUpdate) {
  double oldLower = model->row_lower_[row];
  if (oldLower == newLower) return;
  model->row_lower_[row] = newLower;

  if (!skipRowDualUpdate) {
    if (oldLower == -kHighsInf && newLower != -kHighsInf) {
      // row gained a finite lower bound: dual upper bound loosens from 0
      // to kHighsInf, so implied dual bounds derived through columns in
      // this row may be stale (too tight) and must be reset
      changeRowDualUpper(row, kHighsInf);
      for (const HighsSliceNonzero& nz : getRowVector(row))
        resetRowDualImpliedBoundsDerivedFromCol(nz.index());
    } else if (oldLower != -kHighsInf && newLower == -kHighsInf) {
      // row lost its finite lower bound: dual upper bound tightens from
      // kHighsInf to 0; existing implied dual bounds remain valid and new
      // tightening opportunities are picked up by updateRowDualImpliedBounds
      changeRowDualUpper(row, 0.0);
    }
  }

  // tightening preserves validity of column implied bounds; only reset on
  // loosening
  if (newLower < oldLower) resetColImpliedBoundsDerivedFromRow(row);
  if (mipsolver != nullptr) mipsolver->mipdata_->implications.rowModified(row);
  markChangedRow(row);
}

void HPresolve::changeRowUpper(HighsInt row, double newUpper,
                               bool skipRowDualUpdate) {
  double oldUpper = model->row_upper_[row];
  if (oldUpper == newUpper) return;
  model->row_upper_[row] = newUpper;

  if (!skipRowDualUpdate) {
    if (oldUpper == kHighsInf && newUpper != kHighsInf) {
      // row gained a finite upper bound: dual lower bound loosens from 0
      // to -kHighsInf, so implied dual bounds derived through columns in
      // this row may be stale (too tight) and must be reset
      changeRowDualLower(row, -kHighsInf);
      for (const HighsSliceNonzero& nz : getRowVector(row))
        resetRowDualImpliedBoundsDerivedFromCol(nz.index());
    } else if (oldUpper != kHighsInf && newUpper == kHighsInf) {
      // row lost its finite upper bound: dual lower bound tightens from
      // -kHighsInf to 0; existing implied dual bounds remain valid and new
      // tightening opportunities are picked up by updateRowDualImpliedBounds
      changeRowDualLower(row, 0.0);
    }
  }

  // tightening preserves validity of column implied bounds; only reset on
  // loosening
  if (newUpper > oldUpper) resetColImpliedBoundsDerivedFromRow(row);
  if (mipsolver != nullptr) mipsolver->mipdata_->implications.rowModified(row);
  markChangedRow(row);
}

void HPresolve::addToRowLower(HighsInt row, const HighsCDouble& delta) {
  if (model->row_lower_[row] == -kHighsInf) return;
  changeRowLower(row, static_cast<double>(model->row_lower_[row] + delta));
}

void HPresolve::addToRowUpper(HighsInt row, const HighsCDouble& delta) {
  if (model->row_upper_[row] == kHighsInf) return;
  changeRowUpper(row, static_cast<double>(model->row_upper_[row] + delta));
}

HighsInt HPresolve::findNonzero(HighsInt row, HighsInt col) {
  if (rowroot[row] == -1) return -1;

  auto get_row_left = [&](HighsInt pos) -> HighsInt& { return ARleft[pos]; };
  auto get_row_right = [&](HighsInt pos) -> HighsInt& { return ARright[pos]; };
  auto get_row_key = [&](HighsInt pos) { return Acol[pos]; };
  rowroot[row] =
      highs_splay(col, rowroot[row], get_row_left, get_row_right, get_row_key);

  if (Acol[rowroot[row]] == col) return rowroot[row];

  return -1;
}

void HPresolve::shrinkProblem(HighsPostsolveStack& postsolve_stack) {
  //  printf("HPresolve::shrinkProblem: numDeletedCols = %d; numDeletedRows =
  //  %d\n",
  //         int(numDeletedCols), int(numDeletedRows));
  // The final call to shrinkProblem, or if it's called on return from
  //  if (numDeletedCols == 0 && numDeletedRows == 0) return;
  HighsInt oldNumCol = model->num_col_;
  HighsInt oldNumRow = model->num_row_;
  // If HPresolve::shrinkProblem has been called before setting up the
  // full presolve data structures - implying that presolve has
  // terminated in HPresolveInitialSweep::run, when the model is
  // up-to-date, so no shrinkage is required
  if (!hasPresolveDataStructures()) return;
  assert(colDeleted.size() == static_cast<size_t>(oldNumCol));
  assert(rowDeleted.size() == static_cast<size_t>(oldNumRow));
  model->num_col_ = 0;
  model->num_row_ = 0;
  std::vector<HighsInt> newColIndex(oldNumCol);
  std::vector<HighsInt> newRowIndex(oldNumRow);
  const bool have_col_names = model->col_names_.size() > 0;
  const bool have_row_names = model->row_names_.size() > 0;
  assert(!have_col_names ||
         model->col_names_.size() == static_cast<size_t>(oldNumCol));
  assert(!have_row_names ||
         model->row_names_.size() == static_cast<size_t>(oldNumRow));
  // Shrink the col data
  for (HighsInt i = 0; i != oldNumCol; ++i) {
    if (colDeleted[i])
      newColIndex[i] = -1;
    else {
      newColIndex[i] = model->num_col_++;
      if (newColIndex[i] < i) {
        model->col_cost_[newColIndex[i]] = model->col_cost_[i];
        model->col_lower_[newColIndex[i]] = model->col_lower_[i];
        model->col_upper_[newColIndex[i]] = model->col_upper_[i];
        assert(!std::isnan(model->col_lower_[newColIndex[i]]));
        assert(!std::isnan(model->col_upper_[newColIndex[i]]));
        model->integrality_[newColIndex[i]] = model->integrality_[i];
        implColLower[newColIndex[i]] = implColLower[i];
        implColUpper[newColIndex[i]] = implColUpper[i];
        colLowerSource[newColIndex[i]] = colLowerSource[i];
        colUpperSource[newColIndex[i]] = colUpperSource[i];
        implRowDualSourceByCol[newColIndex[i]] = implRowDualSourceByCol[i];
        colhead[newColIndex[i]] = colhead[i];
        colsize[newColIndex[i]] = colsize[i];
        if (have_col_names)
          model->col_names_[newColIndex[i]] = std::move(model->col_names_[i]);
        changedColFlag[newColIndex[i]] = changedColFlag[i];
      }
    }
  }
  if (model->fme_obj_col_ >= 0)
    model->fme_obj_col_ = newColIndex[model->fme_obj_col_];
  colDeleted.assign(model->num_col_, false);
  model->col_cost_.resize(model->num_col_);
  model->col_lower_.resize(model->num_col_);
  model->col_upper_.resize(model->num_col_);
  model->integrality_.resize(model->num_col_);
  implColLower.resize(model->num_col_);
  implColUpper.resize(model->num_col_);
  colLowerSource.resize(model->num_col_);
  colUpperSource.resize(model->num_col_);
  implRowDualSourceByCol.resize(model->num_col_);
  colhead.resize(model->num_col_);
  colsize.resize(model->num_col_);
  if (have_col_names) model->col_names_.resize(model->num_col_);
  changedColFlag.resize(model->num_col_);
  numDeletedCols = 0;
  // Shrink the row data
  for (HighsInt i = 0; i != oldNumRow; ++i) {
    if (rowDeleted[i])
      newRowIndex[i] = -1;
    else {
      newRowIndex[i] = model->num_row_++;
      if (newRowIndex[i] < i) {
        model->row_lower_[newRowIndex[i]] = model->row_lower_[i];
        model->row_upper_[newRowIndex[i]] = model->row_upper_[i];
        assert(!std::isnan(model->row_lower_[newRowIndex[i]]));
        assert(!std::isnan(model->row_upper_[newRowIndex[i]]));
        rowDualLower[newRowIndex[i]] = rowDualLower[i];
        rowDualUpper[newRowIndex[i]] = rowDualUpper[i];
        implRowDualLower[newRowIndex[i]] = implRowDualLower[i];
        implRowDualUpper[newRowIndex[i]] = implRowDualUpper[i];
        rowDualLowerSource[newRowIndex[i]] = rowDualLowerSource[i];
        rowDualUpperSource[newRowIndex[i]] = rowDualUpperSource[i];
        colImplSourceByRow[newRowIndex[i]] = colImplSourceByRow[i];
        rowroot[newRowIndex[i]] = rowroot[i];
        rowsize[newRowIndex[i]] = rowsize[i];
        rowsizeInteger[newRowIndex[i]] = rowsizeInteger[i];
        rowsizeImplInt[newRowIndex[i]] = rowsizeImplInt[i];
        if (have_row_names)
          model->row_names_[newRowIndex[i]] = std::move(model->row_names_[i]);
        changedRowFlag[newRowIndex[i]] = changedRowFlag[i];
        singleEquationChecked[newRowIndex[i]] = singleEquationChecked[i];
      }
    }
  }
  for (HighsInt i = 0; i != model->num_col_; ++i) {
    if (colLowerSource[i] != -1)
      colLowerSource[i] = newRowIndex[colLowerSource[i]];
    if (colUpperSource[i] != -1)
      colUpperSource[i] = newRowIndex[colUpperSource[i]];
  }

  for (HighsInt i = 0; i != model->num_row_; ++i) {
    if (rowDualLowerSource[i] != -1)
      rowDualLowerSource[i] = newColIndex[rowDualLowerSource[i]];
    if (rowDualUpperSource[i] != -1)
      rowDualUpperSource[i] = newColIndex[rowDualUpperSource[i]];
  }

  for (HighsInt i = 0; i != model->num_col_; ++i) {
    std::set<HighsInt> newSet;
    std::for_each(implRowDualSourceByCol[i].cbegin(),
                  implRowDualSourceByCol[i].cend(), [&](const HighsInt& row) {
                    if (newRowIndex[row] != -1)
                      newSet.emplace(newRowIndex[row]);
                  });
    implRowDualSourceByCol[i] = std::move(newSet);
  }

  for (HighsInt i = 0; i != model->num_row_; ++i) {
    std::set<HighsInt> newSet;
    std::for_each(colImplSourceByRow[i].cbegin(), colImplSourceByRow[i].cend(),
                  [&](const HighsInt& col) {
                    if (newColIndex[col] != -1)
                      newSet.emplace(newColIndex[col]);
                  });
    colImplSourceByRow[i] = std::move(newSet);
  }
  rowDeleted.assign(model->num_row_, false);
  model->row_lower_.resize(model->num_row_);
  model->row_upper_.resize(model->num_row_);
  rowDualLower.resize(model->num_row_);
  rowDualUpper.resize(model->num_row_);
  implRowDualLower.resize(model->num_row_);
  implRowDualUpper.resize(model->num_row_);
  rowDualLowerSource.resize(model->num_row_);
  rowDualUpperSource.resize(model->num_row_);
  colImplSourceByRow.resize(model->num_row_);
  rowroot.resize(model->num_row_);
  rowsize.resize(model->num_row_);
  rowsizeInteger.resize(model->num_row_);
  rowsizeImplInt.resize(model->num_row_);
  if (have_row_names) model->row_names_.resize(model->num_row_);
  changedRowFlag.resize(model->num_row_);
  singleEquationChecked.resize(model->num_row_);

  numDeletedRows = 0;
  postsolve_stack.compressIndexMaps(newRowIndex, newColIndex);
  impliedRowBounds.shrink(newRowIndex, model->num_row_);
  impliedDualRowBounds.shrink(newColIndex, model->num_col_);

  for (size_t i = 0; i != Avalue.size(); ++i) {
    if (Avalue[i] == 0) continue;
    assert(newColIndex[Acol[i]] != -1);
    assert(newRowIndex[Arow[i]] != -1);
    Acol[i] = newColIndex[Acol[i]];
    Arow[i] = newRowIndex[Arow[i]];
  }

  // update index sets
  for (HighsInt& singCol : singletonColumns) singCol = newColIndex[singCol];
  singletonColumns.erase(
      std::remove(singletonColumns.begin(), singletonColumns.end(), -1),
      singletonColumns.end());

  for (HighsInt& chgCol : changedColIndices) chgCol = newColIndex[chgCol];
  changedColIndices.erase(
      std::remove(changedColIndices.begin(), changedColIndices.end(), -1),
      changedColIndices.end());

  for (HighsInt& singRow : singletonRows) singRow = newRowIndex[singRow];
  singletonRows.erase(
      std::remove(singletonRows.begin(), singletonRows.end(), -1),
      singletonRows.end());

  for (HighsInt& chgRow : changedRowIndices) chgRow = newRowIndex[chgRow];
  changedRowIndices.erase(
      std::remove(changedRowIndices.begin(), changedRowIndices.end(), -1),
      changedRowIndices.end());

  for (auto& rowColPair : substitutionOpportunities) {
    // skip deleted elements
    if (rowColPair.first == -1) continue;
    rowColPair.first = newRowIndex[rowColPair.first];
    rowColPair.second = newColIndex[rowColPair.second];
  }
  substitutionOpportunities.erase(
      std::remove_if(substitutionOpportunities.begin(),
                     substitutionOpportunities.end(),
                     [&](const std::pair<HighsInt, HighsInt>& p) {
                       return p.first == -1 || p.second == -1;
                     }),
      substitutionOpportunities.end());

  // todo remove equation set and replace with a vector of doubleton eqs
  equations.clear();
  eqiters.assign(model->num_row_, equations.end());
  for (HighsInt i = 0; i != model->num_row_; ++i) {
    if (isEquation(i)) eqiters[i] = equations.emplace(rowsize[i], i).first;
  }

  if (mipsolver != nullptr) {
    mipsolver->mipdata_->rowMatrixSet = false;
    mipsolver->mipdata_->objectiveFunction = HighsObjectiveFunction(*mipsolver);
    mipsolver->mipdata_->getDomain() = HighsDomain(*mipsolver);
    mipsolver->mipdata_->cliquetable.rebuild(model->num_col_, postsolve_stack,
                                             mipsolver->mipdata_->getDomain(),
                                             newColIndex, newRowIndex);
    presolveCliqueTable.rebuild(mipsolver->mipdata_->cliquetable);
    mipsolver->mipdata_->implications.rebuild(model->num_col_, newColIndex,
                                              newRowIndex);
    mipsolver->mipdata_->getCutPool() =
        HighsCutPool(mipsolver->model_->num_col_,
                     mipsolver->options_mip_->mip_pool_age_limit,
                     mipsolver->options_mip_->mip_pool_soft_limit, 0);
    mipsolver->mipdata_->getConflictPool() =
        HighsConflictPool(5 * mipsolver->options_mip_->mip_pool_age_limit,
                          mipsolver->options_mip_->mip_pool_soft_limit);

    for (HighsInt i = 0; i != oldNumCol; ++i)
      if (newColIndex[i] != -1) numProbes[newColIndex[i]] = numProbes[i];

    mipsolver->mipdata_->debugSolution.shrink(newColIndex);
    numProbes.resize(model->num_col_);
  }
  // Need to set the constraint matrix dimensions
  model->setMatrixDimensions();
  // Need to reset current number of deleted rows and columns in logging
  analysis_.resetNumDeleted();
}

void HPresolve::addToMatrix(const HighsInt row, const HighsInt col,
                            const double val) {
  HighsInt pos = findNonzero(row, col);

  markChangedRow(row);
  markChangedCol(col);

  if (pos == -1) {
    if (freeslots.empty()) {
      pos = static_cast<HighsInt>(Avalue.size());
      Avalue.push_back(val);
      Arow.push_back(row);
      Acol.push_back(col);
      Anext.push_back(-1);
      Aprev.push_back(-1);
      ARleft.push_back(-1);
      ARright.push_back(-1);
    } else {
      pos = freeslots.back();
      freeslots.pop_back();
      Avalue[pos] = val;
      Arow[pos] = row;
      Acol[pos] = col;
      Aprev[pos] = -1;
    }

    link(pos);
    matrixNonZeroChanged(row, col, 0.0, val, false, false);

  } else {
    double sum = Avalue[pos] + val;
    if (std::abs(sum) <= options->small_matrix_value) {
      unlink(pos);
    } else {
      matrixNonZeroChanged(row, col, Avalue[pos], sum, false, false);
      Avalue[pos] = sum;
    }
  }
}

bool HPresolve::addToMatrix(
    HighsPostsolveStack& postsolve_stack, const std::vector<double>& row_lower,
    const std::vector<double>& row_upper,
    const std::vector<std::vector<HighsInt>>& row_indices,
    const std::vector<std::vector<double>>& row_values) {
  // update number of rows
  HighsInt num_rows = static_cast<HighsInt>(row_indices.size());
  if (num_rows == 0) return true;
  HighsInt oldNumRows = model->num_row_;
  model->num_row_ += num_rows;
  model->a_matrix_.num_row_ += num_rows;

  // resize postsolve vectors
  postsolve_stack.appendRowsToModel(num_rows);

  // add row bounds
  model->row_lower_.insert(model->row_lower_.end(), row_lower.begin(),
                           row_lower.end());
  model->row_upper_.insert(model->row_upper_.end(), row_upper.begin(),
                           row_upper.end());

  // initialise row sizes
  if (!okResize(rowroot, model->num_row_, HighsInt{-1})) return false;
  if (!okResize(rowsize, model->num_row_, HighsInt{0})) return false;
  if (!okResize(rowsizeInteger, model->num_row_, HighsInt{0})) return false;
  if (!okResize(rowsizeImplInt, model->num_row_, HighsInt{0})) return false;

  // initialise row duals
  if (!okResize(rowDualLower, model->num_row_, -kHighsInf)) return false;
  if (!okResize(rowDualUpper, model->num_row_, kHighsInf)) return false;
  for (HighsInt i = oldNumRows; i < model->num_row_; i++) {
    if (model->row_lower_[i] == -kHighsInf) rowDualUpper[i] = 0;
    if (model->row_upper_[i] == kHighsInf) rowDualLower[i] = 0;
  }

  // initialise implied row duals
  if (!okResize(implRowDualLower, model->num_row_, -kHighsInf)) return false;
  if (!okResize(implRowDualUpper, model->num_row_, kHighsInf)) return false;
  if (!okResize(rowDualLowerSource, model->num_row_, HighsInt{-1}))
    return false;
  if (!okResize(rowDualUpperSource, model->num_row_, HighsInt{-1}))
    return false;
  if (!okResize(colImplSourceByRow, model->num_row_, std::set<HighsInt>{}))
    return false;

  // initialise flags
  if (!okResize(changedRowFlag, model->num_row_, uint8_t{0})) return false;
  if (!okResize(rowDeleted, model->num_row_, uint8_t{0})) return false;
  if (!okResize(singleEquationChecked, model->num_row_, uint8_t{0}))
    return false;

  // initialise row names
  if (!okResize(model->row_names_, model->num_row_, std::string{}))
    return false;

  // resize vector for equations
  if (!okResize(eqiters, model->num_row_, equations.end())) return false;

  // resize vectors for implied row bounds
  impliedRowBounds.setNumSums(model->num_row_);

  // set bound arrays again (pointers may get invalidated by reallocation)
  impliedDualRowBounds.setBoundArrays(
      rowDualLower.data(), rowDualUpper.data(), implRowDualLower.data(),
      implRowDualUpper.data(), rowDualLowerSource.data(),
      rowDualUpperSource.data());

  if (mipsolver != nullptr)
    mipsolver->mipdata_->implications.resize(mipsolver->numCol(),
                                             mipsolver->numRow());

  for (HighsInt i = 0; i < num_rows; i++) {
    // new row index
    HighsInt row = oldNumRows + i;

    // add non-zeros
    for (size_t j = 0; j < row_indices[i].size(); j++)
      addToMatrix(row, row_indices[i][j], row_values[i][j]);

    // add row singleton
    if (rowsize[row] == 1) singletonRows.push_back(row);

    // add equation
    if (isEquation(row))
      eqiters[row] = equations.emplace(rowsize[row], row).first;
  }

  return true;
}

bool HPresolve::addToMatrix(HighsPostsolveStack& postsolve_stack,
                            double row_lower, double row_upper,
                            const std::vector<HighsInt>& row_indices,
                            const std::vector<double>& row_values) {
  return addToMatrix(postsolve_stack, std::vector<double>{row_lower},
                     std::vector<double>{row_upper},
                     std::vector<std::vector<HighsInt>>{row_indices},
                     std::vector<std::vector<double>>{row_values});
}

HighsTripletListSlice HPresolve::getColumnVector(HighsInt col) const {
  return HighsTripletListSlice(Arow.data(), Avalue.data(), Anext.data(),
                               colhead[col]);
}

HighsTripletTreeSlicePreOrder HPresolve::getRowVector(HighsInt row) const {
  return HighsTripletTreeSlicePreOrder(
      Acol.data(), Avalue.data(), ARleft.data(), ARright.data(), rowroot[row]);
}

HighsTripletTreeSliceInOrder HPresolve::getSortedRowVector(HighsInt row) const {
  return HighsTripletTreeSliceInOrder(Acol.data(), Avalue.data(), ARleft.data(),
                                      ARright.data(), rowroot[row]);
}

void HPresolve::markRowDeleted(HighsInt row) {
  assert(!analysis_.logging_on_);
  assert(!rowDeleted[row]);

  // remove equations from set of equations
  if (isEquation(row) && eqiters[row] != equations.end()) {
    equations.erase(eqiters[row]);
    eqiters[row] = equations.end();
  }

  // prevents row from being added to change vector
  changedRowFlag[row] = true;
  rowDeleted[row] = true;
  ++numDeletedRows;
  resetSingleEquationChecked(row);
}

void HPresolve::markColDeleted(HighsInt col) {
  assert(!analysis_.logging_on_);
  assert(!colDeleted[col]);

  // prevents col from being added to change vector
  changedColFlag[col] = true;
  colDeleted[col] = true;
  ++numDeletedCols;
  if (col == model->fme_obj_col_) model->fme_obj_col_ = -1;
  if (mipsolver != nullptr && mipsolver->mipdata_->cliquesExtracted) {
    presolveCliqueTable.eliminateCol(col);
  }
}

HPresolve::Result HPresolve::changeColUpper(HighsInt col, double newUpper) {
  if (model->integrality_[col] != HighsVarType::kContinuous) {
    newUpper = std::floor(newUpper + primal_feastol);
    if (newUpper == model->col_upper_[col]) return Result::kOk;
  }

  double oldUpper = model->col_upper_[col];
  model->col_upper_[col] = newUpper;

  HPRESOLVE_CHECKED_CALL(checkColBounds(col));

  for (const HighsSliceNonzero& nonzero : getColumnVector(col)) {
    impliedRowBounds.updatedVarUpper(nonzero.index(), col, nonzero.value(),
                                     oldUpper);
    markChangedRow(nonzero.index());
  }
  if (mipsolver != nullptr &&
      model->col_lower_[col] == model->col_upper_[col]) {
    HPRESOLVE_CHECKED_CALL(
        updateCliqueTableFixedCol(col, model->col_lower_[col]));
  }
  return Result::kOk;
}

HPresolve::Result HPresolve::changeColLower(HighsInt col, double newLower) {
  if (model->integrality_[col] != HighsVarType::kContinuous) {
    newLower = std::ceil(newLower - primal_feastol);
    if (newLower == model->col_lower_[col]) return Result::kOk;
  }

  double oldLower = model->col_lower_[col];
  model->col_lower_[col] = newLower;

  HPRESOLVE_CHECKED_CALL(checkColBounds(col));

  for (const HighsSliceNonzero& nonzero : getColumnVector(col)) {
    impliedRowBounds.updatedVarLower(nonzero.index(), col, nonzero.value(),
                                     oldLower);
    markChangedRow(nonzero.index());
  }
  if (mipsolver != nullptr &&
      model->col_lower_[col] == model->col_upper_[col]) {
    HPRESOLVE_CHECKED_CALL(
        updateCliqueTableFixedCol(col, model->col_lower_[col]));
  }
  return Result::kOk;
}

HPresolve::Result HPresolve::changeColBounds(HighsInt col, double newLower,
                                             double newUpper) {
  if (newLower > model->col_upper_[col]) {
    // change upper bound first
    HPRESOLVE_CHECKED_CALL(changeColUpper(col, newUpper));
    HPRESOLVE_CHECKED_CALL(changeColLower(col, newLower));
  } else {
    // change lower bound first
    HPRESOLVE_CHECKED_CALL(changeColLower(col, newLower));
    HPRESOLVE_CHECKED_CALL(changeColUpper(col, newUpper));
  }
  return Result::kOk;
}

HPresolve::Result HPresolve::checkColBounds(HighsInt col, bool* isFixed) {
  assert(!colDeleted[col]);
  double boundDiff = model->col_upper_[col] - model->col_lower_[col];
  if (isFixed != nullptr) *isFixed = false;
  if (boundDiff <= primal_feastol &&
      (boundDiff <= options->small_matrix_value ||
       getMaxAbsColVal(col) * boundDiff <= primal_feastol)) {
    // check for primal infeasibility
    if (boundDiff < -primal_feastol) return Result::kPrimalInfeasible;
    // check for unboundedness
    if (std::abs(model->col_lower_[col]) == kHighsInf)
      return Result::kDualInfeasible;
    // column is fixed
    if (isFixed != nullptr) *isFixed = true;
  }
  return Result::kOk;
}

void HPresolve::changeRowDualUpper(HighsInt row, double newUpper) {
  double oldUpper = rowDualUpper[row];
  rowDualUpper[row] = newUpper;

  // printf("tightening upper bound of column %" HIGHSINT_FORMAT " from %.15g to
  // %.15g\n", col,
  //        oldUpper, newUpper);
  for (const HighsSliceNonzero& nonzero : getRowVector(row)) {
    impliedDualRowBounds.updatedVarUpper(nonzero.index(), row, nonzero.value(),
                                         oldUpper);
    markChangedCol(nonzero.index());
  }
}

void HPresolve::changeRowDualLower(HighsInt row, double newLower) {
  double oldLower = rowDualLower[row];
  rowDualLower[row] = newLower;
  // printf("tightening lower bound of column %" HIGHSINT_FORMAT " from %.15g to
  // %.15g\n", col,
  //        oldLower, newLower);

  for (const HighsSliceNonzero& nonzero : getRowVector(row)) {
    impliedDualRowBounds.updatedVarLower(nonzero.index(), row, nonzero.value(),
                                         oldLower);
    markChangedCol(nonzero.index());
  }
}

void HPresolve::changeImplColUpper(HighsInt col, double newUpper,
                                   HighsInt originRow) {
  double oldImplUpper = implColUpper[col];
  HighsInt oldUpperSource = colUpperSource[col];
  if (oldImplUpper >= model->col_upper_[col] - primal_feastol &&
      newUpper < model->col_upper_[col] - primal_feastol) {
    // the dual constraint can be considered a >= constraint and was free, or a
    // <= constraint before
    markChangedCol(col);
  }
  bool newImpliedFree =
      isLowerImplied(col) &&
      oldImplUpper > model->col_upper_[col] + primal_feastol &&
      newUpper <= model->col_upper_[col] + primal_feastol;

  // remember the source of this upper bound, so that we can correctly identify
  // weak domination
  if (oldUpperSource != originRow) {
    if (oldUpperSource != -1 && oldUpperSource != colLowerSource[col])
      colImplSourceByRow[oldUpperSource].erase(col);
    if (originRow != -1) colImplSourceByRow[originRow].emplace(col);

    colUpperSource[col] = originRow;
  }

  // update implied bound
  implColUpper[col] = newUpper;

  // if the old and the new implied bound are not better than the upper bound,
  // nothing needs to be updated
  if (!newImpliedFree &&
      std::min(oldImplUpper, newUpper) >= model->col_upper_[col])
    return;

  for (const HighsSliceNonzero& nonzero : getColumnVector(col)) {
    impliedRowBounds.updatedImplVarUpper(nonzero.index(), col, nonzero.value(),
                                         oldImplUpper, oldUpperSource);
    if (newImpliedFree && isDualImpliedFree(nonzero.index()))
      substitutionOpportunities.emplace_back(nonzero.index(), col);

    markChangedRow(nonzero.index());
  }
}

void HPresolve::changeImplColLower(HighsInt col, double newLower,
                                   HighsInt originRow) {
  double oldImplLower = implColLower[col];
  HighsInt oldLowerSource = colLowerSource[col];
  if (oldImplLower <= model->col_lower_[col] + primal_feastol &&
      newLower > model->col_lower_[col] + primal_feastol) {
    // the dual constraint can additionally be considered a <= constraint and
    // was free, or a >= constraint before
    markChangedCol(col);
  }
  bool newImpliedFree =
      isUpperImplied(col) &&
      oldImplLower < model->col_lower_[col] - primal_feastol &&
      newLower >= model->col_lower_[col] - primal_feastol;

  // remember the source of this lower bound, so that we can correctly identify
  // weak domination
  if (oldLowerSource != originRow) {
    if (oldLowerSource != -1 && oldLowerSource != colUpperSource[col])
      colImplSourceByRow[oldLowerSource].erase(col);
    if (originRow != -1) colImplSourceByRow[originRow].emplace(col);

    colLowerSource[col] = originRow;
  }

  // update implied bound
  implColLower[col] = newLower;

  // if the old and the new implied bound are not better than the lower bound,
  // nothing needs to be updated
  if (!newImpliedFree &&
      std::max(oldImplLower, newLower) <= model->col_lower_[col])
    return;

  for (const HighsSliceNonzero& nonzero : getColumnVector(col)) {
    impliedRowBounds.updatedImplVarLower(nonzero.index(), col, nonzero.value(),
                                         oldImplLower, oldLowerSource);
    if (newImpliedFree && isDualImpliedFree(nonzero.index()))
      substitutionOpportunities.emplace_back(nonzero.index(), col);

    markChangedRow(nonzero.index());
  }
}

void HPresolve::changeImplRowDualUpper(HighsInt row, double newUpper,
                                       HighsInt originCol) {
  double oldImplUpper = implRowDualUpper[row];
  HighsInt oldUpperSource = rowDualUpperSource[row];

  if (oldImplUpper >= -options->dual_feasibility_tolerance &&
      newUpper < -options->dual_feasibility_tolerance)
    markChangedRow(row);

  bool newDualImplied =
      !isDualImpliedFree(row) &&
      oldImplUpper > rowDualUpper[row] + options->dual_feasibility_tolerance &&
      newUpper <= rowDualUpper[row] + options->dual_feasibility_tolerance;

  // remember the source of this upper bound, so that we can correctly identify
  // weak domination
  if (oldUpperSource != originCol) {
    if (oldUpperSource != -1 && oldUpperSource != rowDualLowerSource[row])
      implRowDualSourceByCol[oldUpperSource].erase(row);
    if (originCol != -1) implRowDualSourceByCol[originCol].emplace(row);

    rowDualUpperSource[row] = originCol;
  }

  // update implied bound
  implRowDualUpper[row] = newUpper;

  // nothing needs to be updated
  if (!newDualImplied && std::min(oldImplUpper, newUpper) >= rowDualUpper[row])
    return;

  for (const HighsSliceNonzero& nonzero : getRowVector(row)) {
    impliedDualRowBounds.updatedImplVarUpper(
        nonzero.index(), row, nonzero.value(), oldImplUpper, oldUpperSource);
    markChangedCol(nonzero.index());

    if (newDualImplied && isImpliedFree(nonzero.index()))
      substitutionOpportunities.emplace_back(row, nonzero.index());
  }
}

void HPresolve::changeImplRowDualLower(HighsInt row, double newLower,
                                       HighsInt originCol) {
  double oldImplLower = implRowDualLower[row];
  HighsInt oldLowerSource = rowDualLowerSource[row];

  if (oldImplLower <= options->dual_feasibility_tolerance &&
      newLower > options->dual_feasibility_tolerance)
    markChangedRow(row);

  bool newDualImplied =
      !isDualImpliedFree(row) &&
      oldImplLower < rowDualLower[row] - options->dual_feasibility_tolerance &&
      newLower >= rowDualLower[row] - options->dual_feasibility_tolerance;

  // remember the source of this lower bound, so that we can correctly identify
  // weak domination
  if (oldLowerSource != originCol) {
    if (oldLowerSource != -1 && oldLowerSource != rowDualUpperSource[row])
      implRowDualSourceByCol[oldLowerSource].erase(row);
    if (originCol != -1) implRowDualSourceByCol[originCol].emplace(row);

    rowDualLowerSource[row] = originCol;
  }

  // update implied bound
  implRowDualLower[row] = newLower;

  // nothing needs to be updated
  if (!newDualImplied && std::max(oldImplLower, newLower) <= rowDualLower[row])
    return;

  for (const HighsSliceNonzero& nonzero : getRowVector(row)) {
    impliedDualRowBounds.updatedImplVarLower(
        nonzero.index(), row, nonzero.value(), oldImplLower, oldLowerSource);
    markChangedCol(nonzero.index());

    if (newDualImplied && isImpliedFree(nonzero.index()))
      substitutionOpportunities.emplace_back(row, nonzero.index());
  }
}

HPresolve::Result HPresolve::scaleMIP(HighsPostsolveStack& postsolve_stack) {
  for (HighsInt i = 0; i < model->num_row_; ++i) {
    if (rowDeleted[i] || rowsize[i] < 1 ||
        rowsizeInteger[i] + rowsizeImplInt[i] == rowsize[i])
      continue;

    storeRow(i);

    double maxAbsVal = 0.0;

    for (size_t j = 0; j < rowpositions.size(); ++j) {
      HighsInt nzPos = rowpositions[j];
      if (model->integrality_[Acol[nzPos]] != HighsVarType::kContinuous)
        continue;

      maxAbsVal = std::max(std::abs(Avalue[nzPos]), maxAbsVal);
    }

    assert(maxAbsVal != 0.0);

    double scale = std::exp2(std::round(-std::log2(maxAbsVal)));
    if (scale == 1.0) continue;

    if (model->row_upper_[i] == kHighsInf) scale = -scale;

    scaleStoredRow(i, scale);
  }

  for (HighsInt i = 0; i < model->num_col_; ++i) {
    if (colDeleted[i] || colsize[i] < 1 ||
        model->integrality_[i] != HighsVarType::kContinuous)
      continue;

    double maxAbsVal = 0;

    for (const HighsSliceNonzero& nonz : getColumnVector(i)) {
      maxAbsVal = std::max(std::abs(nonz.value()), maxAbsVal);
    }

    double scale = std::exp2(std::round(-std::log2(maxAbsVal)));
    if (scale == 1.0) continue;

    HPRESOLVE_CHECKED_CALL(transformColumn(postsolve_stack, i, scale, 0.0));
  }
  return Result::kOk;
}

void HPresolve::getRowPositions(HighsInt row,
                                std::vector<HighsInt>& myrowpositions) const {
  myrowpositions.clear();

  auto rowvector = getSortedRowVector(row);
  for (auto rowiter = rowvector.begin(); rowiter != rowvector.end(); ++rowiter)
    myrowpositions.push_back(rowiter.position());
}

void HPresolve::storeRow(HighsInt row) { getRowPositions(row, rowpositions); }

HighsTripletPositionSlice HPresolve::getStoredRow() const {
  return HighsTripletPositionSlice(Acol.data(), Avalue.data(),
                                   rowpositions.data(), rowpositions.size());
}

bool HPresolve::okFromCSC(const std::vector<double>& Aval,
                          const std::vector<HighsInt>& Aindex,
                          const std::vector<HighsInt>& Astart) {
  Avalue.clear();
  Acol.clear();
  Arow.clear();

  freeslots.clear();
  if (!okAssign(colhead, model->num_col_, HighsInt{-1})) return false;
  if (!okAssign(rowroot, model->num_row_, HighsInt{-1})) return false;
  if (!okAssign(colsize, model->num_col_)) return false;
  if (!okAssign(rowsize, model->num_row_)) return false;
  if (!okAssign(rowsizeInteger, model->num_row_)) return false;
  if (!okAssign(rowsizeImplInt, model->num_row_)) return false;

  impliedRowBounds.setNumSums(0);
  impliedDualRowBounds.setNumSums(0);
  impliedRowBounds.setBoundArrays(
      model->col_lower_.data(), model->col_upper_.data(), implColLower.data(),
      implColUpper.data(), colLowerSource.data(), colUpperSource.data());
  impliedRowBounds.setNumSums(model->num_row_);
  impliedDualRowBounds.setBoundArrays(
      rowDualLower.data(), rowDualUpper.data(), implRowDualLower.data(),
      implRowDualUpper.data(), rowDualLowerSource.data(),
      rowDualUpperSource.data());
  impliedDualRowBounds.setNumSums(model->num_col_);

  HighsInt ncol = static_cast<HighsInt>(Astart.size()) - 1;
  assert(static_cast<size_t>(ncol) == colhead.size());
  HighsInt nnz = static_cast<HighsInt>(Aval.size());

  Avalue = Aval;
  if (!okReserve(Acol, nnz)) return false;
  if (!okReserve(Arow, nnz)) return false;

  for (HighsInt i = 0; i != ncol; ++i) {
    HighsInt collen = Astart[i + 1] - Astart[i];
    Acol.insert(Acol.end(), collen, i);
    Arow.insert(Arow.end(), Aindex.begin() + Astart[i],
                Aindex.begin() + Astart[i + 1]);
  }

  if (!okResize(Anext, nnz)) return false;
  if (!okResize(Aprev, nnz)) return false;
  if (!okResize(ARleft, nnz)) return false;
  if (!okResize(ARright, nnz)) return false;
  // matrixNonZeroChanged cannot be used here: after shrinkProblem,
  // implied bound source sets are non-empty and resets would destroy
  // preserved implied bounds
  for (HighsInt pos = 0; pos != nnz; ++pos) {
    link(pos);
    impliedRowBounds.add(Arow[pos], Acol[pos], Avalue[pos]);
    impliedDualRowBounds.add(Acol[pos], Arow[pos], Avalue[pos]);
  }

  if (equations.empty()) {
    try {
      eqiters.assign(model->num_row_, equations.end());
    } catch (const std::bad_alloc& e) {
      printf("HPresolve::okFromCSC eqiters.assign fails with %s\n", e.what());
      return false;
    }
    for (HighsInt i = 0; i != model->num_row_; ++i) {
      // register equation
      if (isEquation(i)) eqiters[i] = equations.emplace(rowsize[i], i).first;
    }
  }
  return true;
}

bool HPresolve::okFromCSR(const std::vector<double>& ARval,
                          const std::vector<HighsInt>& ARindex,
                          const std::vector<HighsInt>& ARstart) {
  Avalue.clear();
  Acol.clear();
  Arow.clear();

  freeslots.clear();
  if (!okAssign(colhead, model->num_col_, HighsInt{-1})) return false;
  if (!okAssign(rowroot, model->num_row_, HighsInt{-1})) return false;
  if (!okAssign(colsize, model->num_col_)) return false;
  if (!okAssign(rowsize, model->num_row_)) return false;
  if (!okAssign(rowsizeInteger, model->num_row_)) return false;
  if (!okAssign(rowsizeImplInt, model->num_row_)) return false;

  impliedRowBounds.setNumSums(0);
  impliedDualRowBounds.setNumSums(0);
  impliedRowBounds.setBoundArrays(
      model->col_lower_.data(), model->col_upper_.data(), implColLower.data(),
      implColUpper.data(), colLowerSource.data(), colUpperSource.data());
  impliedRowBounds.setNumSums(model->num_row_);
  impliedDualRowBounds.setBoundArrays(
      rowDualLower.data(), rowDualUpper.data(), implRowDualLower.data(),
      implRowDualUpper.data(), rowDualLowerSource.data(),
      rowDualUpperSource.data());
  impliedDualRowBounds.setNumSums(model->num_col_);

  HighsInt nrow = static_cast<HighsInt>(ARstart.size()) - 1;
  assert(static_cast<size_t>(nrow) == rowroot.size());
  HighsInt nnz = static_cast<HighsInt>(ARval.size());

  Avalue = ARval;
  if (!okReserve(Acol, nnz)) return false;
  if (!okReserve(Arow, nnz)) return false;
  //  entries.reserve(nnz);

  for (HighsInt i = 0; i != nrow; ++i) {
    Arow.insert(Arow.end(), ARstart[i + 1] - ARstart[i], i);
    Acol.insert(Acol.end(), ARindex.begin() + ARstart[i],
                ARindex.begin() + ARstart[i + 1]);
  }

  if (!okResize(Anext, nnz)) return false;
  if (!okResize(Aprev, nnz)) return false;
  if (!okResize(ARleft, nnz)) return false;
  if (!okResize(ARright, nnz)) return false;
  // matrixNonZeroChanged cannot be used here: after shrinkProblem,
  // implied bound source sets are non-empty and resets would destroy
  // preserved implied bounds
  for (HighsInt pos = 0; pos != nnz; ++pos) {
    link(pos);
    impliedRowBounds.add(Arow[pos], Acol[pos], Avalue[pos]);
    impliedDualRowBounds.add(Acol[pos], Arow[pos], Avalue[pos]);
  }

  if (equations.empty()) {
    try {
      eqiters.assign(nrow, equations.end());
    } catch (const std::bad_alloc& e) {
      printf("HPresolve::okFromCSR eqiters.assign fails with %s\n", e.what());
      return false;
    }
    for (HighsInt i = 0; i != nrow; ++i) {
      // register equation
      if (isEquation(i)) eqiters[i] = equations.emplace(rowsize[i], i).first;
    }
  }
  return true;
}

void HPresolve::reinsertEquation(HighsInt row) {
  // check if this is an equation row and it now has a different size
  if (isEquation(row) && eqiters[row] != equations.end() &&
      eqiters[row]->first != rowsize[row]) {
    // if that is the case reinsert it into the equation set that is ordered
    // by sparsity
    equations.erase(eqiters[row]);
    eqiters[row] = equations.emplace(rowsize[row], row).first;
  }
}

HPresolve::Result HPresolve::transformColumn(
    HighsPostsolveStack& postsolve_stack, HighsInt col, double scale,
    double constant) {
  // replace column x by x = scale * x' + constant
  if (mipsolver != nullptr)
    mipsolver->mipdata_->implications.columnTransformed(col, scale, constant);

  postsolve_stack.linearTransform(col, scale, constant);

  // new variable x' has the following bounds:
  // scale > 0 --> (lb - constant) / scale <= x' <= (ub - constant) / scale
  // scale < 0 --> (ub - constant) / scale <= x' <= (lb - constant) / scale
  // each matrix coefficient a will be replaced by scale * a and, therefore, the
  // contributions to the constraint activities will change as follows:
  // a * lb --> (a * scale) * (lb - constant) / scale = a * (lb - constant)
  // a * ub --> (a * scale) * (ub - constant) / scale = a * (ub - constant).
  // therefore, for now the scaling can be neglected and the bounds on
  // constraint activities can be updated using the constant term.
  if (constant != 0.0) {
    double oldLower = model->col_lower_[col];
    double oldUpper = model->col_upper_[col];
    model->col_upper_[col] -= constant;
    model->col_lower_[col] -= constant;

    for (const HighsSliceNonzero& nonzero : getColumnVector(col)) {
      impliedRowBounds.updatedVarLower(nonzero.index(), col, nonzero.value(),
                                       oldLower);
      impliedRowBounds.updatedVarUpper(nonzero.index(), col, nonzero.value(),
                                       oldUpper);
    }

    double oldImplLower = implColLower[col];
    double oldImplUpper = implColUpper[col];
    implColLower[col] -= constant;
    implColUpper[col] -= constant;

    for (const HighsSliceNonzero& nonzero : getColumnVector(col)) {
      impliedRowBounds.updatedImplVarLower(nonzero.index(), col,
                                           nonzero.value(), oldImplLower,
                                           colLowerSource[col]);
      impliedRowBounds.updatedImplVarUpper(nonzero.index(), col,
                                           nonzero.value(), oldImplUpper,
                                           colUpperSource[col]);
    }
  }

  // now apply the scaling, which does not change the contributions to the
  // implied row bounds, but requires adjusting the implied bounds of the
  // columns dual constraint
  impliedDualRowBounds.sumScaled(col, scale);

  double boundScale = 1.0 / scale;
  model->col_lower_[col] *= boundScale;
  model->col_upper_[col] *= boundScale;
  implColLower[col] *= boundScale;
  implColUpper[col] *= boundScale;

  if (scale < 0) {
    std::swap(model->col_lower_[col], model->col_upper_[col]);
    std::swap(implColLower[col], implColUpper[col]);
    std::swap(colLowerSource[col], colUpperSource[col]);
  }

  model->offset_ += model->col_cost_[col] * constant;
  model->col_cost_[col] *= scale;

  for (HighsInt coliter = colhead[col]; coliter != -1;
       coliter = Anext[coliter]) {
    double val = Avalue[coliter];
    Avalue[coliter] *= scale;
    HighsInt row = Arow[coliter];
    double rowConstant = val * constant;
    if (model->row_lower_[row] != -kHighsInf)
      model->row_lower_[row] -= rowConstant;
    if (model->row_upper_[row] != kHighsInf)
      model->row_upper_[row] -= rowConstant;
  }

  // finally, use utility methods for rounding scaled bounds of
  // integer-constrained variables and updating bounds on constraint activities
  // accordingly. this should not be done before the preceding bound updates
  // (scaling and swaps) and matrix updates. we rely on the integrality status
  // being already updated to the newly scaled column by the caller, if
  // necessary.
  if (model->integrality_[col] != HighsVarType::kContinuous)
    HPRESOLVE_CHECKED_CALL(
        changeColBounds(col, model->col_lower_[col], model->col_upper_[col]));

  markChangedCol(col);
  return Result::kOk;
}

void HPresolve::scaleRow(HighsInt row, double scale, bool integral) {
  storeRow(row);

  scaleStoredRow(row, scale, integral);
}

void HPresolve::scaleStoredRow(HighsInt row, double scale, bool integral) {
  model->row_upper_[row] *= scale;
  model->row_lower_[row] *= scale;
  implRowDualLower[row] /= scale;
  implRowDualUpper[row] /= scale;
  rowDualLower[row] /= std::copysign(1.0, scale);
  rowDualUpper[row] /= std::copysign(1.0, scale);

  for (size_t j = 0; j < rowpositions.size(); ++j) {
    Avalue[rowpositions[j]] *= scale;
    if (std::abs(Avalue[rowpositions[j]]) <= options->small_matrix_value)
      unlink(rowpositions[j]);
  }

  impliedRowBounds.sumScaled(row, scale);
  if (scale < 0) {
    std::swap(rowDualLower[row], rowDualUpper[row]);
    std::swap(implRowDualLower[row], implRowDualUpper[row]);
    std::swap(rowDualLowerSource[row], rowDualUpperSource[row]);
    std::swap(model->row_lower_[row], model->row_upper_[row]);
  }

  if (integral) {
    if (model->row_upper_[row] != kHighsInf)
      model->row_upper_[row] =
          std::floor(model->row_upper_[row] + primal_feastol);
    if (model->row_lower_[row] != -kHighsInf)
      model->row_lower_[row] =
          std::ceil(model->row_lower_[row] - primal_feastol);
  }
}

void HPresolve::toCSC(std::vector<double>& Aval, std::vector<HighsInt>& Aindex,
                      std::vector<HighsInt>& Astart) {
  // set up the column starts using the column size array
  size_t numcol = colsize.size();
  Astart.resize(numcol + 1);
  HighsInt nnz = 0;
  for (size_t i = 0; i != numcol; ++i) {
    Astart[i] = nnz;
    nnz += colsize[i];
  }
  Astart[numcol] = nnz;

  // now setup the entries of the CSC matrix
  // we reuse the colsize array to count down to zero
  // for determining the position of each nonzero
  Aval.resize(nnz);
  Aindex.resize(nnz);
  size_t numslots = Avalue.size();
  assert(numslots - freeslots.size() == static_cast<size_t>(nnz));
  for (size_t i = 0; i != numslots; ++i) {
    if (Avalue[i] == 0.0) continue;
    assert(Acol[i] >= 0 && Acol[i] < model->num_col_);
    HighsInt pos = Astart[Acol[i] + 1] - colsize[Acol[i]];
    --colsize[Acol[i]];
    assert(colsize[Acol[i]] >= 0);
    Aval[pos] = Avalue[i];
    Aindex[pos] = Arow[i];
  }
}

HPresolve::Result HPresolve::doubletonEq(HighsPostsolveStack& postsolve_stack,
                                         HighsInt row,
                                         HighsPostsolveStack::RowType rowType) {
  assert(this->allow_rule_[kPresolveRuleDoubletonEquation]);
  const bool logging_on = analysis_.logging_on_;
  if (logging_on)
    analysis_.startPresolveRuleLog(kPresolveRuleDoubletonEquation);
  assert(!rowDeleted[row]);
  assert(rowsize[row] == 2);
  assert(isEquation(row));

  // printf("doubleton equation: ");
  // debugPrintRow(row);
  HighsInt nzPos1 = rowroot[row];
  HighsInt nzPos2 = ARright[nzPos1] != -1 ? ARright[nzPos1] : ARleft[nzPos1];

  auto colAtPos1Better = [&]() {
    if (model->integrality_[Acol[nzPos1]] == HighsVarType::kInteger) {
      if (model->integrality_[Acol[nzPos2]] == HighsVarType::kInteger) {
        // both columns integer. For substitution choose smaller absolute
        // coefficient value, or sparser column if values are equal
        if (std::fabs(Avalue[nzPos1]) <
            std::fabs(Avalue[nzPos2]) - options->small_matrix_value) {
          return true;
        } else if (std::fabs(Avalue[nzPos2]) <
                   std::fabs(Avalue[nzPos1]) - options->small_matrix_value) {
          return false;
        } else if (colsize[Acol[nzPos1]] < colsize[Acol[nzPos2]]) {
          return true;
        } else {
          return false;
        }
      } else {
        // one col is integral, substitute the continuous one
        return false;
      }
    } else {
      if (model->integrality_[Acol[nzPos2]] == HighsVarType::kInteger) {
        // one col is integral, substitute the continuous one
        return true;
      } else {
        // both columns continuous the one with a larger absolute coefficient
        // value if the difference is more than factor 2, and otherwise the one
        // with fewer nonzeros if those are equal
        HighsInt col1Size = colsize[Acol[nzPos1]];
        if (col1Size == 1)
          return true;
        else {
          HighsInt col2Size = colsize[Acol[nzPos2]];
          if (col2Size == 1)
            return false;
          else {
            double abs1Val = std::fabs(Avalue[nzPos1]);
            double abs2Val = std::fabs(Avalue[nzPos2]);
            if (col1Size != col2Size &&
                std::max(abs1Val, abs2Val) <= 2.0 * std::min(abs1Val, abs2Val))
              return (col1Size < col2Size);
            else if (abs1Val > abs2Val)
              return true;
            else
              return false;
          }
        }
      }
    }
  };

  HighsInt substcol;
  HighsInt staycol;
  double substcoef;
  double staycoef;

  if (colAtPos1Better()) {
    substcol = Acol[nzPos1];
    staycol = Acol[nzPos2];

    substcoef = Avalue[nzPos1];
    staycoef = Avalue[nzPos2];
  } else {
    substcol = Acol[nzPos2];
    staycol = Acol[nzPos1];

    substcoef = Avalue[nzPos2];
    staycoef = Avalue[nzPos1];
  }

  double rhs = model->row_upper_[row];
  if (model->integrality_[substcol] == HighsVarType::kInteger &&
      model->integrality_[staycol] == HighsVarType::kInteger) {
    // check integrality conditions
    double roundCoef = std::round(staycoef / substcoef) * substcoef;
    if (std::fabs(roundCoef - staycoef) > options->small_matrix_value)
      return Result::kOk;
    staycoef = roundCoef;
    double roundRhs = std::round(rhs / substcoef) * substcoef;
    if (std::fabs(rhs - roundRhs) > primal_feastol)
      return Result::kPrimalInfeasible;
    rhs = roundRhs;
  }

  double oldStayLower = model->col_lower_[staycol];
  double oldStayUpper = model->col_upper_[staycol];
  double substLower = model->col_lower_[substcol];
  double substUpper = model->col_upper_[substcol];

  double stayImplLower;
  double stayImplUpper;
  if (std::signbit(substcoef) != std::signbit(staycoef)) {
    // coefficients have the opposite sign, therefore the implied lower bound of
    // the stay column is computed from the lower bound of the substituted
    // column:
    // staycol * staycoef + substcol * substcoef = rhs
    // staycol = (rhs - substcol * substcoef) / staycoef
    // staycol >= rhs / staycoef + lower(-substcoef/staycoef * substcol)
    // lower(-substcoef/staycoef * substcol) is (-substcoef/staycoef) *
    // substLower if (-substcoef/staycoef) is positive, i.e. if the coefficients
    // have opposite sign
    stayImplLower = substLower == -kHighsInf
                        ? -kHighsInf
                        : static_cast<double>((static_cast<HighsCDouble>(rhs) -
                                               substcoef * substLower) /
                                              staycoef);
    stayImplUpper = substUpper == kHighsInf
                        ? kHighsInf
                        : static_cast<double>((static_cast<HighsCDouble>(rhs) -
                                               substcoef * substUpper) /
                                              staycoef);
  } else {
    stayImplLower = substUpper == kHighsInf
                        ? -kHighsInf
                        : static_cast<double>((static_cast<HighsCDouble>(rhs) -
                                               substcoef * substUpper) /
                                              staycoef);
    stayImplUpper = substLower == -kHighsInf
                        ? kHighsInf
                        : static_cast<double>((static_cast<HighsCDouble>(rhs) -
                                               substcoef * substLower) /
                                              staycoef);
  }

  // possibly tighten bounds of the column that stays
  bool lowerTightened = stayImplLower > oldStayLower + primal_feastol;
  if (lowerTightened)
    HPRESOLVE_CHECKED_CALL(changeColLower(staycol, stayImplLower));

  bool upperTightened = stayImplUpper < oldStayUpper - primal_feastol;
  if (upperTightened)
    HPRESOLVE_CHECKED_CALL(changeColUpper(staycol, stayImplUpper));

  postsolve_stack.doubletonEquation(
      row, substcol, staycol, substcoef, staycoef, rhs, substLower, substUpper,
      model->col_cost_[substcol], lowerTightened, upperTightened, rowType,
      getColumnVector(substcol));

  // finally modify matrix
  HPRESOLVE_CHECKED_CALL(substitute(substcol, staycol, rhs / substcoef,
                                    -staycoef / substcoef, row));

  analysis_.logging_on_ = logging_on;
  if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleDoubletonEquation);

  // since a column was deleted we might have new row singletons which we
  // immediately remove
  HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolve_stack));

  return checkLimits(postsolve_stack);
}

HPresolve::Result HPresolve::singletonRow(HighsPostsolveStack& postsolve_stack,
                                          HighsInt row) {
  const bool logging_on = analysis_.logging_on_;
  if (logging_on) analysis_.startPresolveRuleLog(kPresolveRuleSingletonRow);
  assert(!rowDeleted[row]);
  assert(rowsize[row] == 1);

  // the tree of nonzeros of this row should just contain the single nonzero
  HighsInt nzPos = rowroot[row];
  assert(nzPos != -1);
  // nonzero should have the row in the row array
  assert(Arow[nzPos] == row);
  // tree with one element should not have children
  assert(ARleft[nzPos] == -1);
  assert(ARright[nzPos] == -1);

  HighsInt col = Acol[nzPos];
  double val = Avalue[nzPos];

  // printf("singleton row\n");
  // debugPrintRow(row);
  // delete row singleton nonzero directly, we have all information that we
  // need in local variables
  markRowDeleted(row);
  unlink(nzPos);

  // zeros should not be linked in the matrix
  assert(std::fabs(val) > options->small_matrix_value);

  // check for simple redundancy, compute tightened bounds, and check
  // whether the bounds are equal in tolerances
  double lb, ub;
  bool lowerTightened, upperTightened;
  SingletonRowResult sr = computeSingletonRowBounds(
      val, model->row_lower_[row], model->row_upper_[row],
      model->col_lower_[col], model->col_upper_[col], primal_feastol,
      getMaxAbsColVal(col),
      model->integrality_[col] != HighsVarType::kContinuous, lb, ub,
      lowerTightened, upperTightened);
  if (sr == SingletonRowResult::kRedundant) {
    postsolve_stack.redundantRow(row);
    analysis_.logging_on_ = logging_on;
    if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleSingletonRow);
    return checkLimits(postsolve_stack);
  }
  if (sr == SingletonRowResult::kPrimalInfeasible)
    return Result::kPrimalInfeasible;

  postsolve_stack.singletonRow(row, col, val, lowerTightened, upperTightened);

  // just update bounds (and row activities)
  if (lowerTightened) HPRESOLVE_CHECKED_CALL(changeColLower(col, lb));
  // update bounds, or remove as fixed column directly
  if (ub == lb) {
    postsolve_stack.removedFixedCol(col, lb, model->col_cost_[col],
                                    getColumnVector(col));
    HPRESOLVE_CHECKED_CALL(removeFixedCol(col));
  } else if (upperTightened)
    HPRESOLVE_CHECKED_CALL(changeColUpper(col, ub));

  if (!colDeleted[col] && colsize[col] == 0) {
    Result result = emptyCol(postsolve_stack, col);
    analysis_.logging_on_ = logging_on;
    if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleSingletonRow);
    return result;
  }
  analysis_.logging_on_ = logging_on;
  if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleSingletonRow);
  return checkLimits(postsolve_stack);
}

HPresolve::Result HPresolve::singletonCol(HighsPostsolveStack& postsolve_stack,
                                          HighsInt col, const bool timing) {
  assert(colsize[col] == 1);
  assert(!colDeleted[col]);
  HighsInt nzPos = colhead[col];
  HighsInt row = Arow[nzPos];
  double colCoef = Avalue[nzPos];

  if (rowsize[row] == 1) {
    if (timing)
      analysis_.presolveTimerStart(kPresolveClockSingletonColSingletonRow);
    HPRESOLVE_CHECKED_CALL(singletonRow(postsolve_stack, row););

    if (!colDeleted[col]) {
      assert(colsize[col] == 0);
      HPresolve::Result result = emptyCol(postsolve_stack, col);
      if (timing)
        analysis_.presolveTimerStop(kPresolveClockSingletonColSingletonRow);
      return result;
    }
    if (timing)
      analysis_.presolveTimerStop(kPresolveClockSingletonColSingletonRow);
    return Result::kOk;
  }

  // detect strong / weak domination
  if (timing) analysis_.presolveTimerStart(kPresolveClockSingletonColDominated);
  // Pass handleSingletonRows = false
  HPRESOLVE_CHECKED_CALL(detectDominatedCol(postsolve_stack, col, false));
  if (timing) analysis_.presolveTimerStop(kPresolveClockSingletonColDominated);
  if (colDeleted[col]) return Result::kOk;

  // check if variable is implied integer
  if (mipsolver != nullptr)
    HPRESOLVE_CHECKED_CALL(
        static_cast<Result>(convertImpliedInteger(col, row)));

  // dual fixing
  if (this->allow_rule_[kPresolveRuleDualFixing]) {
    const bool logging_on = analysis_.logging_on_;
    if (logging_on) analysis_.startPresolveRuleLog(kPresolveRuleDualFixing);
    if (timing)
      analysis_.presolveTimerStart(kPresolveClockSingletonColDualFixing);
    HPRESOLVE_CHECKED_CALL(dualFixing(postsolve_stack, col));
    if (timing)
      analysis_.presolveTimerStop(kPresolveClockSingletonColDualFixing);
    analysis_.logging_on_ = logging_on;
    if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleDualFixing);
    if (colDeleted[col]) return Result::kOk;
  }

  // singleton column stuffing
  if (this->allow_rule_[kPresolveRuleColStuffing]) {
    const bool logging_on = analysis_.logging_on_;
    if (logging_on) analysis_.startPresolveRuleLog(kPresolveRuleColStuffing);
    if (timing)
      analysis_.presolveTimerStart(kPresolveClockSingletonColStuffing);
    HPRESOLVE_CHECKED_CALL(singletonColStuffing(postsolve_stack, col));
    if (timing) analysis_.presolveTimerStop(kPresolveClockSingletonColStuffing);
    analysis_.logging_on_ = logging_on;
    if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleColStuffing);
    if (colDeleted[col]) return Result::kOk;
  };

  // update column implied bounds
  if (timing)
    analysis_.presolveTimerStart(kPresolveClockSingletonColImpliedBounds);
  HPRESOLVE_CHECKED_CALL(updateColImpliedBounds(row, col, colCoef));
  if (timing)
    analysis_.presolveTimerStop(kPresolveClockSingletonColImpliedBounds);

  // update row dual implied bounds
  if (model->integrality_[col] != HighsVarType::kInteger) {
    if (timing)
      analysis_.presolveTimerStart(
          kPresolveClockSingletonColRowDualImpliedBounds);
    updateRowDualImpliedBounds(row, col, colCoef);
    if (timing)
      analysis_.presolveTimerStop(
          kPresolveClockSingletonColRowDualImpliedBounds);
  }
  // now check if column is implied free within an equation and substitute the
  // column if that is the case
  if (timing)
    analysis_.presolveTimerStart(kPresolveClockSingletonColDualImpliedFree);
  if (isDualImpliedFree(row) && isImpliedFree(col) &&
      this->allow_rule_[kPresolveRuleFreeColSubstitution]) {
    if (model->integrality_[col] == HighsVarType::kInteger) {
      StatusResult impliedIntegral = isImpliedIntegral(col);
      HPRESOLVE_CHECKED_CALL(static_cast<Result>(impliedIntegral));
      if (!impliedIntegral) {
        if (timing)
          analysis_.presolveTimerStop(
              kPresolveClockSingletonColDualImpliedFree);
        return Result::kOk;
      }
    }
    const bool logging_on = analysis_.logging_on_;

    if (logging_on)
      analysis_.startPresolveRuleLog(kPresolveRuleFreeColSubstitution);

    // todo, store which side of an implied free dual variable needs to be used
    // for substitution
    storeRow(row);

    substituteFreeCol(postsolve_stack, row, col);

    analysis_.logging_on_ = logging_on;
    if (logging_on)
      analysis_.stopPresolveRuleLog(kPresolveRuleFreeColSubstitution);
    if (timing)
      analysis_.presolveTimerStop(kPresolveClockSingletonColDualImpliedFree);
    return checkLimits(postsolve_stack);
  }
  if (timing)
    analysis_.presolveTimerStop(kPresolveClockSingletonColDualImpliedFree);

  if (this->allow_rule_[kPresolveRuleZeroCostSingleton]) {
    // Remove if col is double-sided finite slack
    HPRESOLVE_CHECKED_CALL(zeroCostSingleton(postsolve_stack, col));
  }

  return Result::kOk;
}

HPresolve::Result HPresolve::emptyRow(HighsPostsolveStack& postsolve_stack,
                                      HighsInt row) {
  // Special case of rowPresolve for rows known to be empty
  //
  // Check that the row is feasible
  if (model->row_upper_[row] < -primal_feastol ||
      model->row_lower_[row] > primal_feastol)
    return Result::kPrimalInfeasible;
  postsolve_stack.redundantRow(row);
  return checkLimits(postsolve_stack);
}

HPresolve::Result HPresolve::rowPresolve(HighsPostsolveStack& postsolve_stack,
                                         HighsInt row) {
  assert(!rowDeleted[row]);

  const bool logging_on = analysis_.logging_on_;

  auto checkRowInfeasible = [&](HighsInt row) {
    if (impliedRowBounds.getSumLower(row) >
            model->row_upper_[row] + primal_feastol ||
        impliedRowBounds.getSumUpper(row) <
            model->row_lower_[row] - primal_feastol)
      return Result::kPrimalInfeasible;
    return Result::kOk;
  };

  auto checkRowRedundant = [&](HighsInt row) {
    if (isRedundant(row)) {
      // row is redundant
      int presolveRule =
          rowsize[row] != 0 ? kPresolveRuleRedundantRow : kPresolveRuleEmptyRow;
      if (logging_on) analysis_.startPresolveRuleLog(presolveRule);
      postsolve_stack.redundantRow(row);
      removeRow(row);
      analysis_.logging_on_ = logging_on;
      if (logging_on) analysis_.stopPresolveRuleLog(presolveRule);
      return checkLimits(postsolve_stack);
    }
    return Result::kOk;
  };

  // printf("row presolve: ");
  // debugPrintRow(row);

  // check for infeasibility
  HPRESOLVE_CHECKED_CALL(checkRowInfeasible(row));

  // handle special cases directly via a call to the specialized procedure
  switch (rowsize[row]) {
    default:
      break;
    case 1:
      return singletonRow(postsolve_stack, row);
  }

  // check for redundancy
  HPRESOLVE_CHECKED_CALL(checkRowRedundant(row));
  if (rowDeleted[row]) return Result::kOk;

  if (rowsizeInteger[row] != 0 || rowsizeImplInt[row] != 0) {
    // check if setting variable to opposite bound (with respect to bound used
    // in calculation of bounds on row activity) makes row infeasible.
    // this is a special case of bound tightening, but only performed once for
    // binary variables here (see Suhl, Szymanski: Supernode processing of
    // mixed-integer models. Comput. Optim. Appl. 3(4): 317-331 (1994)).

    for (const HighsSliceNonzero& nonzero : getRowVector(row)) {
      // get column index and coefficient
      HighsInt col = nonzero.index();
      double val = nonzero.value();

      // skip continuous variables and non-binary variables
      if (model->integrality_[col] == HighsVarType::kContinuous ||
          model->col_upper_[col] != model->col_lower_[col] + 1.0)
        continue;

      // lambda for computing offset
      auto computeOffset = [&](HighsInt col, double val) {
        return std::abs(val) *
               (static_cast<HighsCDouble>(model->col_upper_[col]) -
                static_cast<HighsCDouble>(model->col_lower_[col]));
      };

      auto degree1Tests = [&](HighsInt col, double val, HighsInt direction,
                              double rowActivityBound, double rowBound) {
        // return if row is still feasible
        if (direction * rowActivityBound >=
            direction * rowBound - primal_feastol)
          return Result::kOk;

        // tighten bound
        if (direction * val > 0)
          HPRESOLVE_CHECKED_CALL(
              changeColLower(col, model->col_lower_[col] + 1.0));
        else
          HPRESOLVE_CHECKED_CALL(
              changeColUpper(col, model->col_upper_[col] - 1.0));
        return Result::kOk;
      };

      // perform tests
      HPRESOLVE_CHECKED_CALL(degree1Tests(
          col, val, HighsInt{1},
          impliedRowBounds.getSumUpperOrig(row, -computeOffset(col, val)),
          model->row_lower_[row]));
      HPRESOLVE_CHECKED_CALL(degree1Tests(
          col, val, HighsInt{-1},
          impliedRowBounds.getSumLowerOrig(row, computeOffset(col, val)),
          model->row_upper_[row]));
    }
  }

  // check for redundancy again
  HPRESOLVE_CHECKED_CALL(checkRowRedundant(row));
  if (rowDeleted[row]) return Result::kOk;

  auto checkRedundantBounds = [&](HighsInt col, HighsInt row) {
    // check if column singleton has redundant bounds
    assert(model->col_cost_[col] != 0.0);
    if (colsize[col] != 1) return Result::kOk;
    if (model->col_cost_[col] > 0) {
      assert(model->col_lower_[col] == -kHighsInf || !isLowerImplied(col) ||
             colLowerSource[col] == row);
      if (model->col_lower_[col] > implColLower[col] - primal_feastol)
        HPRESOLVE_CHECKED_CALL(changeColLower(col, -kHighsInf));
    } else {
      assert(model->col_upper_[col] == kHighsInf || !isUpperImplied(col) ||
             colUpperSource[col] == row);
      if (model->col_upper_[col] < implColUpper[col] + primal_feastol)
        HPRESOLVE_CHECKED_CALL(changeColUpper(col, kHighsInf));
    }
    return Result::kOk;
  };

  // Store original bounds
  double origRowUpper = model->row_upper_[row];
  double origRowLower = model->row_lower_[row];

  // Convert to equality constraint and record for dual postsolve
  if (!isEquation(row)) {
    if (isImpliedEquationAtLower(row)) {
      // save source before changeRowUpper resets it
      HighsInt dualLowerSourceCol = rowDualLowerSource[row];
      changeRowUpper(row, model->row_lower_[row]);
      postsolve_stack.impliedEquation(row, true, getRowVector(row));
      if (mipsolver == nullptr)
        HPRESOLVE_CHECKED_CALL(checkRedundantBounds(dualLowerSourceCol, row));
    } else if (isImpliedEquationAtUpper(row)) {
      // save source before changeRowLower resets it
      HighsInt dualUpperSourceCol = rowDualUpperSource[row];
      changeRowLower(row, model->row_upper_[row]);
      postsolve_stack.impliedEquation(row, false, getRowVector(row));
      if (mipsolver == nullptr)
        HPRESOLVE_CHECKED_CALL(checkRedundantBounds(dualUpperSourceCol, row));
    }
  }

  // Get row bounds
  double rowUpper = model->row_upper_[row];
  double rowLower = model->row_lower_[row];

  // Handle doubleton equations
  if (rowsize[row] == 2 && rowLower == rowUpper &&
      this->allow_rule_[kPresolveRuleDoubletonEquation]) {
    HighsPostsolveStack::RowType rowType;
    if (origRowLower == origRowUpper) {
      rowType = HighsPostsolveStack::RowType::kEq;
    } else if (origRowUpper != kHighsInf) {
      rowType = HighsPostsolveStack::RowType::kLeq;
    } else {
      assert(origRowLower != -kHighsInf);
      rowType = HighsPostsolveStack::RowType::kGeq;
    }
    return doubletonEq(postsolve_stack, row, rowType);
  }

  // todo: do additional single row presolve for mip here. It may assume a
  // non-redundant and non-infeasible row when considering variable and implied
  // bounds
  if (rowsizeInteger[row] != 0 || rowsizeImplInt[row] != 0) {
    if (rowLower == rowUpper) {
      // equation
      double impliedRowLower = impliedRowBounds.getSumLower(row);
      double impliedRowUpper = impliedRowBounds.getSumUpper(row);
      // see section 3.6 "Simple probing on a single equation", Achterberg et
      // al., Presolve Reductions in Mixed Integer Programming, INFORMS Journal
      // on Computing 32(2):473-506.
      if (impliedRowLower != -kHighsInf && impliedRowUpper != kHighsInf &&
          std::abs(impliedRowLower + impliedRowUpper - 2 * rowUpper) <=
              options->small_matrix_value) {
        // simple probing on equation case
        HighsInt binCol = -1;
        double binCoef = std::abs(impliedRowUpper - rowUpper);
        // store row
        storeRow(row);
        for (const HighsSliceNonzero& nonz : getStoredRow()) {
          if (std::abs(std::abs(nonz.value()) - binCoef) <=
                  options->small_matrix_value &&
              model->integrality_[nonz.index()] == HighsVarType::kInteger &&
              std::abs(model->col_upper_[nonz.index()] -
                       model->col_lower_[nonz.index()] - 1.0) <=
                  primal_feastol) {
            // found a binary variable that implies all other variables to be
            // fixed when it sits at one of its bounds therefore we can
            // substitute all other variables in the row
            binCol = nonz.index();
            // store the binary coefficient with its actual sign
            binCoef = nonz.value();
            break;
          }
        }

        if (binCol != -1) {
          // found binary column for substituting all other columns
          // printf("simple probing case on row of size %" HIGHSINT_FORMAT
          // "\n", rowsize[row]);

          // Snapshot bounds as they may change when removing or substituting
          // a column, which would make future substitutions invalid
          std::vector<std::pair<double, double>> implVarSnapshot;
          implVarSnapshot.reserve(rowpositions.size());
          for (HighsInt rowiter : rowpositions) {
            HighsInt col = Acol[rowiter];
            implVarSnapshot.emplace_back(
                impliedRowBounds.getImplVarLower(row, col),
                impliedRowBounds.getImplVarUpper(row, col));
          }

          // iterate over non-zero positions instead of iterating over the
          // HighsMatrixSlice (provided by HPresolve::getStoredRow) because the
          // latter contains pointers to Acol and Avalue that may be invalidated
          // if these vectors are reallocated (see std::vector::push_back
          // performed in HPresolve::addToMatrix).
          for (size_t i = 0; i < rowpositions.size(); ++i) {
            HighsInt rowiter = rowpositions[i];
            HighsInt col = Acol[rowiter];
            assert(Arow[rowiter] == row);

            // skip binary column
            if (col == binCol) continue;

            // Use eliminateCol rather than substituteCol to avoid
            // cascading clique fixings that could corrupt this loop
            if (mipsolver != nullptr && mipsolver->mipdata_->cliquesExtracted) {
              presolveCliqueTable.eliminateCol(col);
            }

            // get column lower and upper bounds used to compute bounds on row
            // activities
            double col_lower = implVarSnapshot[i].first;
            double col_upper = implVarSnapshot[i].second;
            assert(col_lower != -kHighsInf);
            assert(col_upper != kHighsInf);

            if (model->col_lower_[col] == model->col_upper_[col]) {
              postsolve_stack.removedFixedCol(col, model->col_lower_[col], 0.0,
                                              HighsEmptySlice());
              HPRESOLVE_CHECKED_CALL(removeFixedCol(col));
              continue;
            }

            auto remDoubletonEq = [&](HighsInt col, HighsInt binCol,
                                      HighsInt direction, double lower,
                                      double upper) {
              double bound = direction >= 0 ? upper : lower;
              double scale = direction * (lower - upper);
              double offset = bound - model->col_lower_[binCol] * scale;
              postsolve_stack.doubletonEquation(
                  -1, col, binCol, 1.0, -scale, offset, lower, upper, 0.0,
                  false, false, HighsPostsolveStack::RowType::kEq,
                  HighsEmptySlice());
              return substitute(col, binCol, offset, scale);
            };

            // 1. binary coefficient is positive:
            // setting the binary to its upper bound
            // increases the minimal activity to be equal to the row upper
            // bound and therefore all other variables are fixed to the
            // bound that contributes to the rows minimal activity, i.e. the
            // lower bound for a positive coefficient
            // This case yields the following implications:
            // binCol = ub -> nonzCol = lb
            // binCol = lb -> nonzCol = ub
            // as linear equation:
            // nonzCol = colUb - (colUb - colLb)(binCol - binLb)
            // nonzCol = colUb + binLb * (colUb - colLb) - (colUb - colLb) *
            // binCol
            // 2. the opposite case yields the following implications:
            // binCol = lb -> nonzCol = lb
            // binCol = ub -> nonzCol = ub
            // as linear equation:
            // nonzCol = colLb + (colUb - colLb)(binCol - binLb)
            // nonzCol =
            //    colLb - binLb*(colUb - colLb) + (colUb - colLb)*binCol
            HPRESOLVE_CHECKED_CALL(remDoubletonEq(
                col, binCol,
                std::signbit(binCoef) == std::signbit(Avalue[rowiter])
                    ? HighsInt{1}
                    : HighsInt{-1},
                col_lower, col_upper));
          }

          removeRow(row);
          HPRESOLVE_CHECKED_CALL(checkLimits(postsolve_stack));
          return removeRowSingletons(postsolve_stack);
        }
      }

      if (rowsizeInteger[row] + rowsizeImplInt[row] >= rowsize[row] - 1) {
        HighsInt continuousCol = -1;
        double continuousCoef = 0.0;
        std::vector<double> rowCoefsInt;
        rowCoefsInt.reserve(rowsize[row]);
        storeRow(row);

        for (const HighsSliceNonzero& nonz : getStoredRow()) {
          if (model->integrality_[nonz.index()] == HighsVarType::kContinuous) {
            assert(continuousCoef == 0.0);
            continuousCoef = nonz.value();
            continuousCol = nonz.index();
            continue;
          }

          rowCoefsInt.push_back(nonz.value());
        }

        if (continuousCoef != 0.0) {
          rowCoefsInt.push_back(rowUpper);

          double intScale = HighsIntegers::integralScale(
              rowCoefsInt, options->small_matrix_value,
              options->small_matrix_value);

          if (intScale != 0 && intScale <= 1e3) {
            double scale = 1.0 / std::abs(continuousCoef * intScale);
            if (scale != 1.0) {
              // printf(
              //     "transform continuous column x to implicit integer z with x
              //     "
              //     "= %g * z\n",
              //     scale);
              HPRESOLVE_CHECKED_CALL(
                  transformColumn(postsolve_stack, continuousCol, scale, 0.0));

              HPRESOLVE_CHECKED_CALL(static_cast<Result>(
                  convertImpliedInteger(continuousCol, -1, true)));

              if (intScale != 1.0) scaleStoredRow(row, intScale, true);
            }
          }
        } else {
          double intScale = HighsIntegers::integralScale(
              rowCoefsInt, options->small_matrix_value,
              options->small_matrix_value);

          if (intScale != 0.0 && intScale <= 1e3) {
            double rhs = rowUpper * intScale;
            if (fractionality(rhs) > primal_feastol)
              return Result::kPrimalInfeasible;

            rhs = std::round(rhs);

            HighsInt x1Cand = -1;
            int64_t d = 0;

            for (HighsInt i = 0; i < static_cast<HighsInt>(rowpositions.size());
                 ++i) {
              int64_t newgcd =
                  d == 0 ? int64_t(std::abs(
                               std::round(intScale * Avalue[rowpositions[i]])))
                         : HighsIntegers::gcd(
                               std::abs(std::round(intScale *
                                                   Avalue[rowpositions[i]])),
                               d);
              if (newgcd == 1) {
                // adding this variable would set the gcd to 1, therefore it
                // must be our candidate x1 for substitution. If another
                // candidate already exists no reduction is possible except for
                // scaling the equation
                if (x1Cand != -1) {
                  x1Cand = -1;
                  break;
                }
                x1Cand = i;
              } else {
                d = newgcd;
              }
            }

            if (x1Cand != -1) {
              HighsInt x1Pos = rowpositions[x1Cand];
              HighsInt x1 = Acol[x1Pos];
              double rhs2 = rhs / static_cast<double>(d);
              if (fractionality(rhs2) <= mipsolver->mipdata_->epsilon) {
                // the right hand side is integral, so we can substitute
                // x1 = d * z

                // printf(
                //    "substitute integral column x with integral column z with
                //    " "x = %ld * z\n", d);
                HPRESOLVE_CHECKED_CALL(
                    transformColumn(postsolve_stack, x1, d, 0.0));
              } else {
                // we can substitute x1 = d * z + b, with b = a1^-1 rhs (mod d)

                // first compute the modular multiplicative inverse of a1^-1
                // (mod d) of a1
                int64_t a1 = std::round(intScale * Avalue[x1Pos]);
                a1 = HighsIntegers::mod(a1, d);
                int64_t a1Inverse = HighsIntegers::modularInverse(a1, d);

                // now compute b = a1^-1 rhs (mod d)
                double b =
                    HighsIntegers::mod(a1Inverse * rhs, static_cast<double>(d));

                // printf(
                //     "substitute integral column x with integral column z with
                //     " "x = %ld * z + %g\n", d, b);

                // before we substitute, we check whether the resulting variable
                // z is fixed after rounding its new bounds. If that is the case
                // we directly fix x1 instead of first substituting with d * z +
                // b.
                double zLower = std::ceil((model->col_lower_[x1] - b) /
                                              static_cast<double>(d) -
                                          primal_feastol);
                double zUpper = std::floor((model->col_upper_[x1] - b) /
                                               static_cast<double>(d) +
                                           primal_feastol);

                if (zLower == zUpper) {
                  // Rounded bounds are equal
                  //
                  // Adjust bounds if variable is fixed to a value in between
                  // its bounds
                  double fixVal = zLower * d + b;
                  assert(fixVal > model->col_lower_[x1] - primal_feastol);
                  assert(fixVal < model->col_upper_[x1] + primal_feastol);
                  // change bounds
                  HPRESOLVE_CHECKED_CALL(changeColBounds(x1, fixVal, fixVal));
                  // remove column
                  postsolve_stack.removedFixedCol(
                      x1, fixVal, model->col_cost_[x1], getColumnVector(x1));
                  HPRESOLVE_CHECKED_CALL(removeFixedCol(x1));
                  rowpositions.erase(rowpositions.begin() + x1Cand);
                } else {
                  HPRESOLVE_CHECKED_CALL(
                      transformColumn(postsolve_stack, x1, d, b));
                }
              }

              intScale /= d;
            }

            if (intScale != 1.0) scaleStoredRow(row, intScale, true);
          }
        }
      }
    } else {
      // inequality or ranged row, first store row positions
      storeRow(row);

      if (rowsize[row] == rowsizeInteger[row] + rowsizeImplInt[row]) {
        std::vector<double> rowCoefs;
        std::vector<HighsInt> rowIndex;
        rowCoefs.reserve(rowsize[row]);
        rowIndex.reserve(rowsize[row]);

        for (const HighsSliceNonzero& nonz : getStoredRow()) {
          assert(nonz.value() != 0.0);
          rowCoefs.push_back(nonz.value());
          rowIndex.push_back(nonz.index());
        }

        double intScale = HighsIntegers::integralScale(
            rowCoefs,
            model->row_lower_[row] == -kHighsInf ? primal_feastol
                                                 : options->small_matrix_value,
            model->row_upper_[row] == kHighsInf ? primal_feastol
                                                : options->small_matrix_value);

        auto roundRhs = [&](HighsCDouble rhs, HighsCDouble& roundedRhs,
                            HighsCDouble& fractionRhs, double minRhsTightening,
                            bool& rhsTightened, HighsInt direction) {
          // round rhs
          roundedRhs = direction * floor(direction * rhs + primal_feastol);
          // compute fractional part
          fractionRhs = direction * (rhs - roundedRhs);
          // check if tightened
          rhsTightened =
              fractionRhs >= minRhsTightening - options->small_matrix_value;
        };

        auto checkScaleRow =
            [&](HighsInt row, HighsCDouble& roundedLhs,
                HighsCDouble& roundedRhs, HighsCDouble& fractionLhs,
                HighsCDouble& fractionRhs, double& minLhsTightening,
                double& minRhsTightening, double& maxVal, bool& lhsTightened,
                bool& rhsTightened, bool& isInfeasible, double intScale) {
              HighsCDouble lhs = model->row_lower_[row];
              HighsCDouble rhs = model->row_upper_[row];
              bool lhsFinite = lhs != -kHighsInf;
              bool rhsFinite = rhs != kHighsInf;
              if (lhsFinite) lhs = lhs * intScale;
              if (rhsFinite) rhs = rhs * intScale;
              roundedLhs = -kHighsInf;
              roundedRhs = kHighsInf;
              fractionLhs = 0.0;
              fractionRhs = 0.0;
              minRhsTightening = 0.0;
              minLhsTightening = 0.0;
              maxVal = 0.0;
              lhsTightened = false;
              rhsTightened = false;
              isInfeasible = false;
              for (size_t i = 0; i < rowCoefs.size(); ++i) {
                // computed scaled coefficient
                HighsCDouble scaleCoef =
                    static_cast<HighsCDouble>(rowCoefs[i]) * intScale;
                // round to the nearest integer
                HighsCDouble intCoef = floor(scaleCoef + 0.5);
                // compute difference
                HighsCDouble coefDelta = intCoef - scaleCoef;
                // store integral coefficient and maximum absolute value
                rowCoefs[i] = static_cast<double>(intCoef);
                maxVal = std::max(std::abs(rowCoefs[i]), maxVal);
                // get column upper bound
                double ub = model->col_upper_[rowIndex[i]];
                if (coefDelta < -options->small_matrix_value) {
                  // for the >= side of the constraint a smaller coefficient is
                  // stronger: Therefore we relax the left hand side using the
                  // bound constraint, if the bound is infinite, abort
                  if (lhsFinite) {
                    if (ub == kHighsInf) return false;
                    lhs += ub * coefDelta;
                  }
                  minRhsTightening = std::max(-static_cast<double>(coefDelta),
                                              minRhsTightening);
                } else if (coefDelta > options->small_matrix_value) {
                  if (rhsFinite) {
                    if (ub == kHighsInf) return false;
                    rhs += ub * coefDelta;
                  }
                  // the coefficient was relaxed regarding the rows lower bound.
                  // Therefore the lower bound should be tightened by at least
                  // this amount for the scaled constraint to dominate the
                  // unscaled constraint be rounded by at least this value
                  minLhsTightening = std::max(static_cast<double>(coefDelta),
                                              minLhsTightening);
                }
              }
              // round left-hand and right-hand sides
              if (lhsFinite)
                roundRhs(lhs, roundedLhs, fractionLhs, minLhsTightening,
                         lhsTightened, HighsInt{-1});
              if (rhsFinite)
                roundRhs(rhs, roundedRhs, fractionRhs, minRhsTightening,
                         rhsTightened, HighsInt{1});
              // check for infeasibility
              isInfeasible =
                  lhsFinite && rhsFinite && roundedRhs < roundedLhs - 0.5;
              return true;
            };

        auto scaleRow = [&](HighsInt row, HighsCDouble roundedLhs,
                            HighsCDouble roundedRhs, double scalar,
                            bool checkDelta) {
          // scale the row
          if (roundedLhs != -kHighsInf)
            changeRowLower(row, static_cast<double>(roundedLhs / scalar));
          if (roundedRhs != kHighsInf)
            changeRowUpper(row, static_cast<double>(roundedRhs / scalar));
          for (size_t i = 0; i < rowCoefs.size(); ++i) {
            double delta = static_cast<double>(
                static_cast<HighsCDouble>(rowCoefs[i]) / scalar -
                Avalue[rowpositions[i]]);
            if (!checkDelta || std::fabs(delta) > options->small_matrix_value)
              addToMatrix(row, rowIndex[i], delta);
          }
        };

        auto scaleRowIntVals = [&](HighsInt row, HighsCDouble roundedLhs,
                                   HighsCDouble roundedRhs, double intScale,
                                   double maxVal) {
          // return if maximum value in scaled row or scalar are too large
          if (maxVal > 1000.0 && intScale > 100.0) return false;
          // the scale value is reasonably small, change the row values to be
          // integral
          scaleRow(row, roundedLhs, roundedRhs, 1.0, false);
          return true;
        };

        if (intScale != 0.0) {
          HighsCDouble roundedLhs;
          HighsCDouble roundedRhs;
          HighsCDouble fractionLhs;
          HighsCDouble fractionRhs;
          double minLhsTightening;
          double minRhsTightening;
          double maxVal;
          bool lhsTightened;
          bool rhsTightened;
          bool isInfeasible;
          if (checkScaleRow(row, roundedLhs, roundedRhs, fractionLhs,
                            fractionRhs, minLhsTightening, minRhsTightening,
                            maxVal, lhsTightened, rhsTightened, isInfeasible,
                            intScale)) {
            // check for infeasibility
            if (isInfeasible) return Result::kPrimalInfeasible;
            // only accept row whose sides were tightened
            bool rangedOrEquationRow = lhsTightened && rhsTightened;
            if (rangedOrEquationRow ||
                (lhsTightened && model->row_upper_[row] == kHighsInf) ||
                (rhsTightened && model->row_lower_[row] == -kHighsInf)) {
              // check if constraint can be scaled to integral values
              if (!scaleRowIntVals(row, roundedLhs, roundedRhs, intScale,
                                   maxVal)) {
                if (rangedOrEquationRow) {
                  // ranged or equation row
                  // scale value is large, just tighten the sides
                  roundedLhs /= intScale;
                  roundedRhs /= intScale;
                  if (roundedRhs < model->row_upper_[row] - primal_feastol)
                    changeRowUpper(row, static_cast<double>(roundedRhs));
                  if (roundedLhs > model->row_lower_[row] + primal_feastol)
                    changeRowLower(row, static_cast<double>(roundedLhs));
                } else if ((rhsTightened &&
                            fractionRhs < minRhsTightening - primal_feastol) ||
                           (lhsTightened &&
                            fractionLhs < minLhsTightening - primal_feastol)) {
                  // <= or >= inequality
                  // scale value is large, so we scale back the altered
                  // constraint the scaled back constraint must be stronger
                  // than the original constraint for this to make sense with
                  // is checked with the condition above
                  scaleRow(row, roundedLhs, roundedRhs, intScale, true);
                }
              }
            }
          }
        } else if (!isRanged(row)) {
          // Chvatal-Gomory strengthening
          // See section 3.4 "Chvatal-Gomory strengthening of inequalities",
          // Achterberg et al., Presolve Reductions in Mixed Integer
          // Programming, INFORMS Journal on Computing 32(2):473-506.
          constexpr double maxDynamism = 1e5;
          std::vector<double> roundedRowCoefs;
          roundedRowCoefs.resize(rowsize[row]);
          std::set<double> scalars = {1.0};

          // lambda for reformulating row to only contain variables with lower
          // bounds of zero (if possible)
          auto complementOrShift = [&](HighsInt direction, HighsCDouble& rhs,
                                       double& minAbsCoef, double& maxAbsCoef,
                                       double& dynamism) {
            minAbsCoef = kHighsInf;
            maxAbsCoef = 0.0;
            for (size_t i = 0; i < rowCoefs.size(); ++i) {
              // get column index and (absolute) coefficient
              HighsInt col = rowIndex[i];
              double val = direction * rowCoefs[i];
              double absval = std::abs(val);
              // compute minimum and maximum absolute coefficients along the way
              minAbsCoef = std::min(minAbsCoef, absval);
              maxAbsCoef = std::max(maxAbsCoef, absval);
              if (val < 0.0 && model->col_upper_[col] != kHighsInf) {
                // complement
                rhs -= val * static_cast<HighsCDouble>(model->col_upper_[col]);
              } else if (val > 0.0 && model->col_lower_[col] != -kHighsInf) {
                // shift
                rhs -= val * static_cast<HighsCDouble>(model->col_lower_[col]);
              } else {
                // unbounded variable; cannot shift or complement!
                return false;
              }
            }
            dynamism = maxAbsCoef / minAbsCoef;
            return true;
          };

          // undo shifting and complementation
          auto undoComplementOrShift = [&](HighsInt direction,
                                           HighsCDouble& roundedRhs) {
            for (size_t i = 0; i < rowCoefs.size(); ++i) {
              HighsInt col = rowIndex[i];
              double val = direction * rowCoefs[i];
              if (val < 0.0 && model->col_upper_[col] != kHighsInf) {
                // uncomplement
                roundedRowCoefs[i] = -roundedRowCoefs[i];
                roundedRhs += roundedRowCoefs[i] *
                              static_cast<HighsCDouble>(model->col_upper_[col]);
              } else if (val > 0.0 && model->col_lower_[col] != -kHighsInf) {
                // unshift
                roundedRhs += roundedRowCoefs[i] *
                              static_cast<HighsCDouble>(model->col_lower_[col]);
              }
              // flip coefficient sign for <= inequality
              roundedRowCoefs[i] *= direction;
            }
            // flip rhs sign for <= inequality
            roundedRhs *= direction;
          };

          // round row using given scalar
          auto roundRow = [&](double s, const HighsCDouble& rhs,
                              HighsCDouble& roundedRhs) {
            bool accept = false;
            // round rhs (using feasibility tolerance)
            HighsCDouble scalar = static_cast<HighsCDouble>(s);
            roundedRhs = ceil(rhs * scalar - primal_feastol);
            if (roundedRhs <= primal_feastol) return false;
            HighsCDouble rhsRatio = rhs / roundedRhs;

            for (size_t i = 0; i < rowCoefs.size(); ++i) {
              // coefficient sign has not been flipped for complemented
              // variables; take absolute value of coefficient.
              double absCoef = std::abs(rowCoefs[i]);
              // round coefficient
              roundedRowCoefs[i] =
                  static_cast<double>(ceil(absCoef * scalar - kHighsTiny));
              // compare "normalised" coefficients, i.e. coefficients divided by
              // corresponding rhs.
              double threshold =
                  static_cast<double>(roundedRowCoefs[i] * rhsRatio);
              // return if coefficient is weaker
              if (absCoef < threshold - options->small_matrix_value)
                return false;
              // accept rounding if at least one coefficient is improved
              accept =
                  accept || (absCoef > threshold + options->small_matrix_value);
            }
            return accept;
          };

          // set up scalars suggested by Achterberg et al.
          auto setScalars = [&](double minAbsCoef, double maxAbsCoef) {
            for (HighsInt t = 1; t <= 5; t++) {
              scalars.emplace(t / maxAbsCoef);
              scalars.emplace(t / minAbsCoef);
              scalars.emplace((2 * t - 1) / (2 * minAbsCoef));
            }
          };

          // try different scalars and return if an improving one was found
          auto rowCanBeTightened = [&](const HighsCDouble& rhs,
                                       HighsCDouble& roundedRhs) {
            for (double s : scalars) {
              if (roundRow(s, rhs, roundedRhs)) return true;
            }
            return false;
          };

          // replace the model row by the rounded one
          auto updateRow = [&](HighsInt row, HighsInt direction,
                               const HighsCDouble& roundedRhs) {
            if (direction < 0)
              changeRowUpper(row, static_cast<double>(roundedRhs));
            else
              changeRowLower(row, static_cast<double>(roundedRhs));
            for (size_t i = 0; i < rowCoefs.size(); ++i) {
              double delta = static_cast<double>(
                  static_cast<HighsCDouble>(roundedRowCoefs[i]) - rowCoefs[i]);
              if (std::fabs(delta) > options->small_matrix_value)
                addToMatrix(row, rowIndex[i], delta);
            }
          };

          // convert to >= inequality
          // direction = -1: <= constraint, direction = 1: >= constraint
          HighsInt direction =
              model->row_upper_[row] != kHighsInf ? HighsInt{-1} : HighsInt{1};
          // get rhs
          HighsCDouble rhs =
              direction < 0 ? -model->row_upper_[row] : model->row_lower_[row];

          // initialise
          HighsCDouble roundedRhs = 0.0;
          double minAbsCoef = 0.0;
          double maxAbsCoef = 0.0;
          double dynamism = 0.0;

          // complement or shift variables to have lower bounds of zero
          if (complementOrShift(direction, rhs, minAbsCoef, maxAbsCoef,
                                dynamism)) {
            // identify scalars for row
            if (dynamism <= maxDynamism) setScalars(minAbsCoef, maxAbsCoef);
            // find a scalar that produces improved coefficients
            if (rowCanBeTightened(rhs, roundedRhs)) {
              // undo complementation and shifting
              undoComplementOrShift(direction, roundedRhs);
              // replace row by rounded one
              updateRow(row, direction, roundedRhs);
            }
          }
        }
      }

      // check for redundancy again
      HPRESOLVE_CHECKED_CALL(checkRowRedundant(row));
      if (rowDeleted[row]) return Result::kOk;

      auto strengthenCoefs = [&](HighsCDouble& rhs, HighsInt direction,
                                 HighsCDouble maxAbsCoefValue) {
        // iterate over non-zero positions instead of iterating over the
        // HighsMatrixSlice (provided by HPresolve::getStoredRow) because the
        // latter contains pointers to Acol and Avalue that may be invalidated
        // if these vectors are reallocated (see std::vector::push_back
        // performed in HPresolve::addToMatrix).
        for (HighsInt rowiter : rowpositions) {
          // max. absolute coefficient should not be negative
          assert(maxAbsCoefValue >= 0);

          // get column index and coefficient
          HighsInt col = Acol[rowiter];
          double val = direction * Avalue[rowiter];

          // get lower and upper bounds
          double col_lower = impliedRowBounds.getImplVarLower(row, col);
          double col_upper = impliedRowBounds.getImplVarUpper(row, col);

          // skip continuous variables
          if (model->integrality_[col] == HighsVarType::kContinuous) continue;

          if (val > maxAbsCoefValue + primal_feastol) {
            assert(col_upper != kHighsInf);
            // new matrix coefficient is direction * maxAbsCoefValue; subtract
            // existing matrix coefficient to get delta
            HighsCDouble delta = direction * (maxAbsCoefValue - val);
            addToMatrix(row, col, static_cast<double>(delta));
            rhs += delta * col_upper;
          } else if (val < -maxAbsCoefValue - primal_feastol) {
            assert(col_lower != -kHighsInf);
            // new matrix coefficient is (-direction) * maxAbsCoefValue;
            // subtract existing matrix coefficient to get delta
            HighsCDouble delta = -direction * (maxAbsCoefValue + val);
            addToMatrix(row, col, static_cast<double>(delta));
            rhs += delta * col_lower;
          }
        }
      };

      if (model->row_lower_[row] == -kHighsInf) {
        // <= constraint: try to strengthen coefficients
        double impliedRowUpper = impliedRowBounds.getSumUpper(row);
        if (impliedRowUpper != kHighsInf) {
          HighsCDouble rhs = model->row_upper_[row];
          strengthenCoefs(rhs, HighsInt{1},
                          static_cast<HighsCDouble>(impliedRowUpper) -
                              model->row_upper_[row]);
          changeRowUpper(row, static_cast<double>(rhs));
        }
      }

      if (model->row_upper_[row] == kHighsInf) {
        // >= constraint: try to strengthen coefficients
        double impliedRowLower = impliedRowBounds.getSumLower(row);
        if (impliedRowLower != -kHighsInf) {
          HighsCDouble rhs = model->row_lower_[row];
          strengthenCoefs(rhs, HighsInt{-1},
                          model->row_lower_[row] -
                              static_cast<HighsCDouble>(impliedRowLower));
          changeRowLower(row, static_cast<double>(rhs));
        }
      }
    }
  }  // if (rowsizeInteger[row] != 0 || rowsizeImplInt[row] != 0) {

  // printf("implied bounds without tightenings: [%g,%g]\n", baseiRLower,
  //        baseiRUpper);

  auto checkForcingRow = [&](HighsInt row, HighsInt direction, double rowSide,
                             double impliedRowBound,
                             const HighsCDouble& dynamism,
                             HighsPostsolveStack::RowType rowType) {
    // 1. direction =  1 (>=): forcing row if upper bound on constraint activity
    //                         is equal to row's lower bound
    // 2. direction = -1 (<=): forcing row if lower bound on constraint activity
    //                         is equal to row's upper bound
    // scale tolerance using dynamism (ratio between absolute largest and
    // smallest coefficients)
    if (direction * rowSide == -kHighsInf ||
        direction * impliedRowBound == kHighsInf ||
        abs(static_cast<HighsCDouble>(rowSide) -
            static_cast<HighsCDouble>(impliedRowBound)) >
            std::max(static_cast<double>(primal_feastol / dynamism),
                     std::numeric_limits<double>::epsilon()))
      return Result::kOk;

    // get stored row
    auto rowVector = getStoredRow();

    HighsInt nfixings = 0;
    for (const HighsSliceNonzero& nonzero : rowVector) {
      if (direction * nonzero.value() > 0) {
        if (model->col_upper_[nonzero.index()] <= implColUpper[nonzero.index()])
          ++nfixings;
      } else {
        if (model->col_lower_[nonzero.index()] >= implColLower[nonzero.index()])
          ++nfixings;
      }
    }

    if (nfixings != rowsize[row]) return Result::kOk;

    if (logging_on) analysis_.startPresolveRuleLog(kPresolveRuleForcingRow);
    postsolve_stack.forcingRow(row, rowVector, rowSide, rowType);
    // already mark the row as deleted, since otherwise it would be
    // registered as changed/singleton in the process of fixing and
    // removing the contained columns
    markRowDeleted(row);
    for (const HighsSliceNonzero& nonzero : rowVector) {
      if (direction * nonzero.value() > 0) {
        // the upper bound of the column is as tight as the implied upper
        // bound or comes from this row, which means it is not used in the
        // rows implied bounds. Therefore we can fix the variable at its
        // upper bound.
        HPRESOLVE_CHECKED_CALL(fixColToUpper(postsolve_stack, nonzero.index()));
      } else {
        HPRESOLVE_CHECKED_CALL(fixColToLower(postsolve_stack, nonzero.index()));
      }
    }
    // now the row might be empty, but not necessarily because the implied
    // column bounds might be implied by other rows in which case we
    // cannot fix the column
    postsolve_stack.redundantRow(row);
    // Row removal accounted for above

    // if there are any new row singletons, also remove them immediately
    HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolve_stack));
    analysis_.logging_on_ = logging_on;
    if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleForcingRow);
    return checkLimits(postsolve_stack);
  };

  if (this->allow_rule_[kPresolveRuleForcingRow]) {
    // Allow rule to consider forcing rows

    // store row and compute dynamism
    storeRow(row);
    HighsCDouble dynamism = computeDynamism(getStoredRow());

    // >= inequality
    HPRESOLVE_CHECKED_CALL(
        checkForcingRow(row, HighsInt{1}, model->row_lower_[row],
                        impliedRowBounds.getSumUpperOrig(row), dynamism,
                        HighsPostsolveStack::RowType::kGeq));
    if (rowDeleted[row]) return Result::kOk;

    // <= inequality
    HPRESOLVE_CHECKED_CALL(
        checkForcingRow(row, HighsInt{-1}, model->row_upper_[row],
                        impliedRowBounds.getSumLowerOrig(row), dynamism,
                        HighsPostsolveStack::RowType::kLeq));
    if (rowDeleted[row]) return Result::kOk;
  }

  // update implied bounds of all columns in given row
  HPRESOLVE_CHECKED_CALL(updateColImpliedBounds(row));

  // extract variable bound constraints
  extractVarBounds(row);

  return checkLimits(postsolve_stack);
}

HPresolve::Result HPresolve::emptyCol(HighsPostsolveStack& postsolve_stack,
                                      HighsInt col) {
  const bool logging_on = analysis_.logging_on_;
  if (logging_on) analysis_.startPresolveRuleLog(kPresolveRuleEmptyCol);
  if ((model->col_cost_[col] > 0 && model->col_lower_[col] == -kHighsInf) ||
      (model->col_cost_[col] < 0 && model->col_upper_[col] == kHighsInf)) {
    if (std::abs(model->col_cost_[col]) <= options->dual_feasibility_tolerance)
      model->col_cost_[col] = 0;
    else
      return Result::kDualInfeasible;
  }

  if (model->col_cost_[col] > 0) {
    HPRESOLVE_CHECKED_CALL(fixColToLower(postsolve_stack, col));
  } else if (model->col_cost_[col] < 0 ||
             std::abs(model->col_upper_[col]) <
                 std::abs(model->col_lower_[col])) {
    HPRESOLVE_CHECKED_CALL(fixColToUpper(postsolve_stack, col));
  } else if (model->col_lower_[col] != -kHighsInf) {
    HPRESOLVE_CHECKED_CALL(fixColToLower(postsolve_stack, col));
  } else {
    HPRESOLVE_CHECKED_CALL(fixColToZero(postsolve_stack, col));
  }

  analysis_.logging_on_ = logging_on;
  if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleEmptyCol);
  return checkLimits(postsolve_stack);
}

HPresolve::Result HPresolve::colPresolve(HighsPostsolveStack& postsolve_stack,
                                         HighsInt col, const bool timing) {
  assert(!colDeleted[col]);
  const bool logging_on = analysis_.logging_on_;

  // check bounds
  bool isFixed;
  HPRESOLVE_CHECKED_CALL(checkColBounds(col, &isFixed));
  if (isFixed) {
    // remove fixed column
    if (timing) analysis_.presolveTimerStart(kPresolveClockInitialColIsFixed);
    postsolve_stack.removedFixedCol(col, model->col_lower_[col],
                                    model->col_cost_[col],
                                    getColumnVector(col));
    HPRESOLVE_CHECKED_CALL(removeFixedCol(col));
    if (timing) analysis_.presolveTimerStop(kPresolveClockInitialColIsFixed);
    return checkLimits(postsolve_stack);
  }
  HPresolve::Result result;
  switch (colsize[col]) {
    case 0:
      if (timing) analysis_.presolveTimerStart(kPresolveClockInitialColIsEmpty);
      result = emptyCol(postsolve_stack, col);
      if (timing) analysis_.presolveTimerStop(kPresolveClockInitialColIsEmpty);
      return result;
    case 1:
      if (timing)
        analysis_.presolveTimerStart(kPresolveClockInitialColIsSingleton);
      result = singletonCol(postsolve_stack, col, timing);
      if (timing)
        analysis_.presolveTimerStop(kPresolveClockInitialColIsSingleton);
      return result;
    default:
      break;
  }

  // detect strong / weak domination
  if (timing) analysis_.presolveTimerStart(kPresolveClockInitialColDominated);
  HPRESOLVE_CHECKED_CALL(detectDominatedCol(postsolve_stack, col));
  if (timing) analysis_.presolveTimerStop(kPresolveClockInitialColDominated);
  if (colDeleted[col]) return Result::kOk;

  // column is not (weakly) dominated

  // integer columns cannot be used to tighten bounds on dual multipliers
  if (mipsolver != nullptr) {
    // lambda for changing implied row dual
    auto modifyImpliedRowDualBound = [&](HighsInt col, HighsInt row,
                                         HighsInt direction,
                                         bool isBoundImplied, HighsInt numInf) {
      if (isBoundImplied && row != -1 && numInf == 1 &&
          direction * model->col_cost_[col] >= 0 && !isRanged(row)) {
        HighsInt nzPos = findNonzero(row, col);

        if (model->integrality_[col] != HighsVarType::kInteger ||
            (rowsizeInteger[row] == rowsize[row] &&
             rowCoefficientsIntegral(row, 1.0 / Avalue[nzPos]))) {
          if (direction * Avalue[nzPos] > 0)
            changeImplRowDualLower(row, 0.0, col);
          else
            changeImplRowDualUpper(row, 0.0, col);
        }
      }
    };

    // if there is an infinite or redundant column lower bound, the reduced
    // cost of the column must not be positive (i.e. <= 0; the upper bound on
    // the reduced cost is zero).
    modifyImpliedRowDualBound(col, colLowerSource[col], HighsInt{1},
                              isLowerImplied(col),
                              impliedDualRowBounds.getNumInfSumUpperOrig(col));

    // if there is an infinite or redundant column upper bound, the reduced
    // cost of the column must not be negative (i.e. >= 0; the lower bound on
    // the reduced cost is zero).
    modifyImpliedRowDualBound(col, colUpperSource[col], HighsInt{-1},
                              isUpperImplied(col),
                              impliedDualRowBounds.getNumInfSumLowerOrig(col));

    // check if variable is implied integer
    if (timing)
      analysis_.presolveTimerStart(kPresolveClockInitialColImpliedInteger);
    HPRESOLVE_CHECKED_CALL(static_cast<Result>(convertImpliedInteger(col)));
    if (timing)
      analysis_.presolveTimerStop(kPresolveClockInitialColImpliedInteger);

    // shift "binary" variables to have a lower bound of zero
    if (model->integrality_[col] != HighsVarType::kContinuous &&
        model->col_lower_[col] != 0.0 &&
        (model->col_lower_[col] != -kHighsInf ||
         model->col_upper_[col] != kHighsInf) &&
        model->col_upper_[col] - model->col_lower_[col] > 0.5 &&
        model->col_upper_[col] - model->col_lower_[col] < 1.5) {
      // substitute with the bound that is smaller in magnitude and only
      // substitute if bound is not large for an integer
      if (std::abs(model->col_upper_[col]) > std::abs(model->col_lower_[col])) {
        if (std::abs(model->col_lower_[col]) < 1000.5)
          HPRESOLVE_CHECKED_CALL(transformColumn(postsolve_stack, col, 1.0,
                                                 model->col_lower_[col]));
      } else {
        if (std::abs(model->col_upper_[col]) < 1000.5)
          HPRESOLVE_CHECKED_CALL(transformColumn(postsolve_stack, col, -1.0,
                                                 model->col_upper_[col]));
      }
    }
  }

  // dual fixing
  if (this->allow_rule_[kPresolveRuleDualFixing]) {
    const bool logging_on = analysis_.logging_on_;
    if (logging_on) analysis_.startPresolveRuleLog(kPresolveRuleDualFixing);
    if (timing)
      analysis_.presolveTimerStart(kPresolveClockSingletonColDualFixing);
    HPRESOLVE_CHECKED_CALL(dualFixing(postsolve_stack, col));
    if (timing)
      analysis_.presolveTimerStop(kPresolveClockSingletonColDualFixing);
    analysis_.logging_on_ = logging_on;
    if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleDualFixing);
    if (colDeleted[col]) return Result::kOk;
  }

  // singleton column stuffing
  if (this->allow_rule_[kPresolveRuleColStuffing]) {
    const bool logging_on = analysis_.logging_on_;
    if (logging_on) analysis_.startPresolveRuleLog(kPresolveRuleColStuffing);
    if (timing)
      analysis_.presolveTimerStart(kPresolveClockInitialColSingletonStuffing);
    HPRESOLVE_CHECKED_CALL(singletonColStuffing(postsolve_stack, col));
    if (timing)
      analysis_.presolveTimerStop(kPresolveClockInitialColSingletonStuffing);
    analysis_.logging_on_ = logging_on;
    if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleColStuffing);
    if (colDeleted[col]) return Result::kOk;
  }

  // update dual implied bounds of all rows in given column
  if (model->integrality_[col] != HighsVarType::kInteger)
    updateRowDualImpliedBounds(col);

  return Result::kOk;
}

HPresolve::Result HPresolve::detectDominatedCol(
    HighsPostsolveStack& postsolve_stack, HighsInt col,
    bool handleSingletonRows) {
  assert(!colDeleted[col]);
  // handleSingletonRows is true by default, but set false when
  // calling detectDominatedCol in HPresolve::singletonCol
  //
  // get bounds on column dual
  double colDualUpper =
      -impliedDualRowBounds.getSumLower(col, -model->col_cost_[col]);
  double colDualLower =
      -impliedDualRowBounds.getSumUpper(col, -model->col_cost_[col]);

  // compute dynamism
  HighsCDouble dynamism = computeDynamism(getColumnVector(col));

  const bool logging_on = analysis_.logging_on_;

  auto dominatedCol = [&](HighsInt col, double dualBound, double bound,
                          HighsInt direction) {
    // column is (strongly) dominated if the bounds on the column dual satisfy:
    // 1. lower bound >  dual feasibility tolerance (direction =  1) or
    // 2. upper bound < -dual feasibility tolerance (direction = -1).
    if (direction * dualBound <= options->dual_feasibility_tolerance)
      return Result::kOk;
    // cannot fix to +-infinity -> infeasible
    if (direction * bound == -kHighsInf) return Result::kDualInfeasible;
    if (logging_on) analysis_.startPresolveRuleLog(kPresolveRuleDominatedCol);
    // fix variable
    if (direction > 0)
      HPRESOLVE_CHECKED_CALL(fixColToLower(postsolve_stack, col));
    else
      HPRESOLVE_CHECKED_CALL(fixColToUpper(postsolve_stack, col));
    analysis_.logging_on_ = logging_on;
    if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleDominatedCol);
    // handle row singletons (if requested)
    if (handleSingletonRows)
      HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolve_stack));
    return checkLimits(postsolve_stack);
  };

  auto weaklyDominatedCol = [&](HighsInt col, double dualBound, double bound,
                                double otherBound, const HighsCDouble& dynamism,
                                HighsInt direction) {
    // column is weakly dominated if the bounds on the column dual satisfy:
    // 1. lower bound >= -dual feasibility tolerance (direction =  1) or
    // 2. upper bound <=  dual feasibility tolerance (direction = -1).
    if (direction * dualBound < -options->dual_feasibility_tolerance)
      return Result::kOk;
    if (direction * bound != -kHighsInf) {
      if (logging_on)
        analysis_.startPresolveRuleLog(kPresolveRuleWeaklyDominatedCol);
      // fix variable
      if (direction > 0)
        HPRESOLVE_CHECKED_CALL(fixColToLower(postsolve_stack, col));
      else
        HPRESOLVE_CHECKED_CALL(fixColToUpper(postsolve_stack, col));
      analysis_.logging_on_ = logging_on;
      if (logging_on)
        analysis_.stopPresolveRuleLog(kPresolveRuleWeaklyDominatedCol);
      // handle row singletons (if requested)
      if (handleSingletonRows)
        HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolve_stack));
      return checkLimits(postsolve_stack);
    } else if (this->allow_rule_[kPresolveRuleForcingCol]) {
      // check for forcing column (see Andersen and Andersen, Presolving in
      // linear programming. Math. Program. 71, 221-245, 1995).
      // the column's lower bound is infinite (direction = 1) or its upper
      // bound is infinite (direction = -1).
      // now get lower bound (direction = 1) or upper bound (direction = -1) on
      // column dual using the original bounds on the row duals.
      double boundOnColDual = direction > 0
                                  ? -impliedDualRowBounds.getSumUpperOrig(
                                        col, -model->col_cost_[col])
                                  : -impliedDualRowBounds.getSumLowerOrig(
                                        col, -model->col_cost_[col]);
      if (std::abs(boundOnColDual) <=
          std::max(static_cast<double>(options->dual_feasibility_tolerance /
                                       dynamism),
                   std::numeric_limits<double>::epsilon())) {
        // 1. column dual's upper bound is zero (since the column's lower bound
        // is infinite) and column dual's lower bound is zero as well
        // (direction = 1) or
        // 2. column dual's lower bound is zero (since the column's upper bound
        // is infinite) and column dual's upper bound is zero as well
        // (direction = -1).
        // thus, the column dual is zero, and we can remove the column
        // and all its rows
        if (logging_on) analysis_.startPresolveRuleLog(kPresolveRuleForcingCol);
        postsolve_stack.forcingColumn(
            col, getColumnVector(col), model->col_cost_[col], otherBound,
            direction < 0, model->integrality_[col] == HighsVarType::kInteger);
        markColDeleted(col);
        HighsInt coliter = colhead[col];
        while (coliter != -1) {
          HighsInt row = Arow[coliter];
          double rhs = direction * Avalue[coliter] > 0.0
                           ? model->row_upper_[row]
                           : model->row_lower_[row];
          coliter = Anext[coliter];

          postsolve_stack.forcingColumnRemovedRow(col, row, rhs,
                                                  getRowVector(row));
          removeRow(row);
        }
        analysis_.logging_on_ = logging_on;
        if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleForcingCol);
        return checkLimits(postsolve_stack);
      }
    }
    return Result::kOk;
  };

  // check for dominated column
  HPRESOLVE_CHECKED_CALL(
      dominatedCol(col, colDualLower, model->col_lower_[col], HighsInt{1}));
  if (colDeleted[col]) return Result::kOk;

  HPRESOLVE_CHECKED_CALL(
      dominatedCol(col, colDualUpper, model->col_upper_[col], HighsInt{-1}));
  if (colDeleted[col]) return Result::kOk;

  if (this->allow_rule_[kPresolveRuleWeaklyDominatedCol]) {
    // check for weakly dominated column
    HPRESOLVE_CHECKED_CALL(
        weaklyDominatedCol(col, colDualLower, model->col_lower_[col],
                           model->col_upper_[col], dynamism, HighsInt{1}));
    if (colDeleted[col]) return Result::kOk;

    HPRESOLVE_CHECKED_CALL(
        weaklyDominatedCol(col, colDualUpper, model->col_upper_[col],
                           model->col_lower_[col], dynamism, HighsInt{-1}));
  }
  return Result::kOk;
}

void HPresolve::computeLocks(
    HighsInt col, bool considerObjective,
    std::function<bool(HighsInt, bool, bool)> lockCallback) const {
  // if the callback returns true, we stop examining the locks
  if (considerObjective) {
    // consider objective function
    if (model->col_cost_[col] < 0) {
      // downlock
      if (lockCallback(-1, true, false)) return;
    } else if (model->col_cost_[col] > 0) {
      // uplock
      if (lockCallback(-1, false, true)) return;
    }
  }

  // check coefficients
  for (const auto& nz : getColumnVector(col)) {
    // implied lower bound -> downlock
    bool hasDownLock = yieldsImpliedLowerBound(nz.index(), nz.value());
    // implied upper bound -> uplock
    bool hasUpLock = yieldsImpliedUpperBound(nz.index(), nz.value());
    // callback
    if (hasDownLock || hasUpLock) {
      if (lockCallback(nz.index(), hasDownLock, hasUpLock)) break;
    }
  }
}

double HPresolve::computeImpliedLowerBound(HighsInt col, HighsInt boundCol,
                                           double boundColValue,
                                           HighsInt boundColCoeffPattern) {
  double lowerBound;
  computeColBounds(col, boundCol, boundColValue, boundColCoeffPattern,
                   &lowerBound);
  return lowerBound;
}

double HPresolve::computeImpliedUpperBound(HighsInt col, HighsInt boundCol,
                                           double boundColValue,
                                           HighsInt boundColCoeffPattern) {
  double upperBound;
  computeColBounds(col, boundCol, boundColValue, boundColCoeffPattern, nullptr,
                   &upperBound);
  return upperBound;
}

double HPresolve::computeWorstCaseLowerBound(HighsInt col, HighsInt boundCol,
                                             double boundColValue,
                                             HighsInt boundColCoeffPattern) {
  double lowerBound;
  computeColBounds(col, boundCol, boundColValue, boundColCoeffPattern, nullptr,
                   nullptr, &lowerBound);
  return lowerBound;
}

double HPresolve::computeWorstCaseUpperBound(HighsInt col, HighsInt boundCol,
                                             double boundColValue,
                                             HighsInt boundColCoeffPattern) {
  double upperBound;
  computeColBounds(col, boundCol, boundColValue, boundColCoeffPattern, nullptr,
                   nullptr, nullptr, &upperBound);
  return upperBound;
}

HPresolve::Result HPresolve::checkOriginalModelBounds() {
  // Perform integer rounding of bounds on integer variables and check
  // for trivial bound violations. Only called in HPresolve::presolve,
  // and before the call to HPresolveInitialSweep::run and
  // HPresolve::initialRowAndColPresolve
  const bool is_mip = mipsolver != nullptr;
  assert(!is_mip || model->integrality_.size());
  for (HighsInt iCol = 0; iCol < model->num_col_; iCol++) {
    if (is_mip && model->integrality_[iCol] != HighsVarType::kContinuous) {
      // Perform integer rounding of bounds on integer variables
      model->col_lower_[iCol] =
          std::ceil(model->col_lower_[iCol] - primal_feastol);
      model->col_upper_[iCol] =
          std::floor(model->col_upper_[iCol] + primal_feastol);
    }
    // Check for trivial primal infeasibility
    if (model->col_lower_[iCol] > model->col_upper_[iCol])
      return Result::kPrimalInfeasible;
  }
  // Check for trivial primal infeasibility in rows
  for (HighsInt iRow = 0; iRow < model->num_row_; iRow++) {
    if (model->row_lower_[iRow] > model->row_upper_[iRow])
      return Result::kPrimalInfeasible;
  }
  return Result::kOk;
}

HPresolve::Result HPresolve::initialRowAndColPresolve(
    HighsPostsolveStack& postsolve_stack) {
  // do a full scan over the rows as the singleton arrays and the changed row
  // arrays are not initialized, also unset changedRowFlag so that the row will
  // be added to the changed row vector when it is changed after it was
  // processed
  analysis_.presolveTimerStart(kPresolveClockInitialRow);
  for (HighsInt row = 0; row != model->num_row_; ++row) {
    if (rowDeleted[row]) continue;
    HPRESOLVE_CHECKED_CALL(rowPresolve(postsolve_stack, row));
    changedRowFlag[row] = false;
  }
  analysis_.presolveTimerStop(kPresolveClockInitialRow);

  // same for the columns
  analysis_.presolveTimerStart(kPresolveClockInitialCol);
  const bool timing = analysis_.analyse_presolve_time_;
  for (HighsInt col = 0; col != model->num_col_; ++col) {
    if (colDeleted[col]) continue;
    HPRESOLVE_CHECKED_CALL(colPresolve(postsolve_stack, col, timing));
    changedColFlag[col] = false;
  }
  analysis_.presolveTimerStop(kPresolveClockInitialCol);

  return checkLimits(postsolve_stack);
}

HPresolve::Result HPresolve::fastPresolveLoop(
    HighsPostsolveStack& postsolve_stack) {
  do {
    storeCurrentProblemSize();

    analysis_.presolveTimerStart(kPresolveClockFastLoopRowSingletons);
    HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolve_stack));
    analysis_.presolveTimerStop(kPresolveClockFastLoopRowSingletons);

    analysis_.presolveTimerStart(kPresolveClockFastLoopColSingletons);
    HPRESOLVE_CHECKED_CALL(presolveChangedRows(postsolve_stack));
    analysis_.presolveTimerStop(kPresolveClockFastLoopColSingletons);

    analysis_.presolveTimerStart(kPresolveClockFastLoopDoubletonEquations);
    HPRESOLVE_CHECKED_CALL(removeDoubletonEquations(postsolve_stack));
    analysis_.presolveTimerStop(kPresolveClockFastLoopDoubletonEquations);

    analysis_.presolveTimerStart(kPresolveClockFastLoopChangedRows);
    HPRESOLVE_CHECKED_CALL(presolveColSingletons(postsolve_stack));
    analysis_.presolveTimerStop(kPresolveClockFastLoopChangedRows);

    analysis_.presolveTimerStart(kPresolveClockFastLoopChangedCols);
    HPRESOLVE_CHECKED_CALL(presolveChangedCols(postsolve_stack));
    analysis_.presolveTimerStop(kPresolveClockFastLoopChangedCols);

  } while (problemSizeReduction() > 0.01);

  return Result::kOk;
}

HPresolve::Result HPresolve::presolve(HighsPostsolveStack& postsolve_stack) {
  // For the innermost loop the rounds are nested with an outer loop
  // that layers the newer presolvers
  //
  //    fast presolve loop
  //        - empty, forcing and dominated rows and row singletons immediately
  //        after each forcing row
  //        - doubleton equations and row singletons immediately after each
  //        successful substitution
  //        - col singletons (can this introduce row singletons? If yes then
  //        immediately remove)
  //        - empty, dominated and weakly dominated columns
  //        - row singletons
  //        - if( !has enough changes ) stop
  // main loop
  //    - fast presolve loop
  //    - parallel rows and columns
  //    - if (changes found) fast presolve loop
  //    - aggregator // add limit that catches many substitutions but stops when
  //    many failures, do not run exhaustively as now
  //    - if (changes found) start main loop from beginning
  //    - primal and dual matrix sparsification
  //    - if (changes found) fast presolve loop
  //    - stop
  //

  auto presolveReturn = [&]() {
    if (mipsolver != nullptr) HPRESOLVE_CHECKED_CALL(scaleMIP(postsolve_stack));

    // analysePresolveRuleLog() should return true - no errors
    assert(analysis_.analysePresolveRuleLog());
    // Possibly report presolve log
    analysis_.analysePresolveRuleLog(true);
    return Result::kOk;
  };

  // convert model to minimization problem
  if (model->sense_ == ObjSense::kMaximize) {
    for (HighsInt i = 0; i != model->num_col_; ++i)
      model->col_cost_[i] = -model->col_cost_[i];

    model->offset_ = -model->offset_;
    assert(std::isfinite(model->offset_));
    model->sense_ = ObjSense::kMinimize;
  }

  // Need to check for time-out in checkLimits, so make sure that
  // the timer is well defined, and that its total time clock is
  // running
  assert(this->timer);
  assert(this->timer->running());

  const bool silent = silentLog();
  double report_frequency = 10;
  double current_time = this->timer->read();
  HighsInt last_report_size = model->num_col_ + model->num_row_;
  double last_report_time = current_time;
  if (options->presolve != kHighsOffString && !silent) {
    highsLogUser(options->log_options, HighsLogType::kInfo,
                 "Presolving model\n");
    std::string time_str = highsTimeSecondToString(current_time);
    if (options->timeless_log) time_str = "";
    highsLogUser(options->log_options, HighsLogType::kInfo,
                 "%" HIGHSINT_FORMAT " rows, %" HIGHSINT_FORMAT
                 " cols, %" HIGHSINT_FORMAT " nonzeros %s\n",
                 model->num_row_, model->num_col_, model->numNz(),
                 time_str.c_str());
  }

  // Perform integer rounding of bounds on integer variables and check
  // for trivial bound violations - which yield
  // Result::kPrimalInfeasible
  HPRESOLVE_CHECKED_CALL(checkOriginalModelBounds());

  if (options->presolve != kHighsOffString && mipsolver == nullptr &&
      !options->presolve_rule_test &&
      this->allow_rule_[kPresolveRuleInitialSweep]) {
    numDeletedCols = 0;
    numDeletedRows = 0;
    analysis_.presolveTimerStart(kPresolveClockInitialSweep);
    const bool logging_on = analysis_.logging_on_;
    if (logging_on) analysis_.startPresolveRuleLog(kPresolveRuleInitialSweep);
    HPresolveInitialSweep sweep(*model, *options, allow_rule_, primal_feastol);
    auto sweep_result = sweep.run(postsolve_stack);
    numDeletedCols = sweep.numDeletedCols();
    numDeletedRows = sweep.numDeletedRows();
    analysis_.logging_on_ = logging_on;
    if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleInitialSweep);
    analysis_.presolveTimerStop(kPresolveClockInitialSweep);
    if (sweep_result == HPresolveInitialSweep::Result::kPrimalInfeasible)
      return Result::kPrimalInfeasible;
    if (sweep_result == HPresolveInitialSweep::Result::kDualInfeasible)
      return Result::kDualInfeasible;
    // Check that time and reduction limits have not been reached
    //
    // NB Setting presolve_reduction_limit = 0 ensures that presolve
    // returns after initial sweep
    HPRESOLVE_CHECKED_CALL(checkLimits(postsolve_stack));
  }
  if (!okSetupPresolveDataStructures()) {
    highsLogUser(options->log_options, HighsLogType::kError,
                 "Insufficient memory for presolve data structures\n");
    // Memory allocation error
    return Result::kOutOfMemory;
  }

  // initialize substitution opportunities
  analysis_.presolveTimerStart(kPresolveClockSetupSubstitutionOpportunities);
  setupSubstitutionOpportunities();
  analysis_.presolveTimerStop(kPresolveClockSetupSubstitutionOpportunities);

  if (options->presolve != kHighsOffString) {
    auto report = [&]() {
      if (!silent) {
        HighsInt numCol = model->num_col_ - numDeletedCols;
        HighsInt numRow = model->num_row_ - numDeletedRows;
        HighsInt numNonz =
            static_cast<HighsInt>(Avalue.size() - freeslots.size());
        std::string time_str = highsTimeSecondToString(current_time);
        if (options->timeless_log) time_str = "";
        highsLogUser(options->log_options, HighsLogType::kInfo,
                     "%" HIGHSINT_FORMAT " rows, %" HIGHSINT_FORMAT
                     " cols, %" HIGHSINT_FORMAT " nonzeros %s\n",
                     numRow, numCol, numNonz, time_str.c_str());
      }
    };

    // Need to check for time-out in checkLimits, so make sure that
    // the timer is well defined, and that its total time clock is
    // running
    assert(this->timer);
    assert(this->timer->running());

    // Possibly just perform a bespoke presolve rule test
    if (options->presolve_rule_test) {
      HPRESOLVE_CHECKED_CALL(presolveRuleTest(postsolve_stack));
      return presolveReturn();
    }

    analysis_.presolveTimerStart(kPresolveClockInitial);
    HPRESOLVE_CHECKED_CALL(initialRowAndColPresolve(postsolve_stack));
    analysis_.presolveTimerStop(kPresolveClockInitial);

    HighsInt numParallelRowColCalls = 0;
    // ReductionType::kEqualityRowAddition(s) has no basis postsolve,
    // so can only be used when basis postsolve is not required - when
    // the problem is a MIP, IPM is run without crossover, or when
    // PDLP is used. However, if the LP is reduced to empty, the basis
    // must not be formed in the case of IPM without crossover or PDLP
    bool trySparsify = !may_require_basis_postsolve_;
    bool tryProbing = mipsolver != nullptr;
    bool tryFourierMotzkin = mipsolver != nullptr;
    HighsInt numCliquesBeforeProbing = -1;
    bool domcolAfterProbingCalled = false;
    bool dependentEquationsCalled = mipsolver != nullptr;

    // Start of main presolve loop
    //
    while (true) {
      assert(model->num_col_ >= numDeletedCols);
      assert(model->num_row_ >= numDeletedRows);
      HighsInt current_size =
          model->num_col_ - numDeletedCols + model->num_row_ - numDeletedRows;
      if (options->output_flag) {
        current_time = this->timer->read();
        if (current_size < 0.85 * last_report_size ||
            current_time > last_report_time + report_frequency) {
          last_report_size = current_size;
          last_report_time = current_time;
          report();
        }
      }

      analysis_.presolveTimerStart(kPresolveClockFastLoop);
      HPRESOLVE_CHECKED_CALL(fastPresolveLoop(postsolve_stack));
      analysis_.presolveTimerStop(kPresolveClockFastLoop);

      storeCurrentProblemSize();

      // when presolving after a restart the clique table and implication
      // structure may contain substitutions which we apply directly before
      // running the aggregator as they might lose validity otherwise
      if (mipsolver != nullptr) {
        HighsInt numDelCol = 0;
        HPRESOLVE_CHECKED_CALL(
            applyConflictGraphSubstitutions(postsolve_stack, numDelCol));
      }

      HighsInt numColsEliminatedFourierMotzkin = 0;
      if (tryFourierMotzkin && this->allow_rule_[kPresolveRuleFourierMotzkin]) {
        HPRESOLVE_CHECKED_CALL(
            fourierMotzkin(postsolve_stack, numColsEliminatedFourierMotzkin));
        tryFourierMotzkin = false;
      }

      if (reducedToEmpty()) break;

      if (this->allow_rule_[kPresolveRuleAggregator]) {
        analysis_.presolveTimerStart(kPresolveClockAggregator);
        HPRESOLVE_CHECKED_CALL(aggregator(postsolve_stack));
        analysis_.presolveTimerStop(kPresolveClockAggregator);
      }

      if (problemSizeReduction() > 0.05) continue;

      if (trySparsify && this->allow_rule_[kPresolveRuleSparsify]) {
        HighsInt numNz = numNonzeros();
        analysis_.presolveTimerStart(kPresolveClockSparsify);
        HPRESOLVE_CHECKED_CALL(sparsify(postsolve_stack));
        analysis_.presolveTimerStop(kPresolveClockSparsify);
        double nzReduction =
            100.0 * (1.0 - (numNonzeros() / static_cast<double>(numNz)));

        if (nzReduction > 0) {
          highsLogDev(options->log_options, HighsLogType::kInfo,
                      "Sparsify removed %.1f%% of nonzeros\n", nzReduction);

          analysis_.presolveTimerStart(kPresolveClockFastLoop);
          HPRESOLVE_CHECKED_CALL(fastPresolveLoop(postsolve_stack));
          analysis_.presolveTimerStop(kPresolveClockFastLoop);
        }
        trySparsify = false;
      }

      if (this->allow_rule_[kPresolveRuleParallelRowsAndCols] &&
          numParallelRowColCalls < 5) {
        if (shrinkProblemEnabled && (numDeletedCols >= model->num_col_ / 2 ||
                                     numDeletedRows >= model->num_row_ / 2)) {
          //	analysis_.presolveTimerStart(kPresolveClock@);
          //	analysis_.presolveTimerStop(kPresolveClock@);
          analysis_.presolveTimerStart(kPresolveClockShrinkProblem);
          shrinkProblem(postsolve_stack);
          analysis_.presolveTimerStop(kPresolveClockShrinkProblem);

          toCSC(model->a_matrix_.value_, model->a_matrix_.index_,
                model->a_matrix_.start_);
          okFromCSC(model->a_matrix_.value_, model->a_matrix_.index_,
                    model->a_matrix_.start_);
        }
        storeCurrentProblemSize();
        analysis_.presolveTimerStart(kPresolveClockParallelRowsAndCols);
        HPRESOLVE_CHECKED_CALL(detectParallelRowsAndCols(postsolve_stack));
        analysis_.presolveTimerStop(kPresolveClockParallelRowsAndCols);
        ++numParallelRowColCalls;
        if (problemSizeReduction() > 0.05) continue;
      }

      if (postsolve_stack.numReductions() == 35039) {
      }
      analysis_.presolveTimerStart(kPresolveClockFastLoop);
      HPRESOLVE_CHECKED_CALL(fastPresolveLoop(postsolve_stack));
      analysis_.presolveTimerStop(kPresolveClockFastLoop);

      if (mipsolver != nullptr) {
        HighsInt num_strengthened = -1;
        HPRESOLVE_CHECKED_CALL(
            strengthenInequalities(postsolve_stack, num_strengthened));
        assert(num_strengthened >= 0);
        if (num_strengthened > 0)
          highsLogDev(options->log_options, HighsLogType::kInfo,
                      "Strengthened %" HIGHSINT_FORMAT " coefficients\n",
                      num_strengthened);
      }

      analysis_.presolveTimerStart(kPresolveClockFastLoop);
      HPRESOLVE_CHECKED_CALL(fastPresolveLoop(postsolve_stack));
      analysis_.presolveTimerStop(kPresolveClockFastLoop);

      if (mipsolver != nullptr && numCliquesBeforeProbing == -1) {
        numCliquesBeforeProbing = mipsolver->mipdata_->cliquetable.numCliques();
        storeCurrentProblemSize();
        HPRESOLVE_CHECKED_CALL(dominatedColumns(postsolve_stack));
        if (problemSizeReduction() > 0.0)
          HPRESOLVE_CHECKED_CALL(fastPresolveLoop(postsolve_stack));
        if (problemSizeReduction() > 0.05) continue;
      }

      // enumerate solutions
      if (mipsolver != nullptr && this->allow_rule_[kPresolveRuleEnumeration]) {
        storeCurrentProblemSize();
        HPRESOLVE_CHECKED_CALL(enumerateSolutions(postsolve_stack));
        if (problemSizeReduction() > 0.05) continue;
      }

      if (tryProbing && this->allow_rule_[kPresolveRuleProbing]) {
        HPRESOLVE_CHECKED_CALL(detectImpliedIntegers());
        storeCurrentProblemSize();
        HPRESOLVE_CHECKED_CALL(runProbing(postsolve_stack));
        tryProbing = probingContingent > numProbed &&
                     (problemSizeReduction() > 1.0 || probingEarlyAbort);
        trySparsify = true;
        if (problemSizeReduction() > 0.05 || tryProbing) continue;
        HPRESOLVE_CHECKED_CALL(fastPresolveLoop(postsolve_stack));
      }

      if (!dependentEquationsCalled) {
        if (shrinkProblemEnabled && (numDeletedCols >= model->num_col_ / 2 ||
                                     numDeletedRows >= model->num_row_ / 2)) {
          analysis_.presolveTimerStart(kPresolveClockShrinkProblem);
          shrinkProblem(postsolve_stack);
          analysis_.presolveTimerStop(kPresolveClockShrinkProblem);

          toCSC(model->a_matrix_.value_, model->a_matrix_.index_,
                model->a_matrix_.start_);
          okFromCSC(model->a_matrix_.value_, model->a_matrix_.index_,
                    model->a_matrix_.start_);
        }
        storeCurrentProblemSize();
        if (this->allow_rule_[kPresolveRuleDependentEquations]) {
          analysis_.presolveTimerStart(kPresolveClockDependentEquations);
          HPRESOLVE_CHECKED_CALL(removeDependentEquations(postsolve_stack));
          analysis_.presolveTimerStop(kPresolveClockDependentEquations);
          dependentEquationsCalled = true;
        }
        if (this->allow_rule_[kPresolveRuleDependentFreeCols]) {
          analysis_.presolveTimerStart(kPresolveClockDependentFreeCol);
          HPRESOLVE_CHECKED_CALL(removeDependentFreeCols(postsolve_stack));
          analysis_.presolveTimerStop(kPresolveClockDependentFreeCol);
        }
        if (problemSizeReduction() > 0.05) continue;
      }

      if (mipsolver != nullptr &&
          mipsolver->mipdata_->cliquetable.numCliques() >
              numCliquesBeforeProbing &&
          !domcolAfterProbingCalled) {
        domcolAfterProbingCalled = true;
        storeCurrentProblemSize();
        HPRESOLVE_CHECKED_CALL(dominatedColumns(postsolve_stack));
        if (problemSizeReduction() > 0.0)
          HPRESOLVE_CHECKED_CALL(fastPresolveLoop(postsolve_stack));
        if (problemSizeReduction() > 0.05) continue;
      }

      break;
    }

    if (!reducedToEmpty()) report();
  } else {
    highsLogUser(options->log_options, HighsLogType::kInfo,
                 "\nPresolve is switched off\n");
  }

  return presolveReturn();
}

HPresolve::Result HPresolve::checkTimeLimit() {
  assert(timer);
  if (options->time_limit < kHighsInf && timer->read() >= options->time_limit)
    return Result::kStopped;
  return Result::kOk;
}

HPresolve::Result HPresolve::checkLimits(HighsPostsolveStack& postsolve_stack) {
  size_t numreductions = postsolve_stack.numReductions();

  bool debug_report = false;
  HighsInt check_col = debugGetCheckCol();
  HighsInt check_row = debugGetCheckRow();
  bool col_bound_change = false;
  bool row_bound_change = false;
  if (check_col >= 0 || check_row >= 0) {
    if (check_col >= 0) {
      col_bound_change =
          numreductions == 1 ||
          postsolve_stack.debug_prev_col_lower !=
              model->col_lower_[check_col] ||
          postsolve_stack.debug_prev_col_upper != model->col_upper_[check_col];
      postsolve_stack.debug_prev_col_lower = model->col_lower_[check_col];
      postsolve_stack.debug_prev_col_upper = model->col_upper_[check_col];
    }
    if (check_row >= 0) {
      row_bound_change =
          numreductions == 1 ||
          postsolve_stack.debug_prev_row_lower !=
              model->row_lower_[check_row] ||
          postsolve_stack.debug_prev_row_upper != model->row_upper_[check_row];
      postsolve_stack.debug_prev_row_lower = model->row_lower_[check_row];
      postsolve_stack.debug_prev_row_upper = model->row_upper_[check_row];
    }
    debug_report = numreductions > postsolve_stack.debug_prev_numreductions;
  }
  if (check_col >= 0 && col_bound_change && debug_report) {
    printf("After reduction %4d: col = %4d[%s] has bounds [%11.4g, %11.4g]\n",
           static_cast<int>(numreductions - 1), static_cast<int>(check_col),
           model->col_names_[check_col].c_str(), model->col_lower_[check_col],
           model->col_upper_[check_col]);
    postsolve_stack.debug_prev_numreductions = numreductions;
  }
  if (check_row >= 0 && row_bound_change && debug_report) {
    printf("After reduction %4d: row = %4d[%s] has bounds [%11.4g, %11.4g]\n",
           static_cast<int>(numreductions - 1), static_cast<int>(check_row),
           model->row_names_[check_row].c_str(), model->row_lower_[check_row],
           model->row_upper_[check_row]);
    postsolve_stack.debug_prev_numreductions = numreductions;
  }

  if ((numreductions & 1023u) == 0) HPRESOLVE_CHECKED_CALL(checkTimeLimit());

  // Record the value of numreductions so that the next call with an
  // increase in numreductions can be identified
  this->last_reduction_ = numreductions;
  return numreductions >= this->reductionLimit ? Result::kStopped : Result::kOk;
}

void HPresolve::storeCurrentProblemSize() {
  oldNumCol = model->num_col_ - numDeletedCols;
  oldNumRow = model->num_row_ - numDeletedRows;
}

double HPresolve::problemSizeReduction() const {
  double colReduction =
      100.0 *
      static_cast<double>(oldNumCol - (model->num_col_ - numDeletedCols)) /
      oldNumCol;
  double rowReduction =
      100.0 *
      static_cast<double>(oldNumRow - (model->num_row_ - numDeletedRows)) /
      oldNumRow;

  return std::max(rowReduction, colReduction);
}

void HPresolve::computeColBounds(HighsInt col, HighsInt boundCol,
                                 double boundColValue,
                                 HighsInt boundColCoeffPattern,
                                 double* lowerBound, double* upperBound,
                                 double* worstCaseLowerBound,
                                 double* worstCaseUpperBound) {
  assert(!colDeleted[col]);
  assert(boundCol == -1 || !colDeleted[boundCol]);
  assert(col != boundCol);

  // return if nothing to do
  if (lowerBound == nullptr && upperBound == nullptr &&
      worstCaseLowerBound == nullptr && worstCaseUpperBound == nullptr)
    return;

  // lambda for skipping non-zeros
  auto skipNonZero = [&](HighsInt row, double val) {
    bool lowerRequested =
        lowerBound != nullptr || worstCaseLowerBound != nullptr;
    bool upperRequested =
        upperBound != nullptr || worstCaseUpperBound != nullptr;
    bool hasLower = yieldsImpliedLowerBound(row, val);
    bool hasUpper = yieldsImpliedUpperBound(row, val);
    return ((!upperRequested && lowerRequested && !hasLower) ||
            (!lowerRequested && upperRequested && !hasUpper) ||
            (upperRequested && lowerRequested && !hasLower && !hasUpper));
  };

  // struct for storing non-zeros
  struct nonZeros {
    HighsInt row;
    double jval;
    double kval;
  };

  std::vector<nonZeros> nzs;
  nzs.reserve(colsize[col]);

  auto storeTriplet = [&](HighsInt row, double jval, double kval) {
    // are coefficients required to have the same or opposite signs?
    // otherColCoeffPattern =  0: accept any coefficient pair
    // otherColCoeffPattern =  1: accept coefficients with same signs
    // otherColCoeffPattern = -1: accept coefficients with opposite signs
    if (boundColCoeffPattern != 0 &&
        std::signbit(jval) != std::signbit(boundColCoeffPattern * kval))
      return;
    nzs.push_back({row, jval, kval});
  };

  // store triplets (row, nonzero, nonzero) in a vector to speed up bound
  // computation
  if (boundCol != -1) {
    if (colsize[col] < colsize[boundCol]) {
      for (const auto& colNz : getColumnVector(col)) {
        // skip non-zero if it does not yield requested bounds
        if (skipNonZero(colNz.index(), colNz.value())) continue;
        HighsInt nzPos = findNonzero(colNz.index(), boundCol);
        if (nzPos == -1) continue;
        storeTriplet(colNz.index(), colNz.value(), Avalue[nzPos]);
      }
    } else {
      for (const auto& colNz : getColumnVector(boundCol)) {
        HighsInt nzPos = findNonzero(colNz.index(), col);
        if (nzPos == -1) continue;
        // skip non-zero if it does not yield requested bounds
        if (skipNonZero(colNz.index(), Avalue[nzPos])) continue;
        storeTriplet(colNz.index(), Avalue[nzPos], colNz.value());
      }
    }
  } else {
    for (const auto& colNz : getColumnVector(col)) {
      // skip non-zero if it does not yield requested bounds
      if (skipNonZero(colNz.index(), colNz.value())) continue;
      nzs.push_back({colNz.index(), colNz.value(), -kHighsInf});
    }
  }

  // initialise bounds
  if (lowerBound != nullptr) *lowerBound = -kHighsInf;
  if (upperBound != nullptr) *upperBound = kHighsInf;
  if (worstCaseLowerBound != nullptr) *worstCaseLowerBound = -kHighsInf;
  if (worstCaseUpperBound != nullptr) *worstCaseUpperBound = kHighsInf;

  // lambda for actual bound computation
  auto computeBound = [&](const nonZeros& triplet, double rhs,
                          HighsInt direction, bool isWorstCaseBound) {
    HighsCDouble residual;
    if ((direction > 0 && !isWorstCaseBound) ||
        (direction < 0 && isWorstCaseBound)) {
      residual = isWorstCaseBound ? impliedRowBounds.getResidualSumLower(
                                        triplet.row, col, triplet.jval,
                                        boundCol, triplet.kval, boundColValue)
                                  : impliedRowBounds.getResidualSumLowerOrig(
                                        triplet.row, col, triplet.jval,
                                        boundCol, triplet.kval, boundColValue);
      if (residual == -kHighsInf) return std::copysign(kHighsInf, triplet.jval);
    } else {
      residual = isWorstCaseBound ? impliedRowBounds.getResidualSumUpper(
                                        triplet.row, col, triplet.jval,
                                        boundCol, triplet.kval, boundColValue)
                                  : impliedRowBounds.getResidualSumUpperOrig(
                                        triplet.row, col, triplet.jval,
                                        boundCol, triplet.kval, boundColValue);
      if (residual == kHighsInf) return -std::copysign(kHighsInf, triplet.jval);
    }
    return static_cast<double>((static_cast<HighsCDouble>(rhs) - residual) /
                               triplet.jval);
  };

  // lambda for updating tightest bounds
  auto updateBounds = [&](const nonZeros& triplet, double rhs,
                          HighsInt direction) {
    if (direction * rhs == kHighsInf) return;
    if (direction * triplet.jval < 0) {
      // lower bounds
      if (lowerBound != nullptr)
        *lowerBound =
            std::max(*lowerBound, computeBound(triplet, rhs, direction, false));
      if (worstCaseLowerBound != nullptr && *worstCaseLowerBound != kHighsInf)
        *worstCaseLowerBound = std::max(
            *worstCaseLowerBound, computeBound(triplet, rhs, direction, true));
    } else {
      // upper bounds
      if (upperBound != nullptr)
        *upperBound =
            std::min(*upperBound, computeBound(triplet, rhs, direction, false));
      if (worstCaseUpperBound != nullptr && *worstCaseUpperBound != -kHighsInf)
        *worstCaseUpperBound = std::min(
            *worstCaseUpperBound, computeBound(triplet, rhs, direction, true));
    }
  };

  // compute bounds
  for (const auto& triplet : nzs) {
    // compute bounds using right-hand side (direction = 1) and left-hand side
    // (direction = -1)
    updateBounds(triplet, model->row_upper_[triplet.row], HighsInt{1});
    updateBounds(triplet, model->row_lower_[triplet.row], HighsInt{-1});
  }
}

bool HPresolve::silentLog() const {
  return mipsolver && mipsolver->mipdata_->numRestarts > 0;
}

void HPresolve::moveCutsToPool(HighsPostsolveStack& postsolve_stack) {
  if (mipsolver == nullptr) return;
  std::vector<HighsInt> cutinds;
  std::vector<double> cutvals;
  cutinds.reserve(model->num_col_);
  cutvals.reserve(model->num_col_);
  HighsInt numcuts = 0;
  for (HighsInt i : postsolve_stack.getCutRows()) {
    ++numcuts;
    storeRow(i);
    cutinds.clear();
    cutvals.clear();
    for (HighsInt j : rowpositions) {
      cutinds.push_back(Acol[j]);
      cutvals.push_back(Avalue[j]);
    }

    if (mipsolver != nullptr) {
      mipsolver->mipdata_->getCutPool().addCut(
          *mipsolver, cutinds.data(), cutvals.data(), cutinds.size(),
          model->row_upper_[i],
          rowsizeInteger[i] + rowsizeImplInt[i] == rowsize[i] &&
              rowCoefficientsIntegral(i, 1.0),
          true, false, false);
    }

    markRowDeleted(i);
    for (HighsInt j : rowpositions) unlink(j);
  }

  if (numcuts == 0) return;

  // Compact deleted cut rows. shrinkProblem must not be used here
  // because it replaces the cutpool with a new empty one, destroying
  // the cuts that were just added above.
  HighsInt oldNumRow = model->num_row_;
  std::vector<HighsInt> newRowIndex(oldNumRow);
  HighsInt newNumRow = 0;
  for (HighsInt i = 0; i < oldNumRow; ++i) {
    if (rowDeleted[i])
      newRowIndex[i] = -1;
    else
      newRowIndex[i] = newNumRow++;
  }
  model->num_row_ = newNumRow;

  for (HighsInt i = 0; i < oldNumRow; ++i) {
    if (newRowIndex[i] == -1 || newRowIndex[i] == i) continue;
    model->row_lower_[newRowIndex[i]] = model->row_lower_[i];
    model->row_upper_[newRowIndex[i]] = model->row_upper_[i];
  }
  model->row_lower_.resize(model->num_row_);
  model->row_upper_.resize(model->num_row_);
  model->row_names_.resize(model->num_row_);

  for (size_t i = 0; i < Avalue.size(); ++i) {
    if (Avalue[i] == 0) continue;
    assert(newRowIndex[Arow[i]] != -1);
    Arow[i] = newRowIndex[Arow[i]];
  }

  mipsolver->mipdata_->implications.reindexRows(model->num_row_, newRowIndex);
}

HighsModelStatus HPresolve::run(HighsPostsolveStack& postsolve_stack) {
  presolve_status_ = HighsPresolveStatus::kNotSet;
  shrinkProblemEnabled = true;
  postsolve_stack.debug_prev_numreductions = 0;
  postsolve_stack.debug_prev_col_lower = 0;
  postsolve_stack.debug_prev_col_upper = 0;
  postsolve_stack.debug_prev_row_lower = 0;
  postsolve_stack.debug_prev_row_upper = 0;
  // Presolve should only be called with a model that has a non-empty
  // constraint matrix unless it has no rows
  assert(model->numNz() || model->num_row_ == 0);
  auto reportReductions = [&]() {
    if (options->presolve != kHighsOffString &&
        reductionLimit < kHighsSize_tInf) {
      highsLogUser(options->log_options, HighsLogType::kInfo,
                   "Presolve performed %" PRId64 " of %" PRId64
                   " permitted reductions\n",
                   postsolve_stack.numReductions(), reductionLimit);
    }
  };
  auto reportProfiling = [&]() {
    // Presolve profiling not currently enabled for MIP
    this->analysis_.presolveTimerStop(kPresolveClockPresolve);
    this->analysis_.reportPresolveTimer();
  };

  Result result;
  try {
    result = presolve(postsolve_stack);
  } catch (const std::exception& exception) {
    highsLogDev(options->log_options, HighsLogType::kError,
                "Exception %s in Presolve::presolve\n", exception.what());
    result = handleExceptionIsOom(options->log_options, "presolve", exception)
                 ? Result::kOutOfMemory
                 : Result::kException;
  }
  // Stop any presolve rule logging that is currently running, check
  // the presolve rule logging for errors, and analyse it
  analysis_.stopPresolveRuleLog();
  assert(analysis_.analysePresolveRuleLog());
  analysis_.analysePresolveRuleLog(true);
  switch (result) {
    case Result::kOk:
    case Result::kStopped:
      break;
    case Result::kPrimalInfeasible:
      presolve_status_ = HighsPresolveStatus::kInfeasible;
      reportReductions();
      reportProfiling();
      return HighsModelStatus::kInfeasible;
    case Result::kDualInfeasible:
      presolve_status_ = HighsPresolveStatus::kUnboundedOrInfeasible;
      reportReductions();
      reportProfiling();
      return HighsModelStatus::kUnboundedOrInfeasible;
    case Result::kOutOfMemory:
      presolve_status_ = HighsPresolveStatus::kOutOfMemory;
      return HighsModelStatus::kMemoryLimit;
    case Result::kException:
      presolve_status_ = HighsPresolveStatus::kException;
      return HighsModelStatus::kSolveError;
  }
  assert(result == Result::kOk || result == Result::kStopped);
  // Result::kStopped corresponds to reaching the time or reduction
  // limit, in which case any reductions performed are retained, so
  // complete presolve as if it had run to completion
  reportReductions();
  shrinkProblem(postsolve_stack);

  if (mipsolver != nullptr) {
    mipsolver->mipdata_->cliquetable.setMaxEntries(numNonzeros());
    mipsolver->mipdata_->getDomain().addCutpool(
        mipsolver->mipdata_->getCutPool());
    mipsolver->mipdata_->getDomain().addConflictPool(
        mipsolver->mipdata_->getConflictPool());

    if (mipsolver->mipdata_->numRestarts != 0) moveCutsToPool(postsolve_stack);
  }

  // Possibly populate the model matrix from the presolve matrix data
  // structure
  if (hasPresolveDataStructures())
    toCSC(model->a_matrix_.value_, model->a_matrix_.index_,
          model->a_matrix_.start_);

  reportProfiling();

  if (model->num_col_ == 0) {
    // Reduced to empty
    if (mipsolver) {
      if (model->offset_ > mipsolver->mipdata_->upper_limit) {
        presolve_status_ = HighsPresolveStatus::kInfeasible;
        return HighsModelStatus::kInfeasible;
      }
      mipsolver->mipdata_->lower_bound = 0;
    } else {
      // An LP with no columns must have no rows, unless the reduction
      // limit has been reached. As exposed when studying 2326
      const bool num_row_ok = model->num_row_ == 0 ||
                              postsolve_stack.numReductions() >= reductionLimit;
      assert(num_row_ok);
      if (!num_row_ok) {
        presolve_status_ = HighsPresolveStatus::kNotPresolved;
        return HighsModelStatus::kNotset;
      }
    }
    presolve_status_ = HighsPresolveStatus::kReducedToEmpty;
    // Make sure that zero row activity from the columnless model is
    // consistent with the bounds
    return zeroRowActivityFeasible() ? HighsModelStatus::kOptimal
                                     : HighsModelStatus::kInfeasible;
  } else if (postsolve_stack.numReductions() > 0) {
    // Reductions performed
    presolve_status_ = HighsPresolveStatus::kReduced;
  } else {
    // No reductions performed
    presolve_status_ = HighsPresolveStatus::kNotReduced;
  }

  assert(presolve_status_ != HighsPresolveStatus::kNotSet);
  return HighsModelStatus::kNotset;
}

HPresolve::Result HPresolve::substitute(HighsInt substcol, HighsInt staycol,
                                        double offset, double scale,
                                        HighsInt row) {
  HPRESOLVE_CHECKED_CALL(
      updateCliqueTableSubstituteCol(substcol, staycol, offset, scale));
  markColDeleted(substcol);
  if (row != -1) removeRow(row);

  // Preserve explicit integrality, i.e., upgrade implied integral
  // column if it is substituting an integral column
  if (model->integrality_[substcol] == HighsVarType::kInteger &&
      model->integrality_[staycol] == HighsVarType::kImplicitInteger) {
    model->integrality_[staycol] = HighsVarType::kInteger;
    for (const HighsSliceNonzero& nonzero : getColumnVector(staycol)) {
      ++rowsizeInteger[nonzero.index()];
      --rowsizeImplInt[nonzero.index()];
    }
  }

  // substitute the column in each row where it occurs
  for (HighsInt coliter = colhead[substcol]; coliter != -1;) {
    HighsInt colrow = Arow[coliter];
    double colval = Avalue[coliter];
    // walk to the next position before doing any modifications, because
    // the current position will be deleted in the loop below
    assert(Acol[coliter] == substcol);
    HighsInt colpos = coliter;
    coliter = Anext[coliter];
    assert(!rowDeleted[colrow]);
    unlink(colpos);

    // adjust the sides
    addToRowLower(colrow, -static_cast<HighsCDouble>(colval) * offset);
    addToRowUpper(colrow, -static_cast<HighsCDouble>(colval) * offset);

    addToMatrix(colrow, staycol, scale * colval);
    // printf("after substitution: ");
    // debugPrintRow(colrow);

    // check if this is an equation row and it now has a different size
    reinsertEquation(colrow);
  }

  // substitute column in the objective function
  if (model->col_cost_[substcol] != 0.0) {
    model->offset_ += model->col_cost_[substcol] * offset;
    assert(std::isfinite(model->offset_));

    model->col_cost_[staycol] += scale * model->col_cost_[substcol];

    if (std::abs(model->col_cost_[staycol]) <= options->small_matrix_value)
      model->col_cost_[staycol] = 0.0;
    model->col_cost_[substcol] = 0.0;
  }

  return Result::kOk;
}

HPresolve::Result HPresolve::fixColToLower(HighsPostsolveStack& postsolve_stack,
                                           HighsInt col) {
  double fixval = model->col_lower_[col];
  if (fixval == -kHighsInf) return Result::kDualInfeasible;

  // printf("fixing column %" HIGHSINT_FORMAT " to %.15g\n", col, fixval);

  const bool logging_on = analysis_.logging_on_;
  if (logging_on) analysis_.startPresolveRuleLog(kPresolveRuleFixedCol);
  postsolve_stack.fixedColAtLower(col, fixval, model->col_cost_[col],
                                  getColumnVector(col));
  Result result = removeFixedCol(col, fixval);
  analysis_.logging_on_ = logging_on;
  if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleFixedCol);
  return result;
}

HPresolve::Result HPresolve::fixColToUpper(HighsPostsolveStack& postsolve_stack,
                                           HighsInt col) {
  double fixval = model->col_upper_[col];
  if (fixval == kHighsInf) return Result::kDualInfeasible;

  // printf("fixing column %" HIGHSINT_FORMAT " to %.15g\n", col, fixval);

  const bool logging_on = analysis_.logging_on_;
  if (logging_on) analysis_.startPresolveRuleLog(kPresolveRuleFixedCol);
  postsolve_stack.fixedColAtUpper(col, fixval, model->col_cost_[col],
                                  getColumnVector(col));
  Result result = removeFixedCol(col, fixval);
  analysis_.logging_on_ = logging_on;
  if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleFixedCol);
  return result;
}

HPresolve::Result HPresolve::fixColToZero(HighsPostsolveStack& postsolve_stack,
                                          HighsInt col) {
  const bool logging_on = analysis_.logging_on_;
  if (logging_on) analysis_.startPresolveRuleLog(kPresolveRuleFixedCol);
  postsolve_stack.fixedColAtZero(col, model->col_cost_[col],
                                 getColumnVector(col));
  Result result = removeFixedCol(col, 0.0);
  analysis_.logging_on_ = logging_on;
  if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleFixedCol);
  return result;
}

void HPresolve::unlinkRow(HighsInt row) {
  assert(row >= 0);
  assert(static_cast<size_t>(row) < rowroot.size());
  storeRow(row);
  for (HighsInt rowiter : rowpositions) {
    assert(Arow[rowiter] == row);
    unlink(rowiter);
  }
}

void HPresolve::removeRow(HighsInt row) {
  assert(row >= 0);
  assert(static_cast<size_t>(row) < rowroot.size());
  // first mark the row as logically deleted, so that it is not register as
  // singleton row upon removing its nonzeros
  markRowDeleted(row);
  unlinkRow(row);
}

HPresolve::Result HPresolve::removeFixedCol(HighsInt col) {
  const bool logging_on = analysis_.logging_on_;
  if (logging_on) analysis_.startPresolveRuleLog(kPresolveRuleFixedCol);
  const Result result = removeFixedCol(col, model->col_lower_[col]);
  analysis_.logging_on_ = logging_on;
  if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleFixedCol);
  return result;
}

HPresolve::Result HPresolve::removeFixedCol(HighsInt col, double fixval) {
  HPRESOLVE_CHECKED_CALL(updateCliqueTableFixedCol(col, fixval));
  // mark the column as deleted first so that it is not registered as singleton
  // column upon removing its non-zeros
  markColDeleted(col);

  for (HighsInt coliter = colhead[col]; coliter != -1;) {
    HighsInt colrow = Arow[coliter];
    double colval = Avalue[coliter];
    assert(Acol[coliter] == col);

    HighsInt colpos = coliter;
    coliter = Anext[coliter];

    addToRowLower(colrow, -static_cast<HighsCDouble>(colval) * fixval);
    addToRowUpper(colrow, -static_cast<HighsCDouble>(colval) * fixval);

    unlink(colpos);

    reinsertEquation(colrow);
  }
  model->offset_ += model->col_cost_[col] * fixval;
  assert(std::isfinite(model->offset_));
  model->col_cost_[col] = 0;
  return Result::kOk;
}

HPresolve::Result HPresolve::removeRowSingletons(
    HighsPostsolveStack& postsolve_stack) {
  for (size_t i = 0; i < singletonRows.size(); ++i) {
    HighsInt row = singletonRows[i];
    if (rowDeleted[row] || rowsize[row] > 1) continue;
    // row presolve will delegate to rowSingleton() if the row size is 1
    // if the singleton row has become empty it will also remove the row
    HPRESOLVE_CHECKED_CALL(rowPresolve(postsolve_stack, row));
  }

  singletonRows.clear();

  return Result::kOk;
}

HPresolve::Result HPresolve::presolveColSingletons(
    HighsPostsolveStack& postsolve_stack) {
  for (size_t i = 0; i < singletonColumns.size(); ++i) {
    HighsInt col = singletonColumns[i];
    if (colDeleted[col]) continue;
    HPRESOLVE_CHECKED_CALL(colPresolve(postsolve_stack, col));
  }
  singletonColumns.erase(
      std::remove_if(
          singletonColumns.begin(), singletonColumns.end(),
          [&](HighsInt col) { return colDeleted[col] || colsize[col] > 1; }),
      singletonColumns.end());

  return Result::kOk;
}

HPresolve::Result HPresolve::presolveChangedRows(
    HighsPostsolveStack& postsolve_stack) {
  std::vector<HighsInt> changedRows;
  changedRows.reserve(model->num_row_ - numDeletedRows);
  changedRows.swap(changedRowIndices);
  for (HighsInt row : changedRows) {
    if (rowDeleted[row]) continue;
    HPRESOLVE_CHECKED_CALL(rowPresolve(postsolve_stack, row));
    changedRowFlag[row] = rowDeleted[row];
  }

  return Result::kOk;
}

HPresolve::Result HPresolve::presolveChangedCols(
    HighsPostsolveStack& postsolve_stack) {
  std::vector<HighsInt> changedCols;
  changedCols.reserve(model->num_col_ - numDeletedCols);
  changedCols.swap(changedColIndices);
  for (HighsInt col : changedCols) {
    if (colDeleted[col]) continue;
    HPRESOLVE_CHECKED_CALL(colPresolve(postsolve_stack, col));
    changedColFlag[col] = colDeleted[col];
  }

  return Result::kOk;
}

HPresolve::Result HPresolve::removeDoubletonEquations(
    HighsPostsolveStack& postsolve_stack) {
  auto eq = equations.begin();
  while (eq != equations.end()) {
    HighsInt eqrow = eq->second;
    assert(!rowDeleted[eqrow]);
    assert(eq->first == rowsize[eqrow]);
    assert(isEquation(eqrow));
    if (rowsize[eqrow] > 2) return Result::kOk;
    HPRESOLVE_CHECKED_CALL(rowPresolve(postsolve_stack, eqrow));
    if (rowDeleted[eqrow])
      eq = equations.begin();
    else
      ++eq;
  }

  return Result::kOk;
}

HPresolve::Result HPresolve::strengthenInequalities(
    HighsPostsolveStack& postsolve_stack, HighsInt& num_strengthened) {
  std::vector<int8_t> complementation;
  std::vector<HighsCDouble> reducedcost;
  std::vector<HighsCDouble> upper;
  std::vector<HighsInt> indices;
  std::vector<HighsInt> positions;
  std::vector<HighsInt> stack;
  std::vector<double> coefs;
  std::vector<HighsInt> cover;

  num_strengthened = 0;
  // Check for timeout according to this frequency
  const HighsInt check_time_frequency = 100;

  for (HighsInt row = 0; row != model->num_row_; ++row) {
    if (rowsize[row] <= 1) continue;
    if (isRanged(row)) continue;

    // do not run on very dense rows as this could get expensive
    HighsInt rowsize_limit =
        std::max(HighsInt{1000}, (model->num_col_ - numDeletedCols) / 20);
    if (rowsize[row] > rowsize_limit) continue;

    // printf("strengthening knapsack of %" HIGHSINT_FORMAT " vars\n",
    // rowsize[row]);

    HighsCDouble maxviolation;
    HighsCDouble continuouscontribution = 0.0;
    double scale;

    if (model->row_lower_[row] != -kHighsInf) {
      // ax >= lb, treat as -ax <= -lb --> -ax + lb <= 0
      maxviolation = model->row_lower_[row];
      scale = -1.0;
    } else {
      // ax <= ub --> ax - ub <= 0
      maxviolation = -model->row_upper_[row];
      scale = 1.0;
    }

    complementation.clear();
    reducedcost.clear();
    upper.clear();
    indices.clear();
    positions.clear();
    complementation.reserve(rowsize[row]);
    reducedcost.reserve(rowsize[row]);
    upper.reserve(rowsize[row]);
    indices.reserve(rowsize[row]);
    positions.reserve(rowsize[row]);
    stack.reserve(rowsize[row]);
    stack.push_back(rowroot[row]);

    bool skiprow = false;

    while (!stack.empty()) {
      // pop element from stack
      HighsInt pos = stack.back();
      stack.pop_back();

      // add non-zeros to stack
      if (ARright[pos] != -1) stack.push_back(ARright[pos]);
      if (ARleft[pos] != -1) stack.push_back(ARleft[pos]);

      // get column index
      HighsInt col = Acol[pos];

      // skip row if a column bound is not finite
      skiprow = model->col_lower_[col] == -kHighsInf ||
                model->col_upper_[col] == kHighsInf;
      if (skiprow) break;

      // compute maximum violation
      // scale =  1:  ax <=  ub --> violation =  ax - ub > 0
      // scale = -1: -ax <= -lb --> violation = -ax + lb > 0
      // this means that for scale = 1 we sum up an upper bound on constraint
      // activity, and for scale = -1 we sum up a lower bound on constraint
      // activity.
      int8_t comp;
      double weight = Avalue[pos] * scale;
      HighsCDouble ub = static_cast<HighsCDouble>(model->col_upper_[col]) -
                        static_cast<HighsCDouble>(model->col_lower_[col]);
      if (weight > 0) {
        comp = 1;
        maxviolation +=
            static_cast<HighsCDouble>(model->col_upper_[col]) * weight;
      } else {
        comp = -1;
        maxviolation +=
            static_cast<HighsCDouble>(model->col_lower_[col]) * weight;
        weight = -weight;
      }

      if (ub <= primal_feastol || weight <= primal_feastol) continue;

      if (model->integrality_[col] == HighsVarType::kContinuous) {
        continuouscontribution += weight * ub;
        continue;
      }

      indices.push_back(reducedcost.size());
      positions.push_back(pos);
      reducedcost.push_back(weight);
      complementation.push_back(comp);
      upper.push_back(ub);
    }

    // Check for timeout according to frequency, unless a particularly
    // dense row has just been analysed
    if ((row & check_time_frequency) == 0 || 10 * rowsize[row] > rowsize_limit)
      HPRESOLVE_CHECKED_CALL(checkTimeLimit());

    if (skiprow) {
      stack.clear();
      continue;
    }

    // maxviolation <= 0 implies that the constraint is redundant:
    // scale =  1: upper bound on activity <= model->row_upper_[row]
    // scale = -1: model->row_lower_[row]  <= lower bound on activity
    if (maxviolation <= primal_feastol) {
      HPRESOLVE_CHECKED_CALL(rowPresolve(postsolve_stack, row));
      continue;
    }

    const double smallVal =
        std::max(100 * primal_feastol,
                 primal_feastol * static_cast<double>(maxviolation));
    while (true) {
      if (maxviolation - continuouscontribution <= smallVal || indices.empty())
        break;

      pdqsort(indices.begin(), indices.end(), [&](HighsInt i1, HighsInt i2) {
        return std::make_pair(reducedcost[i1], i1) >
               std::make_pair(reducedcost[i2], i2);
      });

      HighsCDouble lambda = maxviolation - continuouscontribution;

      cover.clear();
      cover.reserve(indices.size());

      for (size_t i = indices.size(); i > 0; --i) {
        HighsInt index = indices[i - 1];
        HighsCDouble delta = upper[index] * reducedcost[index];

        if (upper[index] <= 1000.0 && reducedcost[index] > smallVal &&
            lambda - delta <= smallVal)
          cover.push_back(index);
        else
          lambda -= delta;
      }

      if (cover.empty() || lambda <= smallVal) break;

      HighsInt alpos = *std::min_element(
          cover.begin(), cover.end(), [&](HighsInt i1, HighsInt i2) {
            if (reducedcost[i1] <= 1e-3 || reducedcost[i2] <= 1e-3)
              return reducedcost[i1] > reducedcost[i2];
            return reducedcost[i1] < reducedcost[i2];
          });

      HighsCDouble al = reducedcost[alpos];
      coefs.resize(cover.size());
      HighsCDouble coverrhs = max(ceil(lambda / al - primal_feastol), 1.0);
      HighsCDouble slackupper = -coverrhs;

      HighsCDouble step = kHighsInf;
      for (size_t i = 0; i != cover.size(); ++i) {
        coefs[i] =
            static_cast<double>(ceil(min(reducedcost[cover[i]], lambda) / al -
                                     options->small_matrix_value));
        slackupper += upper[cover[i]] * coefs[i];
        step = min(step, reducedcost[cover[i]] / coefs[i]);
      }
      step = min(step, maxviolation / coverrhs);
      maxviolation -= step * coverrhs;

      HighsInt slackind = reducedcost.size();
      reducedcost.push_back(step);
      upper.push_back(slackupper);

      for (size_t i = 0; i != cover.size(); ++i)
        reducedcost[cover[i]] -= step * coefs[i];

      indices.erase(std::remove_if(indices.begin(), indices.end(),
                                   [&](HighsInt i) {
                                     return reducedcost[i] <= primal_feastol;
                                   }),
                    indices.end());
      indices.push_back(slackind);
    }

    double threshold = static_cast<double>(maxviolation + primal_feastol);

    indices.erase(std::remove_if(indices.begin(), indices.end(),
                                 [&](HighsInt i) {
                                   return static_cast<size_t>(i) >=
                                              positions.size() ||
                                          abs(reducedcost[i]) <= threshold;
                                 }),
                  indices.end());
    if (indices.empty()) continue;

    auto updateNonZeros = [&](HighsInt row, HighsCDouble& rhs,
                              HighsInt direction) {
      for (HighsInt i : indices) {
        assert(Arow[positions[i]] == row);
        HighsCDouble coefdelta = direction * (reducedcost[i] - maxviolation);
        HighsInt col = Acol[positions[i]];

        if (complementation[i] == -1) {
          rhs += coefdelta * model->col_lower_[col];
          addToMatrix(row, col, static_cast<double>(coefdelta));
        } else {
          rhs -= coefdelta * model->col_upper_[col];
          addToMatrix(row, col, static_cast<double>(-coefdelta));
        }
      }
    };

    // update / add non-zeros
    if (scale < 0) {
      HighsCDouble lhs = model->row_lower_[row];
      updateNonZeros(row, lhs, HighsInt{-1});
      changeRowLower(row, static_cast<double>(lhs));
    } else {
      HighsCDouble rhs = model->row_upper_[row];
      updateNonZeros(row, rhs, HighsInt{1});
      changeRowUpper(row, static_cast<double>(rhs));
    }

    num_strengthened += indices.size();
  }

  return Result::kOk;
}

HPresolve::Result HPresolve::detectImpliedIntegers() {
  for (HighsInt col = 0; col != model->num_col_; ++col)
    HPRESOLVE_CHECKED_CALL(static_cast<Result>(convertImpliedInteger(col)));
  return Result::kOk;
}

HPresolve::Result HPresolve::updateCliqueTableFixedCol(const HighsInt col,
                                                       const double val) {
  if (mipsolver == nullptr || !mipsolver->mipdata_->cliquesExtracted ||
      model->integrality_[col] == HighsVarType::kContinuous ||
      (val != 0.0 && val != 1.0)) {
    return Result::kOk;
  }
  std::vector<HighsCliqueTable::CliqueVar> impliedFixings;
  if (!presolveCliqueTable.fixCol(col, static_cast<bool>(val),
                                  impliedFixings)) {
    return Result::kPrimalInfeasible;
  }

  for (const HighsCliqueTable::CliqueVar& fixing : impliedFixings) {
    if (colDeleted[fixing.col]) continue;
    if (fixing.val < model->col_lower_[fixing.col] - primal_feastol ||
        fixing.val > model->col_upper_[fixing.col] + primal_feastol) {
      return Result::kPrimalInfeasible;
    }
    HPRESOLVE_CHECKED_CALL(changeColBounds(fixing.col, fixing.val, fixing.val));
  }
  return Result::kOk;
}

HPresolve::Result HPresolve::updateCliqueTableSubstituteCol(
    const HighsInt substCol, const HighsInt stayCol, const double offset,
    const double scale) {
  if (mipsolver == nullptr || !mipsolver->mipdata_->cliquesExtracted) {
    return Result::kOk;
  }
  bool isBinary =
      model->integrality_[substCol] != HighsVarType::kContinuous &&
      model->integrality_[stayCol] != HighsVarType::kContinuous &&
      model->col_lower_[substCol] == 0.0 && model->col_lower_[stayCol] == 0.0 &&
      model->col_upper_[substCol] == 1.0 && model->col_upper_[stayCol] == 1.0 &&
      ((offset == 0.0 && scale == 1.0) || (offset == 1.0 && scale == -1.0));

  if (!isBinary) return Result::kOk;

  HighsCliqueTable::CliqueVar replacement(
      stayCol, scale == 1.0 ? HighsInt{1} : HighsInt{0});
  std::vector<HighsCliqueTable::CliqueVar> impliedFixings;
  if (!presolveCliqueTable.substituteCol(substCol, replacement, impliedFixings))
    return Result::kPrimalInfeasible;

  for (HighsCliqueTable::CliqueVar& v : impliedFixings) {
    HighsInt col = static_cast<HighsInt>(v.col);
    double val = static_cast<double>(v.val);
    if (colDeleted[col]) continue;
    if (val < model->col_lower_[col] - primal_feastol ||
        val > model->col_upper_[col] + primal_feastol) {
      return Result::kPrimalInfeasible;
    }
    HPRESOLVE_CHECKED_CALL(changeColBounds(col, val, val));
  }
  return Result::kOk;
}

bool HPresolve::zeroRowActivityFeasible() const {
  // Check that zero row activity is feasible - called when reduced model
  // has no columns to assess whether the HighsModelStatus returned is
  // kOptimal or kInfeasible (as was required for 2326)
  for (HighsInt iRow = 0; iRow < model->num_row_; iRow++)
    if (model->row_lower_[iRow] > primal_feastol ||
        model->row_upper_[iRow] < -primal_feastol)
      return false;
  return true;
}

HighsInt HPresolve::debugGetCheckCol() const {
  const std::string check_col_name = "";  // c37";
  HighsInt check_col = -1;
  if (check_col_name == "") return check_col;
  if (model->col_names_.size()) {
    if (model->col_hash_.size() != static_cast<size_t>(model->num_col_))
      model->col_hash_.form(model->col_names_);
    auto search = model->col_hash_.name2index.find(check_col_name);
    if (search != model->col_hash_.name2index.end()) {
      check_col = search->second;
      assert(model->col_names_[check_col] == check_col_name);
    }
  }
  return check_col;
}

HighsInt HPresolve::debugGetCheckRow() const {
  const std::string check_row_name = "";  //"row_ekk_119";
  HighsInt check_row = -1;
  if (check_row_name == "") return check_row;
  if (model->row_names_.size()) {
    if (model->row_hash_.size() != static_cast<size_t>(model->num_row_))
      model->row_hash_.form(model->row_names_);
    auto search = model->row_hash_.name2index.find(check_row_name);
    if (search != model->row_hash_.name2index.end()) {
      check_row = search->second;
      assert(model->row_names_[check_row] == check_row_name);
    }
  }
  return check_row;
}

}  // namespace presolve
