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
#include "lp_data/HighsSolution.h"
#include "mip/HighsCliqueTable.h"
#include "mip/HighsImplications.h"
#include "mip/HighsMipSolverData.h"
#include "mip/HighsObjectiveFunction.h"
#include "mip/MipTimer.h"
#include "presolve/HighsPostsolveStack.h"
#include "test_kkt/DevKkt.h"
#include "util/HFactor.h"
#include "util/HighsCDouble.h"
#include "util/HighsIntegers.h"
#include "util/HighsLinearSumBounds.h"
#include "util/HighsMemoryAllocation.h"
#include "util/HighsSplay.h"
#include "util/HighsUtils.h"

#define ENABLE_SPARSIFY_FOR_LP 0

#define HPRESOLVE_CHECKED_CALL(presolveCall)                           \
  do {                                                                 \
    HPresolve::Result __result = presolveCall;                         \
    if (__result != presolve::HPresolve::Result::kOk) return __result; \
  } while (0)

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

bool HPresolve::okSetInput(HighsLp& model_, const HighsOptions& options_,
                           const HighsInt presolve_reduction_limit,
                           HighsTimer* timer) {
  model = &model_;
  options = &options_;
  this->timer = timer;

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

  if (mipsolver == nullptr) {
    primal_feastol = options->primal_feasibility_tolerance;
    model->integrality_.assign(model->num_col_, HighsVarType::kContinuous);
  } else
    primal_feastol = options->mip_feasibility_tolerance;

  if (model_.a_matrix_.isRowwise()) {
    // Does this even happen?
    assert(model_.a_matrix_.isColwise());
    if (!okFromCSR(model->a_matrix_.value_, model->a_matrix_.index_,
                   model->a_matrix_.start_))
      return false;
  } else {
    if (!okFromCSC(model->a_matrix_.value_, model->a_matrix_.index_,
                   model->a_matrix_.start_))
      return false;
  }

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
  numDeletedCols = 0;
  numDeletedRows = 0;
  // initialize substitution opportunities
  for (HighsInt row = 0; row != model->num_row_; ++row) {
    if (!isDualImpliedFree(row)) continue;
    for (const HighsSliceNonzero& nonzero : getRowVector(row)) {
      if (isImpliedFree(nonzero.index()))
        substitutionOpportunities.emplace_back(row, nonzero.index());
    }
  }
  // Take value passed in as reduction limit, allowing different
  // values to be used for initial presolve, and after restart
  reductionLimit =
      presolve_reduction_limit < 0 ? kHighsSize_tInf : presolve_reduction_limit;
  if (options->presolve != kHighsOffString &&
      reductionLimit < kHighsSize_tInf) {
    highsLogDev(options->log_options, HighsLogType::kInfo,
                "HPresolve::okSetInput reductionLimit = %d\n",
                static_cast<int>(reductionLimit));
  }
  return true;
}

// for MIP presolve
bool HPresolve::okSetInput(HighsMipSolver& mipsolver,
                           const HighsInt presolve_reduction_limit) {
  this->mipsolver = &mipsolver;

  probingContingent = 1000;
  probingNumDelCol = 0;
  numProbed = 0;
  numProbes.assign(mipsolver.numCol(), 0);

  if (mipsolver.model_ != &mipsolver.mipdata_->presolvedModel) {
    mipsolver.mipdata_->presolvedModel = *mipsolver.model_;
    mipsolver.model_ = &mipsolver.mipdata_->presolvedModel;
  } else {
    mipsolver.mipdata_->presolvedModel.col_lower_ =
        mipsolver.mipdata_->domain.col_lower_;
    mipsolver.mipdata_->presolvedModel.col_upper_ =
        mipsolver.mipdata_->domain.col_upper_;
  }

  return okSetInput(mipsolver.mipdata_->presolvedModel, *mipsolver.options_mip_,
                    presolve_reduction_limit, &mipsolver.timer_);
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

bool HPresolve::isUpperImplied(HighsInt col) const {
  return (model->col_upper_[col] == kHighsInf ||
          implColUpper[col] <= model->col_upper_[col] + primal_feastol);
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
          options->small_matrix_value) {
        // update right-hand side and mark row as changed
        model->row_upper_[nz.index()] = rUpper;
        markChangedRow(nz.index());
      }
    }
    if (model->row_lower_[nz.index()] != -kHighsInf) {
      // left-hand side: scale, round up and unscale again
      double rLower =
          std::abs(nz.value()) *
          std::ceil(model->row_lower_[nz.index()] * std::abs(scale) -
                    primal_feastol);
      // check if modification is large enough
      if (std::abs(model->row_lower_[nz.index()] - rLower) >
          options->small_matrix_value) {
        // update left-hand side and mark row as changed
        model->row_lower_[nz.index()] = rLower;
        markChangedRow(nz.index());
      }
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

  // round and update bounds
  changeColLower(col, model->col_lower_[col]);
  changeColUpper(col, model->col_upper_[col]);
  return StatusResult(true);
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

  impliedRowBounds.add(Arow[pos], Acol[pos], Avalue[pos]);
  impliedDualRowBounds.add(Acol[pos], Arow[pos], Avalue[pos]);
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

    impliedDualRowBounds.remove(Acol[pos], Arow[pos], Avalue[pos]);
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
    impliedRowBounds.remove(Arow[pos], Acol[pos], Avalue[pos]);
  }

  // remove implied bounds on row duals that where implied by this column's dual
  // constraint
  resetRowDualImpliedBoundsDerivedFromCol(Acol[pos]);

  // remove implied bounds on columns that where implied by this row
  resetColImpliedBoundsDerivedFromRow(Arow[pos]);

  // modifications to row invalidate lifting opportunities
  clearLiftingOpportunities(Arow[pos]);

  // remove non-zero
  Avalue[pos] = 0;

  freeslots.push_back(pos);
}

void HPresolve::markChangedRow(HighsInt row) {
  if (!changedRowFlag[row]) {
    changedRowIndices.push_back(row);
    changedRowFlag[row] = true;
  }
}

void HPresolve::markChangedCol(HighsInt col) {
  if (!changedColFlag[col]) {
    changedColIndices.push_back(col);
    changedColFlag[col] = true;
  }
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
  // if the column has an infinite lower bound the reduced cost cannot be
  // positive, i.e. the column corresponds to a <= constraint in the dual with
  // right hand side -cost which becomes a >= constraint with side +cost.
  // Furthermore, we can ignore strictly redundant primal
  // column bounds and treat them as if they are infinite
  double impliedMargin = colsize[col] != 1 ? primal_feastol : -primal_feastol;

  double myDualRowLower =
      (model->col_lower_[col] == -kHighsInf) ||
              (implColLower[col] > model->col_lower_[col] + impliedMargin)
          ? model->col_cost_[col]
          : -kHighsInf;

  double myDualRowUpper =
      (model->col_upper_[col] == kHighsInf) ||
              (implColUpper[col] < model->col_upper_[col] - impliedMargin)
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
  // if the column has an infinite lower bound the reduced cost cannot be
  // positive, i.e. the column corresponds to a <= constraint in the dual with
  // right hand side -cost which becomes a >= constraint with side +cost.
  // Furthermore, we can ignore strictly redundant primal
  // column bounds and treat them as if they are infinite
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

bool HPresolve::checkUpdateColImpliedBounds(HighsInt row, double* rowLower,
                                            double* rowUpper) const {
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

void HPresolve::updateColImpliedBounds(HighsInt row, HighsInt col, double val) {
  // propagate implied column bound upper bound if row has an upper bound
  double rowLower, rowUpper;
  if (!checkUpdateColImpliedBounds(row, &rowLower, &rowUpper)) return;

  const double threshold = 1000 * primal_feastol;

  auto checkImpliedBound = [&](HighsInt row, HighsInt col, double val,
                               double rowBnd, double residualAct,
                               HighsInt direction) {
    if (direction * residualAct <= -kHighsInf) return;
    double impliedBound = static_cast<double>(
        (static_cast<HighsCDouble>(rowBnd) - residualAct) / val);

    if (std::abs(impliedBound) * kHighsTiny > primal_feastol) return;

    if (direction * val > 0) {
      // bound is an upper bound
      if (mipsolver != nullptr) {
        // solving a MIP; keep tighter bounds on integer variables
        if (model->integrality_[col] != HighsVarType::kContinuous &&
            impliedBound < model->col_upper_[col] - primal_feastol)
          changeColUpper(col, impliedBound);

        // do not use the implied bound if this a not a model row, since the
        // row can be removed and should not be used, e.g., to identify a
        // column as implied free
        if (mipsolver->mipdata_->postSolveStack.getOrigRowIndex(row) >=
            mipsolver->orig_model_->num_row_) {
          // keep implied bound (as column bound)
          if (impliedBound < model->col_upper_[col] - threshold)
            changeColUpper(col, impliedBound);
          // set to +infinity, so that it is not stored as an implied bound
          impliedBound = kHighsInf;
        }
      }

      // only tighten bound if it is tighter by a wide enough margin
      if (impliedBound < implColUpper[col] - threshold)
        changeImplColUpper(col, impliedBound, row);
    } else {
      // bound is a lower bound
      if (mipsolver != nullptr) {
        // solving a MIP; keep tighter bounds on integer variables
        if (model->integrality_[col] != HighsVarType::kContinuous &&
            impliedBound > model->col_lower_[col] + primal_feastol)
          changeColLower(col, impliedBound);

        // do not use the implied bound if this a not a model row, since the
        // row can be removed and should not be used, e.g., to identify a
        // column as implied free
        if (mipsolver->mipdata_->postSolveStack.getOrigRowIndex(row) >=
            mipsolver->orig_model_->num_row_) {
          // keep implied bound (as column bound)
          if (impliedBound > model->col_lower_[col] + threshold)
            changeColLower(col, impliedBound);
          // set to -infinity, so that it is not stored as an implied bound
          impliedBound = -kHighsInf;
        }
      }

      // only tighten bound if it is tighter by a wide enough margin
      if (impliedBound > implColLower[col] + threshold)
        changeImplColLower(col, impliedBound, row);
    }
  };

  if (rowUpper != kHighsInf)
    checkImpliedBound(row, col, val, rowUpper,
                      impliedRowBounds.getResidualSumLowerOrig(row, col, val),
                      HighsInt{1});

  if (rowLower != -kHighsInf)
    checkImpliedBound(row, col, val, rowLower,
                      impliedRowBounds.getResidualSumUpperOrig(row, col, val),
                      HighsInt{-1});
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
  std::set<HighsInt> affectedCols(colImplSourceByRow[row]);
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
  std::set<HighsInt> affectedRows(implRowDualSourceByCol[col]);
  for (const HighsInt& row : affectedRows) {
    // set implied bounds to infinite values if they were deduced from the
    // given column
    resetRowDualImpliedBounds(row, col);
  }
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
  HighsInt oldNumCol = model->num_col_;
  model->num_col_ = 0;
  std::vector<HighsInt> newColIndex(oldNumCol);
  const bool have_col_names = model->col_names_.size() > 0;
  assert(!have_col_names ||
         model->col_names_.size() == static_cast<size_t>(oldNumCol));
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
  HighsInt oldNumRow = model->num_row_;
  const bool have_row_names = model->row_names_.size() > 0;
  assert(!have_row_names ||
         model->row_names_.size() == static_cast<size_t>(oldNumRow));
  model->num_row_ = 0;
  std::vector<HighsInt> newRowIndex(oldNumRow);
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

  numDeletedRows = 0;
  postsolve_stack.compressIndexMaps(newRowIndex, newColIndex);
  impliedRowBounds.shrink(newRowIndex, model->num_row_);
  impliedDualRowBounds.shrink(newColIndex, model->num_col_);

  HighsInt numNnz = Avalue.size();
  for (HighsInt i = 0; i != numNnz; ++i) {
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
    mipsolver->mipdata_->domain = HighsDomain(*mipsolver);
    mipsolver->mipdata_->cliquetable.rebuild(model->num_col_, postsolve_stack,
                                             mipsolver->mipdata_->domain,
                                             newColIndex, newRowIndex);
    mipsolver->mipdata_->implications.rebuild(model->num_col_, newColIndex,
                                              newRowIndex);
    mipsolver->mipdata_->cutpool =
        HighsCutPool(mipsolver->model_->num_col_,
                     mipsolver->options_mip_->mip_pool_age_limit,
                     mipsolver->options_mip_->mip_pool_soft_limit);
    mipsolver->mipdata_->conflictPool =
        HighsConflictPool(5 * mipsolver->options_mip_->mip_pool_age_limit,
                          mipsolver->options_mip_->mip_pool_soft_limit);

    for (HighsInt i = 0; i != oldNumCol; ++i)
      if (newColIndex[i] != -1) numProbes[newColIndex[i]] = numProbes[i];

    mipsolver->mipdata_->debugSolution.shrink(newColIndex);
    numProbes.resize(model->num_col_);
    // Need to set the constraint matrix dimensions
    model->setMatrixDimensions();
  }
  // Need to set the constraint matrix dimensions
  model->setMatrixDimensions();
  // Need to reset current number of deleted rows and columns in logging
  analysis_.resetNumDeleted();
}

HPresolve::Result HPresolve::dominatedColumns(
    HighsPostsolveStack& postsolve_stack) {
  std::vector<std::pair<uint32_t, uint32_t>> signatures(model->num_col_);

  auto isBinary = [&](HighsInt i) {
    return model->integrality_[i] == HighsVarType::kInteger &&
           model->col_lower_[i] == 0.0 && model->col_upper_[i] == 1.0;
  };

  auto addSignature = [&](HighsInt row, HighsInt col, uint32_t rowLowerFinite,
                          uint32_t rowUpperFinite) {
    HighsInt rowHashedPos = (HighsHashHelpers::hash(row) >> 59);
    assert(rowHashedPos < 32);
    signatures[col].first |= rowLowerFinite << rowHashedPos;
    signatures[col].second |= rowUpperFinite << rowHashedPos;
  };

  auto checkDomination = [&](HighsInt scalj, HighsInt j, HighsInt scalk,
                             HighsInt k) {
    // rule out domination from integers to continuous variables
    if (model->integrality_[j] == HighsVarType::kInteger &&
        model->integrality_[k] != HighsVarType::kInteger)
      return false;

    // check the signatures
    uint32_t sjMinus = signatures[j].first;
    uint32_t sjPlus = signatures[j].second;
    if (scalj == -1) std::swap(sjPlus, sjMinus);

    uint32_t skMinus = signatures[k].first;
    uint32_t skPlus = signatures[k].second;
    if (scalk == -1) std::swap(skPlus, skMinus);

    // the set of rows with a negative coefficient must be a superset of the
    // dominated columns set of rows with a negative coefficient
    if ((~sjMinus & skMinus) != 0) return false;

    // the set of rows with a positive coefficient must be a subset of the
    // dominated columns set of rows with a positive coefficient
    if ((sjPlus & ~skPlus) != 0) return false;

    // next check if the columns cost allows for domination
    double cj = scalj * model->col_cost_[j];
    double ck = scalk * model->col_cost_[k];

    // the dominating columns cost must be smaller or equal to the dominated
    // columns cost
    if (cj > ck + options->small_matrix_value) return false;

    // finally check the column vectors
    for (const HighsSliceNonzero& nonz : getColumnVector(j)) {
      HighsInt row = nonz.index();
      double aj = scalj * nonz.value();

      HighsInt akPos = findNonzero(row, k);
      double ak = scalk * (akPos == -1 ? 0.0 : Avalue[akPos]);

      if (model->row_lower_[row] != -kHighsInf &&
          model->row_upper_[row] != kHighsInf) {
        // the row is an equality or ranged row, therefore the coefficients must
        // be parallel, otherwise one of the inequalities given by the row rules
        // out domination
        if (std::abs(aj - ak) > options->small_matrix_value) return false;
        continue;
      }

      // normalize row to a <= constraint
      if (model->row_upper_[row] == kHighsInf) {
        aj = -aj;
        ak = -ak;
      }

      // the coefficient of the dominating column needs to be smaller than or
      // equal to the coefficient of the dominated column
      if (aj > ak + options->small_matrix_value) return false;
    }

    // check row only occurring in the column vector of k
    for (const HighsSliceNonzero& nonz : getColumnVector(k)) {
      HighsInt row = nonz.index();
      double ak = scalk * nonz.value();

      HighsInt ajPos = findNonzero(row, j);
      // only rows in which aj does not occur are left to check
      if (ajPos != -1) continue;
      double aj = 0.0;

      if (model->row_lower_[row] != -kHighsInf &&
          model->row_upper_[row] != kHighsInf) {
        // the row is an equality or ranged row, therefore the coefficients must
        // be parallel, otherwise one of the inequalities given by the row rules
        // out domination
        if (std::abs(aj - ak) > options->small_matrix_value) return false;
        continue;
      }

      // normalize row to a <= constraint
      if (model->row_upper_[row] == kHighsInf) {
        aj = -aj;
        ak = -ak;
      }

      // the coefficient of the dominating column needs to be smaller than or
      // equal to the coefficient of the dominated column
      if (aj > ak + options->small_matrix_value) return false;
    }

    return true;
  };

  HighsInt numNz = Avalue.size();
  for (HighsInt i = 0; i < numNz; ++i) {
    if (Avalue[i] == 0) continue;

    HighsInt row = Arow[i];
    HighsInt col = Acol[i];
    bool rowLowerFinite = model->row_lower_[row] != -kHighsInf;
    bool rowUpperFinite = model->row_upper_[row] != kHighsInf;

    if (Avalue[i] > 0)
      addSignature(row, col, rowLowerFinite, rowUpperFinite);
    else
      addSignature(row, col, rowUpperFinite, rowLowerFinite);
  }

  HighsInt numFixedCols = 0;
  for (HighsInt j = 0; j < model->num_col_; ++j) {
    if (colDeleted[j]) continue;
    bool upperImplied = isUpperImplied(j);
    bool lowerImplied = isLowerImplied(j);
    bool hasPosCliques = false;
    bool hasNegCliques = false;
    bool colIsBinary = isBinary(j);
    if (colIsBinary) {
      hasPosCliques = mipsolver->mipdata_->cliquetable.numCliques(j, 1) > 0;
      hasNegCliques = mipsolver->mipdata_->cliquetable.numCliques(j, 0) > 0;
    } else if (!upperImplied && !lowerImplied)
      continue;

    HighsInt oldNumFixed = numFixedCols;

    HighsInt bestRowPlus = -1;
    HighsInt bestRowPlusLen = kHighsIInf;
    HighsInt bestRowPlusScale = 0;
    double ajBestRowPlus = 0.0;
    HighsInt bestRowMinus = -1;
    HighsInt bestRowMinusLen = kHighsIInf;
    HighsInt bestRowMinusScale = 0;
    double ajBestRowMinus = 0.0;

    double worstCaseLb = -kHighsInf;
    double worstCaseUb = kHighsInf;

    bool checkPosRow = upperImplied || colIsBinary;
    bool checkNegRow = lowerImplied || colIsBinary;
    for (const HighsSliceNonzero& nonz : getColumnVector(j)) {
      HighsInt row = nonz.index();
      HighsInt scale = model->row_upper_[row] != kHighsInf ? 1 : -1;

      if (colIsBinary) {
        if (model->row_upper_[row] != kHighsInf) {
          if (model->col_cost_[j] >= 0.0 && nonz.value() < 0.0) {
            double maxresact =
                impliedRowBounds.getResidualSumUpper(row, j, nonz.value());
            double wcBound =
                (model->row_upper_[row] - maxresact) / nonz.value();
            worstCaseLb = std::max(wcBound, worstCaseLb);
          } else if (model->col_cost_[j] <= 0.0 && nonz.value() > 0.0) {
            double maxresact =
                impliedRowBounds.getResidualSumUpper(row, j, nonz.value());
            double wcBound =
                (model->row_upper_[row] - maxresact) / nonz.value();
            worstCaseUb = std::min(wcBound, worstCaseUb);
          }
        }

        if (model->row_lower_[row] != -kHighsInf) {
          if (model->col_cost_[j] >= 0.0 && nonz.value() > 0.0) {
            double minresact =
                impliedRowBounds.getResidualSumLower(row, j, nonz.value());
            double wcBound =
                (model->row_lower_[row] - minresact) / nonz.value();
            worstCaseLb = std::max(wcBound, worstCaseLb);
          } else if (model->col_cost_[j] <= 0.0 && nonz.value() < 0.0) {
            double minresact =
                impliedRowBounds.getResidualSumLower(row, j, nonz.value());
            double wcBound =
                (model->row_lower_[row] - minresact) / nonz.value();
            worstCaseUb = std::min(wcBound, worstCaseUb);
          }
        }
      }

      double val = scale * nonz.value();
      if (checkPosRow && val > 0.0 && rowsize[row] < bestRowPlusLen) {
        bestRowPlus = row;
        bestRowPlusLen = rowsize[row];
        bestRowPlusScale = scale;
        ajBestRowPlus = val;
      }

      if (checkNegRow && val < 0.0 && rowsize[row] < bestRowMinusLen) {
        bestRowMinus = row;
        bestRowMinusLen = rowsize[row];
        bestRowMinusScale = scale;
        ajBestRowMinus = val;
      }
    }

    if (colIsBinary) {
      if (model->col_cost_[j] >= 0.0 && worstCaseLb <= 1 + primal_feastol) {
        upperImplied = true;
        if (!lowerImplied && bestRowMinus != -1) {
          storeRow(bestRowMinus);

          bool isEqOrRangedRow =
              model->row_lower_[bestRowMinus] != -kHighsInf &&
              model->row_upper_[bestRowMinus] != kHighsInf;

          for (const HighsSliceNonzero& nonz : getStoredRow()) {
            HighsInt k = nonz.index();
            if (k == j || colDeleted[k]) continue;

            double ak = nonz.value() * bestRowMinusScale;

            if (-ajBestRowMinus <= -ak + options->small_matrix_value &&
                (!isEqOrRangedRow ||
                 -ajBestRowMinus >= -ak - options->small_matrix_value) &&
                checkDomination(-1, j, -1, k)) {
              // case (iii)  lb(x_j) = -inf, -x_j > -x_k: set x_k = ub(x_k)
              ++numFixedCols;
              if (fixColToLowerOrUnbounded(postsolve_stack, j)) {
                // Handle unboundedness
                presolve_status_ = HighsPresolveStatus::kUnboundedOrInfeasible;
                return Result::kDualInfeasible;
              }
              HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolve_stack));
              break;
            } else if (-ajBestRowMinus <= ak + options->small_matrix_value &&
                       (!isEqOrRangedRow ||
                        -ajBestRowMinus >= ak - options->small_matrix_value) &&
                       checkDomination(-1, j, 1, k)) {
              // case (iv)  lb(x_j) = -inf, -x_j > x_k: set x_k = lb(x_k)
              ++numFixedCols;
              if (fixColToLowerOrUnbounded(postsolve_stack, j)) {
                // Handle unboundedness
                presolve_status_ = HighsPresolveStatus::kUnboundedOrInfeasible;
                return Result::kDualInfeasible;
              }
              HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolve_stack));
              break;
            }
          }

          if (colDeleted[j]) {
            HPRESOLVE_CHECKED_CALL(removeDoubletonEquations(postsolve_stack));
            continue;
          }
        }
      }

      if (model->col_cost_[j] <= 0.0 && worstCaseUb >= -primal_feastol) {
        lowerImplied = true;
        if (!upperImplied && bestRowPlus != -1) {
          storeRow(bestRowPlus);
          bool isEqOrRangedRow = model->row_lower_[bestRowPlus] != -kHighsInf &&
                                 model->row_upper_[bestRowPlus] != kHighsInf;
          for (const HighsSliceNonzero& nonz : getStoredRow()) {
            HighsInt k = nonz.index();
            if (k == j || colDeleted[k]) continue;

            double ak = nonz.value() * bestRowPlusScale;

            if (ajBestRowPlus <= ak + options->small_matrix_value &&
                (!isEqOrRangedRow ||
                 ajBestRowPlus >= ak - options->small_matrix_value) &&
                checkDomination(1, j, 1, k)) {
              // case (i)  ub(x_j) = inf, x_j > x_k: set x_k = lb(x_k)
              ++numFixedCols;
              if (fixColToUpperOrUnbounded(postsolve_stack, j)) {
                // Handle unboundedness
                presolve_status_ = HighsPresolveStatus::kUnboundedOrInfeasible;
                return Result::kDualInfeasible;
              }
              HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolve_stack));
              break;
            } else if (ajBestRowPlus <= -ak + options->small_matrix_value &&
                       (!isEqOrRangedRow ||
                        ajBestRowPlus >= -ak - options->small_matrix_value) &&
                       checkDomination(1, j, -1, k)) {
              // case (ii)  ub(x_j) = inf, x_j > -x_k: set x_k = ub(x_k)
              ++numFixedCols;
              if (fixColToUpperOrUnbounded(postsolve_stack, j)) {
                // Handle unboundedness
                presolve_status_ = HighsPresolveStatus::kUnboundedOrInfeasible;
                return Result::kDualInfeasible;
              }
              HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolve_stack));
              break;
            }
          }

          if (colDeleted[j]) {
            HPRESOLVE_CHECKED_CALL(removeDoubletonEquations(postsolve_stack));
            continue;
          }
        }
      }

      if (!upperImplied && !hasPosCliques) bestRowPlus = -1;

      if (!lowerImplied && !hasNegCliques) bestRowMinus = -1;
    }

    if (bestRowPlus != -1) {
      assert(upperImplied || hasPosCliques);
      storeRow(bestRowPlus);
      bool isEqOrRangedRow = model->row_lower_[bestRowPlus] != -kHighsInf &&
                             model->row_upper_[bestRowPlus] != kHighsInf;
      for (const HighsSliceNonzero& nonz : getStoredRow()) {
        HighsInt k = nonz.index();
        if (k == j || colDeleted[k]) continue;

        double ak = nonz.value() * bestRowPlusScale;

        if (model->col_lower_[k] != -kHighsInf &&
            (upperImplied || mipsolver->mipdata_->cliquetable.haveCommonClique(
                                 HighsCliqueTable::CliqueVar(j, 1),
                                 HighsCliqueTable::CliqueVar(k, 1))) &&
            ajBestRowPlus <= ak + options->small_matrix_value &&
            (!isEqOrRangedRow ||
             ajBestRowPlus >= ak - options->small_matrix_value) &&
            checkDomination(1, j, 1, k)) {
          // case (i)  ub(x_j) = inf, x_j > x_k: set x_k = lb(x_k)
          ++numFixedCols;
          if (fixColToLowerOrUnbounded(postsolve_stack, k)) {
            // Handle unboundedness
            presolve_status_ = HighsPresolveStatus::kUnboundedOrInfeasible;
            return Result::kDualInfeasible;
          }
          HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolve_stack));
        } else if (model->col_upper_[k] != kHighsInf &&
                   (upperImplied ||
                    mipsolver->mipdata_->cliquetable.haveCommonClique(
                        HighsCliqueTable::CliqueVar(j, 1),
                        HighsCliqueTable::CliqueVar(k, 0))) &&
                   ajBestRowPlus <= -ak + options->small_matrix_value &&
                   (!isEqOrRangedRow ||
                    ajBestRowPlus >= -ak - options->small_matrix_value) &&
                   checkDomination(1, j, -1, k)) {
          // case (ii)  ub(x_j) = inf, x_j > -x_k: set x_k = ub(x_k)
          ++numFixedCols;
          if (fixColToUpperOrUnbounded(postsolve_stack, k)) {
            // Handle unboundedness
            presolve_status_ = HighsPresolveStatus::kUnboundedOrInfeasible;
            return Result::kDualInfeasible;
          }
          HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolve_stack));
        }
      }
    }

    if (bestRowMinus != -1) {
      assert(lowerImplied || hasNegCliques);
      storeRow(bestRowMinus);

      bool isEqOrRangedRow = model->row_lower_[bestRowMinus] != -kHighsInf &&
                             model->row_upper_[bestRowMinus] != kHighsInf;

      for (const HighsSliceNonzero& nonz : getStoredRow()) {
        HighsInt k = nonz.index();
        if (k == j || colDeleted[k]) continue;

        double ak = nonz.value() * bestRowMinusScale;

        if (model->col_upper_[k] != kHighsInf &&
            (lowerImplied || mipsolver->mipdata_->cliquetable.haveCommonClique(
                                 HighsCliqueTable::CliqueVar(j, 0),
                                 HighsCliqueTable::CliqueVar(k, 0))) &&
            -ajBestRowMinus <= -ak + options->small_matrix_value &&
            (!isEqOrRangedRow ||
             -ajBestRowMinus >= -ak - options->small_matrix_value) &&
            checkDomination(-1, j, -1, k)) {
          // case (iii)  lb(x_j) = -inf, -x_j > -x_k: set x_k = ub(x_k)
          ++numFixedCols;
          if (fixColToUpperOrUnbounded(postsolve_stack, k)) {
            // Handle unboundedness
            presolve_status_ = HighsPresolveStatus::kUnboundedOrInfeasible;
            return Result::kDualInfeasible;
          }
          HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolve_stack));
        } else if (model->col_lower_[k] != -kHighsInf &&
                   (lowerImplied ||
                    mipsolver->mipdata_->cliquetable.haveCommonClique(
                        HighsCliqueTable::CliqueVar(j, 0),
                        HighsCliqueTable::CliqueVar(k, 1))) &&
                   -ajBestRowMinus <= ak + options->small_matrix_value &&
                   (!isEqOrRangedRow ||
                    -ajBestRowMinus >= ak - options->small_matrix_value) &&
                   checkDomination(-1, j, 1, k)) {
          // case (iv)  lb(x_j) = -inf, -x_j > x_k: set x_k = lb(x_k)
          ++numFixedCols;
          if (fixColToLowerOrUnbounded(postsolve_stack, k)) {
            // Handle unboundedness
            presolve_status_ = HighsPresolveStatus::kUnboundedOrInfeasible;
            return Result::kDualInfeasible;
          }
          HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolve_stack));
        }
      }
    }

    if (numFixedCols != oldNumFixed)
      HPRESOLVE_CHECKED_CALL(removeDoubletonEquations(postsolve_stack));
  }

  if (numFixedCols)
    highsLogDev(options->log_options, HighsLogType::kInfo,
                "Fixed %d dominated columns\n", numFixedCols);

  return Result::kOk;
}

HPresolve::Result HPresolve::runProbing(HighsPostsolveStack& postsolve_stack) {
  mipsolver->analysis_.mipTimerStart(kMipClockProbingPresolve);
  probingEarlyAbort = false;
  if (numDeletedCols + numDeletedRows != 0) shrinkProblem(postsolve_stack);

  toCSC(model->a_matrix_.value_, model->a_matrix_.index_,
        model->a_matrix_.start_);
  okFromCSC(model->a_matrix_.value_, model->a_matrix_.index_,
            model->a_matrix_.start_);

  mipsolver->mipdata_->cliquetable.setMaxEntries(numNonzeros());

  // first tighten all bounds if they have an implied bound that is tighter
  // than their column bound before probing this is not done for continuous
  // columns since it may allow stronger dual presolve and more aggregations
  double hugeBound = primal_feastol / kHighsTiny;
  for (HighsInt i = 0; i != model->num_col_; ++i) {
    if (model->col_lower_[i] >= implColLower[i] &&
        model->col_upper_[i] <= implColUpper[i])
      continue;

    if (std::abs(implColLower[i]) <= hugeBound) {
      double newLb = implColLower[i];
      if (newLb > model->col_lower_[i]) changeColLower(i, newLb);
    }

    if (std::abs(implColUpper[i]) <= hugeBound) {
      double newUb = implColUpper[i];
      if (newUb < model->col_upper_[i]) changeColUpper(i, newUb);
    }
  }

  HighsInt oldNumProbed = numProbed;

  mipsolver->mipdata_->setupDomainPropagation();
  HighsDomain& domain = mipsolver->mipdata_->domain;

  domain.propagate();
  if (domain.infeasible()) {
    mipsolver->analysis_.mipTimerStop(kMipClockProbingPresolve);
    return Result::kPrimalInfeasible;
  }
  HighsCliqueTable& cliquetable = mipsolver->mipdata_->cliquetable;
  HighsImplications& implications = mipsolver->mipdata_->implications;
  bool firstCall = !mipsolver->mipdata_->cliquesExtracted;
  mipsolver->mipdata_->cliquesExtracted = true;

  // extract cliques that are part of the formulation every time before probing
  // after the first call we only add cliques that directly correspond to set
  // packing constraints so that the clique merging step can extend/delete them
  if (firstCall) {
    cliquetable.extractCliques(*mipsolver);
    if (domain.infeasible()) {
      mipsolver->analysis_.mipTimerStop(kMipClockProbingPresolve);
      return Result::kPrimalInfeasible;
    }

    // during presolve we keep the objective upper bound without the current
    // offset so we need to update it

    if (mipsolver->mipdata_->upper_limit != kHighsInf) {
      double tmpLimit = mipsolver->mipdata_->upper_limit;
      mipsolver->mipdata_->upper_limit = tmpLimit - model->offset_;
      cliquetable.extractObjCliques(*mipsolver);
      mipsolver->mipdata_->upper_limit = tmpLimit;

      if (domain.infeasible()) {
        mipsolver->analysis_.mipTimerStop(kMipClockProbingPresolve);
        return Result::kPrimalInfeasible;
      }
    }

    domain.propagate();
    if (domain.infeasible()) {
      mipsolver->analysis_.mipTimerStop(kMipClockProbingPresolve);
      return Result::kPrimalInfeasible;
    }
  }

  cliquetable.cleanupFixed(domain);
  if (domain.infeasible()) {
    mipsolver->analysis_.mipTimerStop(kMipClockProbingPresolve);
    return Result::kPrimalInfeasible;
  }

  // store binary variables in vector with their number of implications on
  // other binaries
  std::vector<std::tuple<int64_t, HighsInt, HighsInt, HighsInt>> binaries;

  if (!mipsolver->mipdata_->cliquetable.isFull()) {
    binaries.reserve(model->num_col_);
    HighsRandom random(options->random_seed);
    for (HighsInt i = 0; i != model->num_col_; ++i) {
      if (domain.isBinary(i)) {
        HighsInt implicsUp = cliquetable.getNumImplications(i, 1);
        HighsInt implicsDown = cliquetable.getNumImplications(i, 0);
        binaries.emplace_back(
            -std::min(int64_t{5000}, int64_t(implicsUp) * implicsDown) /
                (1.0 + numProbes[i]),
            -std::min(HighsInt{100}, implicsUp + implicsDown), random.integer(),
            i);
      }
    }
  }
  if (!binaries.empty()) {
    // sort variables with many implications on other binaries first
    pdqsort(binaries.begin(), binaries.end());

    size_t numChangedCols = 0;
    while (domain.getChangedCols().size() != numChangedCols) {
      if (domain.isFixed(domain.getChangedCols()[numChangedCols++]))
        ++probingNumDelCol;
    }

    HighsInt numCliquesStart = cliquetable.numCliques();
    HighsInt numImplicsStart = implications.getNumImplications();
    HighsInt numDelStart = probingNumDelCol;

    HighsInt numDel = probingNumDelCol - numDelStart +
                      implications.substitutions.size() +
                      cliquetable.getSubstitutions().size();
    int64_t splayContingent =
        cliquetable.numNeighbourhoodQueries +
        std::max(mipsolver->submip ? HighsInt{0} : HighsInt{100000},
                 10 * numNonzeros());
    HighsInt numFail = 0;

    // lambda to check if model has enough continuous variables to perform
    // lifting for probing
    auto modelHasPercentageContVars = [&](size_t percentage) {
      size_t num_cols = 0, num_cont_cols = 0;
      for (size_t col = 0; col < colsize.size(); col++) {
        if (colDeleted[col]) continue;
        num_cols++;
        if (model->integrality_[col] == HighsVarType::kContinuous)
          num_cont_cols++;
      }
      return size_t{100} * num_cont_cols >= percentage * num_cols;
    };

    // collect up to 10 lifting opportunities per row
    const size_t maxNumLiftOpps = std::max(
        size_t{100000}, size_t{10} * static_cast<size_t>(model->num_row_));

    // only search for lifting opportunities if at least 2 percent of the
    // variables in the problem are continuous
    size_t numLiftOpps = 0;
    if (mipsolver->options_mip_->mip_lifting_for_probing != -1 &&
        modelHasPercentageContVars(size_t{2})) {
      // store lifting opportunities
      implications.storeLiftingOpportunity = [&](HighsInt row, HighsInt col,
                                                 HighsInt val, double coef) {
        // find lifting opportunities for row
        auto& htree = liftingOpportunities[row];
        // add element
        auto insertresult = htree.insert_or_get(std::make_pair(col, val), coef);
        assert(insertresult.second);
        numLiftOpps++;
      };
    }

    for (const auto& binvar : binaries) {
      HighsInt i = std::get<3>(binvar);

      if (cliquetable.getSubstitution(i) != nullptr || !domain.isBinary(i))
        continue;

      // when a large percentage of columns have been deleted, stop this round
      // of probing
      // if (numDel > std::max(model->num_col_ * 0.2, 1000.)) break;
      probingEarlyAbort =
          numDel >
          std::max(HighsInt{1000}, (model->num_row_ + model->num_col_) / 20);
      if (probingEarlyAbort) break;

      // break in case of too many new implications to not spent ages in
      // probing
      if (cliquetable.isFull() ||
          cliquetable.numCliques() - numCliquesStart >
              std::max(HighsInt{1000000}, 2 * numNonzeros()) ||
          implications.getNumImplications() - numImplicsStart >
              std::max(HighsInt{1000000}, 2 * numNonzeros()))
        break;

      // if (numProbed % 10 == 0)
      //   printf(
      //       "numprobed=%d  numDel=%d  newcliques=%d "
      //       "numNeighbourhoodQueries=%ld  "
      //       "splayContingent=%ld\n",
      //       numProbed, numDel, cliquetable.numCliques() - numCliquesStart,
      //       cliquetable.numNeighbourhoodQueries, splayContingent);
      if (cliquetable.numNeighbourhoodQueries > splayContingent) break;

      if (probingContingent - numProbed < 0) break;

      HighsInt numBoundChgs = 0;
      HighsInt numNewCliques = -cliquetable.numCliques();
      const bool probing_result = implications.runProbing(i, numBoundChgs);
      if (!probing_result) continue;
      probingContingent += numBoundChgs;
      numNewCliques += cliquetable.numCliques();
      numNewCliques = std::max(numNewCliques, HighsInt{0});
      while (domain.getChangedCols().size() != numChangedCols) {
        if (domain.isFixed(domain.getChangedCols()[numChangedCols++]))
          ++probingNumDelCol;
      }
      HighsInt newNumDel = probingNumDelCol - numDelStart +
                           implications.substitutions.size() +
                           cliquetable.getSubstitutions().size();

      if (newNumDel > numDel) {
        probingContingent += numDel;
        if (!mipsolver->submip) {
          splayContingent += 100 * (newNumDel + numDelStart);
          splayContingent += 1000 * numNewCliques;
        }
        numDel = newNumDel;
        numFail = 0;
      } else if (mipsolver->submip || numNewCliques == 0) {
        splayContingent -= 100 * numFail;
        ++numFail;
      } else {
        splayContingent += 1000 * numNewCliques;
        numFail = 0;
      }

      ++numProbed;
      numProbes[i] += 1;

      // Stop collecting lifting opportunities if maximum is reached
      if (numLiftOpps >= maxNumLiftOpps)
        implications.storeLiftingOpportunity = nullptr;

      // printf("nprobed: %" HIGHSINT_FORMAT ", numCliques: %" HIGHSINT_FORMAT
      // "\n", nprobed,
      //       cliquetable.numCliques());
      if (domain.infeasible()) {
        mipsolver->analysis_.mipTimerStop(kMipClockProbingPresolve);
        return Result::kPrimalInfeasible;
      }
    }

    cliquetable.cleanupFixed(domain);

    if (!firstCall) cliquetable.extractCliques(*mipsolver, false);
    cliquetable.runCliqueMerging(domain);

    // apply changes from probing

    // first delete redundant clique inequalities
    for (HighsInt delrow : cliquetable.getDeletedRows())
      if (!rowDeleted[delrow]) removeRow(delrow);
    cliquetable.getDeletedRows().clear();

    // add nonzeros from clique lifting before removing fixed variables, since
    // this might lead to stronger constraint sides
    auto& extensionvars = cliquetable.getCliqueExtensions();
    HighsInt addednnz = extensionvars.size();
    for (const auto& cliqueextension : extensionvars) {
      if (rowDeleted[cliqueextension.first]) {
        --addednnz;
        continue;
      }
      double val;
      if (cliqueextension.second.val == 0) {
        model->row_lower_[cliqueextension.first] -= 1;
        model->row_upper_[cliqueextension.first] -= 1;
        val = -1.0;
      } else
        val = 1.0;
      addToMatrix(cliqueextension.first, cliqueextension.second.col, val);
    }
    extensionvars.clear();

    // now remove fixed columns and tighten domains
    for (HighsInt i = 0; i != model->num_col_; ++i) {
      if (colDeleted[i]) continue;
      if (model->col_lower_[i] < domain.col_lower_[i])
        changeColLower(i, domain.col_lower_[i]);
      if (model->col_upper_[i] > domain.col_upper_[i])
        changeColUpper(i, domain.col_upper_[i]);
      if (domain.isFixed(i)) {
        postsolve_stack.removedFixedCol(i, model->col_lower_[i], 0.0,
                                        HighsEmptySlice());
        removeFixedCol(i);
      }
      HPRESOLVE_CHECKED_CALL(checkLimits(postsolve_stack));
    }

    // finally apply substitutions
    HPRESOLVE_CHECKED_CALL(applyConflictGraphSubstitutions(postsolve_stack));

    highsLogDev(options->log_options, HighsLogType::kInfo,
                "%" HIGHSINT_FORMAT " probing evaluations: %" HIGHSINT_FORMAT
                " deleted rows, %" HIGHSINT_FORMAT
                " deleted "
                "columns, %" HIGHSINT_FORMAT " lifted nonzeros\n",
                numProbed - oldNumProbed, numDeletedRows, numDeletedCols,
                addednnz);

    // lifting for probing
    if (mipsolver->options_mip_->mip_lifting_for_probing != -1) {
      // only perform lifting if probing did not modify the problem so far
      if (numDeletedRows == 0 && numDeletedCols == 0 && addednnz == 0)
        HPRESOLVE_CHECKED_CALL(liftingForProbing(postsolve_stack));
      // clear lifting opportunities
      liftingOpportunities.clear();
      implications.storeLiftingOpportunity = nullptr;
    }
  }

  mipsolver->analysis_.mipTimerStop(kMipClockProbingPresolve);
  return checkLimits(postsolve_stack);
}

HPresolve::Result HPresolve::liftingForProbing(
    HighsPostsolveStack& postsolve_stack) {
  // this method implements lifting for probing as described by Achterberg et
  // al. (2019) Presolve Reductions in Mixed Integer Programming. INFORMS
  // Journal on Computing 32(2):473-506.
  HighsCliqueTable& cliquetable = mipsolver->mipdata_->cliquetable;
  const HighsDomain& domain = mipsolver->mipdata_->domain;

  // collect best lifting opportunity for each row in a vector
  typedef std::pair<HighsCliqueTable::CliqueVar, double> liftingvar;
  typedef std::tuple<HighsInt, std::vector<liftingvar>, double, HighsInt>
      liftingdata;
  std::vector<liftingdata> liftingtable;
  liftingtable.reserve(liftingOpportunities.size());

  // remember overall best score
  double bestscoretotal = -kHighsInf;

  // is lifting allowed to add non-zeros?
  const bool fillallowed = mipsolver->options_mip_->mip_lifting_for_probing > 0;

  // store clique variables and coefficients in a map
  auto comp = [](const HighsCliqueTable::CliqueVar& c1,
                 const HighsCliqueTable::CliqueVar& c2) {
    return c1.col < c2.col || (c1.col == c2.col && c1.val < c2.val);
  };
  std::map<HighsCliqueTable::CliqueVar, std::pair<double, HighsInt>,
           decltype(comp)>
      coefficients(comp);

  // consider lifting opportunities
  size_t numrowsremoved = 0;
  for (const auto& elm : liftingOpportunities) {
    // get row index and skip deleted rows
    HighsInt row = elm.first;
    if (rowDeleted[row]) continue;

    // do not add non-zeros to dense rows
    const bool dense =
        rowsize[row] >
        std::max(HighsInt{1000}, (model->num_col_ - numDeletedCols) / 20);

    // iterate over elements in hash tree
    const auto& htree = elm.second;
    bool isredundant = false;
    htree.for_each([&](const std::pair<HighsInt, HighsInt>& data, double coef) {
      HighsInt col = data.first;
      HighsInt val = data.second;
      HighsInt pos = findNonzero(row, col);
      isredundant = isredundant || htree.contains(std::make_pair(col, 1 - val));
      if (!dense && (fillallowed || pos != -1) && !colDeleted[col] &&
          !domain.isFixed(col))
        coefficients[HighsCliqueTable::CliqueVar{col, val}] = {coef, pos};
    });

    // remove redundant rows
    if (isredundant) {
      numrowsremoved++;
      postsolve_stack.redundantRow(row);
      removeRow(row);
      coefficients.clear();
      HPRESOLVE_CHECKED_CALL(checkLimits(postsolve_stack));
      continue;
    }

    // skip rows with empty coefficient map
    if (coefficients.empty()) continue;

    // vector to hold best clique
    std::vector<liftingvar> bestclique;
    double bestscore = -kHighsInf;
    HighsInt bestnfill = 0;

    // lambda for computing coefficient difference
    auto computeCoeffDiff = [&](double newvalue, HighsInt nzpos) {
      return std::fabs(newvalue - (nzpos == -1 ? 0 : Avalue[nzpos]));
    };

    // store candidates in a vector
    std::vector<HighsCliqueTable::CliqueVar> candidates;
    candidates.reserve(coefficients.size());
    for (const auto& elm : coefficients) {
      candidates.push_back(elm.first);
      // initialize best clique
      double score = computeCoeffDiff(elm.second.first, elm.second.second);
      if (score > bestscore) {
        bestscore = score;
        bestnfill = elm.second.second == -1 ? 1 : 0;
        bestclique = {std::make_pair(elm.first, elm.second.first)};
      }
    }

    if (candidates.size() > 1) {
      // compute cliques
      auto cliques =
          cliquetable.computeMaximalCliques(candidates, primal_feastol);

      // identify clique with highest score
      for (const auto& clique : cliques) {
        HighsCDouble score = 0;
        HighsInt nfill = 0;
        for (const auto& cliquevar : clique) {
          score += computeCoeffDiff(coefficients[cliquevar].first,
                                    coefficients[cliquevar].second);
          if (coefficients[cliquevar].second == -1) nfill++;
        }
        if (score > bestscore) {
          bestscore = static_cast<double>(score);
          bestnfill = nfill;
          bestclique.clear();
          bestclique.reserve(clique.size());
          for (const auto& cliquevar : clique) {
            bestclique.emplace_back(cliquevar, coefficients[cliquevar].first);
          }
        }
      }
    }

    // store best clique
    liftingtable.emplace_back(row, bestclique, bestscore, bestnfill);
    bestscoretotal = std::max(bestscoretotal, bestscore);
    coefficients.clear();
  }

  // lambda for computing score
  auto computeOverallScore = [&](double score, HighsInt numelms,
                                 HighsInt numfillin) {
    const double weight = 0.5;
    return weight * (score / bestscoretotal) +
           (1 - weight) * static_cast<double>(numelms - numfillin) /
               static_cast<double>(numelms);
  };

  // sort according to score
  pdqsort(
      liftingtable.begin(), liftingtable.end(),
      [&](const liftingdata& opp1, const liftingdata& opp2) {
        double score1 = computeOverallScore(
            std::get<2>(opp1), static_cast<HighsInt>(std::get<1>(opp1).size()),
            std::get<3>(opp1));
        double score2 = computeOverallScore(
            std::get<2>(opp2), static_cast<HighsInt>(std::get<1>(opp2).size()),
            std::get<3>(opp2));
        return (score1 == score2 ? std::get<0>(opp1) < std::get<0>(opp2)
                                 : score1 > score2);
      });

  // perform actual lifting
  size_t nfill = 0;
  size_t nmod = 0;
  size_t numrowsmodified = 0;
  const size_t maxnfill = std::max(10 * liftingtable.size(),
                                   static_cast<size_t>(numNonzeros()) / 100);
  for (const auto& lifting : liftingtable) {
    // get clique
    HighsInt row = std::get<0>(lifting);
    const auto& bestclique = std::get<1>(lifting);

    // check against max. fill-in
    size_t newfill = static_cast<size_t>(std::get<3>(lifting));
    if (nfill + newfill > maxnfill) break;
    nfill += newfill;
    nmod += bestclique.size() - newfill;

    // update matrix
    HighsCDouble update = 0.0;
    for (const auto& elm : bestclique) {
      // get data
      const auto& cliquevar = std::get<0>(elm);
      const double& coeff = std::get<1>(elm);
      // add non-zero to matrix
      addToMatrix(row, cliquevar.col, coeff);
      // compute term to update left-hand / right-hand side
      if (cliquevar.val == 0) update += coeff;
    }

    // update left-hand / right-hand sides
    numrowsmodified++;
    if (model->row_lower_[row] != -kHighsInf)
      model->row_lower_[row] += static_cast<double>(update);
    if (model->row_upper_[row] != kHighsInf)
      model->row_upper_[row] += static_cast<double>(update);
  }

  highsLogDev(options->log_options, HighsLogType::kInfo,
              "Lifting for probing removed %d and modified %d row(s), added %d "
              "new and modified %d existing nonzero(s)\n",
              static_cast<int>(numrowsremoved),
              static_cast<int>(numrowsmodified), static_cast<int>(nfill),
              static_cast<int>(nmod));

  return Result::kOk;
}

void HPresolve::addToMatrix(const HighsInt row, const HighsInt col,
                            const double val) {
  HighsInt pos = findNonzero(row, col);

  markChangedRow(row);
  markChangedCol(col);

  if (pos == -1) {
    if (freeslots.empty()) {
      pos = Avalue.size();
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

    // remove implied bounds on row duals that where implied by this column's
    // dual constraint
    resetRowDualImpliedBoundsDerivedFromCol(col);

    // remove implied bounds on columns that where implied by this row
    resetColImpliedBoundsDerivedFromRow(row);

    // modifications to row invalidate lifting opportunities
    clearLiftingOpportunities(row);

  } else {
    double sum = Avalue[pos] + val;
    if (std::abs(sum) <= options->small_matrix_value) {
      unlink(pos);
    } else {
      // remove implied bounds on row duals that where implied by this column's
      // dual constraint
      resetRowDualImpliedBoundsDerivedFromCol(col);

      // remove implied bounds on columns that where implied by this row
      resetColImpliedBoundsDerivedFromRow(row);

      // modifications to row invalidate lifting opportunities
      clearLiftingOpportunities(row);

      // remove the locks and contribution to implied (dual) row bounds, then
      // add then again
      impliedRowBounds.remove(row, col, Avalue[pos]);
      impliedDualRowBounds.remove(col, row, Avalue[pos]);
      Avalue[pos] = sum;
      // value not zero, add new contributions and locks with opposite sign
      impliedRowBounds.add(row, col, Avalue[pos]);
      impliedDualRowBounds.add(col, row, Avalue[pos]);
    }
  }
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
}

void HPresolve::markColDeleted(HighsInt col) {
  assert(!colDeleted[col]);

  // prevents col from being added to change vector
  changedColFlag[col] = true;
  colDeleted[col] = true;
  ++numDeletedCols;
}

void HPresolve::changeColUpper(HighsInt col, double newUpper) {
  if (model->integrality_[col] != HighsVarType::kContinuous) {
    newUpper = std::floor(newUpper + primal_feastol);
    if (newUpper == model->col_upper_[col]) return;
  }

  double oldUpper = model->col_upper_[col];
  model->col_upper_[col] = newUpper;

  for (const HighsSliceNonzero& nonzero : getColumnVector(col)) {
    impliedRowBounds.updatedVarUpper(nonzero.index(), col, nonzero.value(),
                                     oldUpper);
    markChangedRow(nonzero.index());
  }
}

void HPresolve::changeColLower(HighsInt col, double newLower) {
  if (model->integrality_[col] != HighsVarType::kContinuous) {
    newLower = std::ceil(newLower - primal_feastol);
    if (newLower == model->col_lower_[col]) return;
  }

  double oldLower = model->col_lower_[col];
  model->col_lower_[col] = newLower;
  // printf("tightening lower bound of column %" HIGHSINT_FORMAT " from %.15g to
  // %.15g\n", col,
  //        oldLower, newLower);

  for (const HighsSliceNonzero& nonzero : getColumnVector(col)) {
    impliedRowBounds.updatedVarLower(nonzero.index(), col, nonzero.value(),
                                     oldLower);
    markChangedRow(nonzero.index());
  }
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

void HPresolve::scaleMIP(HighsPostsolveStack& postsolve_stack) {
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

    transformColumn(postsolve_stack, i, scale, 0.0);
  }
}

HPresolve::Result HPresolve::applyConflictGraphSubstitutions(
    HighsPostsolveStack& postsolve_stack) {
  HighsCliqueTable& cliquetable = mipsolver->mipdata_->cliquetable;
  HighsImplications& implications = mipsolver->mipdata_->implications;
  for (const auto& substitution : implications.substitutions) {
    if (colDeleted[substitution.substcol] || colDeleted[substitution.staycol])
      continue;

    ++probingNumDelCol;

    postsolve_stack.doubletonEquation(
        -1, substitution.substcol, substitution.staycol, 1.0,
        -substitution.scale, substitution.offset,
        model->col_lower_[substitution.substcol],
        model->col_upper_[substitution.substcol], 0.0, false, false,
        HighsPostsolveStack::RowType::kEq, HighsEmptySlice());
    markColDeleted(substitution.substcol);
    substitute(substitution.substcol, substitution.staycol, substitution.offset,
               substitution.scale);
    HPRESOLVE_CHECKED_CALL(checkLimits(postsolve_stack));
  }

  implications.substitutions.clear();

  for (HighsCliqueTable::Substitution subst : cliquetable.getSubstitutions()) {
    if (colDeleted[subst.substcol] || colDeleted[subst.replace.col]) continue;

    double scale;
    double offset;

    ++probingNumDelCol;

    if (subst.replace.val == 0) {
      scale = -1.0;
      offset = 1.0;
    } else {
      scale = 1.0;
      offset = 0.0;
    }

    postsolve_stack.doubletonEquation(
        -1, subst.substcol, subst.replace.col, 1.0, -scale, offset,
        model->col_lower_[subst.substcol], model->col_upper_[subst.substcol],
        0.0, false, false, HighsPostsolveStack::RowType::kEq,
        HighsEmptySlice());
    markColDeleted(subst.substcol);
    substitute(subst.substcol, subst.replace.col, offset, scale);
    HPRESOLVE_CHECKED_CALL(checkLimits(postsolve_stack));
  }

  cliquetable.getSubstitutions().clear();

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

  HighsInt ncol = Astart.size() - 1;
  assert(ncol == int(colhead.size()));
  HighsInt nnz = Aval.size();

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
  for (HighsInt pos = 0; pos != nnz; ++pos) link(pos);

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

  HighsInt nrow = ARstart.size() - 1;
  assert(nrow == int(rowroot.size()));
  HighsInt nnz = ARval.size();

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
  for (HighsInt pos = 0; pos != nnz; ++pos) link(pos);

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

HighsInt HPresolve::countFillin(HighsInt row) {
  HighsInt fillin = 0;
  for (HighsInt rowiter : rowpositions) {
    if (findNonzero(row, Acol[rowiter]) == -1) fillin += 1;
  }

  return fillin;
}

bool HPresolve::checkFillin(HighsHashTable<HighsInt, HighsInt>& fillinCache,
                            HighsInt row, HighsInt col) {
  // check numerics against markowitz tolerance
  assert(rowpositions.size() == static_cast<size_t>(rowsize[row]));

  // check fillin against max fillin
  HighsInt fillin = -(rowsize[row] + colsize[col] - 1);

#if 1
  // first use fillin for rows where it is already computed
  for (HighsInt coliter = colhead[col]; coliter != -1;
       coliter = Anext[coliter]) {
    if (Arow[coliter] == row) continue;

    auto cachedFillin = fillinCache.find(Arow[coliter]);
    if (cachedFillin == nullptr) continue;

    fillin += (*cachedFillin - 1);
    if (fillin > options->presolve_substitution_maxfillin) return false;
  }

  // iterate over rows of substituted column again to count the fillin for the
  // remaining rows
  for (HighsInt coliter = colhead[col]; coliter != -1;
       coliter = Anext[coliter]) {
    assert(Acol[coliter] == col);

    if (Arow[coliter] == row) continue;

    HighsInt& cachedFillin = fillinCache[Arow[coliter]];

    if (cachedFillin != 0) continue;

    HighsInt rowfillin = countFillin(Arow[coliter]);
    cachedFillin = rowfillin + 1;
    fillin += rowfillin;

    if (fillin > options->presolve_substitution_maxfillin) return false;
    // we count a fillin of 1 if the column is not present in the row and
    // a fillin of zero otherwise. the fillin for the substituted column
    // itself was already counted before the loop so we skip that entry.
  }
#else
  for (HighsInt rowiter : rowpositions) {
    if (rowiter == pos) continue;
    for (coliter = colhead[col]; coliter != -1; coliter = Anext[coliter]) {
      assert(Acol[coliter] == col);

      if (rowiter != coliter &&
          findNonzero(Arow[coliter], Acol[rowiter]) == -1) {
        if (fillin == maxfillin) return false;
        fillin += 1;
      }
    }
  }
#endif

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

void HPresolve::transformColumn(HighsPostsolveStack& postsolve_stack,
                                HighsInt col, double scale, double constant) {
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
  if (model->integrality_[col] != HighsVarType::kContinuous) {
    changeColLower(col, model->col_lower_[col]);
    changeColUpper(col, model->col_upper_[col]);
  }

  markChangedCol(col);
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

  if (integral) {
    if (model->row_upper_[row] != kHighsInf)
      model->row_upper_[row] = std::round(model->row_upper_[row]);
    if (model->row_lower_[row] != -kHighsInf)
      model->row_lower_[row] = std::round(model->row_lower_[row]);
  }

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
}

void HPresolve::substitute(HighsInt row, HighsInt col, double rhs) {
  assert(!rowDeleted[row]);
  assert(!colDeleted[col]);
  HighsInt pos = findNonzero(row, col);
  assert(pos != -1);

  assert(Arow[pos] == row);
  assert(Acol[pos] == col);
  double substrowscale = -1.0 / Avalue[pos];
  assert(isImpliedFree(col));

  markRowDeleted(row);
  markColDeleted(col);

  // substitute the column in each row where it occurs
  for (HighsInt coliter = colhead[col]; coliter != -1;) {
    HighsInt colrow = Arow[coliter];
    double colval = Avalue[coliter];

    // walk to the next position before doing any modifications, because
    // the current position will be deleted in the loop below
    assert(Acol[coliter] == col);
    HighsInt colpos = coliter;
    coliter = Anext[coliter];

    // skip the row that is used for substitution
    if (row == colrow) continue;

    assert(findNonzero(colrow, col) != -1);

    // cancels out and bounds of dual row for this column do not need to be
    // updated
    unlink(colpos);

    // printf("\nbefore substitution: ");
    // debugPrintRow(colrow);

    // determine the scale for the substitution row for addition to this row
    double scale = colval * substrowscale;

    // adjust the sides
    if (model->row_lower_[colrow] != -kHighsInf)
      model->row_lower_[colrow] += scale * rhs;

    if (model->row_upper_[colrow] != kHighsInf)
      model->row_upper_[colrow] += scale * rhs;

    for (HighsInt rowiter : rowpositions) {
      assert(Arow[rowiter] == row);

      if (Acol[rowiter] != col)
        addToMatrix(colrow, Acol[rowiter], scale * Avalue[rowiter]);
    }

    // check if this is an equation row and it now has a different size
    reinsertEquation(colrow);
    // printf("after substitution: ");
    // debugPrintRow(colrow);
  }

  assert(colsize[col] == 1);

  // substitute column in the objective function
  if (model->col_cost_[col] != 0.0) {
    HighsCDouble objscale = model->col_cost_[col] * substrowscale;
    model->offset_ = static_cast<double>(model->offset_ - objscale * rhs);
    assert(std::isfinite(model->offset_));
    for (HighsInt rowiter : rowpositions) {
      // printf("changing col cost to %g = %g + %g * %g\n",
      // double(model->col_cost_[Acol[rowiter]] + objscale * Avalue[rowiter]),
      // model->col_cost_[Acol[rowiter]], double(objscale), Avalue[rowiter]);
      model->col_cost_[Acol[rowiter]] = static_cast<double>(
          model->col_cost_[Acol[rowiter]] + objscale * Avalue[rowiter]);
      if (std::abs(model->col_cost_[Acol[rowiter]]) <=
          options->small_matrix_value)
        model->col_cost_[Acol[rowiter]] = 0.0;
    }
    assert(std::abs(model->col_cost_[col]) <=
           std::max(options->dual_feasibility_tolerance,
                    kHighsTiny * std::abs(static_cast<double>(objscale))));
    model->col_cost_[col] = 0.0;
  }

  // finally remove the entries of the row that was used for substitution
  for (HighsInt rowiter : rowpositions) unlink(rowiter);
}

void HPresolve::toCSC(std::vector<double>& Aval, std::vector<HighsInt>& Aindex,
                      std::vector<HighsInt>& Astart) {
  // set up the column starts using the column size array
  HighsInt numcol = colsize.size();
  Astart.resize(numcol + 1);
  HighsInt nnz = 0;
  for (HighsInt i = 0; i != numcol; ++i) {
    Astart[i] = nnz;
    nnz += colsize[i];
  }
  Astart[numcol] = nnz;

  // now setup the entries of the CSC matrix
  // we reuse the colsize array to count down to zero
  // for determining the position of each nonzero
  Aval.resize(nnz);
  Aindex.resize(nnz);
  HighsInt numslots = Avalue.size();
  assert(numslots - int(freeslots.size()) == nnz);
  for (HighsInt i = 0; i != numslots; ++i) {
    if (Avalue[i] == 0.0) continue;
    assert(Acol[i] >= 0 && Acol[i] < model->num_col_);
    HighsInt pos = Astart[Acol[i] + 1] - colsize[Acol[i]];
    --colsize[Acol[i]];
    assert(colsize[Acol[i]] >= 0);
    Aval[pos] = Avalue[i];
    Aindex[pos] = Arow[i];
  }
}

void HPresolve::toCSR(std::vector<double>& ARval,
                      std::vector<HighsInt>& ARindex,
                      std::vector<HighsInt>& ARstart) {
  // set up the row starts using the row size array
  HighsInt numrow = rowsize.size();
  ARstart.resize(numrow + 1);
  HighsInt nnz = 0;
  for (HighsInt i = 0; i != numrow; ++i) {
    ARstart[i] = nnz;
    nnz += rowsize[i];
  }
  ARstart[numrow] = nnz;

  // now setup the entries of the CSC matrix
  // we reuse the colsize array to count down to zero
  // for determining the position of each nonzero
  ARval.resize(nnz);
  ARindex.resize(nnz);
  for (HighsInt i = 0; i != nnz; ++i) {
    if (Avalue[i] == 0.0) continue;
    HighsInt pos = ARstart[Arow[i] + 1] - rowsize[Arow[i]];
    --rowsize[Arow[i]];
    assert(rowsize[Arow[i]] >= 0);
    ARval[pos] = Avalue[i];
    ARindex[pos] = Acol[i];
  }
}

HPresolve::Result HPresolve::doubletonEq(HighsPostsolveStack& postsolve_stack,
                                         HighsInt row,
                                         HighsPostsolveStack::RowType rowType) {
  assert(analysis_.allow_rule_[kPresolveRuleDoubletonEquation]);
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
  if (lowerTightened) changeColLower(staycol, stayImplLower);

  bool upperTightened = stayImplUpper < oldStayUpper - primal_feastol;
  if (upperTightened) changeColUpper(staycol, stayImplUpper);

  postsolve_stack.doubletonEquation(
      row, substcol, staycol, substcoef, staycoef, rhs, substLower, substUpper,
      model->col_cost_[substcol], lowerTightened, upperTightened, rowType,
      getColumnVector(substcol));

  // finally modify matrix
  markColDeleted(substcol);
  removeRow(row);
  substitute(substcol, staycol, rhs / substcoef, -staycoef / substcoef);

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
  // delete row singleton nonzero directly, we have all information that we need
  // in local variables
  markRowDeleted(row);
  unlink(nzPos);

  // check for simple
  if (val > 0) {
    if (model->col_upper_[col] * val <=
            model->row_upper_[row] + primal_feastol &&
        model->col_lower_[col] * val >=
            model->row_lower_[row] - primal_feastol) {
      postsolve_stack.redundantRow(row);
      analysis_.logging_on_ = logging_on;
      if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleSingletonRow);
      return checkLimits(postsolve_stack);
    }
  } else {
    if (model->col_lower_[col] * val <=
            model->row_upper_[row] + primal_feastol &&
        model->col_upper_[col] * val >=
            model->row_lower_[row] - primal_feastol) {
      postsolve_stack.redundantRow(row);
      analysis_.logging_on_ = logging_on;
      if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleSingletonRow);
      return checkLimits(postsolve_stack);
    }
  }

  // zeros should not be linked in the matrix
  assert(std::fabs(val) > options->small_matrix_value);

  double newColUpper = kHighsInf;
  double newColLower = -kHighsInf;
  if (val > 0) {
    if (model->row_upper_[row] != kHighsInf)
      newColUpper = model->row_upper_[row] / val;
    if (model->row_lower_[row] != -kHighsInf)
      newColLower = model->row_lower_[row] / val;
  } else {
    if (model->row_upper_[row] != kHighsInf)
      newColLower = model->row_upper_[row] / val;
    if (model->row_lower_[row] != -kHighsInf)
      newColUpper = model->row_lower_[row] / val;
  }

  // use either the primal feasibility tolerance for the bound constraint or
  // for the singleton row including scaling, whichever is tighter.
  const double boundTol = primal_feastol / std::max(1.0, std::fabs(val));
  const bool isIntegral = model->integrality_[col] != HighsVarType::kContinuous;

  bool lowerTightened = newColLower > model->col_lower_[col] + boundTol;
  bool upperTightened = newColUpper < model->col_upper_[col] - boundTol;

  double lb, ub;
  if (lowerTightened) {
    if (isIntegral) newColLower = std::ceil(newColLower - boundTol);
    lb = newColLower;
  } else
    lb = model->col_lower_[col];

  if (upperTightened) {
    if (isIntegral) newColUpper = std::floor(newColUpper + boundTol);
    ub = newColUpper;
  } else
    ub = model->col_upper_[col];

  // printf("old bounds [%.15g,%.15g], new bounds [%.15g,%.15g] ... ",
  //        model->col_lower_[col], model->col_upper_[col], lb, ub);
  // check whether the bounds are equal in tolerances
  if (ub <= lb + primal_feastol) {
    // bounds could be infeasible or equal in tolerances, first check infeasible
    if (ub < lb - primal_feastol) return Result::kPrimalInfeasible;

    // bounds are equal in tolerances, if they have a slight infeasibility below
    // those tolerances or they have a slight numerical distance which changes
    // the largest contribution below feasibility tolerance then we can safely
    // set the bound to one of the values. To heuristically get rid of numerical
    // errors we choose the bound that was not tightened, or the midpoint if
    // both where tightened.
    if (ub < lb || (ub > lb && (ub - lb) * std::max(std::fabs(val),
                                                    getMaxAbsColVal(col)) <=
                                   primal_feastol)) {
      if (lowerTightened && upperTightened) {
        ub = 0.5 * (ub + lb);
        lb = ub;
        lowerTightened = lb > model->col_lower_[col];
        upperTightened = ub < model->col_upper_[col];
      } else if (lowerTightened) {
        lb = ub;
        lowerTightened = lb > model->col_lower_[col];
      } else {
        ub = lb;
        upperTightened = ub < model->col_upper_[col];
      }
    }
  }

  // printf("final bounds: [%.15g,%.15g]\n", lb, ub);

  postsolve_stack.singletonRow(row, col, val, lowerTightened, upperTightened);

  // just update bounds (and row activities)
  if (lowerTightened) changeColLower(col, lb);
  // update bounds, or remove as fixed column directly
  if (ub == lb) {
    postsolve_stack.removedFixedCol(col, lb, model->col_cost_[col],
                                    getColumnVector(col));
    removeFixedCol(col);
  } else if (upperTightened)
    changeColUpper(col, ub);

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
                                          HighsInt col) {
  assert(colsize[col] == 1);
  assert(!colDeleted[col]);
  HighsInt nzPos = colhead[col];
  HighsInt row = Arow[nzPos];
  double colCoef = Avalue[nzPos];

  if (rowsize[row] == 1) {
    HPRESOLVE_CHECKED_CALL(singletonRow(postsolve_stack, row););

    if (!colDeleted[col]) {
      assert(colsize[col] == 0);
      return emptyCol(postsolve_stack, col);
    }
    return Result::kOk;
  }

  // detect strong / weak domination
  HPRESOLVE_CHECKED_CALL(detectDominatedCol(postsolve_stack, col, false));
  if (colDeleted[col]) return Result::kOk;

  if (mipsolver != nullptr)
    HPRESOLVE_CHECKED_CALL(
        static_cast<Result>(convertImpliedInteger(col, row)));

  updateColImpliedBounds(row, col, colCoef);

  if (model->integrality_[col] != HighsVarType::kInteger)
    updateRowDualImpliedBounds(row, col, colCoef);

  // now check if column is implied free within an equation and substitute the
  // column if that is the case
  if (isDualImpliedFree(row) && isImpliedFree(col) &&
      analysis_.allow_rule_[kPresolveRuleFreeColSubstitution]) {
    if (model->integrality_[col] == HighsVarType::kInteger) {
      StatusResult impliedIntegral = isImpliedIntegral(col);
      HPRESOLVE_CHECKED_CALL(static_cast<Result>(impliedIntegral));
      if (!impliedIntegral) return Result::kOk;
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
    return checkLimits(postsolve_stack);
  }

  // todo: check for zero cost singleton and remove
  return Result::kOk;
}

void HPresolve::substituteFreeCol(HighsPostsolveStack& postsolve_stack,
                                  HighsInt row, HighsInt col,
                                  bool relaxRowDualBounds) {
  assert(!rowDeleted[row]);
  assert(!colDeleted[col]);
  assert(isDualImpliedFree(row));
  assert(isImpliedFree(col));

  HighsPostsolveStack::RowType rowType;
  double rhs;
  dualImpliedFreeGetRhsAndRowType(row, rhs, rowType, relaxRowDualBounds);

  postsolve_stack.freeColSubstitution(row, col, rhs, model->col_cost_[col],
                                      rowType, getStoredRow(),
                                      getColumnVector(col));

  // todo, check integrality of coefficients and allow this
  substitute(row, col, rhs);
}

HPresolve::Result HPresolve::rowPresolve(HighsPostsolveStack& postsolve_stack,
                                         HighsInt row) {
  assert(!rowDeleted[row]);

  const bool logging_on = analysis_.logging_on_;
  // handle special cases directly via a call to the specialized procedure
  switch (rowsize[row]) {
    default:
      break;
    case 0:
      if (model->row_upper_[row] < -primal_feastol ||
          model->row_lower_[row] > primal_feastol)
        // model infeasible
        return Result::kPrimalInfeasible;
      if (logging_on) analysis_.startPresolveRuleLog(kPresolveRuleEmptyRow);
      postsolve_stack.redundantRow(row);
      markRowDeleted(row);
      analysis_.logging_on_ = logging_on;
      if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleEmptyRow);
      return checkLimits(postsolve_stack);
    case 1:
      return singletonRow(postsolve_stack, row);
  }

  // printf("row presolve: ");
  // debugPrintRow(row);
  double impliedRowUpper = impliedRowBounds.getSumUpper(row);
  double impliedRowLower = impliedRowBounds.getSumLower(row);

  // Allow removal of redundant rows
  if (impliedRowLower > model->row_upper_[row] + primal_feastol ||
      impliedRowUpper < model->row_lower_[row] - primal_feastol) {
    // model infeasible
    return Result::kPrimalInfeasible;
  }

  if (impliedRowLower >= model->row_lower_[row] - primal_feastol &&
      impliedRowUpper <= model->row_upper_[row] + primal_feastol) {
    // row is redundant
    if (logging_on) analysis_.startPresolveRuleLog(kPresolveRuleRedundantRow);
    postsolve_stack.redundantRow(row);
    removeRow(row);
    analysis_.logging_on_ = logging_on;
    if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleRedundantRow);
    return checkLimits(postsolve_stack);
  }

  auto checkRedundantBounds = [&](HighsInt col, HighsInt row) {
    // check if column singleton has redundant bounds
    assert(model->col_cost_[col] != 0.0);
    if (colsize[col] != 1) return;
    if (model->col_cost_[col] > 0) {
      assert(model->col_lower_[col] == -kHighsInf || !isLowerImplied(col) ||
             colLowerSource[col] == row);
      if (model->col_lower_[col] > implColLower[col] - primal_feastol)
        changeColLower(col, -kHighsInf);
    } else {
      assert(model->col_upper_[col] == kHighsInf || !isUpperImplied(col) ||
             colUpperSource[col] == row);
      if (model->col_upper_[col] < implColUpper[col] + primal_feastol)
        changeColUpper(col, kHighsInf);
    }
  };

  // Store original bounds
  double origRowUpper = model->row_upper_[row];
  double origRowLower = model->row_lower_[row];

  if (!isEquation(row)) {
    if (isImpliedEquationAtLower(row)) {
      // Convert to equality constraint (note that currently postsolve will not
      // know about this conversion)
      model->row_upper_[row] = model->row_lower_[row];
      // Since row upper bound is now finite, lower bound on row dual is
      // -kHighsInf
      changeRowDualLower(row, -kHighsInf);
      if (mipsolver == nullptr)
        checkRedundantBounds(rowDualLowerSource[row], row);
    } else if (isImpliedEquationAtUpper(row)) {
      // Convert to equality constraint (note that currently postsolve will not
      // know about this conversion)
      model->row_lower_[row] = model->row_upper_[row];
      // Since row lower bound is now finite, upper bound on row dual is
      // kHighsInf
      changeRowDualUpper(row, kHighsInf);
      if (mipsolver == nullptr)
        checkRedundantBounds(rowDualUpperSource[row], row);
    }
  }

  // Get row bounds
  double rowUpper = model->row_upper_[row];
  double rowLower = model->row_lower_[row];

  // Handle doubleton equations
  if (rowsize[row] == 2 && rowLower == rowUpper &&
      analysis_.allow_rule_[kPresolveRuleDoubletonEquation]) {
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
      if (impliedRowLower != -kHighsInf && impliedRowUpper != kHighsInf &&
          std::abs(impliedRowLower + impliedRowUpper - 2 * rowUpper) <=
              options->small_matrix_value) {
        double binCoef = std::abs(impliedRowUpper - rowUpper);
        // simple probing on equation case
        HighsInt binCol = -1;
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
        // Reduction uses substitution involving range of all columns
        // other than the binary. This is not well defined when any of
        // the columns is not boxed, so look for non-boxed columns
        // Exposed as #1280
        bool all_boxed_column = true;
        for (const HighsSliceNonzero& nonz : getStoredRow()) {
          if (model->col_lower_[nonz.index()] <= -kHighsInf ||
              model->col_upper_[nonz.index()] >= kHighsInf) {
            all_boxed_column = false;
            break;
          }
        }
        if (binCol != -1 && all_boxed_column) {
          // found binary column for substituting all other columns
          // printf("simple probing case on row of size %" HIGHSINT_FORMAT "\n",
          // rowsize[row]);
          for (const HighsSliceNonzero& nonz : getStoredRow()) {
            if (nonz.index() == binCol) continue;

            if (model->col_lower_[nonz.index()] ==
                model->col_upper_[nonz.index()]) {
              postsolve_stack.removedFixedCol(nonz.index(),
                                              model->col_lower_[nonz.index()],
                                              0.0, HighsEmptySlice());
              removeFixedCol(nonz.index());
              continue;
            }

            auto remDoubletonEq = [&](HighsInt col, HighsInt binCol,
                                      HighsInt direction) {
              double bound = direction >= 0 ? model->col_upper_[col]
                                            : model->col_lower_[col];
              double scale =
                  direction * (model->col_lower_[col] - model->col_upper_[col]);
              double offset = bound - model->col_lower_[binCol] * scale;
              postsolve_stack.doubletonEquation(
                  -1, col, binCol, 1.0, -scale, offset, model->col_lower_[col],
                  model->col_upper_[col], 0.0, false, false,
                  HighsPostsolveStack::RowType::kEq, HighsEmptySlice());
              substitute(col, binCol, offset, scale);
            };

            if (std::signbit(binCoef) == std::signbit(nonz.value())) {
              // binary coefficient is positive:
              // setting the binary to its upper bound
              // increases the minimal activity to be equal to the row upper
              // bound and there for all other variables are fixed to the bound
              // that contributes to the rows minimal activity, i.e. the lower
              // bound for a positive coefficient

              // This case yields the following implications:
              // binCol = ub -> nonzCol = lb
              // binCol = lb -> nonzCol = ub
              // as linear equation:
              // nonzCol = colUb - (colUb - colLb)(binCol - binLb)
              // nonzCol = colUb + binLb * (colUb - colLb) - (colUb - colLb) *
              // binCol
              remDoubletonEq(nonz.index(), binCol, HighsInt{1});
            } else {
              // This case yields the following implications:
              // binCol = lb -> nonzCol = lb
              // binCol = ub -> nonzCol = ub
              // as linear equation:
              // nonzCol = colLb + (colUb - colLb)(binCol - binLb)
              // nonzCol =
              //    colLb - binLb*(colUb - colLb) + (colUb - colLb)*binCol
              remDoubletonEq(nonz.index(), binCol, HighsInt{-1});
            }
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
              transformColumn(postsolve_stack, continuousCol, scale, 0.0);

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

            HighsInt rowlen = rowpositions.size();
            HighsInt x1Cand = -1;
            int64_t d = 0;

            for (HighsInt i = 0; i < rowlen; ++i) {
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
                transformColumn(postsolve_stack, x1, d, 0.0);
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
                  if (fixVal > model->col_lower_[x1])
                    changeColLower(x1, fixVal);
                  if (fixVal < model->col_upper_[x1])
                    changeColUpper(x1, fixVal);
                  // Fix variable
                  if (std::abs(model->col_lower_[x1] - fixVal) <=
                      primal_feastol) {
                    if (fixColToLowerOrUnbounded(postsolve_stack, x1)) {
                      // Handle unboundedness
                      presolve_status_ =
                          HighsPresolveStatus::kUnboundedOrInfeasible;
                      return Result::kDualInfeasible;
                    }
                  } else {
                    assert(std::abs(model->col_upper_[x1] - fixVal) <=
                           primal_feastol);
                    if (fixColToUpperOrUnbounded(postsolve_stack, x1)) {
                      // Handle unboundedness
                      presolve_status_ =
                          HighsPresolveStatus::kUnboundedOrInfeasible;
                      return Result::kDualInfeasible;
                    }
                  }
                  rowpositions.erase(rowpositions.begin() + x1Cand);
                } else {
                  transformColumn(postsolve_stack, x1, d, b);
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

        double deltaDown = model->row_lower_[row] == -kHighsInf
                               ? primal_feastol
                               : options->small_matrix_value;
        double deltaUp = model->row_upper_[row] == kHighsInf
                             ? primal_feastol
                             : options->small_matrix_value;

        for (const HighsSliceNonzero& nonz : getStoredRow()) {
          assert(nonz.value() != 0.0);
          rowCoefs.push_back(nonz.value());
          rowIndex.push_back(nonz.index());
        }

        double intScale =
            HighsIntegers::integralScale(rowCoefs, deltaDown, deltaUp);

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
            model->row_lower_[row] = static_cast<double>(roundedLhs / scalar);
          if (roundedRhs != kHighsInf)
            model->row_upper_[row] = static_cast<double>(roundedRhs / scalar);
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
                    model->row_upper_[row] = static_cast<double>(roundedRhs);
                  if (roundedLhs > model->row_lower_[row] + primal_feastol)
                    model->row_lower_[row] = static_cast<double>(roundedLhs);
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
          impliedRowUpper = impliedRowBounds.getSumUpper(row);
          impliedRowLower = impliedRowBounds.getSumLower(row);
        }
      }

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

      if (model->row_lower_[row] == -kHighsInf &&
          impliedRowUpper != kHighsInf) {
        // <= constraint: try to strengthen coefficients
        HighsCDouble rhs = model->row_upper_[row];
        strengthenCoefs(rhs, HighsInt{1},
                        static_cast<HighsCDouble>(impliedRowUpper) -
                            model->row_upper_[row]);
        model->row_upper_[row] = static_cast<double>(rhs);
      }

      if (model->row_upper_[row] == kHighsInf &&
          impliedRowLower != -kHighsInf) {
        // >= constraint: try to strengthen coefficients
        HighsCDouble rhs = model->row_lower_[row];
        strengthenCoefs(rhs, HighsInt{-1},
                        model->row_lower_[row] -
                            static_cast<HighsCDouble>(impliedRowLower));
        model->row_lower_[row] = static_cast<double>(rhs);
      }
    }
  }  // if (rowsizeInteger[row] != 0 || rowsizeImplInt[row] != 0) {

  impliedRowUpper = impliedRowBounds.getSumUpperOrig(row);
  impliedRowLower = impliedRowBounds.getSumLowerOrig(row);

  // printf("implied bounds without tightenings: [%g,%g]\n", baseiRLower,
  //        baseiRUpper);

  auto checkForcingRow = [&](HighsInt row, HighsInt direction, double rowSide,
                             double impliedRowBound, double minAbsCoef,
                             HighsPostsolveStack::RowType rowType) {
    // 1. direction =  1 (>=): forcing row if upper bound on constraint activity
    //                         is equal to row's lower bound
    // 2. direction = -1 (<=): forcing row if lower bound on constraint activity
    //                         is equal to row's upper bound
    // scale tolerance (equivalent to scaling row to have minimum absolute
    // coefficient of 1)
    if (direction * impliedRowBound >
        direction * rowSide + primal_feastol * std::min(1.0, minAbsCoef))
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
        if (model->integrality_[nonzero.index()] != HighsVarType::kContinuous &&
            fractionality(model->col_upper_[nonzero.index()]) >
                mipsolver->options_mip_->mip_feasibility_tolerance) {
          // If a non-continuous variable is fixed at a fractional
          // value then the problem is infeasible
          return Result::kPrimalInfeasible;
        }
        // the upper bound of the column is as tight as the implied upper
        // bound or comes from this row, which means it is not used in the
        // rows implied bounds. Therefore we can fix the variable at its
        // upper bound.
        postsolve_stack.fixedColAtUpper(nonzero.index(),
                                        model->col_upper_[nonzero.index()],
                                        model->col_cost_[nonzero.index()],
                                        getColumnVector(nonzero.index()));
        if (model->col_lower_[nonzero.index()] <
            model->col_upper_[nonzero.index()])
          changeColLower(nonzero.index(), model->col_upper_[nonzero.index()]);
        removeFixedCol(nonzero.index());
      } else {
        if (model->integrality_[nonzero.index()] != HighsVarType::kContinuous &&
            fractionality(model->col_lower_[nonzero.index()]) >
                mipsolver->options_mip_->mip_feasibility_tolerance) {
          // If a non-continuous variable is fixed at a fractional
          // value then the problem is infeasible
          return Result::kPrimalInfeasible;
        }
        postsolve_stack.fixedColAtLower(nonzero.index(),
                                        model->col_lower_[nonzero.index()],
                                        model->col_cost_[nonzero.index()],
                                        getColumnVector(nonzero.index()));
        if (model->col_upper_[nonzero.index()] >
            model->col_lower_[nonzero.index()])
          changeColUpper(nonzero.index(), model->col_lower_[nonzero.index()]);
        removeFixedCol(nonzero.index());
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

  if (analysis_.allow_rule_[kPresolveRuleForcingRow]) {
    // Allow rule to consider forcing rows

    // store row and compute minimum absolute coefficient
    storeRow(row);
    double minAbsCoef = kHighsInf;
    for (const HighsSliceNonzero& nonzero : getStoredRow()) {
      minAbsCoef = std::min(minAbsCoef, std::abs(nonzero.value()));
    }

    // >= inequality
    HPRESOLVE_CHECKED_CALL(checkForcingRow(
        row, HighsInt{1}, model->row_lower_[row], impliedRowUpper, minAbsCoef,
        HighsPostsolveStack::RowType::kGeq));
    if (rowDeleted[row]) return Result::kOk;

    // <= inequality
    HPRESOLVE_CHECKED_CALL(checkForcingRow(
        row, HighsInt{-1}, model->row_upper_[row], impliedRowLower, minAbsCoef,
        HighsPostsolveStack::RowType::kLeq));
    if (rowDeleted[row]) return Result::kOk;
  }

  // implied bounds can only be computed when row bounds are available and
  // bounds on activity contain at most one infinite bound
  if (checkUpdateColImpliedBounds(row)) {
    for (const HighsSliceNonzero& nonzero : getRowVector(row))
      updateColImpliedBounds(row, nonzero.index(), nonzero.value());
  }

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
    if (fixColToLowerOrUnbounded(postsolve_stack, col)) {
      // Handle unboundedness
      presolve_status_ = HighsPresolveStatus::kUnboundedOrInfeasible;
      return Result::kDualInfeasible;
    }
  } else if (model->col_cost_[col] < 0 ||
             std::abs(model->col_upper_[col]) <
                 std::abs(model->col_lower_[col])) {
    if (fixColToUpperOrUnbounded(postsolve_stack, col)) {
      // Handle unboundedness
      presolve_status_ = HighsPresolveStatus::kUnboundedOrInfeasible;
      return Result::kDualInfeasible;
    }
  } else if (model->col_lower_[col] != -kHighsInf) {
    if (fixColToLowerOrUnbounded(postsolve_stack, col)) {
      // Handle unboundedness
      presolve_status_ = HighsPresolveStatus::kUnboundedOrInfeasible;
      return Result::kDualInfeasible;
    }
  } else {
    fixColToZero(postsolve_stack, col);
  }

  analysis_.logging_on_ = logging_on;
  if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleEmptyCol);
  return checkLimits(postsolve_stack);
}

HPresolve::Result HPresolve::colPresolve(HighsPostsolveStack& postsolve_stack,
                                         HighsInt col) {
  assert(!colDeleted[col]);
  const bool logging_on = analysis_.logging_on_;
  double boundDiff = model->col_upper_[col] - model->col_lower_[col];
  if (boundDiff <= primal_feastol) {
    if (boundDiff <= options->small_matrix_value ||
        getMaxAbsColVal(col) * boundDiff <= primal_feastol) {
      if (boundDiff < -primal_feastol) return Result::kPrimalInfeasible;
      postsolve_stack.removedFixedCol(col, model->col_lower_[col],
                                      model->col_cost_[col],
                                      getColumnVector(col));
      removeFixedCol(col);
      return checkLimits(postsolve_stack);
    }
  }

  switch (colsize[col]) {
    case 0:
      return emptyCol(postsolve_stack, col);
    case 1:
      return singletonCol(postsolve_stack, col);
    default:
      break;
  }

  // detect strong / weak domination
  HPRESOLVE_CHECKED_CALL(detectDominatedCol(postsolve_stack, col));
  if (colDeleted[col]) return Result::kOk;

  // column is not (weakly) dominated

  // the associated dual constraint has an upper bound if there is an infinite
  // or redundant column lower bound as then the reduced cost of the column must
  // not be positive i.e. <= 0

  // integer columns cannot be used to tighten bounds on dual multipliers
  if (mipsolver != nullptr) {
    auto modifyImpliedDualRowBound = [&](HighsInt col, HighsInt row,
                                         HighsInt direction,
                                         bool isBoundImplied, HighsInt numInf) {
      if (isBoundImplied && row != -1 && numInf == 1 &&
          direction * model->col_cost_[col] >= 0) {
        if (model->row_lower_[row] == -kHighsInf ||
            model->row_upper_[row] == kHighsInf) {
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
      }
    };

    modifyImpliedDualRowBound(col, colLowerSource[col], HighsInt{1},
                              isLowerImplied(col),
                              impliedDualRowBounds.getNumInfSumUpperOrig(col));

    modifyImpliedDualRowBound(col, colUpperSource[col], HighsInt{-1},
                              isUpperImplied(col),
                              impliedDualRowBounds.getNumInfSumLowerOrig(col));

    HPRESOLVE_CHECKED_CALL(static_cast<Result>(convertImpliedInteger(col)));

    // shift integral variables to have a lower bound of zero
    if (model->integrality_[col] != HighsVarType::kContinuous &&
        model->col_lower_[col] != 0.0 &&
        (model->col_lower_[col] != -kHighsInf ||
         model->col_upper_[col] != kHighsInf) &&
        model->col_upper_[col] - model->col_lower_[col] > 0.5) {
      // substitute with the bound that is smaller in magnitude and only
      // substitute if bound is not large for an integer
      if (std::abs(model->col_upper_[col]) > std::abs(model->col_lower_[col])) {
        if (std::abs(model->col_lower_[col]) < 1000.5)
          transformColumn(postsolve_stack, col, 1.0, model->col_lower_[col]);
      } else {
        if (std::abs(model->col_upper_[col]) < 1000.5)
          transformColumn(postsolve_stack, col, -1.0, model->col_upper_[col]);
      }
    }

    if (model->integrality_[col] == HighsVarType::kInteger) return Result::kOk;
  }

  // now check if we can expect to tighten at least one bound
  if (checkUpdateRowDualImpliedBounds(col)) {
    for (const HighsSliceNonzero& nonzero : getColumnVector(col))
      updateRowDualImpliedBounds(nonzero.index(), col, nonzero.value());
  }

  return Result::kOk;
}

HPresolve::Result HPresolve::detectDominatedCol(
    HighsPostsolveStack& postsolve_stack, HighsInt col,
    bool handleSingletonRows) {
  assert(!colDeleted[col]);

  // get bounds on column dual
  double colDualUpper =
      -impliedDualRowBounds.getSumLower(col, -model->col_cost_[col]);
  double colDualLower =
      -impliedDualRowBounds.getSumUpper(col, -model->col_cost_[col]);

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
    bool unbounded = false;
    if (direction > 0)
      unbounded = fixColToLowerOrUnbounded(postsolve_stack, col);
    else
      unbounded = fixColToUpperOrUnbounded(postsolve_stack, col);
    if (unbounded) {
      // Handle unboundedness
      presolve_status_ = HighsPresolveStatus::kUnboundedOrInfeasible;
      return Result::kDualInfeasible;
    }
    analysis_.logging_on_ = logging_on;
    if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleDominatedCol);
    // handle row singletons (if requested)
    if (handleSingletonRows)
      HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolve_stack));
    return checkLimits(postsolve_stack);
  };

  auto weaklyDominatedCol = [&](HighsInt col, double dualBound, double bound,
                                double otherBound, HighsInt direction) {
    // column is weakly dominated if the bounds on the column dual satisfy:
    // 1. lower bound >= -dual feasibility tolerance (direction =  1) or
    // 2. upper bound <=  dual feasibility tolerance (direction = -1).
    if (direction * dualBound < -options->dual_feasibility_tolerance)
      return Result::kOk;
    if (direction * bound != -kHighsInf) {
      if (logging_on) analysis_.startPresolveRuleLog(kPresolveRuleDominatedCol);
      // fix variable
      bool unbounded = false;
      if (direction > 0)
        unbounded = fixColToLowerOrUnbounded(postsolve_stack, col);
      else
        unbounded = fixColToUpperOrUnbounded(postsolve_stack, col);
      if (unbounded) {
        // Handle unboundedness
        presolve_status_ = HighsPresolveStatus::kUnboundedOrInfeasible;
        return Result::kDualInfeasible;
      }
      analysis_.logging_on_ = logging_on;
      if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleDominatedCol);
      // handle row singletons (if requested)
      if (handleSingletonRows)
        HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolve_stack));
      return checkLimits(postsolve_stack);
    } else if (analysis_.allow_rule_[kPresolveRuleForcingCol]) {
      // get bound on column dual using original bounds on row duals
      double boundOnColDual = direction > 0
                                  ? -impliedDualRowBounds.getSumUpperOrig(
                                        col, -model->col_cost_[col])
                                  : -impliedDualRowBounds.getSumLowerOrig(
                                        col, -model->col_cost_[col]);
      if (boundOnColDual == 0.0) {
        // 1. column's lower bound is infinite (i.e. column dual has upper bound
        // of zero) and column dual's lower bound is zero as well
        // (direction = 1) or
        // 2. column's upper bound is infinite (i.e. column dual has lower bound
        // of zero) and column dual's upper bound is zero as well
        // (direction = -1).
        // thus, the column dual is zero, and we can remove the column and
        // all its rows
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

  // check for weakly dominated column
  HPRESOLVE_CHECKED_CALL(
      weaklyDominatedCol(col, colDualLower, model->col_lower_[col],
                         model->col_upper_[col], HighsInt{1}));
  if (colDeleted[col]) return Result::kOk;

  HPRESOLVE_CHECKED_CALL(
      weaklyDominatedCol(col, colDualUpper, model->col_upper_[col],
                         model->col_lower_[col], HighsInt{-1}));
  return Result::kOk;
}

HPresolve::Result HPresolve::initialRowAndColPresolve(
    HighsPostsolveStack& postsolve_stack) {
  // do a full scan over the rows as the singleton arrays and the changed row
  // arrays are not initialized, also unset changedRowFlag so that the row will
  // be added to the changed row vector when it is changed after it was
  // processed
  for (HighsInt row = 0; row != model->num_row_; ++row) {
    if (rowDeleted[row]) continue;
    HPRESOLVE_CHECKED_CALL(rowPresolve(postsolve_stack, row));
    changedRowFlag[row] = false;
  }

  // same for the columns
  for (HighsInt col = 0; col != model->num_col_; ++col) {
    if (colDeleted[col]) continue;
    if (model->integrality_[col] != HighsVarType::kContinuous) {
      // round and update bounds
      changeColLower(col, model->col_lower_[col]);
      changeColUpper(col, model->col_upper_[col]);
    }
    HPRESOLVE_CHECKED_CALL(colPresolve(postsolve_stack, col));
    changedColFlag[col] = false;
  }

  return checkLimits(postsolve_stack);
}

HPresolve::Result HPresolve::fastPresolveLoop(
    HighsPostsolveStack& postsolve_stack) {
  do {
    storeCurrentProblemSize();

    HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolve_stack));

    HPRESOLVE_CHECKED_CALL(presolveChangedRows(postsolve_stack));

    HPRESOLVE_CHECKED_CALL(removeDoubletonEquations(postsolve_stack));

    HPRESOLVE_CHECKED_CALL(presolveColSingletons(postsolve_stack));

    HPRESOLVE_CHECKED_CALL(presolveChangedCols(postsolve_stack));

  } while (problemSizeReduction() > 0.01);

  return Result::kOk;
}

HPresolve::Result HPresolve::presolve(HighsPostsolveStack& postsolve_stack) {
  // for the inner most loop we take the order roughly from the old presolve
  // but we nest the rounds with a new outer loop which layers the newer
  // presolvers
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

  // convert model to minimization problem
  if (model->sense_ == ObjSense::kMaximize) {
    for (HighsInt i = 0; i != model->num_col_; ++i)
      model->col_cost_[i] = -model->col_cost_[i];

    model->offset_ = -model->offset_;
    assert(std::isfinite(model->offset_));
    model->sense_ = ObjSense::kMinimize;
  }

  // Set up the logic to allow presolve rules, and logging for their
  // effectiveness
  analysis_.setup(this->model, this->options, this->numDeletedRows,
                  this->numDeletedCols);

  if (options->presolve != kHighsOffString) {
    if (mipsolver) mipsolver->mipdata_->cliquetable.setPresolveFlag(true);
    if (!mipsolver || mipsolver->mipdata_->numRestarts == 0)
      highsLogUser(options->log_options, HighsLogType::kInfo,
                   "Presolving model\n");

    auto report = [&]() {
      if (!mipsolver || mipsolver->mipdata_->numRestarts == 0) {
        HighsInt numCol = model->num_col_ - numDeletedCols;
        HighsInt numRow = model->num_row_ - numDeletedRows;
        HighsInt numNonz = Avalue.size() - freeslots.size();
        // Only read the run time if it's to be printed
        const double run_time = options->output_flag ? this->timer->read() : 0;
#ifndef NDEBUG
        std::string time_str = " " + std::to_string(run_time) + "s";
#else
        std::string time_str =
            " " + std::to_string(static_cast<int>(run_time)) + "s";
#endif
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

    HPRESOLVE_CHECKED_CALL(initialRowAndColPresolve(postsolve_stack));

    HighsInt numParallelRowColCalls = 0;
#if ENABLE_SPARSIFY_FOR_LP
    bool trySparsify = true;  // mipsolver != nullptr;
#else
    bool trySparsify =
        mipsolver != nullptr || !options->lp_presolve_requires_basis_postsolve;
#endif
    bool tryProbing =
        mipsolver != nullptr && analysis_.allow_rule_[kPresolveRuleProbing];
    HighsInt numCliquesBeforeProbing = -1;
    bool domcolAfterProbingCalled = false;
    bool dependentEquationsCalled = mipsolver != nullptr;
    HighsInt lastPrintSize = kHighsIInf;

    // Start of main presolve loop
    //
    while (true) {
      HighsInt currSize =
          model->num_col_ - numDeletedCols + model->num_row_ - numDeletedRows;
      if (currSize < 0.85 * lastPrintSize) {
        lastPrintSize = currSize;
        report();
      }

      HPRESOLVE_CHECKED_CALL(fastPresolveLoop(postsolve_stack));

      storeCurrentProblemSize();

      // when presolving after a restart the clique table and implication
      // structure may contain substitutions which we apply directly before
      // running the aggregator as they might lose validity otherwise
      if (mipsolver != nullptr) {
        HPRESOLVE_CHECKED_CALL(
            applyConflictGraphSubstitutions(postsolve_stack));
      }

      if (analysis_.allow_rule_[kPresolveRuleAggregator])
        HPRESOLVE_CHECKED_CALL(aggregator(postsolve_stack));

      if (problemSizeReduction() > 0.05) continue;

      if (trySparsify) {
        HighsInt numNz = numNonzeros();
        HPRESOLVE_CHECKED_CALL(sparsify(postsolve_stack));
        double nzReduction =
            100.0 * (1.0 - (numNonzeros() / static_cast<double>(numNz)));

        if (nzReduction > 0) {
          highsLogDev(options->log_options, HighsLogType::kInfo,
                      "Sparsify removed %.1f%% of nonzeros\n", nzReduction);

          // #1710 exposes that this should not be
          //
          // fastPresolveLoop(postsolve_stack);
          //
          // but
          HPRESOLVE_CHECKED_CALL(fastPresolveLoop(postsolve_stack));
        }
        trySparsify = false;
      }

      if (analysis_.allow_rule_[kPresolveRuleParallelRowsAndCols] &&
          numParallelRowColCalls < 5) {
        if (shrinkProblemEnabled && (numDeletedCols >= model->num_col_ / 2 ||
                                     numDeletedRows >= model->num_row_ / 2)) {
          shrinkProblem(postsolve_stack);

          toCSC(model->a_matrix_.value_, model->a_matrix_.index_,
                model->a_matrix_.start_);
          okFromCSC(model->a_matrix_.value_, model->a_matrix_.index_,
                    model->a_matrix_.start_);
        }
        storeCurrentProblemSize();
        HPRESOLVE_CHECKED_CALL(detectParallelRowsAndCols(postsolve_stack));
        ++numParallelRowColCalls;
        if (problemSizeReduction() > 0.05) continue;
      }

      HPRESOLVE_CHECKED_CALL(fastPresolveLoop(postsolve_stack));

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

      HPRESOLVE_CHECKED_CALL(fastPresolveLoop(postsolve_stack));

      if (mipsolver != nullptr && numCliquesBeforeProbing == -1) {
        numCliquesBeforeProbing = mipsolver->mipdata_->cliquetable.numCliques();
        storeCurrentProblemSize();
        HPRESOLVE_CHECKED_CALL(dominatedColumns(postsolve_stack));
        if (problemSizeReduction() > 0.0)
          HPRESOLVE_CHECKED_CALL(fastPresolveLoop(postsolve_stack));
        if (problemSizeReduction() > 0.05) continue;
      }

      if (tryProbing) {
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
          shrinkProblem(postsolve_stack);

          toCSC(model->a_matrix_.value_, model->a_matrix_.index_,
                model->a_matrix_.start_);
          okFromCSC(model->a_matrix_.value_, model->a_matrix_.index_,
                    model->a_matrix_.start_);
        }
        storeCurrentProblemSize();
        if (analysis_.allow_rule_[kPresolveRuleDependentEquations]) {
          HPRESOLVE_CHECKED_CALL(removeDependentEquations(postsolve_stack));
          dependentEquationsCalled = true;
        }
        if (analysis_.allow_rule_[kPresolveRuleDependentFreeCols])
          HPRESOLVE_CHECKED_CALL(removeDependentFreeCols(postsolve_stack));
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

    // Now consider removing slacks
    if (options->presolve_remove_slacks)
      HPRESOLVE_CHECKED_CALL(removeSlacks(postsolve_stack));

    report();
  } else {
    highsLogUser(options->log_options, HighsLogType::kInfo,
                 "\nPresolve is switched off\n");
  }

  if (mipsolver != nullptr) scaleMIP(postsolve_stack);

  // analysePresolveRuleLog() should return true - no errors
  assert(analysis_.analysePresolveRuleLog());
  // Possibly report presolve log
  analysis_.analysePresolveRuleLog(true);
  return Result::kOk;
}

HPresolve::Result HPresolve::removeSlacks(
    HighsPostsolveStack& postsolve_stack) {
  // SingletonColumns data structure appears not to be retained
  // throughout presolve
  for (HighsInt iCol = 0; iCol != model->num_col_; ++iCol) {
    if (colDeleted[iCol]) continue;
    if (colsize[iCol] != 1) continue;
    if (model->integrality_[iCol] == HighsVarType::kInteger) continue;
    HighsInt coliter = colhead[iCol];
    HighsInt iRow = Arow[coliter];
    assert(Acol[coliter] == iCol);
    assert(!rowDeleted[iRow]);
    if (!isEquation(iRow)) continue;
    double lower = model->col_lower_[iCol];
    double upper = model->col_upper_[iCol];
    double cost = model->col_cost_[iCol];
    double rhs = model->row_lower_[iRow];
    double coeff = Avalue[coliter];
    assert(coeff);
    // Slack is s = (rhs - a^Tx)/coeff
    //
    // Constraint bounds become:
    //
    // For coeff > 0 [rhs - coeff * upper, rhs - coeff * lower]
    //
    // For coeff < 0 [rhs - coeff * lower, rhs - coeff * upper]
    model->row_lower_[iRow] =
        coeff > 0 ? rhs - coeff * upper : rhs - coeff * lower;
    model->row_upper_[iRow] =
        coeff > 0 ? rhs - coeff * lower : rhs - coeff * upper;
    if (cost) {
      // Cost is (cost * rhs / coeff) + (col_cost - (cost/coeff) row_values)^Tx
      double multiplier = cost / coeff;
      for (const HighsSliceNonzero& nonzero : getRowVector(iRow)) {
        HighsInt local_iCol = nonzero.index();
        double local_value = nonzero.value();
        model->col_cost_[local_iCol] -= multiplier * local_value;
      }
      model->offset_ += multiplier * rhs;
    }
    //
    postsolve_stack.slackColSubstitution(iRow, iCol, rhs, getRowVector(iRow));

    markColDeleted(iCol);

    unlink(coliter);
  }
  return Result::kOk;
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

  return numreductions >= reductionLimit ? Result::kStopped : Result::kOk;
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
  assert(model->a_matrix_.numNz() || model->num_row_ == 0);
  auto reportReductions = [&]() {
    if (options->presolve != kHighsOffString &&
        reductionLimit < kHighsSize_tInf) {
      highsLogUser(options->log_options, HighsLogType::kInfo,
                   "Presolve performed %" PRId64 " of %" PRId64
                   " permitted reductions\n",
                   postsolve_stack.numReductions(), reductionLimit);
    }
  };
  switch (presolve(postsolve_stack)) {
    case Result::kStopped:
    case Result::kOk:
      break;
    case Result::kPrimalInfeasible:
      presolve_status_ = HighsPresolveStatus::kInfeasible;
      reportReductions();
      return HighsModelStatus::kInfeasible;
    case Result::kDualInfeasible:
      presolve_status_ = HighsPresolveStatus::kUnboundedOrInfeasible;
      reportReductions();
      return HighsModelStatus::kUnboundedOrInfeasible;
  }
  reportReductions();

  shrinkProblem(postsolve_stack);

  if (mipsolver != nullptr) {
    mipsolver->mipdata_->cliquetable.setPresolveFlag(false);
    mipsolver->mipdata_->cliquetable.setMaxEntries(numNonzeros());
    mipsolver->mipdata_->domain.addCutpool(mipsolver->mipdata_->cutpool);
    mipsolver->mipdata_->domain.addConflictPool(
        mipsolver->mipdata_->conflictPool);

    if (mipsolver->mipdata_->numRestarts != 0) {
      std::vector<HighsInt> cutinds;
      std::vector<double> cutvals;
      cutinds.reserve(model->num_col_);
      cutvals.reserve(model->num_col_);
      HighsInt numcuts = 0;
      for (HighsInt i = model->num_row_ - 1; i >= 0; --i) {
        // check if we already reached the original rows
        if (postsolve_stack.getOrigRowIndex(i) <
            mipsolver->orig_model_->num_row_)
          break;

        // row is a cut, remove it from matrix but add to cutpool
        ++numcuts;
        storeRow(i);
        cutinds.clear();
        cutvals.clear();
        for (HighsInt j : rowpositions) {
          cutinds.push_back(Acol[j]);
          cutvals.push_back(Avalue[j]);
        }

        mipsolver->mipdata_->cutpool.addCut(
            *mipsolver, cutinds.data(), cutvals.data(), cutinds.size(),
            model->row_upper_[i],
            rowsizeInteger[i] + rowsizeImplInt[i] == rowsize[i] &&
                rowCoefficientsIntegral(i, 1.0),
            true, false, false);

        markRowDeleted(i);
        for (HighsInt j : rowpositions) unlink(j);
      }

      model->num_row_ -= numcuts;
      model->row_lower_.resize(model->num_row_);
      model->row_upper_.resize(model->num_row_);
      model->row_names_.resize(model->num_row_);
    }
  }

  toCSC(model->a_matrix_.value_, model->a_matrix_.index_,
        model->a_matrix_.start_);

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
      // limit has been reached
      assert(model->num_row_ == 0 ||
             postsolve_stack.numReductions() >= reductionLimit);
      if (model->num_row_ != 0) {
        presolve_status_ = HighsPresolveStatus::kNotPresolved;
        return HighsModelStatus::kNotset;
      }
    }
    presolve_status_ = HighsPresolveStatus::kReducedToEmpty;
    // Make sure that zero row activity from the column-less model is
    // consistent with the bounds
    return model->num_row_ == 0 || zeroRowActivityFeasible()
               ? HighsModelStatus::kOptimal
               : HighsModelStatus::kInfeasible;
  } else if (postsolve_stack.numReductions() > 0) {
    // Reductions performed
    presolve_status_ = HighsPresolveStatus::kReduced;
  } else {
    // No reductions performed
    presolve_status_ = HighsPresolveStatus::kNotReduced;
  }

  if (!mipsolver && options->use_implied_bounds_from_presolve)
    setRelaxedImpliedBounds();

  assert(presolve_status_ != HighsPresolveStatus::kNotSet);
  return HighsModelStatus::kNotset;
}

void HPresolve::computeIntermediateMatrix(std::vector<HighsInt>& flagRow,
                                          std::vector<HighsInt>& flagCol,
                                          size_t& numreductions) {
  shrinkProblemEnabled = false;
  HighsPostsolveStack stack;
  stack.initializeIndexMaps(flagRow.size(), flagCol.size());
  setReductionLimit(numreductions);
  presolve(stack);
  numreductions = stack.numReductions();

  toCSC(model->a_matrix_.value_, model->a_matrix_.index_,
        model->a_matrix_.start_);

  for (HighsInt i = 0; i != model->num_row_; ++i)
    flagRow[i] = 1 - rowDeleted[i];
  for (HighsInt i = 0; i != model->num_col_; ++i)
    flagCol[i] = 1 - colDeleted[i];
}

HPresolve::Result HPresolve::removeDependentEquations(
    HighsPostsolveStack& postsolve_stack) {
  assert(analysis_.allow_rule_[kPresolveRuleDependentEquations]);
  const bool logging_on = analysis_.logging_on_;
  if (equations.empty()) return Result::kOk;

  if (logging_on)
    analysis_.startPresolveRuleLog(kPresolveRuleDependentEquations);

  HighsSparseMatrix matrix;
  matrix.num_col_ = equations.size();
  matrix.num_row_ = model->num_col_ + 1;
  matrix.start_.resize(matrix.num_col_ + 1);
  matrix.start_[0] = 0;
  const HighsInt maxCapacity = numNonzeros() + matrix.num_col_;
  matrix.value_.reserve(maxCapacity);
  matrix.index_.reserve(maxCapacity);

  std::vector<HighsInt> eqSet(matrix.num_col_);
  HighsInt i = 0;
  for (const std::pair<HighsInt, HighsInt>& p : equations) {
    HighsInt eq = p.second;
    eqSet[i++] = eq;

    // add entries of equation
    for (const HighsSliceNonzero& nonz : getRowVector(eq)) {
      matrix.value_.push_back(nonz.value());
      matrix.index_.push_back(nonz.index());
    }

    // add entry for artificial rhs column
    if (model->row_lower_[eq] != 0.0) {
      matrix.value_.push_back(model->row_lower_[eq]);
      matrix.index_.push_back(model->num_col_);
    }

    matrix.start_[i] = matrix.value_.size();
  }
  std::vector<HighsInt> colSet(matrix.num_col_);
  std::iota(colSet.begin(), colSet.end(), 0);
  HFactor factor;
  factor.setup(matrix, colSet);
  // Set up a time limit to prevent the redundant rows factorization
  // taking forever.
  //
  // Allow no more than 1% of the time limit to be spent on removing
  // dependent equations, but ensure that there is some limit since
  // options->time_limit is infinity by default
  //
  // ToDo: This is strictly non-deterministic, but so conservative
  // that it'll only reap the cases when factor.build never finishes
  const double time_limit =
      std::max(1.0, std::min(0.01 * options->time_limit, 1000.0));
  factor.setTimeLimit(time_limit);
  // Determine rank deficiency of the equations
  highsLogUser(options->log_options, HighsLogType::kInfo,
               "Dependent equations search running on %d equations with time "
               "limit of %.2fs\n",
               static_cast<int>(matrix.num_col_), time_limit);
  double time_taken = -this->timer->read();
  HighsInt build_return = factor.build();
  time_taken += this->timer->read();
  if (build_return == kBuildKernelReturnTimeout) {
    // HFactor::build has timed out, so just return
    highsLogUser(options->log_options, HighsLogType::kInfo,
                 "Dependent equations search terminated after %.3gs due to "
                 "expected time exceeding limit\n",
                 time_taken);
    analysis_.logging_on_ = logging_on;
    if (logging_on)
      analysis_.stopPresolveRuleLog(kPresolveRuleDependentFreeCols);
    return Result::kOk;
  } else {
    double pct_off_timeout =
        1e2 * std::fabs(time_taken - time_limit) / time_limit;
    if (pct_off_timeout < 1.0)
      highsLogUser(options->log_options, HighsLogType::kWarning,
                   "Dependent equations search finished within %.2f%% of limit "
                   "of %.2fs: "
                   "risk of non-deterministic behaviour if solve is repeated\n",
                   pct_off_timeout, time_limit);
  }
  // build_return as rank_deficiency must be valid
  assert(build_return >= 0);
  const HighsInt rank_deficiency = build_return;
  // Analyse what's been removed
  HighsInt num_removed_row = 0;
  HighsInt num_removed_nz = 0;
  HighsInt num_fictitious_rows_skipped = 0;
  for (HighsInt k = 0; k < rank_deficiency; k++) {
    if (factor.var_with_no_pivot[k] >= 0) {
      HighsInt redundant_row = eqSet[factor.var_with_no_pivot[k]];
      num_removed_row++;
      num_removed_nz += rowsize[redundant_row];
      postsolve_stack.redundantRow(redundant_row);
      removeRow(redundant_row);
    } else {
      num_fictitious_rows_skipped++;
    }
  }
  highsLogUser(options->log_options, HighsLogType::kInfo,
               "Dependent equations search removed %d rows and %d nonzeros "
               "in %.2fs (limit = %.2fs)\n",
               static_cast<int>(num_removed_row),
               static_cast<int>(num_removed_nz), time_taken, time_limit);
  if (num_fictitious_rows_skipped)
    highsLogDev(options->log_options, HighsLogType::kInfo,
                ", avoiding %d fictitious rows",
                static_cast<int>(num_fictitious_rows_skipped));
  highsLogDev(options->log_options, HighsLogType::kInfo, "\n");

  analysis_.logging_on_ = logging_on;
  if (logging_on)
    analysis_.stopPresolveRuleLog(kPresolveRuleDependentEquations);
  return Result::kOk;
}

HPresolve::Result HPresolve::removeDependentFreeCols(
    HighsPostsolveStack& postsolve_stack) {
  return Result::kOk;

  // Commented out unreachable code
  //  assert(analysis_.allow_rule_[kPresolveRuleDependentFreeCols]);
  //  const bool logging_on = analysis_.logging_on_;
  //  if (logging_on)
  //    analysis_.startPresolveRuleLog(kPresolveRuleDependentFreeCols);
  //
  //  // todo the postsolve step does not work properly
  //  std::vector<HighsInt> freeCols;
  //  freeCols.reserve(model->num_col_);
  //
  //  for (HighsInt i = 0; i < model->num_col_; ++i) {
  //    if (colDeleted[i]) continue;
  //    if (model->col_lower_[i] == -kHighsInf && model->col_upper_[i] ==
  //    kHighsInf)
  //      freeCols.push_back(i);
  //  }
  //
  //  if (freeCols.empty()) return Result::kOk;
  //
  //  HighsSparseMatrix matrix;
  //  matrix.num_col_ = freeCols.size();
  //  highsLogDev(options->log_options, HighsLogType::kInfo,
  //              "HPresolve::removeDependentFreeCols Got %d free cols, checking
  //              " "for dependent free cols\n", (int)matrix.num_col_);
  //  matrix.num_row_ = model->num_row_ + 1;
  //  matrix.start_.resize(matrix.num_col_ + 1);
  //  matrix.start_[0] = 0;
  //  const HighsInt maxCapacity = numNonzeros() + matrix.num_col_;
  //  matrix.value_.reserve(maxCapacity);
  //  matrix.index_.reserve(maxCapacity);
  //
  //  for (HighsInt i = 0; i < matrix.num_col_; ++i) {
  //    HighsInt col = freeCols[i];
  //    // add entries of free column
  //    for (const HighsSliceNonzero& nonz : getColumnVector(col)) {
  //      matrix.value_.push_back(nonz.value());
  //      matrix.index_.push_back(nonz.index());
  //    }
  //
  //    // add entry for artificial cost row
  //    if (model->col_cost_[col] != 0.0) {
  //      matrix.value_.push_back(model->col_cost_[col]);
  //      matrix.index_.push_back(model->num_row_);
  //    }
  //
  //    matrix.start_[i + 1] = matrix.value_.size();
  //  }
  //  std::vector<HighsInt> colSet(matrix.num_col_);
  //  std::iota(colSet.begin(), colSet.end(), 0);
  //  HFactor factor;
  //  factor.setup(matrix, colSet);
  //  HighsInt rank_deficiency = factor.build();
  //  // Must not have timed out
  //  assert(rank_deficiency >= 0);
  //  highsLogDev(options->log_options, HighsLogType::kInfo,
  //              "HPresolve::removeDependentFreeCols Got %d free cols, checking
  //              " "for dependent free cols\n", (int)matrix.num_col_);
  //  // Analyse what's been removed
  //  HighsInt num_removed_row = 0;
  //  HighsInt num_removed_nz = 0;
  //  HighsInt num_fictitious_cols_skipped = 0;
  //  for (HighsInt k = 0; k < rank_deficiency; k++) {
  //    if (factor.var_with_no_pivot[k] >= 0) {
  //      HighsInt redundant_col = freeCols[factor.var_with_no_pivot[k]];
  //      num_removed_nz += colsize[redundant_col];
  //      fixColToZero(postsolve_stack, redundant_col);
  //    } else {
  //      num_fictitious_cols_skipped++;
  //    }
  //  }
  //  highsLogDev(
  //      options->log_options, HighsLogType::kInfo,
  //      "HPresolve::removeDependentFreeCols Removed %d rows and %d nonzeros",
  //      (int)num_removed_row, (int)num_removed_nz);
  //  if (num_fictitious_cols_skipped)
  //    highsLogDev(options->log_options, HighsLogType::kInfo,
  //                ", avoiding %d fictitious rows",
  //                (int)num_fictitious_cols_skipped);
  //  highsLogDev(options->log_options, HighsLogType::kInfo, "\n");
  //
  //  analysis_.logging_on_ = logging_on;
  //  if (logging_on)
  //  analysis_.stopPresolveRuleLog(kPresolveRuleDependentFreeCols);
  //
  //  return Result::kOk;
}

HPresolve::Result HPresolve::aggregator(HighsPostsolveStack& postsolve_stack) {
  assert(analysis_.allow_rule_[kPresolveRuleAggregator]);
  const bool logging_on = analysis_.logging_on_;
  if (logging_on) analysis_.startPresolveRuleLog(kPresolveRuleAggregator);
  substitutionOpportunities.erase(
      std::remove_if(substitutionOpportunities.begin(),
                     substitutionOpportunities.end(),
                     [&](const std::pair<HighsInt, HighsInt>& p) {
                       HighsInt row = p.first;
                       HighsInt col = p.second;
                       return rowDeleted[row] || colDeleted[col] ||
                              !isImpliedFree(col) || !isDualImpliedFree(row);
                     }),
      substitutionOpportunities.end());

  pdqsort(
      substitutionOpportunities.begin(), substitutionOpportunities.end(),
      [&](const std::pair<HighsInt, HighsInt>& nz1,
          const std::pair<HighsInt, HighsInt>& nz2) {
        HighsInt minLen1 = std::min(rowsize[nz1.first], colsize[nz1.second]);
        HighsInt minLen2 = std::min(rowsize[nz2.first], colsize[nz2.second]);
        if (minLen1 == 2 && minLen2 != 2) return true;
        if (minLen2 == 2 && minLen1 != 2) return false;

        int64_t sizeProd1 = int64_t(rowsize[nz1.first]) * colsize[nz1.second];
        int64_t sizeProd2 = int64_t(rowsize[nz2.first]) * colsize[nz2.second];
        if (sizeProd1 < sizeProd2) return true;
        if (sizeProd2 < sizeProd1) return false;
        if (minLen1 < minLen2) return true;
        if (minLen2 < minLen1) return false;

        return std::make_tuple(HighsHashHelpers::hash(std::make_pair(
                                   uint32_t(nz1.first), uint32_t(nz1.second))),
                               nz1.first, nz1.second) <
               std::make_tuple(HighsHashHelpers::hash(std::make_pair(
                                   uint32_t(nz2.first), uint32_t(nz2.second))),
                               nz2.first, nz2.second);
      });

  HighsInt nfail = 0;
  for (size_t i = 0; i < substitutionOpportunities.size(); ++i) {
    HighsInt row = substitutionOpportunities[i].first;
    HighsInt col = substitutionOpportunities[i].second;

    if (rowDeleted[row] || colDeleted[col] || !isImpliedFree(col) ||
        !isDualImpliedFree(row)) {
      substitutionOpportunities[i].first = -1;
      continue;
    }

    HighsInt nzPos = findNonzero(row, col);
    if (nzPos == -1) {
      substitutionOpportunities[i].first = -1;
      continue;
    }
    if (model->integrality_[col] == HighsVarType::kInteger) {
      StatusResult impliedIntegral = isImpliedIntegral(col);
      HPRESOLVE_CHECKED_CALL(static_cast<Result>(impliedIntegral));
      if (!impliedIntegral) continue;
    }

    // in the case where the row has length two or the column has length two
    // we always do the substitution since the fillin can never be problematic
    if (rowsize[row] == 2 || colsize[col] == 2) {
      storeRow(row);
      substituteFreeCol(postsolve_stack, row, col, true);
      substitutionOpportunities[i].first = -1;
      HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolve_stack));
      HPRESOLVE_CHECKED_CALL(checkLimits(postsolve_stack));
      continue;
    }

    double maxVal = rowsize[row] < colsize[col] ? getMaxAbsRowVal(row)
                                                : getMaxAbsColVal(col);
    if (std::fabs(Avalue[nzPos]) < maxVal * options->presolve_pivot_threshold) {
      maxVal = rowsize[row] < colsize[col] ? getMaxAbsColVal(col)
                                           : getMaxAbsRowVal(row);
      if (std::fabs(Avalue[nzPos]) <
          maxVal * options->presolve_pivot_threshold) {
        substitutionOpportunities[i].first = -1;
        continue;
      }
    }

    storeRow(row);
    HighsInt fillin = -(rowsize[row] + colsize[col] - 1);
    for (const auto& nz : getColumnVector(col)) {
      if (nz.index() == row) continue;
      fillin += countFillin(nz.index());

      if (fillin > options->presolve_substitution_maxfillin) break;
    }

    if (fillin > options->presolve_substitution_maxfillin) {
      ++nfail;
      // if the fill in is too much for multiple tries, then we stop
      // as this indicates that the rows/columns are becoming too dense
      // for substitutions
      if (nfail == 3) break;
      continue;
    }

    nfail = 0;
    substituteFreeCol(postsolve_stack, row, col, true);
    substitutionOpportunities[i].first = -1;
    HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolve_stack));
    HPRESOLVE_CHECKED_CALL(checkLimits(postsolve_stack));
  }

  substitutionOpportunities.erase(
      std::remove_if(
          substitutionOpportunities.begin(), substitutionOpportunities.end(),
          [](const std::pair<HighsInt, HighsInt>& p) { return p.first == -1; }),
      substitutionOpportunities.end());

  analysis_.logging_on_ = logging_on;
  if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleAggregator);
  return Result::kOk;
}

void HPresolve::substitute(HighsInt substcol, HighsInt staycol, double offset,
                           double scale) {
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
    if (model->row_lower_[colrow] != -kHighsInf)
      model->row_lower_[colrow] -= colval * offset;

    if (model->row_upper_[colrow] != kHighsInf)
      model->row_upper_[colrow] -= colval * offset;

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
}

bool HPresolve::fixColToLowerOrUnbounded(HighsPostsolveStack& postsolve_stack,
                                         HighsInt col) {
  double fixval = model->col_lower_[col];
  if (fixval == -kHighsInf) return true;

  const bool logging_on = analysis_.logging_on_;
  if (logging_on) analysis_.startPresolveRuleLog(kPresolveRuleFixedCol);

  // printf("fixing column %" HIGHSINT_FORMAT " to %.15g\n", col, fixval);

  // mark the column as deleted first so that it is not registered as singleton
  // column upon removing its nonzeros
  postsolve_stack.fixedColAtLower(col, fixval, model->col_cost_[col],
                                  getColumnVector(col));
  markColDeleted(col);

  for (HighsInt coliter = colhead[col]; coliter != -1;) {
    HighsInt colrow = Arow[coliter];
    double colval = Avalue[coliter];
    assert(Acol[coliter] == col);

    HighsInt colpos = coliter;
    coliter = Anext[coliter];

    if (model->row_lower_[colrow] != -kHighsInf)
      model->row_lower_[colrow] -= colval * fixval;

    if (model->row_upper_[colrow] != kHighsInf)
      model->row_upper_[colrow] -= colval * fixval;

    unlink(colpos);

    reinsertEquation(colrow);
  }

  model->offset_ += model->col_cost_[col] * fixval;
  assert(std::isfinite(model->offset_));
  model->col_cost_[col] = 0;
  analysis_.logging_on_ = logging_on;
  if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleFixedCol);
  return false;
}

bool HPresolve::fixColToUpperOrUnbounded(HighsPostsolveStack& postsolve_stack,
                                         HighsInt col) {
  double fixval = model->col_upper_[col];
  if (fixval == kHighsInf) return true;

  const bool logging_on = analysis_.logging_on_;
  if (logging_on) analysis_.startPresolveRuleLog(kPresolveRuleFixedCol);

  // printf("fixing column %" HIGHSINT_FORMAT " to %.15g\n", col, fixval);

  // mark the column as deleted first so that it is not registered as singleton
  // column upon removing its nonzeros
  postsolve_stack.fixedColAtUpper(col, fixval, model->col_cost_[col],
                                  getColumnVector(col));
  markColDeleted(col);

  for (HighsInt coliter = colhead[col]; coliter != -1;) {
    HighsInt colrow = Arow[coliter];
    double colval = Avalue[coliter];
    assert(Acol[coliter] == col);

    HighsInt colpos = coliter;
    coliter = Anext[coliter];

    if (model->row_lower_[colrow] != -kHighsInf)
      model->row_lower_[colrow] -= colval * fixval;

    if (model->row_upper_[colrow] != kHighsInf)
      model->row_upper_[colrow] -= colval * fixval;

    unlink(colpos);

    reinsertEquation(colrow);
  }

  model->offset_ += model->col_cost_[col] * fixval;
  assert(std::isfinite(model->offset_));
  model->col_cost_[col] = 0;
  analysis_.logging_on_ = logging_on;
  if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleFixedCol);
  return false;
}

void HPresolve::fixColToZero(HighsPostsolveStack& postsolve_stack,
                             HighsInt col) {
  const bool logging_on = analysis_.logging_on_;
  if (logging_on) analysis_.startPresolveRuleLog(kPresolveRuleFixedCol);
  postsolve_stack.fixedColAtZero(col, model->col_cost_[col],
                                 getColumnVector(col));
  // mark the column as deleted first so that it is not registered as singleton
  // column upon removing its nonzeros
  markColDeleted(col);

  for (HighsInt coliter = colhead[col]; coliter != -1;) {
    HighsInt colrow = Arow[coliter];
    assert(Acol[coliter] == col);

    HighsInt colpos = coliter;
    coliter = Anext[coliter];

    unlink(colpos);

    reinsertEquation(colrow);
  }

  model->col_cost_[col] = 0;
  analysis_.logging_on_ = logging_on;
  if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleFixedCol);
}

void HPresolve::removeRow(HighsInt row) {
  assert(row >= 0);
  assert(static_cast<size_t>(row) < rowroot.size());
  // first mark the row as logically deleted, so that it is not register as
  // singleton row upon removing its nonzeros
  markRowDeleted(row);
  storeRow(row);
  for (HighsInt rowiter : rowpositions) {
    assert(Arow[rowiter] == row);
    unlink(rowiter);
  }
}

void HPresolve::removeFixedCol(HighsInt col) {
  const bool logging_on = analysis_.logging_on_;
  if (logging_on) analysis_.startPresolveRuleLog(kPresolveRuleFixedCol);
  double fixval = model->col_lower_[col];

  markColDeleted(col);

  for (HighsInt coliter = colhead[col]; coliter != -1;) {
    HighsInt colrow = Arow[coliter];
    double colval = Avalue[coliter];
    assert(Acol[coliter] == col);

    HighsInt colpos = coliter;
    coliter = Anext[coliter];

    if (model->row_lower_[colrow] != -kHighsInf)
      model->row_lower_[colrow] -= colval * fixval;

    if (model->row_upper_[colrow] != kHighsInf)
      model->row_upper_[colrow] -= colval * fixval;

    unlink(colpos);

    reinsertEquation(colrow);
  }

  model->offset_ += model->col_cost_[col] * fixval;
  assert(std::isfinite(model->offset_));
  model->col_cost_[col] = 0;
  analysis_.logging_on_ = logging_on;
  if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleFixedCol);
}

HPresolve::Result HPresolve::removeRowSingletons(
    HighsPostsolveStack& postsolve_stack) {
  for (size_t i = 0; i != singletonRows.size(); ++i) {
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
  for (size_t i = 0; i != singletonColumns.size(); ++i) {
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
  std::vector<double> reducedcost;
  std::vector<double> upper;
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
    if (model->row_lower_[row] != -kHighsInf &&
        model->row_upper_[row] != kHighsInf)
      continue;

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
      double ub = model->col_upper_[col] - model->col_lower_[col];
      if (weight > 0) {
        comp = 1;
        maxviolation += model->col_upper_[col] * weight;
      } else {
        comp = -1;
        maxviolation += model->col_lower_[col] * weight;
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
        double delta = upper[index] * reducedcost[index];

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

      double al = reducedcost[alpos];
      coefs.resize(cover.size());
      double coverrhs = std::max(
          std::ceil(static_cast<double>(lambda / al - primal_feastol)), 1.0);
      HighsCDouble slackupper = -coverrhs;

      double step = kHighsInf;
      for (size_t i = 0; i != cover.size(); ++i) {
        coefs[i] = std::ceil(
            std::min(reducedcost[cover[i]], static_cast<double>(lambda)) / al -
            options->small_matrix_value);
        slackupper += upper[cover[i]] * coefs[i];
        step = std::min(step, reducedcost[cover[i]] / coefs[i]);
      }
      step = std::min(step, static_cast<double>(maxviolation / coverrhs));
      maxviolation -= step * coverrhs;

      HighsInt slackind = reducedcost.size();
      reducedcost.push_back(step);
      upper.push_back(static_cast<double>(slackupper));

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
                                          std::abs(reducedcost[i]) <= threshold;
                                 }),
                  indices.end());
    if (indices.empty()) continue;

    auto updateNonZeros = [&](HighsInt row, HighsCDouble& rhs,
                              HighsInt direction) {
      for (HighsInt i : indices) {
        assert(Arow[positions[i]] == row);
        double coefdelta =
            direction * static_cast<double>(reducedcost[i] - maxviolation);
        HighsInt col = Acol[positions[i]];

        if (complementation[i] == -1) {
          rhs += coefdelta * model->col_lower_[col];
          addToMatrix(row, col, coefdelta);
        } else {
          rhs -= coefdelta * model->col_upper_[col];
          addToMatrix(row, col, -coefdelta);
        }
      }
    };

    // update / add non-zeros
    if (scale < 0) {
      HighsCDouble lhs = model->row_lower_[row];
      updateNonZeros(row, lhs, HighsInt{-1});
      model->row_lower_[row] = static_cast<double>(lhs);
    } else {
      HighsCDouble rhs = model->row_upper_[row];
      updateNonZeros(row, rhs, HighsInt{1});
      model->row_upper_[row] = static_cast<double>(rhs);
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

HPresolve::Result HPresolve::detectParallelRowsAndCols(
    HighsPostsolveStack& postsolve_stack) {
  assert(analysis_.allow_rule_[kPresolveRuleParallelRowsAndCols]);
  const bool logging_on = analysis_.logging_on_;
  if (logging_on)
    analysis_.startPresolveRuleLog(kPresolveRuleParallelRowsAndCols);

  std::vector<std::uint64_t> rowHashes;
  std::vector<std::uint64_t> colHashes;
  std::vector<std::pair<double, HighsInt>> rowMax(rowsize.size());
  std::vector<std::pair<double, HighsInt>> colMax(colsize.size());

  HighsHashTable<HighsInt, HighsInt> numRowSingletons;

  HighsInt nnz = Avalue.size();
  rowHashes.assign(rowsize.begin(), rowsize.end());
  colHashes.assign(colsize.begin(), colsize.end());

  // Step 1: Determine scales for rows and columns and remove column singletons
  // from the initial row hashes which are initialized with the row sizes
  for (HighsInt i = 0; i != nnz; ++i) {
    if (Avalue[i] == 0.0) continue;
    assert(!colDeleted[Acol[i]]);
    if (colsize[Acol[i]] == 1) {
      colMax[Acol[i]].first = Avalue[i];
      --rowHashes[Arow[i]];
      numRowSingletons[Arow[i]] += 1;
      continue;
    }
    double absVal = std::abs(Avalue[i]);
    double absRowMax = std::abs(rowMax[Arow[i]].first);

    // among the largest values which are equal in tolerance
    // we use the nonzero with the smallest row/column index for the column/row
    // scale so that we ensure that duplicate rows/columns are scaled to have
    // the same sign
    if (absVal >= absRowMax - options->small_matrix_value) {
      // we are greater or equal with tolerances, check if we are either
      // strictly larger or equal with a smaller index and remember the signed
      // nonzero if one of those things is the case
      if (absVal > absRowMax + options->small_matrix_value ||
          Acol[i] < rowMax[Arow[i]].second) {
        rowMax[Arow[i]].first = Avalue[i];
        rowMax[Arow[i]].second = Acol[i];
      }
    }

    double absColMax = std::abs(colMax[Acol[i]].first);
    if (absVal >= absColMax - options->small_matrix_value) {
      if (absVal > absColMax + options->small_matrix_value ||
          Arow[i] < colMax[Acol[i]].second) {
        colMax[Acol[i]].first = Avalue[i];
        colMax[Acol[i]].second = Arow[i];
      }
    }
  }

  // Step 2: Compute hash values for rows and columns excluding singleton
  // columns
  for (HighsInt i = 0; i != nnz; ++i) {
    if (Avalue[i] == 0.0) continue;
    assert(!rowDeleted[Arow[i]] && !colDeleted[Acol[i]]);
    if (colsize[Acol[i]] == 1) {
      colHashes[Acol[i]] = Arow[i];
    } else {
      HighsHashHelpers::sparse_combine(rowHashes[Arow[i]], Acol[i],
                                       HighsHashHelpers::double_hash_code(
                                           Avalue[i] / rowMax[Arow[i]].first));
      HighsHashHelpers::sparse_combine(colHashes[Acol[i]], Arow[i],
                                       HighsHashHelpers::double_hash_code(
                                           Avalue[i] / colMax[Acol[i]].first));
    }
  }

  // Step 3: Loop over the rows and columns and put them into buckets using the
  // computed hash values. Whenever a bucket already contains a row/column,
  // check if we can apply a (nearly) parallel row reduction or a
  // parallel/dominated column reduction.
  std::unordered_multimap<std::uint64_t, HighsInt> buckets;

  const bool debug_report = false;
  for (HighsInt i = 0; i != model->num_col_; ++i) {
    if (colDeleted[i]) continue;
    if (colsize[i] == 0) {
      HPRESOLVE_CHECKED_CALL(colPresolve(postsolve_stack, i));
      continue;
    }
    auto it = buckets.find(colHashes[i]);
    decltype(it) last = it;

    HighsInt delCol = -1;
    HighsInt parallelColCandidate = -2;

    while (it != buckets.end() && it->first == colHashes[i]) {
      parallelColCandidate = it->second;
      last = it++;

      // we want to check if the columns are parallel, first rule out
      // hash collisions with different size columns
      if (colsize[i] != colsize[parallelColCandidate]) continue;
      // The columns have the same length. Next we determine whether domination
      // is possible in one of the directions, and if it is we designate the
      // dominating column as column 2. The first thing we check is whether the
      // the objective value of one of the (scaled) columns is strictly better
      // then the objective value of the other column which rules out domination
      // in one direction.

      HighsInt col = -1;
      HighsInt duplicateCol = -1;
      double colScale;

      // helpers for checking dominance between parallel columns which is
      // possible for different cases of the variable types: if col can be
      // increased infinitely in which case duplicateCol can be fixed to its
      // lower bound. duplicateCol can be decreased infinitely in which case col
      // can be fixed to its upper bound. for both cases we exploit that the
      // column that remains unfixed can always compensate for the fixed column.
      // This only holds if the compensating column can compensate exactly for
      // feasible value of the fixed column. In the continuous case this
      // trivially holds. In the case where both variables are integer and the
      // scale is +- 1 this also holds trivially. If the scale is > 1 and both
      // variables are integer, this only holds in one direction. We can apply
      // the reduction due to the following reasoning: Applying the scale to
      // col, means we change its meaning and it is not an integer variable
      // anymore, but a variable that moves on multiples of 1/scale. As we have
      // taken care that the scale is >=1 and integral for two integer
      // variables, the scaled column can always exactly compensate for the
      // other column as it can move by 1/k with k being integer. Hence every
      // kth allowed value is integral and no integral value is skipped. If the
      // compensating column is integral
      bool checkColImplBounds = true;
      bool checkDuplicateColImplBounds = true;
      auto isLowerStrictlyImplied = [&](HighsInt col) {
        return (model->col_lower_[col] == -kHighsInf ||
                implColLower[col] > model->col_lower_[col] + primal_feastol);
      };
      auto isUpperStrictlyImplied = [&](HighsInt col) {
        return (model->col_upper_[col] == kHighsInf ||
                implColUpper[col] < model->col_upper_[col] - primal_feastol);
      };
      auto colUpperInf = [&]() {
        if (!checkColImplBounds) return false;
        if (mipsolver == nullptr) {
          // for LP we check strict redundancy of the bounds as otherwise dual
          // postsolve might fail when the bound is used in the optimal solution
          return colScale > 0 ? isUpperStrictlyImplied(col)
                              : isLowerStrictlyImplied(col);
        } else {
          // for MIP we do not need dual postsolve so the reduction is valid if
          // the bound is weakly redundant
          return colScale > 0 ? isUpperImplied(col) : isLowerImplied(col);
        }
      };

      auto colLowerInf = [&]() {
        if (!checkColImplBounds) return false;
        if (mipsolver == nullptr) {
          return colScale > 0 ? isLowerStrictlyImplied(col)
                              : isUpperStrictlyImplied(col);
        } else {
          return colScale > 0 ? isLowerImplied(col) : isUpperImplied(col);
        }
      };

      auto duplicateColUpperInf = [&]() {
        if (!checkDuplicateColImplBounds) return false;
        if (mipsolver == nullptr) {
          return isUpperStrictlyImplied(duplicateCol);
        } else {
          return isUpperImplied(duplicateCol);
        }
      };

      auto duplicateColLowerInf = [&]() {
        if (!checkDuplicateColImplBounds) return false;
        if (mipsolver == nullptr) {
          return isLowerStrictlyImplied(duplicateCol);
        } else {
          return isLowerImplied(duplicateCol);
        }
      };

      // Now check the if the variable types rule out domination in one
      // direction and already skip the column if that rules out domination in
      // both directions due to the previous check on the objective.
      if (model->integrality_[i] == HighsVarType::kInteger &&
          model->integrality_[parallelColCandidate] == HighsVarType::kInteger) {
        // both variables are integral, hence the scale must be integral
        // therefore first choose the smaller colMax value for col2, then check
        // integrality of colMax[col1] / colMax[col2].
        if (std::abs(colMax[i].first) <
            std::abs(colMax[parallelColCandidate].first)) {
          col = i;
          duplicateCol = parallelColCandidate;
        } else {
          col = parallelColCandidate;
          duplicateCol = i;
        }

        double scaleCand = colMax[duplicateCol].first / colMax[col].first;
        if (fractionality(scaleCand, &colScale) > options->small_matrix_value)
          continue;
        assert(std::abs(colScale) >= 1.0);

        // if the scale is larger than 1, duplicate column cannot compensate for
        // all values of scaled col due to integrality as the scaled column
        // moves on a grid of 1/scale.
        //
        // ToDo: Check whether this is too restrictive
        if (colScale != 1.0) checkDuplicateColImplBounds = false;
      } else if (model->integrality_[i] == HighsVarType::kInteger) {
        col = i;
        duplicateCol = parallelColCandidate;
        colScale = colMax[duplicateCol].first / colMax[col].first;

        // as col is integral and duplicateCol is not col cannot compensate for
        // duplicate col
        checkColImplBounds = false;
      } else {
        col = parallelColCandidate;
        duplicateCol = i;
        colScale = colMax[duplicateCol].first / colMax[col].first;

        // as col might be integral and duplicateCol is not integral. In that
        // case col cannot compensate for duplicate col
        checkColImplBounds =
            model->integrality_[parallelColCandidate] != HighsVarType::kInteger;
      }

      double objDiff = static_cast<double>(
          model->col_cost_[col] * static_cast<HighsCDouble>(colScale) -
          model->col_cost_[duplicateCol]);
      // if (std::abs(objDiff) > options->small_matrix_value) continue;
      constexpr HighsInt kMergeParallelCols = 0;
      constexpr HighsInt kDominanceColToUpper = 1;
      constexpr HighsInt kDominanceColToLower = 2;
      constexpr HighsInt kDominanceDuplicateColToUpper = 3;
      constexpr HighsInt kDominanceDuplicateColToLower = 4;

      HighsInt reductionCase = kMergeParallelCols;
      // now do the case distinctions for dominated columns
      // the cases are a lot simpler due to the helper functions
      // for checking the infinite bounds which automatically
      // incorporate the check for the variable types that allow domination.
      if (objDiff < -options->dual_feasibility_tolerance) {
        // scaled col is better than duplicate col
        if (colUpperInf() && model->col_lower_[duplicateCol] != -kHighsInf)
          reductionCase = kDominanceDuplicateColToLower;
        else if (duplicateColLowerInf() &&
                 (colScale < 0 || model->col_upper_[col] != kHighsInf) &&
                 (colScale > 0 || model->col_lower_[col] != -kHighsInf))
          reductionCase =
              colScale > 0 ? kDominanceColToUpper : kDominanceColToLower;
        else
          continue;
      } else if (objDiff > options->dual_feasibility_tolerance) {
        // duplicate col is better than scaled col
        if (colLowerInf() && model->col_upper_[duplicateCol] != kHighsInf)
          reductionCase = kDominanceDuplicateColToUpper;
        else if (duplicateColUpperInf() &&
                 (colScale < 0 || model->col_lower_[col] != -kHighsInf) &&
                 (colScale > 0 || model->col_upper_[col] != kHighsInf))
          reductionCase =
              colScale > 0 ? kDominanceColToLower : kDominanceColToUpper;
        else
          continue;
      } else {
        if (colUpperInf() && model->col_lower_[duplicateCol] != -kHighsInf)
          reductionCase = kDominanceDuplicateColToLower;
        else if (colLowerInf() && model->col_upper_[duplicateCol] != kHighsInf)
          reductionCase = kDominanceDuplicateColToUpper;
        else if (duplicateColUpperInf() &&
                 (colScale < 0 || model->col_lower_[col] != -kHighsInf) &&
                 (colScale > 0 || model->col_upper_[col] != kHighsInf))
          reductionCase =
              colScale > 0 ? kDominanceColToLower : kDominanceColToUpper;
        else if (duplicateColLowerInf() &&
                 (colScale < 0 || model->col_upper_[col] != kHighsInf) &&
                 (colScale > 0 || model->col_lower_[col] != -kHighsInf))
          reductionCase =
              colScale > 0 ? kDominanceColToUpper : kDominanceColToLower;
      }
      if (reductionCase == kMergeParallelCols) {
        const bool x_int = model->integrality_[col] == HighsVarType::kInteger;
        const bool y_int =
            model->integrality_[duplicateCol] == HighsVarType::kInteger;
        bool illegal_scale = true;
        if (x_int) {
          // The only possible reduction if the column parallelism check
          // succeeds is to merge the two columns into one. If one column is
          // integral this means we have restrictions on integers and need to
          // check additional conditions to allow the merging of two integer
          // columns, or a continuous column and an integer.
          if (model->integrality_[duplicateCol] != HighsVarType::kInteger) {
            assert(!y_int);
            // only one column is integral which cannot be duplicateCol due to
            // the way we assign the columns above
            //
            // Scale must not exceed 1/(y_u-y_l) in magnitude
            illegal_scale =
                std::abs(colScale * (model->col_upper_[duplicateCol] -
                                     model->col_lower_[duplicateCol])) <
                1.0 - primal_feastol;
            if (!illegal_scale && debug_report)
              printf(
                  "kMergeParallelCols: T-F is %s legal with scale %.4g and "
                  "duplicateCol = [%.4g, %.4g]\n",
                  illegal_scale ? "not" : "   ", colScale,
                  model->col_lower_[duplicateCol],
                  model->col_upper_[duplicateCol]);
          } else {
            // Both columns integer
            assert(x_int && y_int);
            // Scale must be integer and not exceed (x_u-x_l)+1 in magnitude
            const double scale_limit = model->col_upper_[col] -
                                       model->col_lower_[col] + 1 +
                                       primal_feastol;
            illegal_scale = std::fabs(colScale) > scale_limit;
          }
          if (illegal_scale) continue;
        } else {
          // Neither column integer: no problem with
          assert(!x_int && !y_int);
        }
      }

      bool parallel = true;
      // now check whether the coefficients are actually parallel
      for (const HighsSliceNonzero& colNz : getColumnVector(col)) {
        HighsInt duplicateColRowPos = findNonzero(colNz.index(), duplicateCol);
        parallel = duplicateColRowPos != -1;
        if (!parallel) break;

        parallel = std::abs(static_cast<double>(
                       Avalue[duplicateColRowPos] -
                       static_cast<HighsCDouble>(colScale) * colNz.value())) <=
                   options->small_matrix_value;
        if (!parallel) break;
      }

      if (!parallel) continue;

      switch (reductionCase) {
        case kDominanceDuplicateColToLower:
          delCol = duplicateCol;
          if (colsize[duplicateCol] == 1) {
            HighsInt row = Arow[colhead[duplicateCol]];
            numRowSingletons[row] -= 1;
          }
          if (fixColToLowerOrUnbounded(postsolve_stack, duplicateCol)) {
            // Handle unboundedness
            presolve_status_ = HighsPresolveStatus::kUnboundedOrInfeasible;
            return Result::kDualInfeasible;
          }
          break;
        case kDominanceDuplicateColToUpper:
          delCol = duplicateCol;
          if (colsize[duplicateCol] == 1) {
            HighsInt row = Arow[colhead[duplicateCol]];
            numRowSingletons[row] -= 1;
          }
          if (fixColToUpperOrUnbounded(postsolve_stack, duplicateCol)) {
            // Handle unboundedness
            presolve_status_ = HighsPresolveStatus::kUnboundedOrInfeasible;
            return Result::kDualInfeasible;
          }
          break;
        case kDominanceColToLower:
          delCol = col;
          if (colsize[col] == 1) {
            HighsInt row = Arow[colhead[col]];
            numRowSingletons[row] -= 1;
          }
          if (fixColToLowerOrUnbounded(postsolve_stack, col)) {
            // Handle unboundedness
            presolve_status_ = HighsPresolveStatus::kUnboundedOrInfeasible;
            return Result::kDualInfeasible;
          }
          break;
        case kDominanceColToUpper:
          delCol = col;
          if (colsize[col] == 1) {
            HighsInt row = Arow[colhead[col]];
            numRowSingletons[row] -= 1;
          }
          if (fixColToUpperOrUnbounded(postsolve_stack, col)) {
            // Handle unboundedness
            presolve_status_ = HighsPresolveStatus::kUnboundedOrInfeasible;
            return Result::kDualInfeasible;
          }
          break;
        case kMergeParallelCols:
          const bool ok_merge = postsolve_stack.duplicateColumn(
              colScale, model->col_lower_[col], model->col_upper_[col],
              model->col_lower_[duplicateCol], model->col_upper_[duplicateCol],
              col, duplicateCol,
              model->integrality_[col] == HighsVarType::kInteger,
              model->integrality_[duplicateCol] == HighsVarType::kInteger,
              options->mip_feasibility_tolerance);
          if (!ok_merge && debug_report) {
            printf(
                "HPresolve::detectParallelRowsAndCols Illegal merge "
                "prevented\n");
            break;
          }
          // When merging a continuous variable into an integer
          // variable, the integer will become continuous - since any
          // value in its range can be mapped back to an integer and a
          // continuous variable. Hence the number of integer
          // variables in the rows corresponding to the former integer
          // variable reduces.
          //
          // With the opposite - merging an integer variable into a
          // continuous variable - the retained variable is
          // continuous, so no action is required
          bool rowsizeIntReduction =
              model->integrality_[duplicateCol] != HighsVarType::kInteger &&
              model->integrality_[col] == HighsVarType::kInteger;
          if (rowsizeIntReduction)
            model->integrality_[col] = HighsVarType::kContinuous;

          markChangedCol(col);
          if (colsize[duplicateCol] == 1) {
            HighsInt row = Arow[colhead[duplicateCol]];
            numRowSingletons[row] -= 1;
          }

          // compute bounds of merged variable
          double mergeLower = 0;
          double mergeUpper = 0;
          if (colScale > 0) {
            mergeLower = model->col_lower_[col] +
                         colScale * model->col_lower_[duplicateCol];
            mergeUpper = model->col_upper_[col] +
                         colScale * model->col_upper_[duplicateCol];
          } else {
            mergeLower = model->col_lower_[col] +
                         colScale * model->col_upper_[duplicateCol];
            mergeUpper = model->col_upper_[col] +
                         colScale * model->col_lower_[duplicateCol];
          }

          // change bounds
          changeColLower(col, mergeLower);
          changeColUpper(col, mergeUpper);

          // mark duplicate column as deleted
          markColDeleted(duplicateCol);
          // remove all nonzeros of duplicateCol
          for (HighsInt coliter = colhead[duplicateCol]; coliter != -1;) {
            assert(Acol[coliter] == duplicateCol);

            HighsInt colpos = coliter;
            HighsInt colrow = Arow[coliter];

            // if an integer column was merged into a continuous one make
            // sure to update the integral rowsize
            if (rowsizeIntReduction) rowsizeInteger[colrow] -= 1;

            coliter = Anext[coliter];

            unlink(colpos);

            reinsertEquation(colrow);
          }
          // set cost to zero
          model->col_cost_[duplicateCol] = 0;
          delCol = duplicateCol;

          // remove implied bounds, since they might in general not be valid
          // anymore
          resetColImpliedBounds(col);

          // if an implicit integer and an integer column were merged, check if
          // merged continuous column is implicit integer after merge
          if (rowsizeIntReduction && model->integrality_[duplicateCol] ==
                                         HighsVarType::kImplicitInteger) {
            StatusResult impliedInteger = isImpliedInteger(col);
            HPRESOLVE_CHECKED_CALL(static_cast<Result>(impliedInteger));
            if (impliedInteger) convertImpliedInteger(col, -1, true);
          }
          break;
      }

      break;
    }

    if (delCol != -1) {
      if (delCol != i) buckets.erase(last);

      // we could have new row singletons since a column was removed. Remove
      // those rows immediately
      HPRESOLVE_CHECKED_CALL(checkLimits(postsolve_stack));
      HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolve_stack));
    } else {
      buckets.emplace_hint(last, colHashes[i], i);
    }
  }

  buckets.clear();

  for (HighsInt i = 0; i != model->num_row_; ++i) {
    if (rowDeleted[i]) continue;
    if (rowsize[i] <= 1 || (rowsize[i] == 2 && isEquation(i))) {
      HPRESOLVE_CHECKED_CALL(rowPresolve(postsolve_stack, i));
      continue;
    }
    auto it = buckets.find(rowHashes[i]);
    decltype(it) last = it;

    auto getNumSingletons = [&](HighsInt row) {
      const HighsInt* numSingletonPtr = numRowSingletons.find(row);
      return (numSingletonPtr ? *numSingletonPtr : 0);
    };

    const HighsInt numSingleton = getNumSingletons(i);

#if !ENABLE_SPARSIFY_FOR_LP
    if (mipsolver == nullptr && options->lp_presolve_requires_basis_postsolve &&
        numSingleton != 0)
      continue;
#endif
    HighsInt delRow = -1;
    if (it != buckets.end()) storeRow(i);
    while (it != buckets.end() && it->first == rowHashes[i]) {
      HighsInt parallelRowCand = it->second;
      last = it++;

      const HighsInt numSingletonCandidate = getNumSingletons(parallelRowCand);

#if !ENABLE_SPARSIFY_FOR_LP
      if (mipsolver == nullptr &&
          options->lp_presolve_requires_basis_postsolve &&
          numSingletonCandidate != 0)
        continue;
#endif
      if (rowsize[i] - numSingleton !=
          rowsize[parallelRowCand] - numSingletonCandidate)
        continue;

      if (numSingletonCandidate > 1 || numSingleton > 1) {
        // we only handle the case where the rows have at most one extra
        // singleton except when one row has no extra singleton and is an
        // equation. In that case we sparsify the other row by adding the
        // equation and can subsequently solve it as an individual component as
        // it is a row which only contains singletons
        if ((numSingleton != 0 || !isEquation(i)) &&
            (numSingletonCandidate != 0 || !isEquation(parallelRowCand)))
          continue;
      } else if (numSingletonCandidate != numSingleton) {
        // if only one of the two constraints has an extra singleton,
        // we require at least one of the constraints to be an equation
        // if that is the case we can add that equation to the other row
        // and will make it into either a row singleton or a doubleton equation
        // which is removed afterwards
        if (!isEquation(i) && !isEquation(parallelRowCand)) continue;
      }

      double rowScale = rowMax[parallelRowCand].first / rowMax[i].first;
      // check parallel case
      bool parallel = true;
      for (const HighsSliceNonzero& rowNz : getStoredRow()) {
        if (colsize[rowNz.index()] == 1)  // skip singletons
          continue;
        HighsInt nzPos = findNonzero(parallelRowCand, rowNz.index());
        parallel = nzPos != -1;
        if (!parallel) break;

        parallel = std::abs(static_cast<double>(
                       Avalue[nzPos] -
                       static_cast<HighsCDouble>(rowScale) * rowNz.value())) <=
                   options->small_matrix_value;
        if (!parallel) break;
      }
      if (!parallel) continue;

      if (numSingleton == 0 && numSingletonCandidate == 0) {
        bool rowLowerTightened = false;
        bool rowUpperTightened = false;
        double newUpper;
        double newLower;
        if (rowScale > 0) {
          newUpper = model->row_upper_[i] * rowScale;
          newLower = model->row_lower_[i] * rowScale;
        } else {
          newLower = model->row_upper_[i] * rowScale;
          newUpper = model->row_lower_[i] * rowScale;
        }

        if (newUpper < model->row_upper_[parallelRowCand]) {
          if (newUpper < model->row_lower_[parallelRowCand] - primal_feastol)
            return Result::kPrimalInfeasible;

          if (newUpper <= model->row_lower_[parallelRowCand] + primal_feastol)
            newUpper = model->row_lower_[parallelRowCand];

          if (newUpper < model->row_upper_[parallelRowCand]) {
            rowUpperTightened = true;
            if (rowScale > 0) {
              double tmp = rowDualLower[i] / rowScale;
              rowDualLower[i] = rowDualLower[parallelRowCand] * rowScale;
              rowDualLower[parallelRowCand] = tmp;
            } else {
              double tmp = rowDualUpper[i] / rowScale;
              rowDualUpper[i] = rowDualLower[parallelRowCand] * rowScale;
              rowDualLower[parallelRowCand] = tmp;
            }

            model->row_upper_[parallelRowCand] = newUpper;
          }
        }

        if (newLower > model->row_lower_[parallelRowCand]) {
          if (newLower > model->row_upper_[parallelRowCand] + primal_feastol)
            return Result::kPrimalInfeasible;

          if (newLower >= model->row_upper_[parallelRowCand] - primal_feastol)
            newLower = model->row_upper_[parallelRowCand];

          if (newLower > model->row_lower_[parallelRowCand]) {
            // the rows lower bound is tightened
            // instead of updating the activities of dual constraints, we
            // can simply swap the bounds on the row duals. If the old
            // lower bound on the row dual was finite, the new row dual
            // lower bound is infinite as the new row lower bound must be
            // a finite value. This infinite contribution, was, however,
            // already counted from the parallel row. Therefore by
            // swapping the bounds unlinking the other row will not
            // decrease the infinity counter, but simply remove a bound
            // with zero contribution. For a negative scale we need to
            // swap with the negated upper bound of the row dual of row i.
            rowLowerTightened = true;
            if (rowScale > 0) {
              double tmp = rowDualUpper[i] / rowScale;
              rowDualUpper[i] = rowDualUpper[parallelRowCand] * rowScale;
              rowDualUpper[parallelRowCand] = tmp;
            } else {
              double tmp = rowDualLower[i] / rowScale;
              rowDualLower[i] = rowDualUpper[parallelRowCand] * rowScale;
              rowDualUpper[parallelRowCand] = tmp;
            }

            model->row_lower_[parallelRowCand] = newLower;
          }
        }
        // remove implied bounds, since they might in general not be valid
        // anymore
        resetRowDualImpliedBounds(parallelRowCand);

        postsolve_stack.duplicateRow(parallelRowCand, rowUpperTightened,
                                     rowLowerTightened, i, rowScale);
        delRow = i;
        markRowDeleted(i);
        for (HighsInt rowiter : rowpositions) unlink(rowiter);
        break;
      } else if (isEquation(i)) {
        // row i is equation and parallel (except for singletons)
        // add to the row parallelRowCand
        // printf(
        //    "nearly parallel case with %" HIGHSINT_FORMAT " singletons in eq
        //    row and %" HIGHSINT_FORMAT " " "singletons in other row(eq=%"
        //    HIGHSINT_FORMAT ")\n", numSingleton, numSingletonCandidate,
        //    model->row_lower_[parallelRowCand] ==
        //        model->row_upper_[parallelRowCand]);
        HPRESOLVE_CHECKED_CALL(equalityRowAddition(
            postsolve_stack, i, parallelRowCand, -rowScale, getStoredRow()));
        delRow = parallelRowCand;
      } else if (isEquation(parallelRowCand)) {
        // printf(
        //    "nearly parallel case with %" HIGHSINT_FORMAT " singletons in eq
        //    row and %" HIGHSINT_FORMAT " " "singletons in other inequality
        //    row\n", numSingletonCandidate, numSingleton);
        // the row parallelRowCand is an equation; add it to the other row
        HPRESOLVE_CHECKED_CALL(equalityRowAddition(
            postsolve_stack, parallelRowCand, i,
            -rowMax[i].first / rowMax[parallelRowCand].first,
            getRowVector(parallelRowCand)));
        delRow = i;
      } else {
        assert(numSingleton == 1);
        assert(numSingletonCandidate == 1);

        double rowUpper;
        double rowLower;
        if (rowScale > 0) {
          rowUpper = model->row_upper_[i] * rowScale;
          rowLower = model->row_lower_[i] * rowScale;
        } else {
          rowLower = model->row_upper_[i] * rowScale;
          rowUpper = model->row_lower_[i] * rowScale;
        }
        // todo: two inequalities with one singleton. check whether the rows can
        // be converted to equations by introducing a shared slack variable
        // which is the case if the singletons have similar properties
        // (objective sign, bounds, scaled coefficient) and the scaled right
        // hand sides match. Then the case reduces to adding one equation to the
        // other and substituting one of the singletons due to the resulting
        // doubleton equation.
        //        printf("todo, two inequalities with one additional
        //        singleton\n");
        (void)rowLower;
        (void)rowUpper;
      }
    }

    if (delRow != -1) {
      if (delRow != i) buckets.erase(last);

      HPRESOLVE_CHECKED_CALL(checkLimits(postsolve_stack));
    } else
      buckets.emplace_hint(last, rowHashes[i], i);
  }

  analysis_.logging_on_ = logging_on;
  if (logging_on)
    analysis_.stopPresolveRuleLog(kPresolveRuleParallelRowsAndCols);

  return Result::kOk;
}

template <typename RowStorageFormat>
HPresolve::Result HPresolve::equalityRowAddition(
    HighsPostsolveStack& postsolve_stack, HighsInt stayrow, HighsInt removerow,
    double scale, const HighsMatrixSlice<RowStorageFormat>& rowvector) {
  // extract non-zero positions
  std::vector<HighsInt> stay_rowpositions;
  getRowPositions(stayrow, stay_rowpositions);

  // update postsolve information
  postsolve_stack.equalityRowAddition(removerow, stayrow, scale, rowvector);

  // iterate over non-zero positions instead of iterating over the
  // HighsMatrixSlice because the latter contains pointers to Acol and Avalue
  // that may be invalidated if these vectors are reallocated
  // (see std::vector::push_back performed in HPresolve::addToMatrix).
  for (HighsInt rowiter : stay_rowpositions) {
    HighsInt pos = findNonzero(removerow, Acol[rowiter]);
    if (pos != -1)
      unlink(pos);  // all common nonzeros are cancelled, as the rows are
                    // parallel
    else            // might introduce a singleton
      addToMatrix(removerow, Acol[rowiter], scale * Avalue[rowiter]);
  }

  if (model->row_upper_[removerow] != kHighsInf)
    model->row_upper_[removerow] = static_cast<double>(
        model->row_upper_[removerow] +
        static_cast<HighsCDouble>(scale) * model->row_upper_[stayrow]);
  if (model->row_lower_[removerow] != -kHighsInf)
    model->row_lower_[removerow] = static_cast<double>(
        model->row_lower_[removerow] +
        static_cast<HighsCDouble>(scale) * model->row_upper_[stayrow]);

  // row is now a singleton row, doubleton equation, or a row
  // that contains only singletons and we let the normal row presolve
  // handle the cases
  HPRESOLVE_CHECKED_CALL(rowPresolve(postsolve_stack, removerow));
  return Result::kOk;
}

void HPresolve::setRelaxedImpliedBounds() {
  double hugeBound = primal_feastol / kHighsTiny;
  for (HighsInt i = 0; i != model->num_col_; ++i) {
    if (model->col_lower_[i] >= implColLower[i] &&
        model->col_upper_[i] <= implColUpper[i])
      continue;

    if (std::abs(implColLower[i]) <= hugeBound) {
      // if the bound is derived from a small nonzero value
      // then we want to increase the margin so that we make sure
      // the row it was derived from is violated if the column sits
      // at this relaxed bound in the final solution.
      HighsInt nzPos = findNonzero(colLowerSource[i], i);

      double boundRelax = std::max(1000.0, std::abs(implColLower[i])) *
                          primal_feastol /
                          std::min(1.0, std::abs(Avalue[nzPos]));

      double newLb = implColLower[i] - boundRelax;
      if (newLb > model->col_lower_[i] + boundRelax)
        model->col_lower_[i] = newLb;
    }

    if (std::abs(implColUpper[i]) <= hugeBound) {
      HighsInt nzPos = findNonzero(colUpperSource[i], i);

      double boundRelax = std::max(1000.0, std::abs(implColUpper[i])) *
                          primal_feastol /
                          std::min(1.0, std::abs(Avalue[nzPos]));

      double newUb = implColUpper[i] + boundRelax;
      if (newUb < model->col_upper_[i] - boundRelax)
        model->col_upper_[i] = newUb;
    }
  }
}

// Not currently called
void HPresolve::debug(const HighsLp& lp, const HighsOptions& options) {
  HighsSolution reducedsol;
  HighsBasis reducedbasis;

  HighsSolution sol;
  HighsBasis basis;

  HighsLp model = lp;
  model.integrality_.assign(lp.num_col_, HighsVarType::kContinuous);

  HighsPostsolveStack postsolve_stack;
  postsolve_stack.initializeIndexMaps(lp.num_row_, lp.num_col_);
  {
    HPresolve presolve;
    presolve.okSetInput(model, options, options.presolve_reduction_limit);
    // presolve.setReductionLimit(1622017);
    if (presolve.run(postsolve_stack) != HighsModelStatus::kNotset) return;
    Highs highs;
    highs.passModel(model);
    highs.passOptions(options);
    highs.setOptionValue("presolve", kHighsOffString);
    highs.run();
    if (highs.getModelStatus() != HighsModelStatus::kOptimal) return;
    reducedsol = highs.getSolution();
    reducedbasis = highs.getBasis();
  }
  model = lp;
  sol = reducedsol;
  basis = reducedbasis;
  postsolve_stack.undo(options, sol, basis);
  refineBasis(lp, sol, basis);
  calculateRowValuesQuad(model, sol);
#if 0
  Highs highs;
  highs.passModel(model);
  highs.passOptions(options);
  highs.setSolution(sol);
  basis.debug_origin_name = "HPresolve::debug";
  highs.setBasis(basis);
  highs.run();
  return;
#endif
  std::vector<HighsInt> flagCol(lp.num_col_, 1);
  std::vector<HighsInt> flagRow(lp.num_row_, 1);
  std::vector<HighsInt> Aend;
  std::vector<HighsInt> ARstart;
  std::vector<HighsInt> ARindex;
  std::vector<double> ARvalue;
  dev_kkt_check::KktInfo kktinfo = dev_kkt_check::initInfo();
  Aend.assign(model.a_matrix_.start_.begin() + 1, model.a_matrix_.start_.end());
  highsSparseTranspose(model.num_row_, model.num_col_, model.a_matrix_.start_,
                       model.a_matrix_.index_, model.a_matrix_.value_, ARstart,
                       ARindex, ARvalue);
  dev_kkt_check::State state(
      model.num_col_, model.num_row_, model.a_matrix_.start_, Aend,
      model.a_matrix_.index_, model.a_matrix_.value_, ARstart, ARindex, ARvalue,
      model.col_cost_, model.col_lower_, model.col_upper_, model.row_lower_,
      model.row_upper_, flagCol, flagRow, sol.col_value, sol.col_dual,
      sol.row_value, sol.row_dual, basis.col_status, basis.row_status);
  bool checkResult = dev_kkt_check::checkKkt(state, kktinfo);
  if (checkResult && kktinfo.pass_bfs) {
    printf("kkt check of postsolved solution and basis passed\n");
    return;
  }
  size_t good = postsolve_stack.numReductions();
  size_t bad = 0;
  size_t reductionLim = (good + bad) / 2;

  // good = 1734357, bad = 1734289;
  // good = 1050606, bad = 1050605;
  // good = 1811527, bad = 1811526;
  // reductionLim = bad;
  do {
    model = lp;
    model.integrality_.assign(lp.num_col_, HighsVarType::kContinuous);

    {
      HPresolve presolve;
      presolve.okSetInput(model, options, options.presolve_reduction_limit);
      presolve.computeIntermediateMatrix(flagRow, flagCol, reductionLim);
    }
#if 1
    model = lp;
    model.integrality_.assign(lp.num_col_, HighsVarType::kContinuous);
    HPresolve presolve;
    presolve.okSetInput(model, options, options.presolve_reduction_limit);
    HighsPostsolveStack tmp;
    tmp.initializeIndexMaps(model.num_row_, model.num_col_);
    presolve.setReductionLimit(reductionLim);
    presolve.run(tmp);

    sol = reducedsol;
    basis = reducedbasis;
    postsolve_stack.undoUntil(options, sol, basis, tmp.numReductions());

    HighsBasis temp_basis;
    HighsSolution temp_sol;
    temp_basis.col_status.resize(model.num_col_);
    temp_sol.col_dual.resize(model.num_col_);
    temp_sol.col_value.resize(model.num_col_);
    for (HighsInt i = 0; i != model.num_col_; ++i) {
      temp_sol.col_dual[i] = sol.col_dual[tmp.getOrigColIndex(i)];
      temp_sol.col_value[i] = sol.col_value[tmp.getOrigColIndex(i)];
      temp_basis.col_status[i] = basis.col_status[tmp.getOrigColIndex(i)];
    }

    temp_basis.row_status.resize(model.num_row_);
    temp_sol.row_dual.resize(model.num_row_);
    for (HighsInt i = 0; i != model.num_row_; ++i) {
      temp_sol.row_dual[i] = sol.row_dual[tmp.getOrigRowIndex(i)];
      temp_basis.row_status[i] = basis.row_status[tmp.getOrigRowIndex(i)];
    }
    temp_sol.row_value.resize(model.num_row_);
    calculateRowValuesQuad(model, sol);
    temp_basis.valid = true;
    temp_basis.useful = true;
    refineBasis(model, temp_sol, temp_basis);
    Highs highs;
    highs.passOptions(options);
    highs.passModel(model);
    temp_basis.debug_origin_name = "HPresolve::debug";
    highs.setBasis(temp_basis);
    // highs.writeModel("model.mps");
    // highs.writeBasis("bad.bas");
    highs.run();
    printf("simplex iterations with postsolved basis: %" HIGHSINT_FORMAT "\n",
           highs.getInfo().simplex_iteration_count);
    checkResult = highs.getInfo().simplex_iteration_count == 0;
#else

    if (reductionLim == good) break;

    Aend.assign(model.a_matrix_.start_.begin() + 1,
                model.a_matrix_.start_.end());
    highsSparseTranspose(model.num_row_, model.num_col_, model.a_matrix_.start_,
                         model.a_matrix_.index_, model.a_matrix_.value_,
                         ARstart, ARindex, ARvalue);
    sol = reducedsol;
    basis = reducedbasis;
    postsolve_stack.undoUntil(options, sol, basis, reductionLim);

    calculateRowValuesQuad(model, sol);
    kktinfo = dev_kkt_check::initInfo();
    checkResult = dev_kkt_check::checkKkt(state, kktinfo);
    checkResult = checkResult && kktinfo.pass_bfs;
#endif
    if (bad == good - 1) break;

    if (checkResult) {
      good = reductionLim;
    } else {
      bad = reductionLim;
    }
    reductionLim = (bad + good) / 2;
    printf("binary search ongoing: good=%zu, bad=%zu\n", good, bad);
  } while (true);

  printf("binary search finished: good=%zu, bad=%zu\n", good, bad);
  assert(false);
}

HPresolve::Result HPresolve::sparsify(HighsPostsolveStack& postsolve_stack) {
  std::vector<HighsPostsolveStack::Nonzero> sparsifyRows;
  HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolve_stack));
  HPRESOLVE_CHECKED_CALL(removeDoubletonEquations(postsolve_stack));
  std::vector<HighsInt> tmpEquations;
  tmpEquations.reserve(equations.size());

  const double minNonzeroVal = std::sqrt(primal_feastol);

  for (const auto& eq : equations) tmpEquations.emplace_back(eq.second);
  for (HighsInt eqrow : tmpEquations) {
    if (rowDeleted[eqrow]) continue;

    assert(!rowDeleted[eqrow]);
    assert(isEquation(eqrow));

    storeRow(eqrow);

    HighsInt sparsestCol = -1;
    HighsInt secondSparsestCol = -1;
    HighsInt sparsestColLen = kHighsIInf;
    HighsInt secondSparsestColLen = kHighsIInf;
    for (HighsInt nzPos : rowpositions) {
      HighsInt col = Acol[nzPos];
      if (colsize[col] < sparsestColLen) {
        secondSparsestCol = sparsestCol;
        secondSparsestColLen = sparsestColLen;
        sparsestCol = col;
        sparsestColLen = colsize[col];
      } else if (colsize[col] < secondSparsestColLen) {
        secondSparsestCol = col;
        secondSparsestColLen = colsize[col];
      }
    }

    assert(sparsestCol != -1 && secondSparsestCol != -1);
    assert(colsize[sparsestCol] <= colsize[secondSparsestCol]);

    std::map<double, HighsInt> possibleScales;
    sparsifyRows.clear();

    for (const HighsSliceNonzero& colNz : getColumnVector(sparsestCol)) {
      HighsInt candRow = colNz.index();
      if (candRow == eqrow) continue;

      possibleScales.clear();

      HighsInt misses = 0;
      // allow no fillin if a completely continuous row is used to cancel a row
      // that has integers as there are instances where this leads to a huge
      // deterioration of cut performance
      HighsInt maxMisses = 1;
      if (rowsizeInteger[eqrow] == 0 && rowsizeInteger[candRow] != 0)
        --maxMisses;
      for (const HighsSliceNonzero& nonzero : getStoredRow()) {
        double candRowVal;
        if (nonzero.index() == sparsestCol) {
          candRowVal = colNz.value();
        } else {
          HighsInt nzPos = findNonzero(candRow, nonzero.index());
          if (nzPos == -1) {
            if (model->integrality_[nonzero.index()] ==
                    HighsVarType::kInteger &&
                model->col_upper_[nonzero.index()] -
                        model->col_lower_[nonzero.index()] >
                    1.5) {
              // do not allow fillin of general integers
              misses = 2;
              break;
            }
            ++misses;
            if (misses > maxMisses) break;
            continue;
          }
          candRowVal = Avalue[nzPos];
        }

        double scale = -candRowVal / nonzero.value();
        if (std::abs(scale) > 1e3) continue;

        double scaleTolerance = minNonzeroVal / std::abs(nonzero.value());
        auto it = possibleScales.lower_bound(scale - scaleTolerance);
        if (it != possibleScales.end() &&
            std::abs(it->first - scale) <= scaleTolerance) {
          // there already is a scale that is very close and could produce
          // a matrix value for this nonzero that is below the allowed
          // threshold. Therefore we check if the matrix value is small enough
          // for this nonzero to be deleted, in which case the number of
          // deleted nonzeros for the other scale is increased. If it is not
          // small enough we do not use this scale or the other one because
          // such small matrix values may lead to numerical troubles.

          // scale is already marked to be numerically bad
          if (it->second == -1) continue;

          if (std::abs(it->first * nonzero.value() + candRowVal) <=
              options->small_matrix_value)
            it->second += 1;
          else
            it->second = -1;
        } else
          possibleScales.emplace(scale, 1);
      }

      if (misses > maxMisses || possibleScales.empty()) continue;

      HighsInt numCancel = 0;
      double scale = 0.0;

      for (const auto& s : possibleScales) {
        if (s.second <= misses) continue;

        if (s.second > numCancel ||
            (s.second == numCancel && std::abs(s.first) < std::abs(scale))) {
          scale = s.first;
          numCancel = s.second;
        }
      }

      assert(scale != 0.0 || numCancel == 0);

      // cancels at least one nonzero if the scale cancels more than there is
      // fillin
      if (numCancel > misses) sparsifyRows.emplace_back(candRow, scale);
    }

    if (model->integrality_[sparsestCol] != HighsVarType::kInteger ||
        (model->col_upper_[sparsestCol] - model->col_lower_[sparsestCol]) <
            1.5) {
      // now check for rows which do not contain the sparsest column but all
      // other columns by scanning the second sparsest column
      for (const HighsSliceNonzero& colNz :
           getColumnVector(secondSparsestCol)) {
        HighsInt candRow = colNz.index();
        if (candRow == eqrow) continue;

        if (rowsizeInteger[eqrow] == 0 && rowsizeInteger[candRow] != 0)
          continue;

        HighsInt sparsestColPos = findNonzero(candRow, sparsestCol);

        // if the row has a nonzero for the sparsest column we have already
        // checked it
        if (sparsestColPos != -1) continue;

        possibleScales.clear();
        bool skip = false;
        for (const HighsSliceNonzero& nonzero : getStoredRow()) {
          double candRowVal;
          if (nonzero.index() == secondSparsestCol) {
            candRowVal = colNz.value();
          } else {
            HighsInt nzPos = findNonzero(candRow, nonzero.index());
            // we already have a miss for the sparsest column, so with another
            // one we want to skip the row
            skip = nzPos == -1;
            if (skip) break;

            candRowVal = Avalue[nzPos];
          }

          double scale = -candRowVal / nonzero.value();
          if (std::abs(scale) > 1e3) continue;

          double scaleTolerance = minNonzeroVal / std::abs(nonzero.value());
          auto it = possibleScales.lower_bound(scale - scaleTolerance);
          if (it != possibleScales.end() &&
              std::abs(it->first - scale) <= scaleTolerance) {
            // there already is a scale that is very close and could produce
            // a matrix value for this nonzero that is below the allowed
            // threshold. Therefore we check if the matrix value is small enough
            // for this nonzero to be deleted, in which case the number of
            // deleted nonzeros for the other scale is increased. If it is not
            // small enough we do not use this scale or the other one because
            // such small matrix values may lead to numerical troubles.

            // scale is already marked to be numerically bad
            if (it->second == -1) continue;

            if (std::abs(it->first * nonzero.value() + candRowVal) <=
                options->small_matrix_value) {
              it->second += 1;
            } else {
              // mark scale to be numerically bad
              it->second = -1;
              continue;
            }
          } else
            possibleScales.emplace(scale, 1);
        }

        if (skip || possibleScales.empty()) continue;

        HighsInt numCancel = 0;
        double scale = 0.0;

        for (const auto& s : possibleScales) {
          if (s.second <= 1) continue;
          if (s.second > numCancel ||
              (s.second == numCancel && std::abs(s.first) < std::abs(scale))) {
            scale = s.first;
            numCancel = s.second;
          }
        }

        assert(scale != 0.0 || numCancel == 0);

        // cancels at least one nonzero if the scale cancels more than there is
        // fillin
        if (numCancel > 1) sparsifyRows.emplace_back(candRow, scale);
      }
    }

    if (sparsifyRows.empty()) continue;

    postsolve_stack.equalityRowAdditions(eqrow, getStoredRow(), sparsifyRows);
    double rhs = model->row_lower_[eqrow];
    for (const auto& sparsifyRow : sparsifyRows) {
      HighsInt row = sparsifyRow.index;
      double scale = sparsifyRow.value;

      if (model->row_lower_[row] != -kHighsInf)
        model->row_lower_[row] += scale * rhs;

      if (model->row_upper_[row] != kHighsInf)
        model->row_upper_[row] += scale * rhs;

      for (HighsInt pos : rowpositions)
        addToMatrix(row, Acol[pos], scale * Avalue[pos]);

      reinsertEquation(row);
    }

    HPRESOLVE_CHECKED_CALL(checkLimits(postsolve_stack));
    HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolve_stack));
    HPRESOLVE_CHECKED_CALL(removeDoubletonEquations(postsolve_stack));
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
    if (model->col_hash_.name2index.size() !=
        static_cast<size_t>(model->num_col_))
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
    if (model->row_hash_.name2index.size() !=
        static_cast<size_t>(model->num_row_))
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
