#include "presolve/HPresolve.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>

#include "Highs.h"
#include "io/HighsIO.h"
#include "lp_data/HConst.h"
#include "lp_data/HStruct.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsSolution.h"
#include "mip/HighsCliqueTable.h"
#include "mip/HighsImplications.h"
#include "mip/HighsMipSolverData.h"
#include "presolve/HighsPostsolveStack.h"
#include "test/DevKkt.h"
#include "util/HighsCDouble.h"
#include "util/HighsIntegers.h"
#include "util/HighsLinearSumBounds.h"
#include "util/HighsSplay.h"
#include "util/HighsUtils.h"

#define ENABLE_SPARSIFY_FOR_LP 0

#define HPRESOLVE_CHECKED_CALL(presolveCall)                          \
  do {                                                                \
    HPresolve::Result __result = presolveCall;                        \
    if (__result != presolve::HPresolve::Result::Ok) return __result; \
  } while (0)

namespace presolve {

#ifndef NDEBUG
void HPresolve::debugPrintRow(HighsPostsolveStack& postSolveStack,
                              HighsInt row) {
  printf("(row %" HIGHSINT_FORMAT ") %.15g (impl: %.15g) <= ",
         postSolveStack.getOrigRowIndex(row), model->rowLower_[row],
         impliedRowBounds.getSumLower(row));

  for (const HighsSliceNonzero& nonzero : getSortedRowVector(row)) {
    // for (HighsInt rowiter = rowhead[row]; rowiter != -1; rowiter =
    // ARnext[rowiter]) {
    char colchar = model->integrality_[nonzero.index()] == HighsVarType::INTEGER
                       ? 'y'
                       : 'x';
    char signchar = nonzero.value() < 0 ? '-' : '+';
    printf("%c%g %c%" HIGHSINT_FORMAT " ", signchar, std::abs(nonzero.value()),
           colchar, postSolveStack.getOrigColIndex(nonzero.index()));
  }

  printf("<= %.15g (impl: %.15g)\n", model->rowUpper_[row],
         impliedRowBounds.getSumUpper(row));
}
#endif

void HPresolve::setInput(HighsLp& model_, const HighsOptions& options_) {
  model = &model_;
  options = &options_;

  colLowerSource.resize(model->numCol_, -1);
  colUpperSource.resize(model->numCol_, -1);
  implColLower.resize(model->numCol_, -HIGHS_CONST_INF);
  implColUpper.resize(model->numCol_, HIGHS_CONST_INF);

  rowDualLower.resize(model->numRow_, -HIGHS_CONST_INF);
  rowDualUpper.resize(model->numRow_, HIGHS_CONST_INF);
  implRowDualLower.resize(model->numRow_, -HIGHS_CONST_INF);
  implRowDualUpper.resize(model->numRow_, HIGHS_CONST_INF);
  rowDualUpperSource.resize(model->numRow_, -1);
  rowDualLowerSource.resize(model->numRow_, -1);

  for (HighsInt i = 0; i != model->numRow_; ++i) {
    if (model->rowLower_[i] == -HIGHS_CONST_INF) rowDualLower[i] = 0;
    if (model->rowUpper_[i] == HIGHS_CONST_INF) rowDualUpper[i] = 0;
  }

  if (mipsolver == nullptr)
    model->integrality_.assign(model->numCol_, HighsVarType::CONTINUOUS);

  if (model_.orientation_ == MatrixOrientation::ROWWISE)
    fromCSR(model->Avalue_, model->Aindex_, model->Astart_);
  else
    fromCSC(model->Avalue_, model->Aindex_, model->Astart_);

  // initialize everything as changed, but do not add all indices
  // since the first thing presolve will do is a scan for easy reductions
  // of each row and column and set the flag of processed columns to false
  // from then on they are added to the vector whenever there are changes
  changedRowFlag.resize(model->numRow_, true);
  rowDeleted.resize(model->numRow_, false);
  changedRowIndices.reserve(model->numRow_);
  changedColFlag.resize(model->numCol_, true);
  colDeleted.resize(model->numCol_, false);
  changedColIndices.reserve(model->numCol_);
  numDeletedCols = 0;
  numDeletedRows = 0;
  reductionLimit = std::numeric_limits<size_t>::max();
}

// for MIP presolve
void HPresolve::setInput(HighsMipSolver& mipsolver) {
  this->mipsolver = &mipsolver;

  probingContingent = 1000;
  probingNumDelCol = 0;
  numProbed = 0;

  if (mipsolver.model_ != &mipsolver.mipdata_->presolvedModel) {
    mipsolver.mipdata_->presolvedModel = *mipsolver.model_;
    mipsolver.model_ = &mipsolver.mipdata_->presolvedModel;
  } else {
    mipsolver.mipdata_->presolvedModel.colLower_ =
        mipsolver.mipdata_->domain.colLower_;
    mipsolver.mipdata_->presolvedModel.colUpper_ =
        mipsolver.mipdata_->domain.colUpper_;
  }

  setInput(mipsolver.mipdata_->presolvedModel, *mipsolver.options_mip_);
}

bool HPresolve::rowCoefficientsIntegral(HighsInt row, double scale) const {
  for (const HighsSliceNonzero& nz : getRowVector(row)) {
    double val = nz.value() * scale;
    if (std::abs(val - std::round(val)) > options->mip_epsilon) return false;
  }

  return true;
}

bool HPresolve::isLowerImplied(HighsInt col) const {
  return (model->colLower_[col] == -HIGHS_CONST_INF ||
          implColLower[col] >=
              model->colLower_[col] - options->primal_feasibility_tolerance);
}

bool HPresolve::isUpperImplied(HighsInt col) const {
  return (model->colUpper_[col] == HIGHS_CONST_INF ||
          implColUpper[col] <=
              model->colUpper_[col] + options->primal_feasibility_tolerance);
}

bool HPresolve::isImpliedFree(HighsInt col) const {
  return (model->colLower_[col] == -HIGHS_CONST_INF ||
          implColLower[col] >=
              model->colLower_[col] - options->primal_feasibility_tolerance) &&
         (model->colUpper_[col] == HIGHS_CONST_INF ||
          implColUpper[col] <=
              model->colUpper_[col] + options->primal_feasibility_tolerance);
}

bool HPresolve::isDualImpliedFree(HighsInt row) const {
  return model->rowLower_[row] == model->rowUpper_[row] ||
         (model->rowUpper_[row] != HIGHS_CONST_INF &&
          implRowDualLower[row] >= -options->dual_feasibility_tolerance) ||
         (model->rowLower_[row] != -HIGHS_CONST_INF &&
          implRowDualUpper[row] <= options->dual_feasibility_tolerance);
}

bool HPresolve::isImpliedIntegral(HighsInt col) {
  bool runDualDetection = true;

  assert(model->integrality_[col] == HighsVarType::INTEGER);

  for (const HighsSliceNonzero& nz : getColumnVector(col)) {
    // if not all other columns are integer, skip row and also do not try the
    // dual detection in the second loop as it must hold for all rows
    if (rowsizeInteger[nz.index()] < rowsize[nz.index()]) {
      runDualDetection = false;
      continue;
    }

    double rowLower =
        implRowDualLower[nz.index()] > options->dual_feasibility_tolerance
            ? model->rowUpper_[nz.index()]
            : model->rowLower_[nz.index()];

    double rowUpper =
        implRowDualUpper[nz.index()] < -options->dual_feasibility_tolerance
            ? model->rowLower_[nz.index()]
            : model->rowUpper_[nz.index()];

    if (rowUpper == rowLower) {
      // if there is an equation the dual detection does not need to be tried
      runDualDetection = false;
      double scale = 1.0 / nz.value();
      if (!rowCoefficientsIntegral(nz.index(), scale)) continue;

      double rhs = model->rowLower_[nz.index()] * scale;

      if (std::abs(rhs - std::round(rhs)) >
          options->mip_feasibility_tolerance) {
        // todo infeasible
      }

      return true;
    }
  }

  if (!runDualDetection) return false;

  for (const HighsSliceNonzero& nz : getColumnVector(col)) {
    double scale = 1.0 / nz.value();
    if (!rowCoefficientsIntegral(nz.index(), scale)) return false;
    if (model->rowUpper_[nz.index()] != HIGHS_CONST_INF) {
      double rUpper =
          std::abs(nz.value()) *
          std::floor(model->rowUpper_[nz.index()] * std::abs(scale) +
                     options->mip_feasibility_tolerance);
      if (std::abs(model->rowUpper_[nz.index()] - rUpper) >
          options->mip_epsilon) {
        model->rowUpper_[nz.index()] = rUpper;
        markChangedRow(nz.index());
      }
    } else {
      assert(model->rowLower_[nz.index()] != -HIGHS_CONST_INF);
      double rLower = std::abs(nz.value()) *
                      std::ceil(model->rowUpper_[nz.index()] * std::abs(scale) -
                                options->mip_feasibility_tolerance);
      if (std::abs(model->rowLower_[nz.index()] - rLower) >
          options->mip_epsilon) {
        model->rowUpper_[nz.index()] = rLower;
        markChangedRow(nz.index());
      }
    }
  }

  return true;
}

bool HPresolve::isImpliedInteger(HighsInt col) {
  bool runDualDetection = true;

  assert(model->integrality_[col] == HighsVarType::CONTINUOUS);

  for (const HighsSliceNonzero& nz : getColumnVector(col)) {
    // if not all other columns are integer, skip row and also do not try the
    // dual detection in the second loop as it must hold for all rows
    if (rowsizeInteger[nz.index()] + rowsizeImplInt[nz.index()] <
        rowsize[nz.index()] - 1) {
      runDualDetection = false;
      continue;
    }

    double rowLower =
        implRowDualLower[nz.index()] > options->dual_feasibility_tolerance
            ? model->rowUpper_[nz.index()]
            : model->rowLower_[nz.index()];

    double rowUpper =
        implRowDualUpper[nz.index()] < -options->dual_feasibility_tolerance
            ? model->rowLower_[nz.index()]
            : model->rowUpper_[nz.index()];

    if (rowUpper == rowLower) {
      // if there is an equation the dual detection does not need to be tried
      runDualDetection = false;
      double scale = 1.0 / nz.value();
      double rhs = model->rowLower_[nz.index()] * scale;

      if (std::abs(rhs - std::round(rhs)) >
          options->mip_feasibility_tolerance) {
        continue;
      }

      if (!rowCoefficientsIntegral(nz.index(), scale)) continue;

      return true;
    }
  }

  if (!runDualDetection) return false;

  if ((model->colLower_[col] != -HIGHS_CONST_INF &&
       std::abs(std::round(model->colLower_[col]) - model->colLower_[col]) >
           options->mip_epsilon) ||
      (model->colUpper_[col] != -HIGHS_CONST_INF &&
       std::abs(std::round(model->colUpper_[col]) - model->colUpper_[col]) >
           options->mip_epsilon))
    return false;

  for (const HighsSliceNonzero& nz : getColumnVector(col)) {
    double scale = 1.0 / nz.value();
    if (model->rowUpper_[nz.index()] != HIGHS_CONST_INF) {
      double rhs = model->rowUpper_[nz.index()];
      if (std::abs(rhs - std::round(rhs)) > options->mip_feasibility_tolerance)
        return false;
    }

    if (model->rowLower_[nz.index()] != -HIGHS_CONST_INF) {
      double rhs = model->rowLower_[nz.index()];
      if (std::abs(rhs - std::round(rhs)) > options->mip_feasibility_tolerance)
        return false;
    }

    if (!rowCoefficientsIntegral(nz.index(), scale)) return false;
  }

  return true;
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
  if (model->integrality_[Acol[pos]] == HighsVarType::INTEGER)
    ++rowsizeInteger[Arow[pos]];
  else if (model->integrality_[Acol[pos]] == HighsVarType::IMPLICIT_INTEGER)
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
    if (colsize[Acol[pos]] <= 1)
      singletonColumns.push_back(Acol[pos]);
    else
      markChangedCol(Acol[pos]);

    impliedDualRowBounds.remove(Acol[pos], Arow[pos], Avalue[pos]);
    if (colUpperSource[Acol[pos]] == Arow[pos])
      changeImplColUpper(Acol[pos], HIGHS_CONST_INF, -1);

    if (colLowerSource[Acol[pos]] == Arow[pos])
      changeImplColLower(Acol[pos], -HIGHS_CONST_INF, -1);
  }

  auto get_row_left = [&](HighsInt pos) -> HighsInt& { return ARleft[pos]; };
  auto get_row_right = [&](HighsInt pos) -> HighsInt& { return ARright[pos]; };
  auto get_row_key = [&](HighsInt pos) { return Acol[pos]; };
  highs_splay_unlink(pos, rowroot[Arow[pos]], get_row_left, get_row_right,
                     get_row_key);
  --rowsize[Arow[pos]];
  if (model->integrality_[Acol[pos]] == HighsVarType::INTEGER)
    --rowsizeInteger[Arow[pos]];
  else if (model->integrality_[Acol[pos]] == HighsVarType::IMPLICIT_INTEGER)
    --rowsizeImplInt[Arow[pos]];

  if (!rowDeleted[Arow[pos]]) {
    if (rowsize[Arow[pos]] <= 1)
      singletonRows.push_back(Arow[pos]);
    else
      markChangedRow(Arow[pos]);
    impliedRowBounds.remove(Arow[pos], Acol[pos], Avalue[pos]);

    if (rowDualUpperSource[Arow[pos]] == Acol[pos])
      changeImplRowDualUpper(Arow[pos], HIGHS_CONST_INF, -1);

    if (rowDualLowerSource[Arow[pos]] == Acol[pos])
      changeImplRowDualLower(Arow[pos], -HIGHS_CONST_INF, -1);
  }

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

void HPresolve::updateRowDualImpliedBounds(HighsInt row, HighsInt col,
                                           double val) {
  // propagate implied row dual bound bound
  // if the column has an infinite lower bound the reduced cost cannot be
  // positive, i.e. the column corresponds to a <= constraint in the dual with
  // right hand side -cost furthermore, we can ignore strictly redundant primal
  // column bounds and treat them as if they are infinite
  double dualRowUpper =
      (model->colLower_[col] == -HIGHS_CONST_INF) ||
              (implColLower[col] >
               model->colLower_[col] + options->primal_feasibility_tolerance)
          ? -model->colCost_[col]
          : HIGHS_CONST_INF;

  double dualRowLower =
      (model->colUpper_[col] == HIGHS_CONST_INF) ||
              (implColUpper[col] <
               model->colUpper_[col] - options->primal_feasibility_tolerance)
          ? -model->colCost_[col]
          : -HIGHS_CONST_INF;

  if (dualRowUpper != HIGHS_CONST_INF) {
    // get minimal value of other row duals in the column
    double residualMinAct =
        impliedDualRowBounds.getResidualSumLowerOrig(col, row, val);
    if (residualMinAct != -HIGHS_CONST_INF) {
      double impliedBound =
          double((HighsCDouble(dualRowUpper) - residualMinAct) / val);

      if (std::abs(impliedBound) * HIGHS_CONST_TINY <=
          options->dual_feasibility_tolerance) {
        if (val > 0) {
          // only tighten bound if it is tighter by a wide enough margin
          if (impliedBound < implRowDualUpper[row] -
                                 1000 * options->dual_feasibility_tolerance)
            changeImplRowDualUpper(row, impliedBound, col);
        } else {
          if (impliedBound > implRowDualLower[row] +
                                 1000 * options->primal_feasibility_tolerance)
            changeImplRowDualLower(row, impliedBound, col);
        }
      }
    }
  }

  if (dualRowLower != -HIGHS_CONST_INF) {
    // get maximal value of other columns in the row
    double residualMaxAct =
        impliedDualRowBounds.getResidualSumUpperOrig(col, row, val);
    if (residualMaxAct != HIGHS_CONST_INF) {
      double impliedBound =
          double((HighsCDouble(dualRowLower) - residualMaxAct) / val);

      if (std::abs(impliedBound) * HIGHS_CONST_TINY <=
          options->dual_feasibility_tolerance) {
        if (val > 0) {
          // only tighten bound if it is tighter by a wide enough margin
          if (impliedBound > implRowDualLower[row] +
                                 1000 * options->primal_feasibility_tolerance)
            changeImplRowDualLower(row, impliedBound, col);
        } else {
          if (impliedBound < implRowDualUpper[row] -
                                 1000 * options->dual_feasibility_tolerance)
            changeImplRowDualUpper(row, impliedBound, col);
        }
      }
    }
  }
}

void HPresolve::updateColImpliedBounds(HighsInt row, HighsInt col, double val) {
  // propagate implied column bound upper bound if row has an upper bound
  double rowUpper = implRowDualUpper[row] < -options->dual_feasibility_tolerance
                        ? model->rowLower_[row]
                        : model->rowUpper_[row];
  double rowLower = implRowDualLower[row] > options->dual_feasibility_tolerance
                        ? model->rowUpper_[row]
                        : model->rowLower_[row];

  assert(rowLower != HIGHS_CONST_INF);
  assert(rowUpper != -HIGHS_CONST_INF);

  if (rowUpper != HIGHS_CONST_INF) {
    // get minimal value of other columns in the row
    double residualMinAct =
        impliedRowBounds.getResidualSumLowerOrig(row, col, val);
    if (residualMinAct != -HIGHS_CONST_INF) {
      double impliedBound =
          double((HighsCDouble(rowUpper) - residualMinAct) / val);

      if (std::abs(impliedBound) * HIGHS_CONST_TINY <=
          options->primal_feasibility_tolerance) {
        if (val > 0) {
          // bound is an upper bound
          // check if we may round the bound due to integrality restrictions
          if (mipsolver != nullptr) {
            if (model->integrality_[col] != HighsVarType::CONTINUOUS) {
              double roundedBound =
                  std::floor(impliedBound + options->mip_feasibility_tolerance);

              if (roundedBound < model->colUpper_[col])
                changeColUpper(col, roundedBound);
            }

            if (mipsolver->mipdata_->postSolveStack.getOrigRowIndex(row) >=
                mipsolver->orig_model_->numRow_) {
              if (impliedBound <
                  model->colUpper_[col] -
                      1000 * options->primal_feasibility_tolerance)
                changeColUpper(col, impliedBound);

              impliedBound = HIGHS_CONST_INF;
            }
          }

          // only tighten bound if it is tighter by a wide enough margin
          if (impliedBound <
              implColUpper[col] - 1000 * options->primal_feasibility_tolerance)
            changeImplColUpper(col, impliedBound, row);
        } else {
          // bound is a lower bound
          // check if we may round the bound due to integrality restrictions
          if (mipsolver != nullptr) {
            if (model->integrality_[col] != HighsVarType::CONTINUOUS) {
              double roundedBound =
                  std::ceil(impliedBound - options->mip_feasibility_tolerance);

              if (roundedBound > model->colLower_[col])
                changeColLower(col, roundedBound);
            }

            // do not use the implied bound if this a not a model row, since the
            // row can be removed and should not be used, e.g., to identify a
            // column as implied free
            if (mipsolver->mipdata_->postSolveStack.getOrigRowIndex(row) >=
                mipsolver->orig_model_->numRow_) {
              if (impliedBound >
                  model->colLower_[col] +
                      1000 * options->primal_feasibility_tolerance)
                changeColLower(col, impliedBound);

              impliedBound = -HIGHS_CONST_INF;
            }
          }

          // only tighten bound if it is tighter by a wide enough margin
          if (impliedBound >
              implColLower[col] + 1000 * options->primal_feasibility_tolerance)
            changeImplColLower(col, impliedBound, row);
        }
      }
    }
  }

  if (rowLower != -HIGHS_CONST_INF) {
    // get maximal value of other columns in the row
    double residualMaxAct =
        impliedRowBounds.getResidualSumUpperOrig(row, col, val);
    if (residualMaxAct != HIGHS_CONST_INF) {
      double impliedBound =
          double((HighsCDouble(rowLower) - residualMaxAct) / val);

      if (std::abs(impliedBound) * HIGHS_CONST_TINY <=
          options->primal_feasibility_tolerance) {
        if (val > 0) {
          // bound is a lower bound
          // check if we may round the bound due to integrality restrictions
          if (mipsolver != nullptr) {
            if (model->integrality_[col] != HighsVarType::CONTINUOUS) {
              double roundedBound =
                  std::ceil(impliedBound - options->mip_feasibility_tolerance);

              // change bounds of integers immediately
              if (roundedBound > model->colLower_[col])
                changeColLower(col, roundedBound);
            }

            if (mipsolver->mipdata_->postSolveStack.getOrigRowIndex(row) >=
                mipsolver->orig_model_->numRow_) {
              if (impliedBound >
                  model->colLower_[col] +
                      1000 * options->primal_feasibility_tolerance)
                changeColLower(col, impliedBound);

              impliedBound = -HIGHS_CONST_INF;
            }
          }

          // only tighten bound if it is tighter by a wide enough margin
          if (impliedBound >
              implColLower[col] + 1000 * options->primal_feasibility_tolerance)
            changeImplColLower(col, impliedBound, row);
        } else {
          // bound is an upper bound
          // check if we may round the bound due to integrality restrictions
          if (mipsolver != nullptr) {
            if (model->integrality_[col] != HighsVarType::CONTINUOUS) {
              double roundedBound =
                  std::floor(impliedBound + options->mip_feasibility_tolerance);

              // change bounds of integers immediately
              if (roundedBound < model->colUpper_[col])
                changeColUpper(col, roundedBound);
            }

            if (mipsolver->mipdata_->postSolveStack.getOrigRowIndex(row) >=
                mipsolver->orig_model_->numRow_) {
              if (impliedBound <
                  model->colUpper_[col] -
                      1000 * options->primal_feasibility_tolerance)
                changeColUpper(col, impliedBound);

              impliedBound = HIGHS_CONST_INF;
            }
          }

          // only tighten bound if it is tighter by a wide enough margin
          if (impliedBound <
              implColUpper[col] - 1000 * options->primal_feasibility_tolerance)
            changeImplColUpper(col, impliedBound, row);
        }
      }
    }
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

void HPresolve::shrinkProblem(HighsPostsolveStack& postSolveStack) {
  HighsInt oldNumCol = model->numCol_;
  model->numCol_ = 0;
  std::vector<HighsInt> newColIndex(oldNumCol);
  for (HighsInt i = 0; i != oldNumCol; ++i) {
    if (colDeleted[i])
      newColIndex[i] = -1;
    else {
      newColIndex[i] = model->numCol_++;
      model->colCost_[newColIndex[i]] = model->colCost_[i];
      model->colLower_[newColIndex[i]] = model->colLower_[i];
      model->colUpper_[newColIndex[i]] = model->colUpper_[i];
      assert(!std::isnan(model->colLower_[newColIndex[i]]));
      assert(!std::isnan(model->colUpper_[newColIndex[i]]));
      model->integrality_[newColIndex[i]] = model->integrality_[i];
      implColLower[newColIndex[i]] = implColLower[i];
      implColUpper[newColIndex[i]] = implColUpper[i];
      colLowerSource[newColIndex[i]] = colLowerSource[i];
      colUpperSource[newColIndex[i]] = colUpperSource[i];
      colhead[newColIndex[i]] = colhead[i];
      colsize[newColIndex[i]] = colsize[i];
      if ((HighsInt)model->col_names_.size() > 0)
        model->col_names_[newColIndex[i]] = std::move(model->col_names_[i]);
      changedColFlag[newColIndex[i]] = changedColFlag[i];
    }
  }
  colDeleted.assign(model->numCol_, false);
  model->colCost_.resize(model->numCol_);
  model->colLower_.resize(model->numCol_);
  model->colUpper_.resize(model->numCol_);
  model->integrality_.resize(model->numCol_);
  implColLower.resize(model->numCol_);
  implColUpper.resize(model->numCol_);
  colLowerSource.resize(model->numCol_);
  colUpperSource.resize(model->numCol_);
  colhead.resize(model->numCol_);
  colsize.resize(model->numCol_);
  if ((HighsInt)model->col_names_.size() > 0)
    model->col_names_.resize(model->numCol_);
  changedColFlag.resize(model->numCol_);
  numDeletedCols = 0;
  HighsInt oldNumRow = model->numRow_;
  model->numRow_ = 0;
  std::vector<HighsInt> newRowIndex(oldNumRow);
  for (HighsInt i = 0; i != oldNumRow; ++i) {
    if (rowDeleted[i])
      newRowIndex[i] = -1;
    else {
      newRowIndex[i] = model->numRow_++;
      model->rowLower_[newRowIndex[i]] = model->rowLower_[i];
      model->rowUpper_[newRowIndex[i]] = model->rowUpper_[i];
      assert(!std::isnan(model->rowLower_[newRowIndex[i]]));
      assert(!std::isnan(model->rowUpper_[newRowIndex[i]]));
      rowDualLower[newRowIndex[i]] = rowDualLower[i];
      rowDualUpper[newRowIndex[i]] = rowDualUpper[i];
      implRowDualLower[newRowIndex[i]] = implRowDualLower[i];
      implRowDualUpper[newRowIndex[i]] = implRowDualUpper[i];
      rowDualLowerSource[newRowIndex[i]] = rowDualLowerSource[i];
      rowDualUpperSource[newRowIndex[i]] = rowDualUpperSource[i];
      rowroot[newRowIndex[i]] = rowroot[i];
      rowsize[newRowIndex[i]] = rowsize[i];
      rowsizeInteger[newRowIndex[i]] = rowsizeInteger[i];
      rowsizeImplInt[newRowIndex[i]] = rowsizeImplInt[i];
      if ((HighsInt)model->row_names_.size() > 0)
        model->row_names_[newRowIndex[i]] = std::move(model->row_names_[i]);
      changedRowFlag[newRowIndex[i]] = changedRowFlag[i];
    }
  }

  for (HighsInt i = 0; i != model->numCol_; ++i) {
    if (colLowerSource[i] != -1)
      colLowerSource[i] = newRowIndex[colLowerSource[i]];
    if (colUpperSource[i] != -1)
      colUpperSource[i] = newRowIndex[colUpperSource[i]];
  }

  for (HighsInt i = 0; i != model->numRow_; ++i) {
    if (rowDualLowerSource[i] != -1)
      rowDualLowerSource[i] = newColIndex[rowDualLowerSource[i]];
    if (rowDualUpperSource[i] != -1)
      rowDualUpperSource[i] = newColIndex[rowDualUpperSource[i]];
  }

  rowDeleted.assign(model->numRow_, false);
  model->rowLower_.resize(model->numRow_);
  model->rowUpper_.resize(model->numRow_);
  rowDualLower.resize(model->numRow_);
  rowDualUpper.resize(model->numRow_);
  implRowDualLower.resize(model->numRow_);
  implRowDualUpper.resize(model->numRow_);
  rowDualLowerSource.resize(model->numRow_);
  rowDualUpperSource.resize(model->numRow_);
  rowroot.resize(model->numRow_);
  rowsize.resize(model->numRow_);
  rowsizeInteger.resize(model->numRow_);
  rowsizeImplInt.resize(model->numRow_);
  if ((HighsInt)model->row_names_.size() > 0)
    model->row_names_.resize(model->numRow_);
  changedRowFlag.resize(model->numRow_);

  numDeletedRows = 0;
  postSolveStack.compressIndexMaps(newRowIndex, newColIndex);
  impliedRowBounds.shrink(newRowIndex, model->numRow_);
  impliedDualRowBounds.shrink(newColIndex, model->numCol_);

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
  eqiters.assign(model->numRow_, equations.end());
  for (HighsInt i = 0; i != model->numRow_; ++i) {
    if (model->rowLower_[i] == model->rowUpper_[i])
      eqiters[i] = equations.emplace(rowsize[i], i).first;
  }

  if (mipsolver != nullptr) {
    mipsolver->mipdata_->rowMatrixSet = false;
    mipsolver->mipdata_->domain = HighsDomain(*mipsolver);
    mipsolver->mipdata_->cliquetable.rebuild(
        model->numCol_, mipsolver->mipdata_->domain, newColIndex, newRowIndex);
    mipsolver->mipdata_->implications.rebuild(model->numCol_, newColIndex,
                                              newRowIndex);
    mipsolver->mipdata_->cutpool = HighsCutPool(
        mipsolver->model_->numCol_, mipsolver->options_mip_->mip_pool_age_limit,
        mipsolver->options_mip_->mip_pool_soft_limit);
  }
}

HPresolve::Result HPresolve::runProbing(HighsPostsolveStack& postSolveStack) {
  if (numDeletedCols + numDeletedRows != 0) shrinkProblem(postSolveStack);

  toCSC(model->Avalue_, model->Aindex_, model->Astart_);
  fromCSC(model->Avalue_, model->Aindex_, model->Astart_);

  // first tighten all bounds if they have an implied bound that is tighter
  // thatn their column bound before probing this is not done for continuous
  // columns since it may allow stronger dual presolve and more aggregations
  double hugeBound = options->primal_feasibility_tolerance / HIGHS_CONST_TINY;
  for (HighsInt i = 0; i != model->numCol_; ++i) {
    if (model->colLower_[i] >= implColLower[i] &&
        model->colUpper_[i] <= implColUpper[i])
      continue;

    if (std::abs(implColLower[i]) <= hugeBound) {
      double newLb = implColLower[i];
      if (newLb > model->colLower_[i]) changeColLower(i, newLb);
    }

    if (std::abs(implColUpper[i]) <= hugeBound) {
      double newUb = implColUpper[i];
      if (newUb < model->colUpper_[i]) changeColUpper(i, newUb);
    }
  }

  HighsInt oldNumProbed = numProbed;

  mipsolver->mipdata_->setupDomainPropagation();
  HighsDomain& domain = mipsolver->mipdata_->domain;

  domain.propagate();
  if (domain.infeasible()) return Result::PrimalInfeasible;
  HighsCliqueTable& cliquetable = mipsolver->mipdata_->cliquetable;
  HighsImplications& implications = mipsolver->mipdata_->implications;
  bool firstCall = !mipsolver->mipdata_->cliquesExtracted;
  mipsolver->mipdata_->cliquesExtracted = true;

  // extract cliques that are part of the formulation every time before probing
  // after the first call we only add cliques that directly correspond to set
  // packing constraints so that the clique merging step can extend/delete them
  if (firstCall) {
    cliquetable.extractCliques(*mipsolver);
    if (domain.infeasible()) return Result::PrimalInfeasible;

    // during presolve we keep the objective upper bound without the current
    // offset so we need to update it

    if (mipsolver->mipdata_->upper_limit != HIGHS_CONST_INF) {
      double tmpLimit = mipsolver->mipdata_->upper_limit;
      mipsolver->mipdata_->upper_limit = tmpLimit - model->offset_;
      cliquetable.extractObjCliques(*mipsolver);
      mipsolver->mipdata_->upper_limit = tmpLimit;

      if (domain.infeasible()) return Result::PrimalInfeasible;
    }

    domain.propagate();
    if (domain.infeasible()) return Result::PrimalInfeasible;
  }

  cliquetable.cleanupFixed(domain);
  if (domain.infeasible()) return Result::PrimalInfeasible;

  // store binary variables in vector with their number of implications on
  // other binaries
  std::vector<std::tuple<int64_t, HighsInt, HighsInt, HighsInt>> binaries;
  binaries.reserve(model->numCol_);
  HighsRandom random(options->highs_random_seed);
  for (HighsInt i = 0; i != model->numCol_; ++i) {
    if (domain.isBinary(i)) {
      HighsInt implicsUp = cliquetable.getNumImplications(i, 1);
      HighsInt implicsDown = cliquetable.getNumImplications(i, 0);
      binaries.emplace_back(-int64_t(implicsUp) * implicsDown,
                            -std::min(HighsInt{100}, implicsUp + implicsDown),
                            random.integer(), i);
    }
  }
  if (!binaries.empty()) {
    // sort variables with many implications on other binaries first
    std::sort(binaries.begin(), binaries.end());

    size_t numChangedCols = 0;
    while (domain.getChangedCols().size() != numChangedCols) {
      if (domain.isFixed(domain.getChangedCols()[numChangedCols++]))
        ++probingNumDelCol;
    }

    HighsInt numCliquesStart = cliquetable.numCliques();
    HighsInt numDelStart = probingNumDelCol;

    HighsInt numDel = probingNumDelCol - numDelStart +
                      implications.substitutions.size() +
                      cliquetable.getSubstitutions().size();

    // printf("start probing wit %" HIGHSINT_FORMAT " cliques\n");
    for (const std::tuple<int64_t, HighsInt, HighsInt, HighsInt>& binvar :
         binaries) {
      HighsInt i = std::get<3>(binvar);

      if (cliquetable.getSubstitution(i) != nullptr) continue;

      if (domain.isBinary(i)) {
        // break in case of too many new implications to not spent ages in
        // probing
        if (cliquetable.numCliques() - numCliquesStart >
            std::max(HighsInt{1000000}, numNonzeros()))
          break;

        // when a large percentage of columns have been deleted, stop this round
        // of probing
        // if (numDel > std::max(model->numCol_ * 0.2, 1000.)) break;
        if (numDel > (model->numRow_ + model->numCol_) * 0.05) break;
        if (probingContingent - numProbed < 0) break;

        HighsInt numBoundChgs = 0;

        if (!implications.runProbing(i, numBoundChgs)) continue;
        probingContingent += numBoundChgs;

        while (domain.getChangedCols().size() != numChangedCols) {
          if (domain.isFixed(domain.getChangedCols()[numChangedCols++]))
            ++probingNumDelCol;
        }
        HighsInt newNumDel = probingNumDelCol - numDelStart +
                             implications.substitutions.size() +
                             cliquetable.getSubstitutions().size();

        if (newNumDel > numDel) {
          probingContingent += numDel;
          numDel = newNumDel;
        }

        ++numProbed;

        // printf("nprobed: %" HIGHSINT_FORMAT ", numCliques: %" HIGHSINT_FORMAT
        // "\n", nprobed,
        //       cliquetable.numCliques());
        if (domain.infeasible()) {
          return Result::PrimalInfeasible;
        }
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

    // add nonzeros from clique lifting before removign fixed variables, since
    // this might lead to stronger constraint sides
    std::vector<std::pair<HighsInt, HighsCliqueTable::CliqueVar>>&
        extensionvars = cliquetable.getCliqueExtensions();
    HighsInt addednnz = extensionvars.size();
    for (std::pair<HighsInt, HighsCliqueTable::CliqueVar> cliqueextension :
         extensionvars) {
      if (rowDeleted[cliqueextension.first]) {
        --addednnz;
        continue;
      }
      assert(findNonzero(cliqueextension.first, cliqueextension.second.col) ==
             -1);
      double val;
      if (cliqueextension.second.val == 0) {
        model->rowLower_[cliqueextension.first] -= 1;
        model->rowUpper_[cliqueextension.first] -= 1;
        val = -1.0;
      } else
        val = 1.0;
      addToMatrix(cliqueextension.first, cliqueextension.second.col, val);
    }
    extensionvars.clear();

    // now remove fixed columns and tighten domains
    for (HighsInt i = 0; i != model->numCol_; ++i) {
      if (colDeleted[i]) continue;
      if (model->colLower_[i] < domain.colLower_[i])
        changeColLower(i, domain.colLower_[i]);
      if (model->colUpper_[i] > domain.colUpper_[i])
        changeColUpper(i, domain.colUpper_[i]);
      if (domain.isFixed(i)) {
        postSolveStack.removedFixedCol(i, model->colLower_[i], 0.0,
                                       HighsEmptySlice());
        removeFixedCol(i);
      }
      HPRESOLVE_CHECKED_CALL(checkLimits(postSolveStack));
    }

    // finally apply substitutions
    HPRESOLVE_CHECKED_CALL(applyConflictGraphSubstitutions(postSolveStack));

    highsLogUser(options->log_options, HighsLogType::INFO,
                 "%" HIGHSINT_FORMAT " probing evaluations: %" HIGHSINT_FORMAT
                 " deleted rows, %" HIGHSINT_FORMAT
                 " deleted "
                 "columns, %" HIGHSINT_FORMAT " lifted nonzeros\n",
                 numProbed - oldNumProbed, numDeletedRows, numDeletedCols,
                 addednnz);
  }

  return checkLimits(postSolveStack);
}

void HPresolve::addToMatrix(HighsInt row, HighsInt col, double val) {
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
  } else {
    double sum = Avalue[pos] + val;
    if (std::abs(sum) <= options->small_matrix_value) {
      unlink(pos);
    } else {
      // remove implied bounds on the row dual that where implied by this
      // columns dual constraint
      if (rowDualUpperSource[row] == col)
        changeImplRowDualUpper(row, HIGHS_CONST_INF, -1);

      if (rowDualLowerSource[row] == col)
        changeImplRowDualLower(row, -HIGHS_CONST_INF, -1);

      // remove implied bounds on the column that where implied by this row
      if (colUpperSource[col] == row)
        changeImplColUpper(col, HIGHS_CONST_INF, -1);

      if (colLowerSource[col] == row)
        changeImplColLower(col, -HIGHS_CONST_INF, -1);

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
  if (model->rowLower_[row] == model->rowUpper_[row] &&
      eqiters[row] != equations.end()) {
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
  if (model->integrality_[col] != HighsVarType::CONTINUOUS) {
    newUpper = std::floor(newUpper + options->mip_feasibility_tolerance);
    if (newUpper == model->colUpper_[col]) return;
  }

  double oldUpper = model->colUpper_[col];
  model->colUpper_[col] = newUpper;

  for (const HighsSliceNonzero& nonzero : getColumnVector(col)) {
    impliedRowBounds.updatedVarUpper(nonzero.index(), col, nonzero.value(),
                                     oldUpper);
    markChangedRow(nonzero.index());
  }
}

void HPresolve::changeColLower(HighsInt col, double newLower) {
  if (model->integrality_[col] != HighsVarType::CONTINUOUS) {
    newLower = std::ceil(newLower - options->mip_feasibility_tolerance);
    if (newLower == model->colLower_[col]) return;
  }

  double oldLower = model->colLower_[col];
  model->colLower_[col] = newLower;
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
  if (oldImplUpper >=
          model->colUpper_[col] - options->primal_feasibility_tolerance &&
      newUpper <
          model->colUpper_[col] - options->primal_feasibility_tolerance) {
    // the dual constraint can be considered a >= constraint and was free, or a
    // <= constraint before
    markChangedCol(col);
  }
  bool newImpliedFree =
      isLowerImplied(col) &&
      oldImplUpper >
          model->colUpper_[col] + options->primal_feasibility_tolerance &&
      newUpper <= model->colUpper_[col] + options->primal_feasibility_tolerance;

  // remember the source of this lower bound, so that we can correctly identify
  // weak domination
  colUpperSource[col] = originRow;
  implColUpper[col] = newUpper;

  // if the old and the new implied bound not better than the lower bound
  // nothing
  // needs to be updated
  if (!newImpliedFree &&
      std::min(oldImplUpper, newUpper) >= model->colUpper_[col])
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
  if (oldImplLower <=
          model->colLower_[col] + options->primal_feasibility_tolerance &&
      newLower >
          model->colLower_[col] + options->primal_feasibility_tolerance) {
    // the dual constraint can additionally be considered a <= constraint and
    // was free, or a
    // >= constraint before
    markChangedCol(col);
  }
  bool newImpliedFree =
      isUpperImplied(col) &&
      oldImplLower <
          model->colLower_[col] - options->primal_feasibility_tolerance &&
      newLower >= model->colLower_[col] - options->primal_feasibility_tolerance;

  // remember the source of this lower bound, so that we can correctly identify
  // weak domination
  colLowerSource[col] = originRow;
  implColLower[col] = newLower;

  // if the old and the new implied bound not better than the lower bound
  // nothing needs to be updated
  if (!newImpliedFree &&
      std::max(oldImplLower, newLower) <= model->colLower_[col])
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

  // remember the source of this lower bound, so that we can correctly identify
  // weakdomination
  rowDualUpperSource[row] = originCol;
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

  if (oldImplLower >= -options->dual_feasibility_tolerance &&
      newLower < -options->dual_feasibility_tolerance)
    markChangedRow(row);

  bool newDualImplied =
      !isDualImpliedFree(row) &&
      oldImplLower < rowDualLower[row] - options->dual_feasibility_tolerance &&
      newLower >= rowDualLower[row] - options->dual_feasibility_tolerance;

  // remember the source of this lower bound, so that we can correctly identify
  // a weakly domination
  rowDualLowerSource[row] = originCol;
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

HPresolve::Result HPresolve::applyConflictGraphSubstitutions(
    HighsPostsolveStack& postSolveStack) {
  HighsCliqueTable& cliquetable = mipsolver->mipdata_->cliquetable;
  HighsImplications& implications = mipsolver->mipdata_->implications;
  for (const auto& substitution : implications.substitutions) {
    if (colDeleted[substitution.substcol] || colDeleted[substitution.staycol])
      continue;

    ++probingNumDelCol;

    postSolveStack.doubletonEquation(-1, substitution.substcol,
                                     substitution.staycol, 1.0,
                                     -substitution.scale, substitution.offset,
                                     model->colLower_[substitution.substcol],
                                     model->colUpper_[substitution.substcol],
                                     0.0, false, false, HighsEmptySlice());
    markColDeleted(substitution.substcol);
    substitute(substitution.substcol, substitution.staycol, substitution.offset,
               substitution.scale);
    HPRESOLVE_CHECKED_CALL(checkLimits(postSolveStack));
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

    postSolveStack.doubletonEquation(
        -1, subst.substcol, subst.replace.col, 1.0, -scale, offset,
        model->colLower_[subst.substcol], model->colUpper_[subst.substcol], 0.0,
        false, false, HighsEmptySlice());
    markColDeleted(subst.substcol);
    substitute(subst.substcol, subst.replace.col, offset, scale);
    HPRESOLVE_CHECKED_CALL(checkLimits(postSolveStack));
  }

  cliquetable.getSubstitutions().clear();

  return Result::Ok;
}

void HPresolve::storeRow(HighsInt row) {
  rowpositions.clear();

  auto rowVec = getSortedRowVector(row);
  auto rowVecEnd = rowVec.end();
  for (auto iter = rowVec.begin(); iter != rowVecEnd; ++iter)
    rowpositions.push_back(iter.position());
}

HighsTripletPositionSlice HPresolve::getStoredRow() const {
  return HighsTripletPositionSlice(Acol.data(), Avalue.data(),
                                   rowpositions.data(), rowpositions.size());
}

void HPresolve::fromCSC(const std::vector<double>& Aval,
                        const std::vector<HighsInt>& Aindex,
                        const std::vector<HighsInt>& Astart) {
  Avalue.clear();
  Acol.clear();
  Arow.clear();

  freeslots.clear();
  colhead.assign(model->numCol_, -1);
  rowroot.assign(model->numRow_, -1);
  colsize.assign(model->numCol_, 0);
  rowsize.assign(model->numRow_, 0);
  rowsizeInteger.assign(model->numRow_, 0);
  rowsizeImplInt.assign(model->numRow_, 0);

  impliedRowBounds.setNumSums(0);
  impliedDualRowBounds.setNumSums(0);
  impliedRowBounds.setBoundArrays(
      model->colLower_.data(), model->colUpper_.data(), implColLower.data(),
      implColUpper.data(), colLowerSource.data(), colUpperSource.data());
  impliedRowBounds.setNumSums(model->numRow_);
  impliedDualRowBounds.setBoundArrays(
      rowDualLower.data(), rowDualUpper.data(), implRowDualLower.data(),
      implRowDualUpper.data(), rowDualLowerSource.data(),
      rowDualUpperSource.data());
  impliedDualRowBounds.setNumSums(model->numCol_);

  HighsInt ncol = Astart.size() - 1;
  assert(ncol == int(colhead.size()));
  HighsInt nnz = Aval.size();

  Avalue = Aval;
  Acol.reserve(nnz);
  Arow.reserve(nnz);

  for (HighsInt i = 0; i != ncol; ++i) {
    HighsInt collen = Astart[i + 1] - Astart[i];
    Acol.insert(Acol.end(), collen, i);
    Arow.insert(Arow.end(), Aindex.begin() + Astart[i],
                Aindex.begin() + Astart[i + 1]);
  }

  Anext.resize(nnz);
  Aprev.resize(nnz);
  ARleft.resize(nnz);
  ARright.resize(nnz);
  for (HighsInt pos = 0; pos != nnz; ++pos) link(pos);

  if (equations.empty()) {
    eqiters.assign(model->numRow_, equations.end());
    for (HighsInt i = 0; i != model->numRow_; ++i) {
      // register equation
      if (model->rowLower_[i] == model->rowUpper_[i])
        eqiters[i] = equations.emplace(rowsize[i], i).first;
    }
  }
}

void HPresolve::fromCSR(const std::vector<double>& ARval,
                        const std::vector<HighsInt>& ARindex,
                        const std::vector<HighsInt>& ARstart) {
  Avalue.clear();
  Acol.clear();
  Arow.clear();

  freeslots.clear();
  colhead.assign(model->numCol_, -1);
  rowroot.assign(model->numRow_, -1);
  colsize.assign(model->numCol_, 0);
  rowsize.assign(model->numRow_, 0);
  rowsizeInteger.assign(model->numRow_, 0);
  rowsizeImplInt.assign(model->numRow_, 0);

  impliedRowBounds.setNumSums(0);
  impliedDualRowBounds.setNumSums(0);
  impliedRowBounds.setBoundArrays(
      model->colLower_.data(), model->colUpper_.data(), implColLower.data(),
      implColUpper.data(), colLowerSource.data(), colUpperSource.data());
  impliedRowBounds.setNumSums(model->numRow_);
  impliedDualRowBounds.setBoundArrays(
      rowDualLower.data(), rowDualUpper.data(), implRowDualLower.data(),
      implRowDualUpper.data(), rowDualLowerSource.data(),
      rowDualUpperSource.data());
  impliedDualRowBounds.setNumSums(model->numCol_);

  HighsInt nrow = ARstart.size() - 1;
  assert(nrow == int(rowroot.size()));
  HighsInt nnz = ARval.size();

  Avalue = ARval;
  Acol.reserve(nnz);
  Arow.reserve(nnz);
  //  entries.reserve(nnz);

  for (HighsInt i = 0; i != nrow; ++i) {
    HighsInt rowlen = ARstart[i + 1] - ARstart[i];
    Arow.insert(Arow.end(), rowlen, i);
    Acol.insert(Acol.end(), ARindex.begin() + ARstart[i],
                ARindex.begin() + ARstart[i + 1]);
  }

  Anext.resize(nnz);
  Aprev.resize(nnz);
  ARleft.resize(nnz);
  ARright.resize(nnz);
  for (HighsInt pos = 0; pos != nnz; ++pos) link(pos);

  if (equations.empty()) {
    eqiters.assign(nrow, equations.end());
    for (HighsInt i = 0; i != nrow; ++i) {
      // register equation
      if (model->rowLower_[i] == model->rowUpper_[i])
        eqiters[i] = equations.emplace(rowsize[i], i).first;
    }
  }
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
  assert(int(rowpositions.size()) == rowsize[row]);

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
    if (model->rowLower_[colrow] != -HIGHS_CONST_INF)
      model->rowLower_[colrow] += scale * rhs;

    if (model->rowUpper_[colrow] != HIGHS_CONST_INF)
      model->rowUpper_[colrow] += scale * rhs;

    for (HighsInt rowiter : rowpositions) {
      assert(Arow[rowiter] == row);

      if (Acol[rowiter] != col)
        addToMatrix(colrow, Acol[rowiter], scale * Avalue[rowiter]);
    }

    // check if this is an equation row and it now has a different size
    if (model->rowLower_[colrow] == model->rowUpper_[colrow] &&
        eqiters[colrow] != equations.end() &&
        eqiters[colrow]->first != rowsize[colrow]) {
      // if that is the case reinsert it into the equation set that is ordered
      // by sparsity
      equations.erase(eqiters[colrow]);
      eqiters[colrow] = equations.emplace(rowsize[colrow], colrow).first;
    }
    // printf("after substitution: ");
    // debugPrintRow(colrow);
  }

  assert(colsize[col] == 1);

  // substitute column in the objective function
  if (model->colCost_[col] != 0.0) {
    HighsCDouble objscale = model->colCost_[col] * substrowscale;
    model->offset_ = double(model->offset_ - objscale * rhs);
    assert(std::isfinite(model->offset_));
    for (HighsInt rowiter : rowpositions) {
      // printf("changing col cost to %g = %g + %g * %g\n",
      // double(model->colCost_[Acol[rowiter]] + objscale * Avalue[rowiter]),
      // model->colCost_[Acol[rowiter]], double(objscale), Avalue[rowiter]);
      model->colCost_[Acol[rowiter]] =
          double(model->colCost_[Acol[rowiter]] + objscale * Avalue[rowiter]);
      if (std::abs(model->colCost_[Acol[rowiter]]) <=
          options->small_matrix_value)
        model->colCost_[Acol[rowiter]] = 0.0;
    }
    assert(model->colCost_[col] == 0);
    model->colCost_[col] = 0.0;
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
    assert(Acol[i] >= 0 && Acol[i] < model->numCol_);
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

HPresolve::Result HPresolve::doubletonEq(HighsPostsolveStack& postSolveStack,
                                         HighsInt row) {
  assert(!rowDeleted[row]);
  assert(rowsize[row] == 2);
  assert(model->rowLower_[row] == model->rowUpper_[row]);
  // printf("doubleton equation: ");
  // debugPrintRow(row);
  HighsInt nzPos1 = rowroot[row];
  HighsInt nzPos2 = ARright[nzPos1] != -1 ? ARright[nzPos1] : ARleft[nzPos1];

  HighsInt substcol;
  HighsInt staycol;
  double substcoef;
  double staycoef;
  double rhs = model->rowUpper_[row];
  if (model->integrality_[Acol[nzPos1]] == HighsVarType::INTEGER) {
    if (model->integrality_[Acol[nzPos2]] == HighsVarType::INTEGER) {
      // both columns integer. For substitution choose smaller absolute
      // coefficient value, or sparser column if values are equal
      if (std::abs(Avalue[nzPos1]) <
          std::abs(Avalue[nzPos2]) - options->mip_epsilon) {
        substcol = Acol[nzPos1];
        staycol = Acol[nzPos2];

        substcoef = Avalue[nzPos1];
        staycoef = Avalue[nzPos2];
      } else if (std::abs(Avalue[nzPos2]) <
                 std::abs(Avalue[nzPos1]) - options->mip_epsilon) {
        substcol = Acol[nzPos2];
        staycol = Acol[nzPos1];

        substcoef = Avalue[nzPos2];
        staycoef = Avalue[nzPos1];
      } else if (colsize[Acol[nzPos1]] < colsize[Acol[nzPos2]]) {
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

      // check integrality conditions
      double roundCoef = std::round(staycoef / substcoef) * substcoef;
      if (std::abs(roundCoef - staycoef) > options->mip_epsilon)
        return Result::Ok;
      staycoef = roundCoef;
      double roundRhs = std::round(rhs / substcoef) * substcoef;
      if (std::abs(rhs - roundRhs) > options->mip_feasibility_tolerance)
        return Result::PrimalInfeasible;
      rhs = roundRhs;
    } else {
      // one col is integral, substitute the continuous one
      substcol = Acol[nzPos2];
      staycol = Acol[nzPos1];

      substcoef = Avalue[nzPos2];
      staycoef = Avalue[nzPos1];
    }
  } else {
    if (model->integrality_[Acol[nzPos2]] == HighsVarType::INTEGER) {
      // one col is integral, substitute the continuous one
      substcol = Acol[nzPos1];
      staycol = Acol[nzPos2];

      substcoef = Avalue[nzPos1];
      staycoef = Avalue[nzPos2];
    } else {
      // both columns continuous the one with a larger absolute coefficient
      // value if the difference is more than factor 2, and otherwise the one
      // with fewer nonzeros if those are equal

      double abs1Val = std::abs(Avalue[nzPos1]);
      double abs2Val = std::abs(Avalue[nzPos2]);
      bool colAtPos1Better;
      if (abs1Val > 0.5 * abs2Val)
        colAtPos1Better = true;
      else if (abs2Val > 0.5 * abs1Val)
        colAtPos1Better = false;
      else
        colAtPos1Better = colsize[Acol[nzPos1]] < colsize[Acol[nzPos2]];

      if (colAtPos1Better) {
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
    }
  }

  double oldStayLower = model->colLower_[staycol];
  double oldStayUpper = model->colUpper_[staycol];
  double substLower = model->colLower_[substcol];
  double substUpper = model->colUpper_[substcol];

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
    stayImplLower =
        substLower == -HIGHS_CONST_INF
            ? -HIGHS_CONST_INF
            : double((HighsCDouble(rhs) - substcoef * substLower) / staycoef);
    stayImplUpper =
        substUpper == HIGHS_CONST_INF
            ? HIGHS_CONST_INF
            : double((HighsCDouble(rhs) - substcoef * substUpper) / staycoef);
  } else {
    stayImplLower =
        substUpper == HIGHS_CONST_INF
            ? -HIGHS_CONST_INF
            : double((HighsCDouble(rhs) - substcoef * substUpper) / staycoef);
    stayImplUpper =
        substLower == -HIGHS_CONST_INF
            ? HIGHS_CONST_INF
            : double((HighsCDouble(rhs) - substcoef * substLower) / staycoef);
  }

  // possibly tighten bounds of the column that stays
  bool lowerTightened = false;
  bool upperTightened = false;
  if (stayImplLower > oldStayLower + options->primal_feasibility_tolerance) {
    lowerTightened = true;
    changeColLower(staycol, stayImplLower);
  }

  if (stayImplUpper < oldStayUpper - options->primal_feasibility_tolerance) {
    upperTightened = true;
    changeColUpper(staycol, stayImplUpper);
  }

  postSolveStack.doubletonEquation(row, substcol, staycol, substcoef, staycoef,
                                   rhs, substLower, substUpper,
                                   model->colCost_[substcol], lowerTightened,
                                   upperTightened, getColumnVector(substcol));

  // finally modify matrix
  markColDeleted(substcol);
  removeRow(row);
  substitute(substcol, staycol, rhs / substcoef, -staycoef / substcoef);

  // since a column was deleted we might have new row singletons which we
  // immediately remove
  HPRESOLVE_CHECKED_CALL(removeRowSingletons(postSolveStack));

  return checkLimits(postSolveStack);
}

HPresolve::Result HPresolve::singletonRow(HighsPostsolveStack& postSolveStack,
                                          HighsInt row) {
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
    if (model->colUpper_[col] * val <=
            model->rowUpper_[row] + options->primal_feasibility_tolerance &&
        model->colLower_[col] * val >=
            model->rowLower_[row] - options->primal_feasibility_tolerance) {
      postSolveStack.redundantRow(row);
      return checkLimits(postSolveStack);
    }
  } else {
    if (model->colLower_[col] * val <=
            model->rowUpper_[row] + options->primal_feasibility_tolerance &&
        model->colUpper_[col] * val >=
            model->rowLower_[row] - options->primal_feasibility_tolerance) {
      postSolveStack.redundantRow(row);
      return checkLimits(postSolveStack);
    }
  }

  // zeros should not be linked in the matrix
  assert(std::abs(val) > options->small_matrix_value);

  double newColUpper = HIGHS_CONST_INF;
  double newColLower = -HIGHS_CONST_INF;
  if (val > 0) {
    if (model->rowUpper_[row] != HIGHS_CONST_INF)
      newColUpper = model->rowUpper_[row] / val;
    if (model->rowLower_[row] != -HIGHS_CONST_INF)
      newColLower = model->rowLower_[row] / val;
  } else {
    if (model->rowUpper_[row] != HIGHS_CONST_INF)
      newColLower = model->rowUpper_[row] / val;
    if (model->rowLower_[row] != -HIGHS_CONST_INF)
      newColUpper = model->rowLower_[row] / val;
  }

  bool lowerTightened = newColLower > model->colLower_[col] +
                                          options->primal_feasibility_tolerance;
  bool upperTightened = newColUpper < model->colUpper_[col] -
                                          options->primal_feasibility_tolerance;
  double lb = lowerTightened ? newColLower : model->colLower_[col];
  double ub = upperTightened ? newColUpper : model->colUpper_[col];

  // printf("old bounds [%.15g,%.15g], new bounds [%.15g,%.15g] ... ",
  //        model->colLower_[col], model->colUpper_[col], lb, ub);
  // check whether the bounds are equal in tolerances
  if (ub <= lb + options->primal_feasibility_tolerance) {
    // bounds could be infeasible or equal in tolerances, first check infeasible
    if (ub < lb - options->primal_feasibility_tolerance)
      return Result::PrimalInfeasible;

    // bounds are equal in tolerances, if they have a slight infeasibility below
    // those tolerances or they have a slight numerical distance which changes
    // the largest contribution below feasibility tolerance then we can safely
    // set the bound to one of the values. To heuristically get rid of numerical
    // errors we choose the bound that was not tightened, or the midpoint if
    // both where tightened.
    if (ub < lb || (ub > lb && (ub - lb) * getMaxAbsColVal(col) <=
                                   options->primal_feasibility_tolerance)) {
      if (lowerTightened && upperTightened) {
        ub = 0.5 * (ub + lb);
        lb = ub;
        lowerTightened = lb > model->colLower_[col];
        upperTightened = ub < model->colUpper_[col];
      } else if (lowerTightened) {
        lb = ub;
        lowerTightened = lb > model->colLower_[col];
      } else {
        ub = lb;
        upperTightened = ub < model->colUpper_[col];
      }
    }
  }

  // printf("final bounds: [%.15g,%.15g]\n", lb, ub);

  postSolveStack.singletonRow(row, col, val, lowerTightened, upperTightened);

  // just update bounds (and row activities)
  if (lowerTightened) changeColLower(col, lb);
  // update bounds, or remove as fixed column directly
  if (ub == lb) {
    postSolveStack.removedFixedCol(col, lb, model->colCost_[col],
                                   getColumnVector(col));
    removeFixedCol(col);
  } else if (upperTightened)
    changeColUpper(col, ub);

  if (!colDeleted[col] && colsize[col] == 0)
    return emptyCol(postSolveStack, col);

  return checkLimits(postSolveStack);
}

HPresolve::Result HPresolve::singletonCol(HighsPostsolveStack& postSolveStack,
                                          HighsInt col) {
  assert(colsize[col] == 1);
  assert(!colDeleted[col]);
  HighsInt nzPos = colhead[col];
  HighsInt row = Arow[nzPos];
  double colCoef = Avalue[nzPos];

  double colDualUpper =
      impliedDualRowBounds.getSumUpper(col, model->colCost_[col]);
  double colDualLower =
      impliedDualRowBounds.getSumLower(col, model->colCost_[col]);

  // check for dominated column
  if (colDualLower > options->dual_feasibility_tolerance) {
    if (model->colLower_[col] == -HIGHS_CONST_INF)
      return Result::DualInfeasible;
    else
      fixColToLower(postSolveStack, col);
    return checkLimits(postSolveStack);
  }

  if (colDualUpper < -options->dual_feasibility_tolerance) {
    if (model->colUpper_[col] == HIGHS_CONST_INF)
      return Result::DualInfeasible;
    else
      fixColToUpper(postSolveStack, col);
    return checkLimits(postSolveStack);
  }

  // check for weakly dominated column
  if (colDualUpper <= options->dual_feasibility_tolerance) {
    if (model->colUpper_[col] != HIGHS_CONST_INF)
      fixColToUpper(postSolveStack, col);
    else if (impliedDualRowBounds.getSumUpperOrig(col) == 0.0) {
      // todo: forcing column, since this implies colDual >= 0 and we
      // already checked that colDual <= 0 and since the cost are 0.0
      // all the rows are at a dual multiplier of zero and we can determine
      // one nonbasic row in postsolve, and make the other rows and the column
      // basic. The columns primal value is computed from the non-basic row
      // which is chosen such that the values of all rows are primal feasible
      // printf("removing forcing column of size %" HIGHSINT_FORMAT "\n",
      // colsize[col]);
      postSolveStack.forcingColumn(col, getColumnVector(col),
                                   model->colCost_[col], model->colLower_[col],
                                   true);
      markColDeleted(col);
      HighsInt coliter = colhead[col];
      while (coliter != -1) {
        HighsInt row = Arow[coliter];
        double rhs = Avalue[coliter] > 0.0 ? model->rowLower_[row]
                                           : model->rowUpper_[row];
        coliter = Anext[coliter];

        postSolveStack.forcingColumnRemovedRow(col, row, rhs,
                                               getRowVector(row));
        removeRow(row);
      }
    }
    return checkLimits(postSolveStack);
  }
  if (colDualLower >= -options->dual_feasibility_tolerance) {
    if (model->colLower_[col] != -HIGHS_CONST_INF)
      fixColToLower(postSolveStack, col);
    else if (impliedDualRowBounds.getSumLowerOrig(col) == 0.0) {
      // forcing column, since this implies colDual <= 0 and we already checked
      // that colDual >= 0
      // printf("removing forcing column of size %" HIGHSINT_FORMAT "\n",
      // colsize[col]);
      postSolveStack.forcingColumn(col, getColumnVector(col),
                                   model->colCost_[col], model->colUpper_[col],
                                   false);
      markColDeleted(col);
      HighsInt coliter = colhead[col];
      while (coliter != -1) {
        HighsInt row = Arow[coliter];
        double rhs = Avalue[coliter] > 0.0 ? model->rowUpper_[row]
                                           : model->rowLower_[row];
        coliter = Anext[coliter];

        postSolveStack.forcingColumnRemovedRow(col, row, rhs,
                                               getRowVector(row));
        removeRow(row);
      }
    }
    return checkLimits(postSolveStack);
  }

  if (mipsolver != nullptr &&
      model->integrality_[col] == HighsVarType::CONTINUOUS &&
      isImpliedInteger(col)) {
    model->integrality_[col] = HighsVarType::IMPLICIT_INTEGER;
    ++rowsizeImplInt[row];
    double ceilLower =
        std::ceil(model->colLower_[col] - options->mip_feasibility_tolerance);
    double floorUpper =
        std::floor(model->colUpper_[col] + options->mip_feasibility_tolerance);

    if (ceilLower > model->colLower_[col]) changeColLower(col, ceilLower);
    if (floorUpper < model->colUpper_[col]) changeColUpper(col, floorUpper);
  }

  updateColImpliedBounds(row, col, colCoef);

  if (model->integrality_[col] != HighsVarType::INTEGER)
    updateRowDualImpliedBounds(row, col, colCoef);

  // now check if column is implied free within an equation and substitute the
  // column if that is the case
  if (isDualImpliedFree(row) && isImpliedFree(col)) {
    if (model->integrality_[col] == HighsVarType::INTEGER &&
        !isImpliedIntegral(col))
      return Result::Ok;
    // todo, store which side of an implied free dual variable needs to be used
    // for substitution
    storeRow(row);

    HighsPostsolveStack::RowType rowType = HighsPostsolveStack::RowType::Eq;
    double rhs;
    if (model->rowLower_[row] == model->rowUpper_[row]) {
      rhs = model->rowUpper_[row];
      rowType = HighsPostsolveStack::RowType::Eq;
    } else if ((model->rowUpper_[row] != HIGHS_CONST_INF &&
                implRowDualLower[row] >=
                    -options->dual_feasibility_tolerance)) {
      rhs = model->rowUpper_[row];
      rowType = HighsPostsolveStack::RowType::Leq;
    } else {
      rhs = model->rowLower_[row];
      rowType = HighsPostsolveStack::RowType::Geq;
    }

    postSolveStack.freeColSubstitution(row, col, rhs, model->colCost_[col],
                                       rowType, getStoredRow(),
                                       getColumnVector(col));
    // todo, check integrality of coefficients and allow this
    substitute(row, col, rhs);

    return checkLimits(postSolveStack);
  }

  // todo: check for zero cost singleton and remove
  return Result::Ok;
}

HPresolve::Result HPresolve::rowPresolve(HighsPostsolveStack& postSolveStack,
                                         HighsInt row) {
  assert(!rowDeleted[row]);

  // handle special cases directly via a call to the specialized procedure
  switch (rowsize[row]) {
    default:
      break;
    case 0:
      if (model->rowUpper_[row] < -options->primal_feasibility_tolerance ||
          model->rowLower_[row] > options->primal_feasibility_tolerance)
        // model infeasible
        return Result::PrimalInfeasible;
      postSolveStack.redundantRow(row);
      markRowDeleted(row);
      return checkLimits(postSolveStack);
    case 1:
      return singletonRow(postSolveStack, row);
  }

  // printf("row presolve: ");
  // debugPrintRow(row);
  double impliedRowUpper = impliedRowBounds.getSumUpper(row);
  double impliedRowLower = impliedRowBounds.getSumLower(row);

  if (impliedRowLower >
          model->rowUpper_[row] + options->primal_feasibility_tolerance ||
      impliedRowUpper <
          model->rowLower_[row] - options->primal_feasibility_tolerance) {
    // model infeasible
    return Result::PrimalInfeasible;
  }

  if (impliedRowLower >=
          model->rowLower_[row] - options->primal_feasibility_tolerance &&
      impliedRowUpper <=
          model->rowUpper_[row] + options->primal_feasibility_tolerance) {
    // row is redundant
    postSolveStack.redundantRow(row);
    removeRow(row);
    return checkLimits(postSolveStack);
  }

  // todo: do additional single row presolve for mip here. It may assume a
  // non-redundant and non-infeasible row when considering variable and implied
  // bounds
  if (rowsizeInteger[row] != 0 || rowsizeImplInt[row] != 0) {
    if (model->rowLower_[row] == model->rowUpper_[row]) {
      if (rowsize[row] == 2) return doubletonEq(postSolveStack, row);
      // equation
      if (impliedRowLower != -HIGHS_CONST_INF &&
          impliedRowUpper != HIGHS_CONST_INF &&
          std::abs(impliedRowLower + impliedRowUpper -
                   2 * model->rowUpper_[row]) <= options->mip_epsilon) {
        double binCoef = std::abs(impliedRowUpper - model->rowUpper_[row]);
        // simple probing on equation case
        HighsInt binCol = -1;
        storeRow(row);
        for (const HighsSliceNonzero& nonz : getStoredRow()) {
          if (std::abs(std::abs(nonz.value()) - binCoef) <=
                  options->mip_epsilon &&
              model->integrality_[nonz.index()] == HighsVarType::INTEGER &&
              std::abs(model->colUpper_[nonz.index()] -
                       model->colLower_[nonz.index()] - 1.0) <=
                  options->mip_feasibility_tolerance) {
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
          // printf("simple probing case on row of size %" HIGHSINT_FORMAT "\n",
          // rowsize[row]);
          for (const HighsSliceNonzero& nonz : getStoredRow()) {
            if (nonz.index() == binCol) continue;

            if (model->colLower_[nonz.index()] ==
                model->colUpper_[nonz.index()]) {
              postSolveStack.removedFixedCol(nonz.index(),
                                             model->colLower_[nonz.index()],
                                             0.0, HighsEmptySlice());
              removeFixedCol(nonz.index());
              continue;
            }

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
              double scale = model->colLower_[nonz.index()] -
                             model->colUpper_[nonz.index()];
              double offset = model->colUpper_[nonz.index()] -
                              model->colLower_[binCol] * scale;
              postSolveStack.doubletonEquation(
                  -1, nonz.index(), binCol, 1.0, -scale, offset,
                  model->colLower_[nonz.index()],
                  model->colUpper_[nonz.index()], 0.0, false, false,
                  HighsEmptySlice());
              substitute(nonz.index(), binCol, offset, scale);
            } else {
              // This case yields the following implications:
              // binCol = lb -> nonzCol = lb
              // binCol = ub -> nonzCol = ub
              // as linear equation:
              // nonzCol = colLb + (colUb - colLb)(binCol - binLb)
              // nonzCol =
              //    colLb - binLb*(colUb - colLb) + (colUb - colLb)*binCol
              double scale = model->colUpper_[nonz.index()] -
                             model->colLower_[nonz.index()];
              double offset = model->colLower_[nonz.index()] -
                              model->colLower_[binCol] * scale;
              postSolveStack.doubletonEquation(
                  -1, nonz.index(), binCol, 1.0, -scale, offset,
                  model->colLower_[nonz.index()],
                  model->colUpper_[nonz.index()], 0.0, false, false,
                  HighsEmptySlice());
              substitute(nonz.index(), binCol, offset, scale);
            }
          }

          removeRow(row);
          HPRESOLVE_CHECKED_CALL(checkLimits(postSolveStack));
          return removeRowSingletons(postSolveStack);
        }
      }
    } else {
      // inequality or ranged row

      if (rowsize[row] == rowsizeInteger[row] + rowsizeImplInt[row]) {
        std::vector<double> rowCoefs;
        std::vector<HighsInt> rowIndex;
        rowCoefs.reserve(rowsize[row]);
        rowIndex.reserve(rowsize[row]);

        double deltaDown = model->rowLower_[row] == -HIGHS_CONST_INF
                               ? options->mip_feasibility_tolerance
                               : options->mip_epsilon;
        double deltaUp = model->rowUpper_[row] == HIGHS_CONST_INF
                             ? options->mip_feasibility_tolerance
                             : options->mip_epsilon;

        storeRow(row);
        for (const HighsSliceNonzero& nonz : getStoredRow()) {
          rowCoefs.push_back(nonz.value());
          rowIndex.push_back(nonz.index());
        }

        double intScale =
            HighsIntegers::integralScale(rowCoefs, deltaDown, deltaUp);

        if (intScale != 0.0 &&
            std::abs(intScale - 1.0) > options->mip_epsilon) {
          if (model->rowLower_[row] == -HIGHS_CONST_INF) {
            // <= inequality
            HighsCDouble rhs = model->rowUpper_[row] * intScale;
            bool success = true;
            double minRhsTightening = 0.0;
            for (HighsInt i = 0; i != rowsize[row]; ++i) {
              double coef = rowCoefs[i];
              HighsCDouble scaleCoef = HighsCDouble(coef) * intScale;
              HighsCDouble intCoef = floor(scaleCoef + 0.5);
              HighsCDouble coefDelta = intCoef - scaleCoef;
              rowCoefs[i] = double(intCoef);
              if (coefDelta < -options->mip_epsilon) {
                minRhsTightening =
                    std::max(-double(coefDelta), minRhsTightening);
              } else if (coefDelta > options->mip_epsilon) {
                if (model->colUpper_[rowIndex[i]] == HIGHS_CONST_INF) {
                  success = false;
                  break;
                }

                rhs += model->colUpper_[rowIndex[i]] * coefDelta;
              }
            }

            if (success) {
              HighsCDouble roundRhs =
                  floor(rhs + options->mip_feasibility_tolerance);
              if (rhs - roundRhs >= minRhsTightening - options->mip_epsilon) {
                // scaled and rounded is not weaker than the original constraint
                if (intScale < 100.0) {
                  // printf(
                  //     "scaling constraint to integral values with scale %g, "
                  //     "rounded scaled side from %g to %g\n",
                  //     intScale, double(rhs), double(roundRhs));
                  // the scale value is reasonably small, change the row values
                  // to be integral
                  model->rowUpper_[row] = double(roundRhs);
                  for (HighsInt i = 0; i != rowsize[row]; ++i)
                    addToMatrix(row, rowIndex[i],
                                rowCoefs[i] - Avalue[rowpositions[i]]);
                } else if (rhs - roundRhs <
                           minRhsTightening -
                               options->mip_feasibility_tolerance) {
                  // printf(
                  //     "tightening right hand side from %g to %g due to "
                  //     "rounding with integral scale %g\n",
                  //     model->rowUpper_[row], double(roundRhs / intScale),
                  //     intScale);
                  // scale value is large, so we scale back the altered
                  // constraint the scaled back constraint must be stronger than
                  // the original constraint for this to make sense with is
                  // checked with the condition above
                  model->rowUpper_[row] = double(roundRhs / intScale);
                  for (HighsInt i = 0; i != rowsize[row]; ++i) {
                    double delta = double(HighsCDouble(rowCoefs[i]) / intScale -
                                          Avalue[rowpositions[i]]);
                    if (std::abs(delta) > options->mip_epsilon)
                      addToMatrix(row, rowIndex[i], delta);
                  }
                }
              }
            }
          } else if (model->rowUpper_[row] == HIGHS_CONST_INF) {
            // >= inequality
            HighsCDouble rhs = model->rowLower_[row] * intScale;
            bool success = true;
            double minRhsTightening = 0.0;
            for (HighsInt i = 0; i != rowsize[row]; ++i) {
              double coef = rowCoefs[i];
              HighsCDouble scaleCoef = HighsCDouble(coef) * intScale;
              HighsCDouble intCoef = floor(scaleCoef + 0.5);
              HighsCDouble coefDelta = intCoef - scaleCoef;
              rowCoefs[i] = double(intCoef);
              if (coefDelta < -options->mip_epsilon) {
                if (model->colUpper_[rowIndex[i]] == HIGHS_CONST_INF) {
                  success = false;
                  break;
                }

                rhs += model->colUpper_[rowIndex[i]] * coefDelta;
              } else if (coefDelta > options->mip_epsilon) {
                minRhsTightening =
                    std::max(-double(coefDelta), minRhsTightening);
              }
            }

            if (success) {
              HighsCDouble roundRhs =
                  ceil(rhs - options->mip_feasibility_tolerance);
              if (rhs - roundRhs <= minRhsTightening + options->mip_epsilon) {
                // scaled and rounded is not weaker than the original constraint
                if (intScale < 100.0) {
                  // printf(
                  //     "scaling constraint to integral values with scale %g, "
                  //     "rounded scaled side from %g to %g\n",
                  //     intScale, double(rhs), double(roundRhs));
                  // the scale value is reasonably small, change the row values
                  // to be integral
                  model->rowLower_[row] = double(roundRhs);
                  for (HighsInt i = 0; i != rowsize[row]; ++i)
                    addToMatrix(row, rowIndex[i],
                                rowCoefs[i] - Avalue[rowpositions[i]]);
                } else if (rhs - roundRhs >
                           minRhsTightening +
                               options->mip_feasibility_tolerance) {
                  // scale value is large, so we scale back the altered
                  // constraint the scaled back constraint must be stronger than
                  // the original constraint for this to make sense with is
                  // checked with the condition above
                  // printf(
                  //     "tightening left hand side from %g to %g due to
                  //     rounding " "with integral scale %g\n",
                  //     model->rowLower_[row], double(roundRhs / intScale),
                  //     intScale);
                  model->rowLower_[row] = double(roundRhs / intScale);
                  for (HighsInt i = 0; i != rowsize[row]; ++i) {
                    double delta = double(HighsCDouble(rowCoefs[i]) / intScale -
                                          Avalue[rowpositions[i]]);
                    if (std::abs(delta) > options->mip_epsilon)
                      addToMatrix(row, rowIndex[i], delta);
                  }
                }
              }
            }
          } else {
            // ranged row or equation, can maybe tighten sides and
            HighsCDouble lhs = model->rowLower_[row] * intScale;
            HighsCDouble rhs = model->rowUpper_[row] * intScale;
            bool success = true;
            double minRhsTightening = 0.0;
            double minLhsTightening = 0.0;
            for (HighsInt i = 0; i != rowsize[row]; ++i) {
              double coef = rowCoefs[i];
              HighsCDouble scaleCoef = HighsCDouble(coef) * intScale;
              HighsCDouble intCoef = floor(scaleCoef + 0.5);
              HighsCDouble coefDelta = intCoef - scaleCoef;
              rowCoefs[i] = double(intCoef);
              if (coefDelta < -options->mip_epsilon) {
                // for the >= side of the constraint a smaller coefficient is
                // stronger: Therefore we relax the left hand side using the
                // bound constraint, if the bound is infinite, abort
                if (model->colUpper_[rowIndex[i]] == HIGHS_CONST_INF) {
                  success = false;
                  break;
                }

                lhs += model->colUpper_[rowIndex[i]] * coefDelta;
                minRhsTightening =
                    std::max(-double(coefDelta), minRhsTightening);
              } else if (coefDelta > options->mip_epsilon) {
                if (model->colUpper_[rowIndex[i]] == HIGHS_CONST_INF) {
                  success = false;
                  break;
                }

                rhs += model->colUpper_[rowIndex[i]] * coefDelta;

                // the coefficient was relaxed regarding the rows lower bound.
                // Therefore the lower bound should be tightened by at least
                // this amount for the scaled constraint to dominate the
                // unscaled constraint be rounded by at least this value
                minLhsTightening =
                    std::max(double(coefDelta), minLhsTightening);
              }
            }

            if (success) {
              HighsCDouble roundLhs =
                  ceil(lhs - options->mip_feasibility_tolerance);
              HighsCDouble roundRhs =
                  floor(rhs + options->mip_feasibility_tolerance);

              // rounded row proves infeasibility regardless of coefficient
              // values
              if (roundRhs - roundLhs < -0.5) return Result::PrimalInfeasible;

              if (roundLhs >= intScale * model->rowLower_[row] +
                                  minLhsTightening - options->mip_epsilon &&
                  roundRhs <= intScale * model->rowUpper_[row] -
                                  minRhsTightening + options->mip_epsilon) {
                // scaled row with adjusted coefficients and sides is not weaker
                // than the original row
                if (intScale < 100.0) {
                  // printf(
                  //     "scaling constraint to integral values with scale %g, "
                  //     "rounded scaled sides from %g to %g and %g to %g\n",
                  //     intScale, double(rhs), double(roundRhs), double(lhs),
                  //     double(roundLhs));
                  // the scale value is reasonably small, change the row values
                  // to be integral
                  model->rowLower_[row] = double(roundLhs);
                  model->rowUpper_[row] = double(roundRhs);
                  for (HighsInt i = 0; i != rowsize[row]; ++i)
                    addToMatrix(row, rowIndex[i],
                                rowCoefs[i] - Avalue[rowpositions[i]]);
                } else {
                  // scale value is large, just tighten the sides
                  roundLhs /= intScale;
                  roundRhs /= intScale;
                  if (roundRhs < model->rowUpper_[row] -
                                     options->mip_feasibility_tolerance)
                    model->rowUpper_[row] = double(roundRhs);
                  if (roundLhs > model->rowLower_[row] +
                                     options->mip_feasibility_tolerance)
                    model->rowLower_[row] = double(roundLhs);
                }
              }
            }
          }

          impliedRowUpper = impliedRowBounds.getSumUpper(row);
          impliedRowLower = impliedRowBounds.getSumLower(row);
        }
      }

      if (model->rowLower_[row] == -HIGHS_CONST_INF &&
          impliedRowUpper != HIGHS_CONST_INF) {
        HighsInt numTightened = 0;
        double maxCoefValue = impliedRowUpper - model->rowUpper_[row];
        HighsCDouble rhs = model->rowUpper_[row];
        for (const HighsSliceNonzero& nonz : getRowVector(row)) {
          if (model->integrality_[nonz.index()] == HighsVarType::CONTINUOUS)
            continue;

          if (nonz.value() >
              maxCoefValue + options->mip_feasibility_tolerance) {
            // <= contraint, we decrease the coefficient value and the right
            // hand side
            double delta = maxCoefValue - nonz.value();
            addToMatrix(row, nonz.index(), delta);
            rhs += delta * model->colUpper_[nonz.index()];
            ++numTightened;
          } else if (nonz.value() <
                     -maxCoefValue - options->mip_feasibility_tolerance) {
            double delta = -maxCoefValue - nonz.value();
            addToMatrix(row, nonz.index(), delta);
            rhs += delta * model->colLower_[nonz.index()];
            ++numTightened;
          }
        }

        model->rowUpper_[row] = double(rhs);
      }

      if (model->rowUpper_[row] == HIGHS_CONST_INF &&
          impliedRowLower != -HIGHS_CONST_INF) {
        HighsInt numTightened = 0;
        double maxCoefValue = model->rowLower_[row] - impliedRowLower;
        HighsCDouble rhs = model->rowLower_[row];
        for (const HighsSliceNonzero& nonz : getRowVector(row)) {
          if (model->integrality_[nonz.index()] == HighsVarType::CONTINUOUS)
            continue;

          if (nonz.value() >
              maxCoefValue + options->mip_feasibility_tolerance) {
            double delta = maxCoefValue - nonz.value();
            addToMatrix(row, nonz.index(), delta);
            rhs += delta * model->colLower_[nonz.index()];
            ++numTightened;
          } else if (nonz.value() <
                     -maxCoefValue - options->mip_feasibility_tolerance) {
            double delta = -maxCoefValue - nonz.value();
            addToMatrix(row, nonz.index(), delta);
            rhs += delta * model->colUpper_[nonz.index()];
            ++numTightened;
          }
        }

        model->rowLower_[row] = double(rhs);
      }
    }
  }

  impliedRowUpper = impliedRowBounds.getSumUpperOrig(row);
  impliedRowLower = impliedRowBounds.getSumLowerOrig(row);

  // printf("implied bounds without tightenings: [%g,%g]\n", baseiRLower,
  //        baseiRUpper);

  if (impliedRowUpper <=  // check for forcing row on the row lower bound
      model->rowLower_[row] + options->primal_feasibility_tolerance) {
    // the row upper bound that is implied by the column bounds is equal to
    // the row lower bound there for we can fix all columns at their bound
    // as this is the only feasible assignment for this row and then find a
    // suitable dual multiplier in postsolve. First we store the row on the
    // postsolve stack (forcingRow() call) afterwards we store each column
    // fixing on the postsolve stack. As the postsolve goes over the stack
    // in reverse, it will first restore the column primal and dual values
    // as the dual values are required to find the proper dual multiplier for
    // the row and the column that we put in the basis.
    storeRow(row);
    auto rowVector = getStoredRow();

    HighsInt nfixings = 0;
    for (const HighsSliceNonzero& nonzero : rowVector) {
      if (nonzero.value() > 0) {
        if (model->colUpper_[nonzero.index()] <= implColUpper[nonzero.index()])
          ++nfixings;
      } else {
        if (model->colLower_[nonzero.index()] >= implColLower[nonzero.index()])
          ++nfixings;
      }
    }

    if (nfixings == rowsize[row]) {
      postSolveStack.forcingRow(row, rowVector, model->rowLower_[row],
                                HighsPostsolveStack::RowType::Geq);
      // already mark the row as deleted, since otherwise it would be registered
      // as changed/singleton in the process of fixing and removing the
      // contained columns

      markRowDeleted(row);
      for (const HighsSliceNonzero& nonzero : rowVector) {
        if (nonzero.value() > 0) {
          // the upper bound of the column is as tight as the implied upper
          // bound or comes from this row, which means it is not used in the
          // rows implied bounds. Therefore we can fix the variable at its
          // upper bound.
          postSolveStack.fixedColAtUpper(nonzero.index(),
                                         model->colUpper_[nonzero.index()],
                                         model->colCost_[nonzero.index()],
                                         getColumnVector(nonzero.index()));
          if (model->colLower_[nonzero.index()] <
              model->colUpper_[nonzero.index()])
            changeColLower(nonzero.index(), model->colUpper_[nonzero.index()]);
          removeFixedCol(nonzero.index());
        } else {
          postSolveStack.fixedColAtLower(nonzero.index(),
                                         model->colLower_[nonzero.index()],
                                         model->colCost_[nonzero.index()],
                                         getColumnVector(nonzero.index()));

          if (model->colUpper_[nonzero.index()] >
              model->colLower_[nonzero.index()])
            changeColUpper(nonzero.index(), model->colLower_[nonzero.index()]);
          removeFixedCol(nonzero.index());
        }
      }
      // now the row might be empty, but not necessarily because the implied
      // column bounds might be implied by other rows in which case we cannot
      // fix the column
      postSolveStack.redundantRow(row);

      HPRESOLVE_CHECKED_CALL(removeRowSingletons(postSolveStack));
      return checkLimits(postSolveStack);
    }
    // if there are any new row singletons, also remove them immediately
  } else if (impliedRowLower >=
             model->rowUpper_[row] - options->primal_feasibility_tolerance) {
    // forcing row in the other direction
    storeRow(row);
    auto rowVector = getStoredRow();

    HighsInt nfixings = 0;
    for (const HighsSliceNonzero& nonzero : rowVector) {
      if (nonzero.value() < 0) {
        if (model->colUpper_[nonzero.index()] <= implColUpper[nonzero.index()])
          ++nfixings;
      } else {
        if (model->colLower_[nonzero.index()] >= implColLower[nonzero.index()])
          ++nfixings;
      }
    }
    if (nfixings == rowsize[row]) {
      postSolveStack.forcingRow(row, rowVector, model->rowUpper_[row],
                                HighsPostsolveStack::RowType::Leq);
      markRowDeleted(row);
      for (const HighsSliceNonzero& nonzero : rowVector) {
        if (nonzero.value() < 0) {
          postSolveStack.fixedColAtUpper(nonzero.index(),
                                         model->colUpper_[nonzero.index()],
                                         model->colCost_[nonzero.index()],
                                         getColumnVector(nonzero.index()));
          if (model->colLower_[nonzero.index()] <
              model->colUpper_[nonzero.index()])
            changeColLower(nonzero.index(), model->colUpper_[nonzero.index()]);

          removeFixedCol(nonzero.index());
        } else {
          postSolveStack.fixedColAtLower(nonzero.index(),
                                         model->colLower_[nonzero.index()],
                                         model->colCost_[nonzero.index()],
                                         getColumnVector(nonzero.index()));
          if (model->colUpper_[nonzero.index()] >
              model->colLower_[nonzero.index()])
            changeColUpper(nonzero.index(), model->colLower_[nonzero.index()]);

          removeFixedCol(nonzero.index());
        }
      }

      postSolveStack.redundantRow(row);

      HPRESOLVE_CHECKED_CALL(removeRowSingletons(postSolveStack));
      return checkLimits(postSolveStack);
    }
  }

  if (rowsize[row] == 2 && model->rowLower_[row] == model->rowUpper_[row])
    return doubletonEq(postSolveStack, row);

  bool hasRowUpper =
      model->rowUpper_[row] != HIGHS_CONST_INF ||
      implRowDualLower[row] > options->dual_feasibility_tolerance;
  bool hasRowLower =
      model->rowLower_[row] != HIGHS_CONST_INF ||
      implRowDualUpper[row] < -options->dual_feasibility_tolerance;

  if ((hasRowUpper && impliedRowBounds.getNumInfSumLowerOrig(row) <= 1) ||
      (hasRowLower && impliedRowBounds.getNumInfSumUpperOrig(row) <= 1)) {
    for (const HighsSliceNonzero& nonzero : getRowVector(row))
      updateColImpliedBounds(row, nonzero.index(), nonzero.value());
  }

  return checkLimits(postSolveStack);
}

HPresolve::Result HPresolve::emptyCol(HighsPostsolveStack& postSolveStack,
                                      HighsInt col) {
  if ((model->colCost_[col] > 0 && model->colLower_[col] == -HIGHS_CONST_INF) ||
      (model->colCost_[col] < 0 && model->colUpper_[col] == HIGHS_CONST_INF)) {
    if (std::abs(model->colCost_[col]) <= options->dual_feasibility_tolerance)
      model->colCost_[col] = 0;
    else
      return Result::DualInfeasible;
  }

  if (model->colCost_[col] > 0)
    fixColToLower(postSolveStack, col);
  else if (model->colCost_[col] < 0 ||
           std::abs(model->colUpper_[col]) < std::abs(model->colLower_[col]))
    fixColToUpper(postSolveStack, col);
  else if (model->colLower_[col] != -HIGHS_CONST_INF)
    fixColToLower(postSolveStack, col);
  else
    fixColToZero(postSolveStack, col);

  return checkLimits(postSolveStack);
}

HPresolve::Result HPresolve::colPresolve(HighsPostsolveStack& postSolveStack,
                                         HighsInt col) {
  assert(!colDeleted[col]);

  double boundDiff = model->colUpper_[col] - model->colLower_[col];
  if (boundDiff <= options->primal_feasibility_tolerance) {
    if (boundDiff <= options->small_matrix_value ||
        getMaxAbsColVal(col) * boundDiff <=
            options->primal_feasibility_tolerance) {
      if (boundDiff < -options->primal_feasibility_tolerance)
        return Result::PrimalInfeasible;
      postSolveStack.removedFixedCol(col, model->colLower_[col],
                                     model->colCost_[col],
                                     getColumnVector(col));
      removeFixedCol(col);
      return checkLimits(postSolveStack);
    }
  }

  switch (colsize[col]) {
    case 0:
      return emptyCol(postSolveStack, col);
    case 1:
      return singletonCol(postSolveStack, col);
    default:
      break;
  }

  double colDualUpper =
      impliedDualRowBounds.getSumUpper(col, model->colCost_[col]);
  double colDualLower =
      impliedDualRowBounds.getSumLower(col, model->colCost_[col]);

  // check for dominated column
  if (colDualLower > options->dual_feasibility_tolerance) {
    if (model->colLower_[col] == -HIGHS_CONST_INF)
      return Result::DualInfeasible;
    else {
      fixColToLower(postSolveStack, col);
      HPRESOLVE_CHECKED_CALL(removeRowSingletons(postSolveStack));
    }
    return checkLimits(postSolveStack);
  }

  if (colDualUpper < -options->dual_feasibility_tolerance) {
    if (model->colUpper_[col] == HIGHS_CONST_INF)
      return Result::DualInfeasible;
    else {
      fixColToUpper(postSolveStack, col);
      HPRESOLVE_CHECKED_CALL(removeRowSingletons(postSolveStack));
    }
    return checkLimits(postSolveStack);
  }

  // check for weakly dominated column
  if (colDualUpper <= options->dual_feasibility_tolerance) {
    if (model->colUpper_[col] != HIGHS_CONST_INF) {
      fixColToUpper(postSolveStack, col);
      HPRESOLVE_CHECKED_CALL(removeRowSingletons(postSolveStack));
      return checkLimits(postSolveStack);
    } else if (impliedDualRowBounds.getSumUpperOrig(col) == 0.0) {
      postSolveStack.forcingColumn(col, getColumnVector(col),
                                   model->colCost_[col], model->colLower_[col],
                                   true);
      markColDeleted(col);
      HighsInt coliter = colhead[col];
      while (coliter != -1) {
        HighsInt row = Arow[coliter];
        double rhs = Avalue[coliter] > 0.0 ? model->rowLower_[row]
                                           : model->rowUpper_[row];
        coliter = Anext[coliter];
        postSolveStack.forcingColumnRemovedRow(col, row, rhs,
                                               getRowVector(row));
        removeRow(row);
      }
    }
  } else if (colDualLower >= -options->dual_feasibility_tolerance) {
    // symmetric case for fixing to the lower bound
    if (model->colLower_[col] != -HIGHS_CONST_INF) {
      fixColToLower(postSolveStack, col);
      HPRESOLVE_CHECKED_CALL(removeRowSingletons(postSolveStack));
      return checkLimits(postSolveStack);
    } else if (impliedDualRowBounds.getSumLowerOrig(col) == 0.0) {
      postSolveStack.forcingColumn(col, getColumnVector(col),
                                   model->colCost_[col], model->colUpper_[col],
                                   false);
      markColDeleted(col);
      HighsInt coliter = colhead[col];
      while (coliter != -1) {
        HighsInt row = Arow[coliter];
        double rhs = Avalue[coliter] > 0.0 ? model->rowUpper_[row]
                                           : model->rowLower_[row];
        coliter = Anext[coliter];
        postSolveStack.forcingColumnRemovedRow(col, row, rhs,
                                               getRowVector(row));
        removeRow(row);
      }
    }
  }

  // column is not (weakly) dominated

  // integer columns cannot be used to tighten bounds on dual multipliers
  if (mipsolver != nullptr) {
    if (model->integrality_[col] == HighsVarType::INTEGER)
      return Result::Ok;
    else if (model->integrality_[col] == HighsVarType::CONTINUOUS &&
             isImpliedInteger(col)) {
      model->integrality_[col] = HighsVarType::IMPLICIT_INTEGER;
      for (const HighsSliceNonzero& nonzero : getColumnVector(col))
        ++rowsizeImplInt[nonzero.index()];
      double ceilLower =
          std::ceil(model->colLower_[col] - options->mip_feasibility_tolerance);
      double floorUpper = std::floor(model->colUpper_[col] +
                                     options->mip_feasibility_tolerance);

      if (ceilLower > model->colLower_[col]) changeColLower(col, ceilLower);
      if (floorUpper < model->colUpper_[col]) changeColUpper(col, floorUpper);
    }
  }

  // the associated dual constraint has an upper bound if there is an infinite
  // or redundant column lower bound as then the reduced cost of the column must
  // not be positive i.e. <= 0
  bool dualConsHasUpper = isLowerImplied(col);
  bool dualConsHasLower = isUpperImplied(col);

  // now check if we can expect to tighten at least one bound
  if ((dualConsHasLower && impliedDualRowBounds.getNumInfSumUpper(col) <= 1) ||
      (dualConsHasUpper && impliedDualRowBounds.getNumInfSumLower(col) <= 1)) {
    for (const HighsSliceNonzero& nonzero : getColumnVector(col))
      updateRowDualImpliedBounds(nonzero.index(), col, nonzero.value());
  }

  return Result::Ok;
}

HPresolve::Result HPresolve::initialRowAndColPresolve(
    HighsPostsolveStack& postSolveStack) {
  // do a full scan over the rows as the singleton arrays and the changed row
  // arrays are not initialized, also unset changedRowFlag so that the row will
  // be added to the changed row vector when it is changed after it was
  // processed
  for (HighsInt row = 0; row != model->numRow_; ++row) {
    if (rowDeleted[row]) continue;
    HPRESOLVE_CHECKED_CALL(rowPresolve(postSolveStack, row));
    changedRowFlag[row] = false;
  }

  // same for the columns
  for (HighsInt col = 0; col != model->numCol_; ++col) {
    if (colDeleted[col]) continue;
    HPRESOLVE_CHECKED_CALL(colPresolve(postSolveStack, col));
    changedColFlag[col] = false;
  }

  return checkLimits(postSolveStack);
}

HPresolve::Result HPresolve::fastPresolveLoop(
    HighsPostsolveStack& postSolveStack) {
  do {
    storeCurrentProblemSize();

    HPRESOLVE_CHECKED_CALL(presolveChangedRows(postSolveStack));

    HPRESOLVE_CHECKED_CALL(removeDoubletonEquations(postSolveStack));

    HPRESOLVE_CHECKED_CALL(presolveColSingletons(postSolveStack));

    HPRESOLVE_CHECKED_CALL(presolveChangedCols(postSolveStack));

  } while (problemSizeReduction() > 0.01);

  return Result::Ok;
}

HPresolve::Result HPresolve::presolve(HighsPostsolveStack& postSolveStack) {
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
  //    - aggregator // add limit that catches many subsitutions but stops when
  //    many failures, do not run exhaustively as now
  //    - if (changes found) start main loop from beginning
  //    - primal and dual matrix sparsification
  //    - if (changes found) fast presolve loop
  //    - stop
  //

  // convert model to minimization problem
  if (model->sense_ == ObjSense::MAXIMIZE) {
    for (HighsInt i = 0; i != model->numCol_; ++i)
      model->colCost_[i] = -model->colCost_[i];

    model->offset_ = -model->offset_;
    assert(std::isfinite(model->offset_));
    model->sense_ = ObjSense::MINIMIZE;
  }

  if (options->presolve != "off") {
    highsLogUser(options->log_options, HighsLogType::INFO,
                 "\nPresolving model\n");

    auto report = [&]() {
      HighsInt numCol = model->numCol_ - numDeletedCols;
      HighsInt numRow = model->numRow_ - numDeletedRows;
      HighsInt numNonz = Avalue.size() - freeslots.size();
      highsLogUser(options->log_options, HighsLogType::INFO,
                   "%" HIGHSINT_FORMAT " rows, %" HIGHSINT_FORMAT
                   " cols, %" HIGHSINT_FORMAT " nonzeros\n",
                   numRow, numCol, numNonz);
    };

    HPRESOLVE_CHECKED_CALL(initialRowAndColPresolve(postSolveStack));

    HighsInt numParallelRowColCalls = 0;
#if ENABLE_SPARSIFY_FOR_LP
    bool trySparsify = true;  // mipsolver != nullptr;
#else
    bool trySparsify = mipsolver != nullptr;
#endif
    bool tryProbing = mipsolver != nullptr;
    while (true) {
      report();

      HPRESOLVE_CHECKED_CALL(fastPresolveLoop(postSolveStack));

      storeCurrentProblemSize();

      // when presolving after a restart the clique table and implication
      // structure may contain substitutions which we apply directly before
      // running the aggregator as they might loose validity otherwise
      if (mipsolver != nullptr) {
        HPRESOLVE_CHECKED_CALL(applyConflictGraphSubstitutions(postSolveStack));
      }

      HPRESOLVE_CHECKED_CALL(aggregator(postSolveStack));

      if (problemSizeReduction() > 0.05) continue;

      if (trySparsify) {
        HighsInt numNz = numNonzeros();
        HPRESOLVE_CHECKED_CALL(sparsify(postSolveStack));
        double nzReduction = 100.0 * (1.0 - (numNonzeros() / (double)numNz));

        if (nzReduction > 0) {
          highsLogUser(options->log_options, HighsLogType::INFO,
                       "Sparsify removed %.1f%% of nonzeros\n", nzReduction);

          fastPresolveLoop(postSolveStack);
        }
        trySparsify = false;
      }

      if (numParallelRowColCalls < 5) {
        if (shrinkProblemEnabled && (numDeletedCols >= 0.5 * model->numCol_ ||
                                     numDeletedRows >= 0.5 * model->numRow_)) {
          shrinkProblem(postSolveStack);

          toCSC(model->Avalue_, model->Aindex_, model->Astart_);
          fromCSC(model->Avalue_, model->Aindex_, model->Astart_);
        }
        storeCurrentProblemSize();
        HPRESOLVE_CHECKED_CALL(detectParallelRowsAndCols(postSolveStack));
        ++numParallelRowColCalls;
        if (problemSizeReduction() > 0.05) continue;
      }

      HPRESOLVE_CHECKED_CALL(fastPresolveLoop(postSolveStack));

      strengthenInequalities();

      HPRESOLVE_CHECKED_CALL(fastPresolveLoop(postSolveStack));

      if (tryProbing) {
        detectImpliedIntegers();
        storeCurrentProblemSize();
        HPRESOLVE_CHECKED_CALL(runProbing(postSolveStack));
        tryProbing =
            probingContingent > numProbed && problemSizeReduction() > 1.0;
        trySparsify = true;
        if (problemSizeReduction() > 0.05) continue;
        HPRESOLVE_CHECKED_CALL(fastPresolveLoop(postSolveStack));
      }

      break;
    }

    report();
  } else {
    highsLogUser(options->log_options, HighsLogType::INFO,
                 "\nPresolve is switched off\n");
  }

  return Result::Ok;
}

HPresolve::Result HPresolve::checkLimits(HighsPostsolveStack& postSolveStack) {
  // todo: check timelimit
#if 0
  for (HighsInt row = 0; row != model->numRow_; ++row) {
    row = 515251;
    // if (rowDeleted[row]) continue;

    if (model->rowLower_[row] == model->rowUpper_[row]) {
      assert(eqiters[row] != equations.end());
      assert(eqiters[row]->first == rowsize[row]);
      assert(eqiters[row]->second == row);
    }

    // debugPrintRow(row);

    double iRUpper = 0.0;
    double iRLower = 0.0;
    HighsInt rowlen = 0;
    for (const HighsSliceNonzero& nonz : getRowVector(row)) {
      double lb = colLowerSource[nonz.index()] == row
                      ? model->colLower_[nonz.index()]
                      : std::max(implColLower[nonz.index()],
                                 model->colLower_[nonz.index()]);
      double ub = colUpperSource[nonz.index()] == row
                      ? model->colUpper_[nonz.index()]
                      : std::min(implColUpper[nonz.index()],
                                 model->colUpper_[nonz.index()]);
      ++rowlen;
      if (nonz.value() > 0) {
        iRUpper += ub * nonz.value();
        iRLower += lb * nonz.value();
      } else {
        iRUpper += lb * nonz.value();
        iRLower += ub * nonz.value();
      }
    }

    assert(rowsize[row] == rowlen);

    double impliedRowLower = impliedRowBounds.getSumLower(row);
    double impliedRowUpper = impliedRowBounds.getSumUpper(row);
    assert(iRUpper == impliedRowUpper ||
           std::abs(iRUpper - impliedRowUpper) <=
               options->primal_feasibility_tolerance);
    assert(iRLower == impliedRowLower ||
           std::abs(iRLower - impliedRowLower) <=
               options->primal_feasibility_tolerance);
    break;
  }
//#else
  for (HighsInt col = 0; col != model->numCol_; ++col) {
    if (colDeleted[col]) continue;
    double iDRUpper = 0.0;
    double iDRLower = 0.0;
    HighsInt collen = 0;
    for (const HighsSliceNonzero& nonz : getColumnVector(col)) {
      double rdUpper = rowDualUpperSource[nonz.index()] != col
                           ? std::min(rowDualUpper[nonz.index()],
                                      implRowDualUpper[nonz.index()])
                           : rowDualUpper[nonz.index()];
      double rdLower = rowDualLowerSource[nonz.index()] != col
                           ? std::max(rowDualLower[nonz.index()],
                                      implRowDualLower[nonz.index()])
                           : rowDualLower[nonz.index()];
      ++collen;
      if (nonz.value() > 0) {
        iDRUpper += rdUpper * nonz.value();
        iDRLower += rdLower * nonz.value();
      } else {
        iDRUpper += rdLower * nonz.value();
        iDRLower += rdUpper * nonz.value();
      }
    }

    assert(colsize[col] == collen);

    double impliedDualRowLower = impliedDualRowBounds.getSumLower(col);
    double impliedDualRowUpper = impliedDualRowBounds.getSumUpper(col);

    assert(iDRUpper == impliedDualRowUpper ||
           std::abs(iDRUpper - impliedDualRowUpper) <=
               options->primal_feasibility_tolerance);
    assert(iDRLower == impliedDualRowLower ||
           std::abs(iDRLower - impliedDualRowLower) <=
               options->primal_feasibility_tolerance);
  }
#endif

  return postSolveStack.numReductions() >= reductionLimit ? Result::Stopped
                                                          : Result::Ok;
}

void HPresolve::storeCurrentProblemSize() {
  oldNumCol = model->numCol_ - numDeletedCols;
  oldNumRow = model->numRow_ - numDeletedRows;
}

double HPresolve::problemSizeReduction() {
  double colReduction =
      100.0 * double(oldNumCol - (model->numCol_ - numDeletedCols)) / oldNumCol;
  double rowReduction =
      100.0 * double(oldNumRow - (model->numRow_ - numDeletedRows)) / oldNumRow;

  return std::max(rowReduction, colReduction);
}

HighsModelStatus HPresolve::run(HighsPostsolveStack& postSolveStack) {
  shrinkProblemEnabled = true;
  switch (presolve(postSolveStack)) {
    case Result::Stopped:
    case Result::Ok:
      break;
    case Result::PrimalInfeasible:
      return HighsModelStatus::PRIMAL_INFEASIBLE;
    case Result::DualInfeasible:
      return HighsModelStatus::DUAL_INFEASIBLE;
  }

  shrinkProblem(postSolveStack);

  if (mipsolver != nullptr) {
    mipsolver->mipdata_->domain.addCutpool(mipsolver->mipdata_->cutpool);

    if (mipsolver->mipdata_->numRestarts != 0) {
      std::vector<HighsInt> cutinds;
      std::vector<double> cutvals;
      cutinds.reserve(model->numCol_);
      cutvals.reserve(model->numCol_);
      HighsInt numcuts = 0;
      for (HighsInt i = model->numRow_ - 1; i >= 0; --i) {
        // check if we already reached the original rows
        if (postSolveStack.getOrigRowIndex(i) < mipsolver->orig_model_->numRow_)
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
            model->rowUpper_[i],
            rowsizeInteger[i] + rowsizeImplInt[i] == rowsize[i] &&
                rowCoefficientsIntegral(i, 1.0));

        markRowDeleted(i);
        for (HighsInt j : rowpositions) unlink(j);
      }

      model->numRow_ -= numcuts;
      model->rowLower_.resize(model->numRow_);
      model->rowUpper_.resize(model->numRow_);
      model->row_names_.resize(model->numRow_);
    }
  }

  toCSC(model->Avalue_, model->Aindex_, model->Astart_);

  if (model->numCol_ == 0) {
    if (mipsolver) {
      mipsolver->mipdata_->upper_bound = 0;
      mipsolver->mipdata_->lower_bound = 0;
    }
    return HighsModelStatus::OPTIMAL;
  }

  if (!mipsolver) setRelaxedImpliedBounds();

  return HighsModelStatus::NOTSET;
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

  toCSC(model->Avalue_, model->Aindex_, model->Astart_);

  for (HighsInt i = 0; i != model->numRow_; ++i) flagRow[i] = 1 - rowDeleted[i];
  for (HighsInt i = 0; i != model->numCol_; ++i) flagCol[i] = 1 - colDeleted[i];
}

HPresolve::Result HPresolve::aggregator(HighsPostsolveStack& postSolveStack) {
  HighsInt numsubst = 0;
  HighsInt numsubstint = 0;
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

  std::sort(
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
    if (model->integrality_[col] == HighsVarType::INTEGER &&
        !isImpliedIntegral(col))
      continue;

    // in the case where the row has length two or the column has length two
    // we always do the substitution since the fillin can never be problematic
    if (rowsize[row] == 2 || colsize[col] == 2) {
      double rhs;
      HighsPostsolveStack::RowType rowType;
      if (model->rowLower_[row] == model->rowUpper_[row]) {
        rowType = HighsPostsolveStack::RowType::Eq;
        rhs = model->rowUpper_[row];
      } else if ((model->rowUpper_[row] != HIGHS_CONST_INF &&
                  implRowDualLower[row] >=
                      -options->dual_feasibility_tolerance)) {
        rowType = HighsPostsolveStack::RowType::Leq;
        rhs = model->rowUpper_[row];
        changeRowDualLower(row, -HIGHS_CONST_INF);
      } else {
        rowType = HighsPostsolveStack::RowType::Geq;
        rhs = model->rowLower_[row];
        changeRowDualUpper(row, HIGHS_CONST_INF);
      }

      ++numsubst;
      if (model->integrality_[col] == HighsVarType::INTEGER) ++numsubstint;
      storeRow(row);

      postSolveStack.freeColSubstitution(row, col, rhs, model->colCost_[col],
                                         rowType, getStoredRow(),
                                         getColumnVector(col));
      substitutionOpportunities[i].first = -1;

      substitute(row, col, rhs);
      HPRESOLVE_CHECKED_CALL(removeRowSingletons(postSolveStack));
      HPRESOLVE_CHECKED_CALL(checkLimits(postSolveStack));
      continue;
    }

    if (rowsize[row] < colsize[col]) {
      double maxVal = getMaxAbsRowVal(row);
      if (std::abs(Avalue[nzPos]) <
          maxVal * options->presolve_pivot_threshold) {
        maxVal = getMaxAbsColVal(col);
        if (std::abs(Avalue[nzPos]) <
            maxVal * options->presolve_pivot_threshold)
          substitutionOpportunities[i].first = -1;
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
    ++numsubst;
    if (model->integrality_[col] == HighsVarType::INTEGER) ++numsubstint;
    double rhs;
    HighsPostsolveStack::RowType rowType;
    if (model->rowLower_[row] == model->rowUpper_[row]) {
      rowType = HighsPostsolveStack::RowType::Eq;
      rhs = model->rowUpper_[row];
    } else if ((model->rowUpper_[row] != HIGHS_CONST_INF &&
                implRowDualLower[row] >=
                    -options->dual_feasibility_tolerance)) {
      rowType = HighsPostsolveStack::RowType::Leq;
      rhs = model->rowUpper_[row];
      changeRowDualLower(row, -HIGHS_CONST_INF);
    } else {
      rowType = HighsPostsolveStack::RowType::Geq;
      rhs = model->rowLower_[row];
      changeRowDualUpper(row, HIGHS_CONST_INF);
    }

    postSolveStack.freeColSubstitution(row, col, rhs, model->colCost_[col],
                                       rowType, getStoredRow(),
                                       getColumnVector(col));
    substitutionOpportunities[i].first = -1;
    substitute(row, col, rhs);
    HPRESOLVE_CHECKED_CALL(removeRowSingletons(postSolveStack));
    HPRESOLVE_CHECKED_CALL(checkLimits(postSolveStack));
  }

  substitutionOpportunities.erase(
      std::remove_if(
          substitutionOpportunities.begin(), substitutionOpportunities.end(),
          [](const std::pair<HighsInt, HighsInt>& p) { return p.first == -1; }),
      substitutionOpportunities.end());

  return Result::Ok;
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
    if (model->rowLower_[colrow] != -HIGHS_CONST_INF)
      model->rowLower_[colrow] -= colval * offset;

    if (model->rowUpper_[colrow] != HIGHS_CONST_INF)
      model->rowUpper_[colrow] -= colval * offset;

    addToMatrix(colrow, staycol, scale * colval);
    // printf("after substitution: ");
    // debugPrintRow(colrow);

    // check if this is an equation row and it now has a different size
    if (model->rowLower_[colrow] == model->rowUpper_[colrow] &&
        eqiters[colrow] != equations.end() &&
        eqiters[colrow]->first != rowsize[colrow]) {
      // if that is the case reinsert it into the equation set that is ordered
      // by sparsity
      equations.erase(eqiters[colrow]);
      eqiters[colrow] = equations.emplace(rowsize[colrow], colrow).first;
    }
  }

  // substitute column in the objective function
  if (model->colCost_[substcol] != 0.0) {
    model->offset_ += model->colCost_[substcol] * offset;
    assert(std::isfinite(model->offset_));

    model->colCost_[staycol] += scale * model->colCost_[substcol];

    if (std::abs(model->colCost_[staycol]) <= options->small_matrix_value)
      model->colCost_[staycol] = 0.0;
    model->colCost_[substcol] = 0.0;
  }
}

void HPresolve::fixColToLower(HighsPostsolveStack& postSolveStack,
                              HighsInt col) {
  double fixval = model->colLower_[col];

  // printf("fixing column %" HIGHSINT_FORMAT " to %.15g\n", col, fixval);

  // mark the column as deleted first so that it is not registered as singleton
  // column upon removing its nonzeros
  postSolveStack.fixedColAtLower(col, fixval, model->colCost_[col],
                                 getColumnVector(col));
  markColDeleted(col);

  for (HighsInt coliter = colhead[col]; coliter != -1;) {
    HighsInt colrow = Arow[coliter];
    double colval = Avalue[coliter];
    assert(Acol[coliter] == col);

    HighsInt colpos = coliter;
    coliter = Anext[coliter];

    if (model->rowLower_[colrow] != -HIGHS_CONST_INF)
      model->rowLower_[colrow] -= colval * fixval;

    if (model->rowUpper_[colrow] != HIGHS_CONST_INF)
      model->rowUpper_[colrow] -= colval * fixval;

    unlink(colpos);

    if (model->rowLower_[colrow] == model->rowUpper_[colrow] &&
        eqiters[colrow] != equations.end() &&
        eqiters[colrow]->first != rowsize[colrow]) {
      // if that is the case reinsert it into the equation set that is ordered
      // by sparsity
      equations.erase(eqiters[colrow]);
      eqiters[colrow] = equations.emplace(rowsize[colrow], colrow).first;
    }
  }

  model->offset_ += model->colCost_[col] * fixval;
  assert(std::isfinite(model->offset_));
  model->colCost_[col] = 0;
}

void HPresolve::fixColToUpper(HighsPostsolveStack& postSolveStack,
                              HighsInt col) {
  double fixval = model->colUpper_[col];

  // printf("fixing column %" HIGHSINT_FORMAT " to %.15g\n", col, fixval);

  // mark the column as deleted first so that it is not registered as singleton
  // column upon removing its nonzeros
  postSolveStack.fixedColAtUpper(col, fixval, model->colCost_[col],
                                 getColumnVector(col));
  markColDeleted(col);

  for (HighsInt coliter = colhead[col]; coliter != -1;) {
    HighsInt colrow = Arow[coliter];
    double colval = Avalue[coliter];
    assert(Acol[coliter] == col);

    HighsInt colpos = coliter;
    coliter = Anext[coliter];

    if (model->rowLower_[colrow] != -HIGHS_CONST_INF)
      model->rowLower_[colrow] -= colval * fixval;

    if (model->rowUpper_[colrow] != HIGHS_CONST_INF)
      model->rowUpper_[colrow] -= colval * fixval;

    unlink(colpos);

    if (model->rowLower_[colrow] == model->rowUpper_[colrow] &&
        eqiters[colrow] != equations.end() &&
        eqiters[colrow]->first != rowsize[colrow]) {
      // if that is the case reinsert it into the equation set that is ordered
      // by sparsity
      equations.erase(eqiters[colrow]);
      eqiters[colrow] = equations.emplace(rowsize[colrow], colrow).first;
    }
  }

  model->offset_ += model->colCost_[col] * fixval;
  assert(std::isfinite(model->offset_));
  model->colCost_[col] = 0;
}

void HPresolve::fixColToZero(HighsPostsolveStack& postSolveStack,
                             HighsInt col) {
  postSolveStack.fixedColAtZero(col, model->colCost_[col],
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

    if (model->rowLower_[colrow] == model->rowUpper_[colrow] &&
        eqiters[colrow] != equations.end() &&
        eqiters[colrow]->first != rowsize[colrow]) {
      // if that is the case reinsert it into the equation set that is ordered
      // by sparsity
      equations.erase(eqiters[colrow]);
      eqiters[colrow] = equations.emplace(rowsize[colrow], colrow).first;
    }
  }

  model->colCost_[col] = 0;
}

void HPresolve::removeRow(HighsInt row) {
  assert(row < int(rowroot.size()));
  assert(row >= 0);
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
  double fixval = model->colLower_[col];

  markColDeleted(col);

  for (HighsInt coliter = colhead[col]; coliter != -1;) {
    HighsInt colrow = Arow[coliter];
    double colval = Avalue[coliter];
    assert(Acol[coliter] == col);

    HighsInt colpos = coliter;
    coliter = Anext[coliter];

    if (model->rowLower_[colrow] != -HIGHS_CONST_INF)
      model->rowLower_[colrow] -= colval * fixval;

    if (model->rowUpper_[colrow] != HIGHS_CONST_INF)
      model->rowUpper_[colrow] -= colval * fixval;

    unlink(colpos);

    if (model->rowLower_[colrow] == model->rowUpper_[colrow] &&
        eqiters[colrow] != equations.end() &&
        eqiters[colrow]->first != rowsize[colrow]) {
      // if that is the case reinsert it into the equation set that is ordered
      // by sparsity
      equations.erase(eqiters[colrow]);
      eqiters[colrow] = equations.emplace(rowsize[colrow], colrow).first;
    }
  }

  model->offset_ += model->colCost_[col] * fixval;
  assert(std::isfinite(model->offset_));
  model->colCost_[col] = 0;
}

HPresolve::Result HPresolve::removeRowSingletons(
    HighsPostsolveStack& postSolveStack) {
  for (size_t i = 0; i != singletonRows.size(); ++i) {
    HighsInt row = singletonRows[i];
    if (rowDeleted[row] || rowsize[row] > 1) continue;
    // row presolve will delegate to rowSingleton() if the row size is 1
    // if the singleton row has become empty it will also remove the row
    HPRESOLVE_CHECKED_CALL(rowPresolve(postSolveStack, row));
  }

  singletonRows.clear();

  return Result::Ok;
}

HPresolve::Result HPresolve::presolveColSingletons(
    HighsPostsolveStack& postSolveStack) {
  for (size_t i = 0; i != singletonColumns.size(); ++i) {
    HighsInt col = singletonColumns[i];
    if (colDeleted[col]) continue;
    HPRESOLVE_CHECKED_CALL(colPresolve(postSolveStack, col));
  }
  singletonColumns.erase(
      std::remove_if(
          singletonColumns.begin(), singletonColumns.end(),
          [&](HighsInt col) { return colDeleted[col] || colsize[col] > 1; }),
      singletonColumns.end());

  return Result::Ok;
}

HPresolve::Result HPresolve::presolveChangedRows(
    HighsPostsolveStack& postSolveStack) {
  std::vector<HighsInt> changedRows;
  changedRows.reserve(model->numRow_ - numDeletedRows);
  changedRows.swap(changedRowIndices);
  for (HighsInt row : changedRows) {
    if (rowDeleted[row]) continue;
    HPRESOLVE_CHECKED_CALL(rowPresolve(postSolveStack, row));
    changedRowFlag[row] = rowDeleted[row];
  }

  return Result::Ok;
}

HPresolve::Result HPresolve::presolveChangedCols(
    HighsPostsolveStack& postSolveStack) {
  std::vector<HighsInt> changedCols;
  changedCols.reserve(model->numCol_ - numDeletedCols);
  changedCols.swap(changedColIndices);
  for (HighsInt col : changedCols) {
    if (colDeleted[col]) continue;
    HPRESOLVE_CHECKED_CALL(colPresolve(postSolveStack, col));
    changedColFlag[col] = colDeleted[col];
  }

  return Result::Ok;
}

HPresolve::Result HPresolve::removeDoubletonEquations(
    HighsPostsolveStack& postSolveStack) {
  auto eq = equations.begin();
  while (eq != equations.end()) {
    HighsInt eqrow = eq->second;
    assert(!rowDeleted[eqrow]);
    assert(eq->first == rowsize[eqrow]);
    assert(model->rowLower_[eqrow] == model->rowUpper_[eqrow]);
    if (rowsize[eqrow] > 2) return Result::Ok;
    HPRESOLVE_CHECKED_CALL(rowPresolve(postSolveStack, eqrow));
    if (rowDeleted[eqrow])
      eq = equations.begin();
    else
      ++eq;
  }

  return Result::Ok;
}

HighsInt HPresolve::strengthenInequalities() {
  std::vector<int8_t> complementation;
  std::vector<double> reducedcost;
  std::vector<double> upper;
  std::vector<HighsInt> indices;
  std::vector<HighsInt> positions;
  std::vector<HighsInt> stack;
  std::vector<double> coefs;
  std::vector<HighsInt> cover;

  HighsInt numstrenghtened = 0;

  for (HighsInt row = 0; row != model->numRow_; ++row) {
    if (rowsize[row] <= 1) continue;
    if (model->rowLower_[row] != -HIGHS_CONST_INF &&
        model->rowUpper_[row] != HIGHS_CONST_INF)
      continue;

    // do not run on very dense rows as this could get expensive
    if (rowsize[row] >
        std::max(HighsInt{1000},
                 HighsInt(0.05 * (model->numCol_ - numDeletedCols))))
      continue;

    // printf("strengthening knapsack of %" HIGHSINT_FORMAT " vars\n",
    // rowsize[row]);

    HighsCDouble maxviolation;
    HighsCDouble continuouscontribution = 0.0;
    double scale;

    if (model->rowLower_[row] != -HIGHS_CONST_INF) {
      maxviolation = model->rowLower_[row];
      scale = -1.0;
    } else {
      maxviolation = -model->rowUpper_[row];
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
    stack.reserve(rowsize[row]);
    stack.push_back(rowroot[row]);

    bool skiprow = false;

    while (!stack.empty()) {
      HighsInt pos = stack.back();
      stack.pop_back();

      if (ARright[pos] != -1) stack.push_back(ARright[pos]);
      if (ARleft[pos] != -1) stack.push_back(ARleft[pos]);

      int8_t comp;
      double weight;
      double ub;
      weight = Avalue[pos] * scale;
      HighsInt col = Acol[pos];
      ub = model->colUpper_[col] - model->colLower_[col];

      if (ub == HIGHS_CONST_INF) {
        skiprow = true;
        break;
      }

      if (weight > 0) {
        if (model->colUpper_[col] == HIGHS_CONST_INF) {
          skiprow = true;
          break;
        }

        comp = 1;
        maxviolation += model->colUpper_[col] * weight;
      } else {
        if (model->colLower_[col] == -HIGHS_CONST_INF) {
          skiprow = true;
          break;
        }
        comp = -1;
        maxviolation += model->colLower_[col] * weight;
        weight = -weight;
      }

      if (ub <= options->mip_feasibility_tolerance ||
          weight <= options->mip_feasibility_tolerance)
        continue;

      if (model->integrality_[col] == HighsVarType::CONTINUOUS) {
        continuouscontribution += weight * ub;
        continue;
      }

      indices.push_back(reducedcost.size());
      positions.push_back(pos);
      reducedcost.push_back(weight);
      complementation.push_back(comp);
      upper.push_back(ub);
    }

    if (skiprow) {
      stack.clear();
      continue;
    }

    const double smallVal =
        std::max(10 * options->mip_feasibility_tolerance,
                 options->mip_feasibility_tolerance * double(maxviolation));
    while (true) {
      if (maxviolation - continuouscontribution <= smallVal || indices.empty())
        break;

      std::sort(indices.begin(), indices.end(), [&](HighsInt i1, HighsInt i2) {
        return std::make_pair(reducedcost[i1], i1) >
               std::make_pair(reducedcost[i2], i2);
      });

      HighsCDouble lambda = maxviolation - continuouscontribution;

      cover.clear();
      cover.reserve(indices.size());

      for (HighsInt i = indices.size() - 1; i >= 0; --i) {
        double delta = upper[indices[i]] * reducedcost[indices[i]];

        if (reducedcost[indices[i]] > smallVal && lambda - delta <= smallVal)
          cover.push_back(indices[i]);
        else
          lambda -= delta;
      }

      if (cover.empty() || lambda <= smallVal) break;

      HighsInt alpos = *std::min_element(
          cover.begin(), cover.end(), [&](HighsInt i1, HighsInt i2) {
            return reducedcost[i1] < reducedcost[i2];
          });

      HighsInt coverend = cover.size();

      double al = reducedcost[alpos];
      coefs.resize(coverend);
      double coverrhs = std::max(
          std::ceil(double(lambda / al - options->mip_feasibility_tolerance)),
          1.0);
      HighsCDouble slackupper = -coverrhs;

      double step = HIGHS_CONST_INF;
      for (HighsInt i = 0; i != coverend; ++i) {
        coefs[i] =
            std::ceil(std::min(reducedcost[cover[i]], double(lambda)) / al -
                      options->small_matrix_value);
        slackupper += upper[cover[i]] * coefs[i];
        step = std::min(step, reducedcost[cover[i]] / coefs[i]);
      }
      step = std::min(step, double(maxviolation / coverrhs));
      maxviolation -= step * coverrhs;

      HighsInt slackind = reducedcost.size();
      reducedcost.push_back(step);
      upper.push_back(double(slackupper));

      for (HighsInt i = 0; i != coverend; ++i)
        reducedcost[cover[i]] -= step * coefs[i];

      indices.erase(std::remove_if(indices.begin(), indices.end(),
                                   [&](HighsInt i) {
                                     return reducedcost[i] <=
                                            options->mip_feasibility_tolerance;
                                   }),
                    indices.end());
      indices.push_back(slackind);
    }

    double threshold =
        double(maxviolation + options->mip_feasibility_tolerance);

    indices.erase(std::remove_if(indices.begin(), indices.end(),
                                 [&](HighsInt i) {
                                   return i >= (HighsInt)positions.size() ||
                                          std::abs(reducedcost[i]) <= threshold;
                                 }),
                  indices.end());
    if (indices.empty()) continue;

    if (scale == -1.0) {
      HighsCDouble lhs = model->rowLower_[row];
      for (HighsInt i : indices) {
        double coefdelta = double(reducedcost[i] - maxviolation);
        HighsInt pos = positions[i];

        if (complementation[i] == -1) {
          lhs -= coefdelta * model->colLower_[Acol[pos]];
          addToMatrix(row, Acol[pos], -coefdelta);
        } else {
          lhs += coefdelta * model->colUpper_[Acol[pos]];
          addToMatrix(row, Acol[pos], coefdelta);
        }
      }

      model->rowLower_[row] = double(lhs);
    } else {
      HighsCDouble rhs = model->rowUpper_[row];
      for (HighsInt i : indices) {
        double coefdelta = double(reducedcost[i] - maxviolation);
        HighsInt pos = positions[i];

        if (complementation[i] == -1) {
          rhs += coefdelta * model->colLower_[Acol[pos]];
          addToMatrix(row, Acol[pos], coefdelta);
        } else {
          rhs -= coefdelta * model->colUpper_[Acol[pos]];
          addToMatrix(row, Acol[pos], -coefdelta);
        }
      }

      model->rowUpper_[row] = double(rhs);
    }

    numstrenghtened += indices.size();
  }

  return numstrenghtened;
}

HighsInt HPresolve::detectImpliedIntegers() {
  HighsInt numImplInt = 0;

  for (HighsInt col = 0; col != model->numCol_; ++col) {
    if (colDeleted[col]) continue;
    if (model->integrality_[col] == HighsVarType::CONTINUOUS &&
        isImpliedInteger(col)) {
      ++numImplInt;
      model->integrality_[col] = HighsVarType::IMPLICIT_INTEGER;

      for (const HighsSliceNonzero& nonzero : getColumnVector(col))
        ++rowsizeImplInt[nonzero.index()];

      double ceilLower =
          std::ceil(model->colLower_[col] - options->mip_feasibility_tolerance);
      double floorUpper = std::floor(model->colUpper_[col] +
                                     options->mip_feasibility_tolerance);

      if (ceilLower > model->colLower_[col]) changeColLower(col, ceilLower);
      if (floorUpper < model->colUpper_[col]) changeColUpper(col, floorUpper);
    }
  }

  return numImplInt;
}

HPresolve::Result HPresolve::detectParallelRowsAndCols(
    HighsPostsolveStack& postSolveStack) {
  std::vector<std::uint64_t> rowHashes;
  std::vector<std::uint64_t> colHashes;
  std::vector<std::pair<double, HighsInt>> rowMax(rowsize.size());
  std::vector<std::pair<double, HighsInt>> colMax(colsize.size());

  HighsHashTable<HighsInt, HighsInt> numRowSingletons;

  HighsInt nnz = Avalue.size();
  rowHashes.assign(rowsize.begin(), rowsize.end());
  colHashes.assign(colsize.begin(), colsize.end());

  // Step 1: Determine scales for rows and columns and remove column singletons
  // from the intial row hashes which are initialized with the row sizes
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
    // we use the nonzero with the smalles row/column index for the column/row
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
                                           Avalue[i] / colMax[Acol[i]].second));
    }
  }

  // Step 3: Loop over the rows and columns and put them into buckets using the
  // computed hash values. Whenever a bucket already contains a row/column,
  // check if we can apply a (nearly) parallel row reduction or a
  // parallel/dominated column reduction.
  HighsInt numRowBuckets = 0;
  HighsInt numColBuckets = 0;

  std::unordered_multimap<std::uint64_t, HighsInt> buckets;

  for (HighsInt i = 0; i != model->numCol_; ++i) {
    if (colDeleted[i]) continue;
    if (colsize[i] == 0) {
      HPRESOLVE_CHECKED_CALL(colPresolve(postSolveStack, i));
      continue;
    }
    auto it = buckets.find(colHashes[i]);
    decltype(it) last = it;

    HighsInt delCol = -1;
    HighsInt parallelColCandidate = -2;

    if (it == buckets.end()) ++numColBuckets;
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
      auto colUpperInf = [&]() {
        if (!checkColImplBounds) return false;
        if (mipsolver == nullptr) {
          // for LP we check strict reduncancy of the bounds as otherwise dual
          // postsolve might fail when the bound is used in the optimal solution
          return colScale > 0
                     ? model->colUpper_[col] == HIGHS_CONST_INF ||
                           implColUpper[col] <
                               model->colUpper_[col] -
                                   options->primal_feasibility_tolerance
                     : model->colLower_[col] == -HIGHS_CONST_INF ||
                           implColLower[col] >
                               model->colLower_[col] +
                                   options->primal_feasibility_tolerance;
        } else {
          // for MIP we do not need dual postsolve so the reduction is valid if
          // the bound is weakly redundant
          return colScale > 0 ? model->colUpper_[col] == HIGHS_CONST_INF ||
                                    implColUpper[col] <=
                                        model->colUpper_[col] +
                                            options->mip_feasibility_tolerance
                              : model->colLower_[col] == -HIGHS_CONST_INF ||
                                    implColLower[col] >=
                                        model->colLower_[col] -
                                            options->mip_feasibility_tolerance;
        }
      };

      auto colLowerInf = [&]() {
        if (!checkColImplBounds) return false;
        if (mipsolver == nullptr) {
          return colScale > 0
                     ? model->colLower_[col] == -HIGHS_CONST_INF ||
                           implColLower[col] >
                               model->colLower_[col] +
                                   options->primal_feasibility_tolerance
                     : model->colUpper_[col] == HIGHS_CONST_INF ||
                           implColUpper[col] <
                               model->colUpper_[col] -
                                   options->primal_feasibility_tolerance;
        } else {
          return colScale > 0 ? model->colLower_[col] == -HIGHS_CONST_INF ||
                                    implColLower[col] >=
                                        model->colLower_[col] -
                                            options->mip_feasibility_tolerance
                              : model->colUpper_[col] == HIGHS_CONST_INF ||
                                    implColUpper[col] <=
                                        model->colUpper_[col] +
                                            options->mip_feasibility_tolerance;
        }
      };

      auto duplicateColUpperInf = [&]() {
        if (!checkDuplicateColImplBounds) return false;
        if (mipsolver == nullptr) {
          return model->colUpper_[duplicateCol] == HIGHS_CONST_INF ||
                 implColUpper[duplicateCol] <
                     model->colUpper_[duplicateCol] -
                         options->primal_feasibility_tolerance;
        } else {
          return model->colUpper_[duplicateCol] == HIGHS_CONST_INF ||
                 implColUpper[duplicateCol] <=
                     model->colUpper_[duplicateCol] +
                         options->mip_feasibility_tolerance;
        }
      };

      auto duplicateColLowerInf = [&]() {
        if (!checkDuplicateColImplBounds) return false;
        if (mipsolver == nullptr) {
          return model->colLower_[duplicateCol] == -HIGHS_CONST_INF ||
                 implColLower[duplicateCol] >
                     model->colLower_[duplicateCol] +
                         options->primal_feasibility_tolerance;
        } else {
          return model->colLower_[duplicateCol] == -HIGHS_CONST_INF ||
                 implColLower[duplicateCol] >=
                     model->colLower_[duplicateCol] -
                         options->mip_feasibility_tolerance;
        }
      };

      // Now check the if the variable types rule out domination in one
      // direction and already skip the column if that rules out domination in
      // both directions due to the previous check on the objective.
      if (model->integrality_[i] == HighsVarType::INTEGER &&
          model->integrality_[parallelColCandidate] == HighsVarType::INTEGER) {
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
        colScale = std::round(scaleCand);
        assert(std::abs(colScale) >= 1.0);
        if (std::abs(colScale - scaleCand) > options->mip_epsilon) continue;

        // if the scale is larger than 1, duplicate column cannot compensate for
        // all values of scaled col due to integrality as the scaled column
        // moves on a grid of 1/scale.
        if (colScale != 1.0) checkDuplicateColImplBounds = false;
      } else if (model->integrality_[i] == HighsVarType::INTEGER) {
        col = i;
        duplicateCol = parallelColCandidate;
        colScale = colMax[duplicateCol].first / colMax[col].first;

        // as col is integral and dulicateCol is not col cannot compensate for
        // duplicate col
        checkColImplBounds = false;
      } else {
        col = parallelColCandidate;
        duplicateCol = i;
        colScale = colMax[duplicateCol].first / colMax[col].first;

        // as col might be integral and dulicateCol is not integral. In that
        // case col cannot compensate for duplicate col
        checkColImplBounds =
            model->integrality_[parallelColCandidate] != HighsVarType::INTEGER;
      }

      double objDiff = double(model->colCost_[col] * HighsCDouble(colScale) -
                              model->colCost_[duplicateCol]);
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
        if (colUpperInf() && model->colLower_[duplicateCol] != HIGHS_CONST_INF)
          reductionCase = kDominanceDuplicateColToLower;
        else if (duplicateColLowerInf() &&
                 (colScale < 0 || model->colUpper_[col] != HIGHS_CONST_INF) &&
                 (colScale > 0 || model->colLower_[col] != -HIGHS_CONST_INF))
          reductionCase =
              colScale > 0 ? kDominanceColToUpper : kDominanceColToLower;
        else
          continue;
      } else if (objDiff > options->dual_feasibility_tolerance) {
        // duplicate col is better than scaled col
        if (colLowerInf() && model->colUpper_[duplicateCol] != HIGHS_CONST_INF)
          reductionCase = kDominanceDuplicateColToUpper;
        else if (duplicateColUpperInf() &&
                 (colScale < 0 || model->colLower_[col] != -HIGHS_CONST_INF) &&
                 (colScale > 0 || model->colUpper_[col] != HIGHS_CONST_INF))
          reductionCase =
              colScale > 0 ? kDominanceColToLower : kDominanceColToUpper;
        else
          continue;
      } else {
        if (colUpperInf() && model->colLower_[duplicateCol] != -HIGHS_CONST_INF)
          reductionCase = kDominanceDuplicateColToLower;
        else if (colLowerInf() &&
                 model->colUpper_[duplicateCol] != HIGHS_CONST_INF)
          reductionCase = kDominanceDuplicateColToUpper;
        else if (duplicateColUpperInf() &&
                 (colScale < 0 || model->colLower_[col] != -HIGHS_CONST_INF) &&
                 (colScale > 0 || model->colUpper_[col] != HIGHS_CONST_INF))
          reductionCase =
              colScale > 0 ? kDominanceColToLower : kDominanceColToUpper;
        else if (duplicateColLowerInf() &&
                 (colScale < 0 || model->colUpper_[col] != HIGHS_CONST_INF) &&
                 (colScale > 0 || model->colLower_[col] != -HIGHS_CONST_INF))
          reductionCase =
              colScale > 0 ? kDominanceColToUpper : kDominanceColToLower;
      }
      double mergeLower = 0;
      double mergeUpper = 0;
      if (reductionCase == kMergeParallelCols) {
        if (colScale > 0) {
          mergeLower =
              model->colLower_[col] + colScale * model->colLower_[duplicateCol];
          mergeUpper =
              model->colUpper_[col] + colScale * model->colUpper_[duplicateCol];
        } else {
          mergeLower =
              model->colLower_[col] + colScale * model->colUpper_[duplicateCol];
          mergeUpper =
              model->colUpper_[col] + colScale * model->colLower_[duplicateCol];
        }
        if (model->integrality_[col] == HighsVarType::INTEGER) {
          // the only possible reduction if the column parallelism check
          // succeeds is to merge the two columns into one. If one column is
          // integral this means we have restrictions on integers and need to
          // check additional conditions to allow the merging of two integer
          // columns, or a continuous column and an integer.
          if (model->integrality_[duplicateCol] != HighsVarType::INTEGER) {
            // only one column is integral which cannot be duplicateCol due to
            // the way we assign the columns above
            if (std::abs(colScale * (model->colUpper_[duplicateCol] -
                                     model->colLower_[duplicateCol])) <
                1.0 - options->primal_feasibility_tolerance)
              continue;
          } else if (colScale > 1.0) {
            // round bounds to exact integer values to make sure they are not
            // wrongly truncated in conversions happening below
            mergeLower = std::round(mergeLower);
            mergeUpper = std::round(mergeUpper);

            // this should not happen, since this would allow domination and
            // would have been caught by the cases above
            assert(mergeLower != -HIGHS_CONST_INF);
            assert(mergeUpper != HIGHS_CONST_INF);

            HighsInt kMax = mergeUpper;
            bool representable = true;
            for (HighsInt k = mergeLower; k <= kMax; ++k) {
              // we loop over the domain of the merged variable to check whether
              // there exists a value for col and duplicateCol so that both are
              // within their bounds. since the merged column y is defined as y
              // = col + colScale * duplicateCol, we know that the value of col
              // can be computed as col = y - colScale * duplicateCol. Hence we
              // loop over the domain of col2 until we verify that a suitable
              // value of column 1 exists to yield the desired value for y.
              double mergeVal = mergeLower + k;
              HighsInt k2Max = model->colUpper_[duplicateCol];
              assert(k2Max == model->colUpper_[duplicateCol]);
              representable = false;
              for (HighsInt k2 = model->colLower_[duplicateCol]; k2 <= k2Max;
                   ++k2) {
                double colVal = mergeVal - colScale * k2;
                if (colVal >= model->colLower_[col] -
                                  options->primal_feasibility_tolerance &&
                    colVal <= model->colUpper_[col] +
                                  options->primal_feasibility_tolerance) {
                  representable = true;
                  break;
                }
              }

              if (!representable) break;
            }

            if (!representable) continue;
          }
        }
      }

      bool parallel = true;
      // now check whether the coefficients are actually parallel
      for (const HighsSliceNonzero& colNz : getColumnVector(col)) {
        HighsInt duplicateColRowPos = findNonzero(colNz.index(), duplicateCol);
        if (duplicateColRowPos == -1) {
          parallel = false;
          break;
        }

        double difference = std::abs(
            double(Avalue[duplicateColRowPos] - colScale * colNz.value()));
        if (difference > options->small_matrix_value) {
          parallel = false;
          break;
        }
      }

      if (!parallel) continue;

      switch (reductionCase) {
        case kDominanceDuplicateColToLower:
          delCol = duplicateCol;
          if (colsize[duplicateCol] == 1) {
            HighsInt row = Arow[colhead[duplicateCol]];
            numRowSingletons[row] -= 1;
          }
          fixColToLower(postSolveStack, duplicateCol);
          break;
        case kDominanceDuplicateColToUpper:
          delCol = duplicateCol;
          if (colsize[duplicateCol] == 1) {
            HighsInt row = Arow[colhead[duplicateCol]];
            numRowSingletons[row] -= 1;
          }
          fixColToUpper(postSolveStack, duplicateCol);
          break;
        case kDominanceColToLower:
          delCol = col;
          if (colsize[col] == 1) {
            HighsInt row = Arow[colhead[col]];
            numRowSingletons[row] -= 1;
          }
          fixColToLower(postSolveStack, col);
          break;
        case kDominanceColToUpper:
          delCol = col;
          if (colsize[col] == 1) {
            HighsInt row = Arow[colhead[col]];
            numRowSingletons[row] -= 1;
          }
          fixColToUpper(postSolveStack, col);
          break;
        case kMergeParallelCols:
          postSolveStack.duplicateColumn(
              colScale, model->colLower_[col], model->colUpper_[col],
              model->colLower_[duplicateCol], model->colUpper_[duplicateCol],
              col, duplicateCol,
              model->integrality_[col] == HighsVarType::INTEGER,
              model->integrality_[duplicateCol] == HighsVarType::INTEGER);

          markChangedCol(col);
          if (colsize[duplicateCol] == 1) {
            HighsInt row = Arow[colhead[duplicateCol]];
            numRowSingletons[row] -= 1;
          }

          // by updating the bounds properly, the unlink calls will update the
          // implied row upper bounds to the correct values. For finite bounds
          // simply setting the bounds of duplicate col to zero suffices. For
          // infinite bounds we need to make sure the counters for the number of
          // infinite bounds that contribute to the implied row bounds are
          // updated correctly and that all finite contributions are removed.
          if (colScale > 0) {
            if (mergeUpper == HIGHS_CONST_INF &&
                model->colUpper_[col] != HIGHS_CONST_INF)
              model->colUpper_[duplicateCol] = model->colUpper_[col] / colScale;
            else
              model->colUpper_[duplicateCol] = 0;

            if (mergeLower == -HIGHS_CONST_INF &&
                model->colLower_[col] != -HIGHS_CONST_INF)
              // make sure that upon removal of the duplicate column the finite
              // contribution of col's lower bound is removed and the infinite
              // contribution of duplicateCol is retained
              model->colLower_[duplicateCol] = model->colLower_[col] / colScale;
            else
              model->colLower_[duplicateCol] = 0;
          } else {
            if (mergeUpper == HIGHS_CONST_INF &&
                model->colUpper_[col] != HIGHS_CONST_INF)
              model->colLower_[duplicateCol] = model->colUpper_[col] / colScale;
            else
              model->colLower_[duplicateCol] = 0;

            if (mergeLower == -HIGHS_CONST_INF &&
                model->colLower_[col] != -HIGHS_CONST_INF)
              // make sure that upon removal of the duplicate column the finite
              // contribution of col's lower bound is removed and the infinite
              // contribution of duplicateCol is retained
              model->colUpper_[duplicateCol] = model->colLower_[col] / colScale;
            else
              model->colUpper_[duplicateCol] = 0;
          }

          model->colLower_[col] = mergeLower;
          model->colUpper_[col] = mergeUpper;

          // mark duplicate column as deleted
          markColDeleted(duplicateCol);
          // remove all nonzeros of duplicateCol
          for (HighsInt coliter = colhead[duplicateCol]; coliter != -1;) {
            assert(Acol[coliter] == duplicateCol);

            HighsInt colpos = coliter;
            HighsInt colrow = Arow[coliter];
            coliter = Anext[coliter];

            unlink(colpos);

            if (model->rowLower_[colrow] == model->rowUpper_[colrow] &&
                eqiters[colrow] != equations.end() &&
                eqiters[colrow]->first != rowsize[colrow]) {
              // if that is the case reinsert it into the equation set that is
              // ordered by sparsity
              equations.erase(eqiters[colrow]);
              eqiters[colrow] =
                  equations.emplace(rowsize[colrow], colrow).first;
            }
          }
          // set cost to zero
          model->colCost_[duplicateCol] = 0;
          delCol = duplicateCol;

          // remove implied bounds, since they might in general not be valid
          // anymore
          if (colLowerSource[col] != -1)
            changeImplColLower(col, -HIGHS_CONST_INF, -1);

          if (colUpperSource[col] != -1)
            changeImplColUpper(col, HIGHS_CONST_INF, -1);

          break;
      }

      break;
    }

    if (delCol != -1) {
      if (delCol != i) buckets.erase(last);

      // we could have new row singletons since a column was removed. Remove
      // those rows immediately
      HPRESOLVE_CHECKED_CALL(checkLimits(postSolveStack));
      HPRESOLVE_CHECKED_CALL(removeRowSingletons(postSolveStack));
    } else {
      buckets.emplace_hint(last, colHashes[i], i);
    }
  }

  buckets.clear();

  for (HighsInt i = 0; i != model->numRow_; ++i) {
    if (rowDeleted[i]) continue;
    if (rowsize[i] <= 1 ||
        (rowsize[i] == 2 && model->rowLower_[i] == model->rowUpper_[i])) {
      HPRESOLVE_CHECKED_CALL(rowPresolve(postSolveStack, i));
      ++numRowBuckets;
      continue;
    }
    auto it = buckets.find(rowHashes[i]);
    decltype(it) last = it;

    const HighsInt* numSingletonPtr = numRowSingletons.find(i);
    HighsInt numSingleton = numSingletonPtr ? *numSingletonPtr : 0;

#if !ENABLE_SPARSIFY_FOR_LP
    if (mipsolver == nullptr && numSingleton != 0) continue;
#endif
    HighsInt delRow = -1;
    if (it == buckets.end())
      ++numRowBuckets;
    else
      storeRow(i);
    while (it != buckets.end() && it->first == rowHashes[i]) {
      HighsInt parallelRowCand = it->second;
      last = it++;

      numSingletonPtr = numRowSingletons.find(parallelRowCand);
      const HighsInt numSingletonCandidate =
          numSingletonPtr ? *numSingletonPtr : 0;
#if !ENABLE_SPARSIFY_FOR_LP
      if (mipsolver == nullptr && numSingletonCandidate != 0) continue;
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
        if ((numSingleton != 0 || model->rowLower_[i] != model->rowUpper_[i]) &&
            (numSingletonCandidate != 0 ||
             model->rowLower_[parallelRowCand] !=
                 model->rowUpper_[parallelRowCand]))
          continue;
      } else if (numSingletonCandidate != numSingleton) {
        // if only one of the two constraints has an extra singleton,
        // we require at least one of the constraints to be an equation
        // if that is the case we can add that equation to the other row
        // and will make it into either a row singleton or a doubleton equation
        // which is removed afterwards
        if (model->rowLower_[i] != model->rowUpper_[i] &&
            model->rowLower_[parallelRowCand] !=
                model->rowUpper_[parallelRowCand])
          continue;
      }

      double rowScale = rowMax[parallelRowCand].first / rowMax[i].first;
      // check parallel case
      bool parallel = true;
      for (const HighsSliceNonzero& rowNz : getStoredRow()) {
        if (colsize[rowNz.index()] == 1)  // skip singletons
          continue;
        HighsInt nzPos = findNonzero(parallelRowCand, rowNz.index());
        if (nzPos == -1) {
          parallel = false;
          break;
        }

        if (std::abs(double(Avalue[nzPos] -
                            HighsCDouble(rowScale) * rowNz.value())) >
            options->small_matrix_value) {
          parallel = false;
          break;
        }
      }
      if (!parallel) continue;

      if (numSingleton == 0 && numSingletonCandidate == 0) {
        bool rowLowerTightened = false;
        bool rowUpperTightened = false;
        double newUpper;
        double newLower;
        if (rowScale > 0) {
          newUpper = model->rowUpper_[i] * rowScale;
          newLower = model->rowLower_[i] * rowScale;
        } else {
          newLower = model->rowUpper_[i] * rowScale;
          newUpper = model->rowLower_[i] * rowScale;
        }

        if (newUpper < model->rowUpper_[parallelRowCand]) {
          if (newUpper < model->rowLower_[parallelRowCand] -
                             options->primal_feasibility_tolerance)
            return Result::PrimalInfeasible;

          if (newUpper <= model->rowLower_[parallelRowCand] +
                              options->primal_feasibility_tolerance)
            newUpper = model->rowLower_[parallelRowCand];

          if (newUpper < model->rowUpper_[parallelRowCand]) {
            rowUpperTightened = true;
            if (rowScale > 0) {
              double tmp = rowDualUpper[i] / rowScale;
              rowDualUpper[i] = rowDualUpper[parallelRowCand] * rowScale;
              rowDualUpper[parallelRowCand] = tmp;
            } else {
              double tmp = rowDualLower[i] / rowScale;
              rowDualLower[i] = rowDualUpper[parallelRowCand] * rowScale;
              rowDualUpper[parallelRowCand] = tmp;
            }

            model->rowUpper_[parallelRowCand] = newUpper;
          }
        }

        if (newLower > model->rowLower_[parallelRowCand]) {
          if (newLower > model->rowUpper_[parallelRowCand] +
                             options->primal_feasibility_tolerance)
            return Result::PrimalInfeasible;

          if (newLower >= model->rowUpper_[parallelRowCand] -
                              options->primal_feasibility_tolerance)
            newLower = model->rowUpper_[parallelRowCand];

          if (newLower > model->rowLower_[parallelRowCand]) {
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
              double tmp = rowDualLower[i] / rowScale;
              rowDualLower[i] = rowDualLower[parallelRowCand] * rowScale;
              rowDualLower[parallelRowCand] = tmp;
            } else {
              double tmp = rowDualUpper[i] / rowScale;
              rowDualUpper[i] = rowDualLower[parallelRowCand] * rowScale;
              rowDualLower[parallelRowCand] = tmp;
            }

            model->rowLower_[parallelRowCand] = newLower;
          }
        }
        if (rowDualLowerSource[parallelRowCand] != -1)
          changeImplRowDualLower(parallelRowCand, -HIGHS_CONST_INF, -1);
        if (rowDualUpperSource[parallelRowCand] != -1)
          changeImplRowDualUpper(parallelRowCand, HIGHS_CONST_INF, -1);

        postSolveStack.duplicateRow(parallelRowCand, rowUpperTightened,
                                    rowLowerTightened, i, rowScale);
        delRow = i;
        markRowDeleted(i);
        for (HighsInt rowiter : rowpositions) unlink(rowiter);
        break;
      } else if (model->rowLower_[i] == model->rowUpper_[i]) {
        // row i is equation and parallel (except for singletons)
        // add to the row parallelRowCand
        // printf(
        //    "nearly parallel case with %" HIGHSINT_FORMAT " singletons in eq
        //    row and %" HIGHSINT_FORMAT " " "singletons in other row(eq=%"
        //    HIGHSINT_FORMAT ")\n", numSingleton, numSingletonCandidate,
        //    model->rowLower_[parallelRowCand] ==
        //        model->rowUpper_[parallelRowCand]);
        postSolveStack.equalityRowAddition(parallelRowCand, i, -rowScale,
                                           getStoredRow());
        for (const HighsSliceNonzero& rowNz : getStoredRow()) {
          HighsInt pos = findNonzero(parallelRowCand, rowNz.index());
          if (pos != -1)
            unlink(pos);  // all common nonzeros are cancelled, as the rows are
                          // parallel
          else            // might introduce a singleton
            addToMatrix(parallelRowCand, rowNz.index(),
                        -rowScale * rowNz.value());
        }

        if (model->rowUpper_[parallelRowCand] != HIGHS_CONST_INF)
          model->rowUpper_[parallelRowCand] =
              double(model->rowUpper_[parallelRowCand] -
                     HighsCDouble(rowScale) * model->rowUpper_[i]);
        if (model->rowLower_[parallelRowCand] != -HIGHS_CONST_INF)
          model->rowLower_[parallelRowCand] =
              double(model->rowLower_[parallelRowCand] -
                     HighsCDouble(rowScale) * model->rowUpper_[i]);

        // parallelRowCand is now a singleton row, doubleton equation, or a row
        // that contains only singletons and we let the normal row presolve
        // handle the cases
        HPRESOLVE_CHECKED_CALL(rowPresolve(postSolveStack, parallelRowCand));
        delRow = parallelRowCand;
      } else if (model->rowLower_[parallelRowCand] ==
                 model->rowUpper_[parallelRowCand]) {
        // printf(
        //    "nearly parallel case with %" HIGHSINT_FORMAT " singletons in eq
        //    row and %" HIGHSINT_FORMAT " " "singletons in other inequality
        //    row\n", numSingletonCandidate, numSingleton);
        // the row parallelRowCand is an equation; add it to the other row
        double scale = -rowMax[i].first / rowMax[parallelRowCand].first;
        postSolveStack.equalityRowAddition(i, parallelRowCand, scale,
                                           getRowVector(parallelRowCand));
        for (const HighsSliceNonzero& rowNz : getRowVector(parallelRowCand)) {
          HighsInt pos = findNonzero(i, rowNz.index());
          if (pos != -1)
            unlink(pos);  // all common nonzeros are cancelled, as the rows are
                          // parallel
          else            // might introduce a singleton
            addToMatrix(i, rowNz.index(), scale * rowNz.value());
        }

        if (model->rowUpper_[i] != HIGHS_CONST_INF)
          model->rowUpper_[i] =
              double(model->rowUpper_[i] +
                     HighsCDouble(scale) * model->rowUpper_[parallelRowCand]);
        if (model->rowLower_[i] != -HIGHS_CONST_INF)
          model->rowLower_[i] =
              double(model->rowLower_[i] +
                     HighsCDouble(scale) * model->rowUpper_[parallelRowCand]);

        HPRESOLVE_CHECKED_CALL(rowPresolve(postSolveStack, i));
        delRow = i;
      } else {
        assert(numSingleton == 1);
        assert(numSingletonCandidate == 1);

        double rowUpper;
        double rowLower;
        if (rowScale > 0) {
          rowUpper = model->rowUpper_[i] * rowScale;
          rowLower = model->rowLower_[i] * rowScale;
        } else {
          rowLower = model->rowUpper_[i] * rowScale;
          rowUpper = model->rowLower_[i] * rowScale;
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

      HPRESOLVE_CHECKED_CALL(checkLimits(postSolveStack));
    } else
      buckets.emplace_hint(last, rowHashes[i], i);
  }

  return Result::Ok;
}

void HPresolve::setRelaxedImpliedBounds() {
  double hugeBound = options->primal_feasibility_tolerance / HIGHS_CONST_TINY;
  for (HighsInt i = 0; i != model->numCol_; ++i) {
    if (model->colLower_[i] >= implColLower[i] &&
        model->colUpper_[i] <= implColUpper[i])
      continue;

    if (std::abs(implColLower[i]) <= hugeBound) {
      // if the bound is derived from a small nonzero value
      // then we want to increase the margin so that we make sure
      // the row it was derived from is violated if the column sits
      // at this relaxed bound in the final solution.
      HighsInt nzPos = findNonzero(colLowerSource[i], i);
      double boundRelax = 128.0 * options->primal_feasibility_tolerance /
                          std::min(1.0, std::abs(Avalue[nzPos]));
      double newLb =
          implColLower[i] -
          std::max(boundRelax, options->primal_feasibility_tolerance *
                                   std::abs(implColLower[i]));
      if (newLb > model->colLower_[i]) model->colLower_[i] = newLb;
    }

    if (std::abs(implColUpper[i]) <= hugeBound) {
      HighsInt nzPos = findNonzero(colUpperSource[i], i);
      double boundRelax = 128.0 * options->primal_feasibility_tolerance /
                          std::min(1.0, std::abs(Avalue[nzPos]));
      double newUb =
          implColUpper[i] +
          std::max(boundRelax, options->primal_feasibility_tolerance *
                                   std::abs(implColUpper[i]));
      if (newUb < model->colUpper_[i]) model->colUpper_[i] = newUb;
    }
  }
}

void HPresolve::debug(const HighsLp& lp, const HighsOptions& options) {
  HighsSolution reducedsol;
  HighsBasis reducedbasis;

  HighsSolution sol;
  HighsBasis basis;

  HighsLp model = lp;
  model.integrality_.assign(lp.numCol_, HighsVarType::CONTINUOUS);

  HighsPostsolveStack postSolveStack;
  postSolveStack.initializeIndexMaps(lp.numRow_, lp.numCol_);
  {
    HPresolve presolve;
    presolve.setInput(model, options);
    // presolve.setReductionLimit(1622017);
    if (presolve.run(postSolveStack) != HighsModelStatus::NOTSET) return;
    Highs highs;
    highs.passModel(model);
    highs.passHighsOptions(options);
    highs.setHighsOptionValue("presolve", "off");
    highs.run();
    if (highs.getModelStatus(true) != HighsModelStatus::OPTIMAL) return;
    reducedsol = highs.getSolution();
    reducedbasis = highs.getBasis();
  }
  model = lp;
  sol = reducedsol;
  basis = reducedbasis;
  postSolveStack.undo(options, sol, basis);
  refineBasis(lp, sol, basis);
  calculateRowValues(model, sol);
#if 0
  Highs highs;
  highs.passModel(model);
  highs.passHighsOptions(options);
  highs.setSolution(sol);
  highs.setBasis(basis);
  highs.run();
  return;
#endif
  std::vector<HighsInt> flagCol(lp.numCol_, 1);
  std::vector<HighsInt> flagRow(lp.numRow_, 1);
  std::vector<HighsInt> Aend;
  std::vector<HighsInt> ARstart;
  std::vector<HighsInt> ARindex;
  std::vector<double> ARvalue;
  dev_kkt_check::KktInfo kktinfo = dev_kkt_check::initInfo();
  Aend.assign(model.Astart_.begin() + 1, model.Astart_.end());
  highsSparseTranspose(model.numRow_, model.numCol_, model.Astart_,
                       model.Aindex_, model.Avalue_, ARstart, ARindex, ARvalue);
  dev_kkt_check::State state(
      model.numCol_, model.numRow_, model.Astart_, Aend, model.Aindex_,
      model.Avalue_, ARstart, ARindex, ARvalue, model.colCost_, model.colLower_,
      model.colUpper_, model.rowLower_, model.rowUpper_, flagCol, flagRow,
      sol.col_value, sol.col_dual, sol.row_value, sol.row_dual,
      basis.col_status, basis.row_status);
  bool checkResult = dev_kkt_check::checkKkt(state, kktinfo);
  if (checkResult && kktinfo.pass_bfs) {
    printf("kkt check of postsolved solution and basis passed\n");
    return;
  }
  size_t good = postSolveStack.numReductions();
  size_t bad = 0;
  size_t reductionLim = (good + bad) / 2;

  // good = 1734357, bad = 1734289;
  // good = 1050606, bad = 1050605;
  // good = 1811527, bad = 1811526;
  // reductionLim = bad;
  do {
    model = lp;
    model.integrality_.assign(lp.numCol_, HighsVarType::CONTINUOUS);

    {
      HPresolve presolve;
      presolve.setInput(model, options);
      presolve.computeIntermediateMatrix(flagRow, flagCol, reductionLim);
    }
#if 1
    model = lp;
    model.integrality_.assign(lp.numCol_, HighsVarType::CONTINUOUS);
    HPresolve presolve;
    presolve.setInput(model, options);
    HighsPostsolveStack tmp;
    tmp.initializeIndexMaps(model.numRow_, model.numCol_);
    presolve.setReductionLimit(reductionLim);
    presolve.run(tmp);

    sol = reducedsol;
    basis = reducedbasis;
    postSolveStack.undoUntil(options, flagRow, flagCol, sol, basis,
                             tmp.numReductions());

    HighsBasis tmpBasis;
    HighsSolution tmpSol;
    tmpBasis.col_status.resize(model.numCol_);
    tmpSol.col_dual.resize(model.numCol_);
    tmpSol.col_value.resize(model.numCol_);
    for (HighsInt i = 0; i != model.numCol_; ++i) {
      tmpSol.col_dual[i] = sol.col_dual[tmp.getOrigColIndex(i)];
      tmpSol.col_value[i] = sol.col_value[tmp.getOrigColIndex(i)];
      tmpBasis.col_status[i] = basis.col_status[tmp.getOrigColIndex(i)];
    }

    tmpBasis.row_status.resize(model.numRow_);
    tmpSol.row_dual.resize(model.numRow_);
    for (HighsInt i = 0; i != model.numRow_; ++i) {
      tmpSol.row_dual[i] = sol.row_dual[tmp.getOrigRowIndex(i)];
      tmpBasis.row_status[i] = basis.row_status[tmp.getOrigRowIndex(i)];
    }
    tmpSol.row_value.resize(model.numRow_);
    calculateRowValues(model, sol);
    tmpBasis.valid_ = true;
    refineBasis(model, tmpSol, tmpBasis);
    Highs highs;
    highs.passHighsOptions(options);
    highs.passModel(model);
    highs.setBasis(tmpBasis);
    // highs.writeModel("model.mps");
    // highs.writeBasis("bad.bas");
    highs.run();
    printf("simplex iterations with postsolved basis: %" HIGHSINT_FORMAT "\n",
           highs.getSimplexIterationCount());
    checkResult = highs.getSimplexIterationCount() == 0;
#else

    if (reductionLim == good) break;

    Aend.assign(model.Astart_.begin() + 1, model.Astart_.end());
    highsSparseTranspose(model.numRow_, model.numCol_, model.Astart_,
                         model.Aindex_, model.Avalue_, ARstart, ARindex,
                         ARvalue);
    sol = reducedsol;
    basis = reducedbasis;
    postSolveStack.undoUntil(options, flagRow, flagCol, sol, basis,
                             reductionLim);

    calculateRowValues(model, sol);
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

HPresolve::Result HPresolve::sparsify(HighsPostsolveStack& postSolveStack) {
  std::vector<HighsPostsolveStack::Nonzero> sparsifyRows;
  HPRESOLVE_CHECKED_CALL(removeRowSingletons(postSolveStack));
  HPRESOLVE_CHECKED_CALL(removeDoubletonEquations(postSolveStack));
  std::vector<HighsInt> tmpEquations;
  tmpEquations.reserve(equations.size());

  const double minNonzeroVal = std::sqrt(options->primal_feasibility_tolerance);

  for (const auto& eq : equations) tmpEquations.emplace_back(eq.second);
  for (HighsInt eqrow : tmpEquations) {
    if (rowDeleted[eqrow]) continue;

    assert(!rowDeleted[eqrow]);
    assert(model->rowLower_[eqrow] == model->rowUpper_[eqrow]);

    storeRow(eqrow);

    HighsInt secondSparsestColumn = -1;
    HighsInt sparsestCol = Acol[rowpositions[0]];
    HighsInt sparsestColLen = HIGHS_CONST_I_INF;
    for (size_t i = 1; i < rowpositions.size(); ++i) {
      HighsInt col = Acol[rowpositions[i]];
      if (colsize[col] < sparsestColLen) {
        sparsestColLen = colsize[col];
        secondSparsestColumn = sparsestCol;
        sparsestCol = col;
      }
    }

    if (colsize[secondSparsestColumn] < colsize[sparsestCol])
      std::swap(sparsestCol, secondSparsestColumn);

    assert(sparsestCol != -1 && secondSparsestColumn != -1);

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
            if (model->integrality_[nonzero.index()] == HighsVarType::INTEGER &&
                model->colUpper_[nonzero.index()] -
                        model->colLower_[nonzero.index()] >
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
          // there already is a scale that is very close and could produces
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

    if (model->integrality_[sparsestCol] != HighsVarType::INTEGER ||
        (model->colUpper_[sparsestCol] - model->colLower_[sparsestCol]) < 1.5) {
      // now check for rows which do not contain the sparsest column but all
      // other columns by scanning the second sparsest column
      for (const HighsSliceNonzero& colNz :
           getColumnVector(secondSparsestColumn)) {
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
          if (nonzero.index() == secondSparsestColumn) {
            candRowVal = colNz.value();
          } else {
            HighsInt nzPos = findNonzero(candRow, nonzero.index());
            if (nzPos == -1) {
              // we already have a miss for the sparsest column, so with another
              // one we want to skip the row
              skip = true;
              break;
            }

            candRowVal = Avalue[nzPos];
          }

          double scale = -candRowVal / nonzero.value();
          if (std::abs(scale) > 1e3) continue;

          double scaleTolerance = minNonzeroVal / std::abs(nonzero.value());
          auto it = possibleScales.lower_bound(scale - scaleTolerance);
          if (it != possibleScales.end() &&
              std::abs(it->first - scale) <= scaleTolerance) {
            // there already is a scale that is very close and could produces
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

    postSolveStack.equalityRowAdditions(eqrow, getStoredRow(), sparsifyRows);
    double rhs = model->rowLower_[eqrow];
    for (const auto& sparsifyRow : sparsifyRows) {
      HighsInt row = sparsifyRow.index;
      double scale = sparsifyRow.value;

      if (model->rowLower_[row] != -HIGHS_CONST_INF)
        model->rowLower_[row] += scale * rhs;

      if (model->rowUpper_[row] != HIGHS_CONST_INF)
        model->rowUpper_[row] += scale * rhs;

      for (HighsInt pos : rowpositions)
        addToMatrix(row, Acol[pos], scale * Avalue[pos]);

      if (model->rowLower_[row] == model->rowUpper_[row] &&
          eqiters[row] != equations.end() &&
          eqiters[row]->first != rowsize[row]) {
        // if that is the case reinsert it into the equation set that is ordered
        // by sparsity
        equations.erase(eqiters[row]);
        eqiters[row] = equations.emplace(rowsize[row], row).first;
      }
    }

    HPRESOLVE_CHECKED_CALL(checkLimits(postSolveStack));
    HPRESOLVE_CHECKED_CALL(removeRowSingletons(postSolveStack));
    HPRESOLVE_CHECKED_CALL(removeDoubletonEquations(postSolveStack));
  }

  return Result::Ok;
}

}  // namespace presolve
