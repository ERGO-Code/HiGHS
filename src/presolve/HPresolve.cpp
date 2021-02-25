#include "presolve/HPresolve.h"

#include <algorithm>
#include <atomic>

#include "presolve/HighsPostsolveStack.h"
#include "util/HighsSplay.h"

#define HPRESOLVE_CHECKED_CALL(presolveCall)                          \
  do {                                                                \
    HPresolve::Result __result = presolveCall;                        \
    if (__result != presolve::HPresolve::Result::Ok) return __result; \
  } while (0)

namespace presolve {

#ifndef NDEBUG
void HPresolve::debugPrintRow(int row) {
  printf("(row %d) %.15g (impl: %.15g) <= ", row, model->rowLower_[row],
         impliedRowBounds.getSumLower(row));

  for (const HighsSliceNonzero& nonzero : getRowVector(row)) {
    // for (int rowiter = rowhead[row]; rowiter != -1; rowiter =
    // ARnext[rowiter]) {
    char colchar = model->integrality_[nonzero.index()] == HighsVarType::INTEGER
                       ? 'y'
                       : 'x';
    char signchar = nonzero.value() < 0 ? '-' : '+';
    printf("%c%g %c%d ", signchar, std::abs(nonzero.value()), colchar,
           nonzero.index());
  }

  printf("<= %.15g (impl: %.15g)\n", model->rowUpper_[row],
         impliedRowBounds.getSumUpper(row));
}
#endif
#if 0
void HPresolve::debugPrintSubMatrix(int row, int col) {
  printf("submatrix for col %d and row %d:\n", col, row);
  debugPrintRow(row);
  for (int coliter = colhead[col]; coliter != -1; coliter = Anext[coliter]) {
    int r = Arow[coliter];

    if (r == row) continue;

    printf("(row %d) %g <= ... ", r, rowLower[r]);

    loopRow(row, [&](int rowiter) {
      // for (int rowiter = rowhead[row]; rowiter != -1; rowiter =
      // ARnext[rowiter]) {
      auto it = entries.find(std::make_pair(r, Acol[rowiter]));
      if (it != entries.end()) {
        assert(Acol[it->second] == Acol[rowiter]);
        char colchar =
            integrality[Acol[it->second]] == HighsVarType::INTEGER ? 'y' : 'x';
        char signchar = Avalue[it->second] < 0 ? '-' : '+';
        printf("%c%g %c%d ", signchar, std::abs(Avalue[it->second]), colchar,
               Acol[it->second]);
      }

      return false;
    });

    printf(" ... <= %g\n", rowUpper[r]);
  }
}
#endif

void HPresolve::setInput(HighsLp& model_, const HighsOptions& options_) {
  model = &model_;
  options = &options_;

  colhead.resize(model->numCol_, -1);
  colsize.resize(model->numCol_);
  col_numerics_threshold.resize(model->numCol_);
  colLowerSource.resize(model->numCol_, -1);
  colUpperSource.resize(model->numCol_, -1);
  implColLower.resize(model->numCol_, -HIGHS_CONST_INF);
  implColUpper.resize(model->numCol_, HIGHS_CONST_INF);

  rowroot.resize(model->numRow_, -1);
  rowsize.resize(model->numRow_);

  impliedRowBounds.setBoundArrays(
      model->colLower_.data(), model->colUpper_.data(), implColLower.data(),
      implColUpper.data(), colLowerSource.data(), colUpperSource.data());
  impliedRowBounds.setNumSums(model->numRow_);

  rowDualLower.resize(model->numRow_, -HIGHS_CONST_INF);
  rowDualUpper.resize(model->numRow_, HIGHS_CONST_INF);
  implRowDualLower.resize(model->numRow_, -HIGHS_CONST_INF);
  implRowDualUpper.resize(model->numRow_, HIGHS_CONST_INF);
  rowDualUpperSource.resize(model->numRow_, -1);
  rowDualLowerSource.resize(model->numRow_, -1);

  for (int i = 0; i != model->numRow_; ++i) {
    if (model->rowLower_[i] == -HIGHS_CONST_INF) rowDualLower[i] = 0;
    if (model->rowUpper_[i] == HIGHS_CONST_INF) rowDualUpper[i] = 0;
  }

  impliedDualRowBounds.setBoundArrays(
      rowDualLower.data(), rowDualUpper.data(), implRowDualLower.data(),
      implRowDualUpper.data(), rowDualLowerSource.data(),
      rowDualUpperSource.data());
  impliedDualRowBounds.setNumSums(model->numCol_);

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

bool HPresolve::isRowInteger(int row, double scale) const {
  for (const HighsSliceNonzero& nz : getRowVector(row)) {
    // first check coefficient, as this requires no additional memory access
    // into the model->integrality array
    double val = nz.value() * scale;
    if (std::abs(val - std::round(val)) > options->mip_epsilon) return false;

    if (model->integrality_[nz.index()] != HighsVarType::INTEGER) return false;
  }

  return true;
}

bool HPresolve::isRowIntegral(int row, double scale) const {
  for (const HighsSliceNonzero& nz : getRowVector(row)) {
    // first check coefficient, as this requires no additional memory access
    // into the model->integrality array
    double val = nz.value() * scale;
    if (std::abs(val - std::round(val)) > options->mip_epsilon) return false;

    if (model->integrality_[nz.index()] == HighsVarType::CONTINUOUS)
      return false;
  }

  return true;
}

bool HPresolve::isLowerImplied(int col) const {
  return (model->colLower_[col] == -HIGHS_CONST_INF ||
          implColLower[col] >=
              model->colLower_[col] - options->primal_feasibility_tolerance);
}

bool HPresolve::isUpperImplied(int col) const {
  return (model->colUpper_[col] == HIGHS_CONST_INF ||
          implColUpper[col] <=
              model->colUpper_[col] + options->primal_feasibility_tolerance);
}

bool HPresolve::isImpliedFree(int col) const {
  return (model->colLower_[col] == -HIGHS_CONST_INF ||
          implColLower[col] >=
              model->colLower_[col] - options->primal_feasibility_tolerance) &&
         (model->colUpper_[col] == HIGHS_CONST_INF ||
          implColUpper[col] <=
              model->colUpper_[col] + options->primal_feasibility_tolerance);
}

bool HPresolve::isDualImpliedFree(int row) const {
  return (rowDualLower[row] == -HIGHS_CONST_INF ||
          implRowDualLower[row] >=
              rowDualLower[row] - options->dual_feasibility_tolerance) &&
         (rowDualUpper[row] == HIGHS_CONST_INF ||
          implRowDualUpper[row] <=
              rowDualUpper[row] + options->dual_feasibility_tolerance);
}

bool HPresolve::isImpliedIntegral(int col) {
  bool haveEquation = false;
  for (const HighsSliceNonzero& nz : getColumnVector(col)) {
    double rowLower =
        implRowDualLower[nz.index()] > options->dual_feasibility_tolerance
            ? model->rowUpper_[nz.index()]
            : model->rowLower_[nz.index()];

    double rowUpper =
        implRowDualUpper[nz.index()] < -options->dual_feasibility_tolerance
            ? model->rowLower_[nz.index()]
            : model->rowUpper_[nz.index()];

    if (rowUpper == rowLower) {
      haveEquation = true;
      double scale = 1.0 / nz.value();
      if (!isRowInteger(nz.index(), scale)) continue;

      double rhs = model->rowLower_[nz.index()] * scale;

      if (std::abs(rhs - std::round(rhs)) >
          options->mip_feasibility_tolerance) {
        // todo infeasible
      }

      return true;
    }
  }

  if (haveEquation) return false;

  for (const HighsSliceNonzero& nz : getColumnVector(col)) {
    double scale = 1.0 / nz.value();
    if (!isRowInteger(nz.index(), scale)) return false;
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

void HPresolve::link(int pos) {
  Anext[pos] = colhead[Acol[pos]];
  Aprev[pos] = -1;
  colhead[Acol[pos]] = pos;
  if (Anext[pos] != -1) Aprev[Anext[pos]] = pos;

  ++colsize[Acol[pos]];
  col_numerics_threshold[Acol[pos]] =
      std::max(options->presolve_pivot_threshold * std::abs(Avalue[pos]),
               col_numerics_threshold[Acol[pos]]);

  auto get_row_left = [&](int pos) -> int& { return ARleft[pos]; };
  auto get_row_right = [&](int pos) -> int& { return ARright[pos]; };
  auto get_row_key = [&](int pos) { return Acol[pos]; };
  highs_splay_link(pos, rowroot[Arow[pos]], get_row_left, get_row_right,
                   get_row_key);

  impliedRowBounds.add(Arow[pos], Acol[pos], Avalue[pos]);
  impliedDualRowBounds.add(Acol[pos], Arow[pos], Avalue[pos]);
  ++rowsize[Arow[pos]];
}

void HPresolve::unlink(int pos) {
  int next = Anext[pos];
  int prev = Aprev[pos];

  if (next != -1) Aprev[next] = prev;

  if (prev != -1)
    Anext[prev] = next;
  else
    colhead[Acol[pos]] = next;
  --colsize[Acol[pos]];

  if (!colDeleted[Acol[pos]]) {
    if (colsize[Acol[pos]] == 1) singletonColumns.push_back(Acol[pos]);

    impliedDualRowBounds.remove(Acol[pos], Arow[pos], Avalue[pos]);
    if (colUpperSource[Acol[pos]] == Arow[pos])
      changeImplColUpper(Acol[pos], HIGHS_CONST_INF, -1);

    if (colLowerSource[Acol[pos]] == Arow[pos])
      changeImplColLower(Acol[pos], -HIGHS_CONST_INF, -1);
  }

  auto get_row_left = [&](int pos) -> int& { return ARleft[pos]; };
  auto get_row_right = [&](int pos) -> int& { return ARright[pos]; };
  auto get_row_key = [&](int pos) { return Acol[pos]; };
  highs_splay_unlink(pos, rowroot[Arow[pos]], get_row_left, get_row_right,
                     get_row_key);
  --rowsize[Arow[pos]];

  if (!rowDeleted[Arow[pos]]) {
    if (rowsize[Arow[pos]] == 1) singletonRows.push_back(Arow[pos]);
    impliedRowBounds.remove(Arow[pos], Acol[pos], Avalue[pos]);

    if (rowDualUpperSource[Arow[pos]] == Acol[pos])
      changeImplRowDualUpper(Arow[pos], HIGHS_CONST_INF, -1);

    if (rowDualLowerSource[Arow[pos]] == Acol[pos])
      changeImplRowDualLower(Arow[pos], -HIGHS_CONST_INF, -1);
  }

  Avalue[pos] = 0;

  freeslots.push_back(pos);
}

void HPresolve::markChangedRow(int row) {
  if (!changedRowFlag[row]) {
    changedRowIndices.push_back(row);
    changedRowFlag[row] = true;
  }
}

void HPresolve::markChangedCol(int col) {
  if (!changedColFlag[col]) {
    changedColIndices.push_back(col);
    changedColFlag[col] = true;
  }
}

double HPresolve::getMaxAbsColVal(int col) const {
  double maxVal = 0.0;

  for (const auto& nz : getColumnVector(col))
    maxVal = std::max(std::abs(nz.value()), maxVal);

  return maxVal;
}

double HPresolve::getMaxAbsRowVal(int row) const {
  double maxVal = 0.0;

  for (const auto& nz : getRowVector(row))
    maxVal = std::max(std::abs(nz.value()), maxVal);

  return maxVal;
}

void HPresolve::updateRowDualImpliedBounds(int row, int col, double val) {
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
        impliedDualRowBounds.getResidualSumLower(col, row, val);
    if (residualMinAct != -HIGHS_CONST_INF) {
      double impliedBound =
          double((HighsCDouble(dualRowUpper) - residualMinAct) / val);

      if (val > 0) {
        // only tighten bound if it is tighter by a wide enough margin
        if (impliedBound <
            implRowDualUpper[row] - 1000 * options->dual_feasibility_tolerance)
          changeImplRowDualUpper(row, impliedBound, col);
      } else {
        if (impliedBound > implRowDualLower[row] +
                               1000 * options->primal_feasibility_tolerance)
          changeImplRowDualLower(row, impliedBound, col);
      }
    }
  }

  if (dualRowLower != -HIGHS_CONST_INF) {
    // get maximal value of other columns in the row
    double residualMaxAct =
        impliedDualRowBounds.getResidualSumUpper(col, row, val);
    if (residualMaxAct != HIGHS_CONST_INF) {
      double impliedBound =
          double((HighsCDouble(dualRowLower) - residualMaxAct) / val);

      if (val > 0) {
        // only tighten bound if it is tighter by a wide enough margin
        if (impliedBound > implRowDualLower[row] +
                               1000 * options->primal_feasibility_tolerance)
          changeImplRowDualLower(row, impliedBound, col);
      } else {
        if (impliedBound <
            implRowDualUpper[row] - 1000 * options->dual_feasibility_tolerance)
          changeImplRowDualUpper(row, impliedBound, col);
      }
    }
  }
}

void HPresolve::updateColImpliedBounds(int row, int col, double val) {
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
    double residualMinAct = impliedRowBounds.getResidualSumLower(row, col, val);
    if (residualMinAct != -HIGHS_CONST_INF) {
      double impliedBound =
          double((HighsCDouble(rowUpper) - residualMinAct) / val);

      if (val > 0) {
        // bound is an upper bound
        // check if we may round the bound due to integrality restrictions
        if (model->integrality_[col] != HighsVarType::CONTINUOUS) {
          impliedBound =
              std::floor(impliedBound + options->mip_feasibility_tolerance);

          if (model->integrality_[col] == HighsVarType::INTEGER &&
              impliedBound < model->colUpper_[col])
            changeColUpper(col, impliedBound);
        }

        // only tighten bound if it is tighter by a wide enough margin
        if (impliedBound <
            implColUpper[col] - 1000 * options->primal_feasibility_tolerance)
          changeImplColUpper(col, impliedBound, row);
      } else {
        // bound is a lower bound
        // check if we may round the bound due to integrality restrictions
        if (model->integrality_[col] != HighsVarType::CONTINUOUS) {
          impliedBound =
              std::ceil(impliedBound - options->mip_feasibility_tolerance);

          if (model->integrality_[col] == HighsVarType::INTEGER &&
              impliedBound > model->colLower_[col])
            changeColLower(col, impliedBound);
        }

        // only tighten bound if it is tighter by a wide enough margin
        if (impliedBound >
            implColLower[col] + 1000 * options->primal_feasibility_tolerance)
          changeImplColLower(col, impliedBound, row);
      }
    }
  }

  if (rowLower != -HIGHS_CONST_INF) {
    // get maximal value of other columns in the row
    double residualMaxAct = impliedRowBounds.getResidualSumUpper(row, col, val);
    if (residualMaxAct != HIGHS_CONST_INF) {
      double impliedBound =
          double((HighsCDouble(rowLower) - residualMaxAct) / val);
      if (val > 0) {
        // bound is a lower bound
        // check if we may round the bound due to integrality restrictions
        if (model->integrality_[col] != HighsVarType::CONTINUOUS) {
          impliedBound =
              std::ceil(impliedBound - options->mip_feasibility_tolerance);

          // change bounds of integers immediately
          if (model->integrality_[col] == HighsVarType::INTEGER &&
              impliedBound > model->colLower_[col])
            changeColLower(col, impliedBound);
        }

        // only tighten bound if it is tighter by a wide enough margin
        if (impliedBound >
            implColLower[col] + 1000 * options->primal_feasibility_tolerance)
          changeImplColLower(col, impliedBound, row);
      } else {
        // bound is an upper bound
        // check if we may round the bound due to integrality restrictions
        if (model->integrality_[col] != HighsVarType::CONTINUOUS) {
          impliedBound =
              std::floor(impliedBound + options->mip_feasibility_tolerance);

          // change bounds of integers immediately
          if (model->integrality_[col] == HighsVarType::INTEGER &&
              impliedBound < model->colUpper_[col])
            changeColUpper(col, impliedBound);
        }

        // only tighten bound if it is tighter by a wide enough margin
        if (impliedBound <
            implColUpper[col] - 1000 * options->primal_feasibility_tolerance)
          changeImplColUpper(col, impliedBound, row);
      }
    }
  }
}

int HPresolve::findNonzero(int row, int col) {
  if (rowroot[row] == -1) return -1;

  auto get_row_left = [&](int pos) -> int& { return ARleft[pos]; };
  auto get_row_right = [&](int pos) -> int& { return ARright[pos]; };
  auto get_row_key = [&](int pos) { return Acol[pos]; };
  rowroot[row] =
      highs_splay(col, rowroot[row], get_row_left, get_row_right, get_row_key);

  if (Acol[rowroot[row]] == col) return rowroot[row];

  return -1;
}

void HPresolve::addToMatrix(int row, int col, double val) {
  int pos = findNonzero(row, col);

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

HighsTripletListSlice HPresolve::getColumnVector(int col) const {
  return HighsTripletListSlice(Arow.data(), Avalue.data(), Anext.data(),
                               colhead[col]);
}

HighsTripletTreeSlicePreOrder HPresolve::getRowVector(int row) const {
  return HighsTripletTreeSlicePreOrder(
      Acol.data(), Avalue.data(), ARleft.data(), ARright.data(), rowroot[row]);
}

HighsTripletTreeSliceInOrder HPresolve::getSortedRowVector(int row) const {
  return HighsTripletTreeSliceInOrder(Acol.data(), Avalue.data(), ARleft.data(),
                                      ARright.data(), rowroot[row]);
}

void HPresolve::markRowDeleted(int row) {
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

void HPresolve::markColDeleted(int col) {
  assert(!colDeleted[col]);
  // prevents col from being added to change vector
  changedColFlag[col] = true;
  colDeleted[col] = true;
  ++numDeletedCols;
}

void HPresolve::changeColUpper(int col, double newUpper) {
  if (model->integrality_[col] != HighsVarType::CONTINUOUS) {
    newUpper = std::floor(newUpper + options->mip_feasibility_tolerance);
    if (newUpper == model->colUpper_[col]) return;
  }

  double oldUpper = model->colUpper_[col];
  model->colUpper_[col] = newUpper;

  // printf("tightening upper bound of column %d from %.15g to %.15g\n", col,
  //        oldUpper, newUpper);

  if (implColUpper[col] <= std::min(newUpper, oldUpper)) {
    // only source row of the implied bound needs updating
    // all other rows already use the tighter bound
    int nz = findNonzero(colUpperSource[col], col);
    impliedRowBounds.updatedVarUpper(colUpperSource[col], col, Avalue[nz],
                                     oldUpper);
    return;
  }

  for (const HighsSliceNonzero& nonzero : getColumnVector(col)) {
    impliedRowBounds.updatedVarUpper(nonzero.index(), col, nonzero.value(),
                                     oldUpper);
    markChangedRow(nonzero.index());
  }
}

void HPresolve::changeColLower(int col, double newLower) {
  if (model->integrality_[col] != HighsVarType::CONTINUOUS) {
    newLower = std::ceil(newLower - options->mip_feasibility_tolerance);
    if (newLower == model->colLower_[col]) return;
  }

  double oldLower = model->colLower_[col];
  model->colLower_[col] = newLower;
  // printf("tightening lower bound of column %d from %.15g to %.15g\n", col,
  //        oldLower, newLower);

  if (implColLower[col] >= std::max(newLower, oldLower)) {
    // only source row of the implied bound needs updating
    // all other rows already used the tighter implied bound
    int nz = findNonzero(colLowerSource[col], col);
    impliedRowBounds.updatedVarLower(colLowerSource[col], col, Avalue[nz],
                                     oldLower);
    return;
  }

  for (const HighsSliceNonzero& nonzero : getColumnVector(col)) {
    impliedRowBounds.updatedVarLower(nonzero.index(), col, nonzero.value(),
                                     oldLower);
    markChangedRow(nonzero.index());
  }
}

void HPresolve::changeImplColUpper(int col, double newUpper, int originRow) {
  double oldImplUpper = implColUpper[col];
  int oldUpperSource = colUpperSource[col];
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

void HPresolve::changeImplColLower(int col, double newLower, int originRow) {
  double oldImplLower = implColLower[col];
  int oldLowerSource = colLowerSource[col];
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

void HPresolve::changeImplRowDualUpper(int row, double newUpper,
                                       int originCol) {
  double oldImplUpper = implRowDualUpper[row];
  int oldUpperSource = rowDualUpperSource[row];

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

void HPresolve::changeImplRowDualLower(int row, double newLower,
                                       int originCol) {
  double oldImplLower = implRowDualLower[row];
  int oldLowerSource = rowDualLowerSource[row];

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

void HPresolve::storeRow(int row) {
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
                        const std::vector<int>& Aindex,
                        const std::vector<int>& Astart) {
  Avalue.clear();
  Acol.clear();
  Arow.clear();

  int ncol = Astart.size() - 1;
  assert(ncol == int(colhead.size()));
  int nnz = Aval.size();

  Avalue = Aval;
  Acol.reserve(nnz);
  Arow.reserve(nnz);

  for (int i = 0; i != ncol; ++i) {
    int collen = Astart[i + 1] - Astart[i];
    Acol.insert(Acol.end(), collen, i);
    Arow.insert(Arow.end(), Aindex.begin() + Astart[i],
                Aindex.begin() + Astart[i + 1]);
  }

  Anext.resize(nnz);
  Aprev.resize(nnz);
  ARleft.resize(nnz);
  ARright.resize(nnz);
  for (int pos = 0; pos != nnz; ++pos) link(pos);

  eqiters.assign(model->numRow_, equations.end());
  for (int i = 0; i != model->numRow_; ++i) {
    // register equation
    if (model->rowLower_[i] == model->rowUpper_[i])
      eqiters[i] = equations.emplace(rowsize[i], i).first;
  }
}

void HPresolve::fromCSR(const std::vector<double>& ARval,
                        const std::vector<int>& ARindex,
                        const std::vector<int>& ARstart) {
  Avalue.clear();
  Acol.clear();
  Arow.clear();

  int nrow = ARstart.size() - 1;
  assert(nrow == int(rowroot.size()));
  int nnz = ARval.size();

  Avalue = ARval;
  Acol.reserve(nnz);
  Arow.reserve(nnz);
  //  entries.reserve(nnz);

  for (int i = 0; i != nrow; ++i) {
    int rowlen = ARstart[i + 1] - ARstart[i];
    Arow.insert(Arow.end(), rowlen, i);
    Acol.insert(Acol.end(), ARindex.begin() + ARstart[i],
                ARindex.begin() + ARstart[i + 1]);
  }

  for (int pos = 0; pos != nnz; ++pos) link(pos);

  eqiters.assign(nrow, equations.end());
  for (int i = 0; i != nrow; ++i) {
    // register equation
    if (model->rowLower_[i] == model->rowUpper_[i])
      eqiters[i] = equations.emplace(rowsize[i], i).first;
  }
}

int HPresolve::countFillin(int row) {
  int fillin = 0;
  for (int rowiter : rowpositions) {
    if (findNonzero(row, Acol[rowiter]) == -1) fillin += 1;
  }

  return fillin;
}

bool HPresolve::checkFillin(HighsHashTable<int, int>& fillinCache, int row,
                            int col) {
  // check numerics against markowitz tolerance
  assert(int(rowpositions.size()) == rowsize[row]);

  // check fillin against max fillin
  int fillin = -(rowsize[row] + colsize[col] - 1);

#if 1
  // first use fillin for rows where it is already computed
  for (int coliter = colhead[col]; coliter != -1; coliter = Anext[coliter]) {
    if (Arow[coliter] == row) continue;

    auto cachedFillin = fillinCache.find(Arow[coliter]);
    if (cachedFillin == nullptr) continue;

    fillin += (*cachedFillin - 1);
    if (fillin > options->presolve_substitution_maxfillin) return false;
  }

  // iterate over rows of substituted column again to count the fillin for the
  // remaining rows
  for (int coliter = colhead[col]; coliter != -1; coliter = Anext[coliter]) {
    assert(Acol[coliter] == col);

    if (Arow[coliter] == row) continue;

    int& cachedFillin = fillinCache[Arow[coliter]];

    if (cachedFillin != 0) continue;

    int rowfillin = countFillin(Arow[coliter]);
    cachedFillin = rowfillin + 1;
    fillin += rowfillin;

    if (fillin > options->presolve_substitution_maxfillin) return false;
    // we count a fillin of 1 if the column is not present in the row and
    // a fillin of zero otherwise. the fillin for the substituted column
    // itself was already counted before the loop so we skip that entry.
  }
#else
  for (int rowiter : rowpositions) {
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

void HPresolve::substitute(int row, int col, double rhs) {
  assert(!rowDeleted[row]);
  assert(!colDeleted[col]);
  int pos = findNonzero(row, col);
  assert(pos != -1);

  assert(Arow[pos] == row);
  assert(Acol[pos] == col);
  double substrowscale = -1.0 / Avalue[pos];
  assert(isImpliedFree(col));

  markRowDeleted(row);
  markColDeleted(col);

  // substitute the column in each row where it occurs
  for (int coliter = colhead[col]; coliter != -1;) {
    int colrow = Arow[coliter];
    double colval = Avalue[coliter];

    // walk to the next position before doing any modifications, because
    // the current position will be deleted in the loop below
    assert(Acol[coliter] == col);
    int colpos = coliter;
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

    for (int rowiter : rowpositions) {
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
    double objscale = model->colCost_[col] * substrowscale;
    model->offset_ -= objscale * rhs;
    for (int rowiter : rowpositions) {
      model->colCost_[Acol[rowiter]] += objscale * Avalue[rowiter];
      if (std::abs(model->colCost_[Acol[rowiter]]) <=
          options->small_matrix_value)
        model->colCost_[Acol[rowiter]] = 0.0;
    }
    assert(model->colCost_[col] == 0);
    model->colCost_[col] = 0.0;
  }

  // finally remove the entries of the row that was used for substitution
  for (int rowiter : rowpositions) unlink(rowiter);
}

void HPresolve::toCSC(std::vector<double>& Aval, std::vector<int>& Aindex,
                      std::vector<int>& Astart) {
  // set up the column starts using the column size array
  int numcol = colsize.size();
  Astart.resize(numcol + 1);
  int nnz = 0;
  for (int i = 0; i != numcol; ++i) {
    Astart[i] = nnz;
    nnz += colsize[i];
  }
  Astart[numcol] = nnz;

  // now setup the entries of the CSC matrix
  // we reuse the colsize array to count down to zero
  // for determining the position of each nonzero
  Aval.resize(nnz);
  Aindex.resize(nnz);
  int numslots = Avalue.size();
  assert(numslots - int(freeslots.size()) == nnz);
  for (int i = 0; i != numslots; ++i) {
    if (Avalue[i] == 0.0) continue;
    int pos = Astart[Acol[i] + 1] - colsize[Acol[i]];
    --colsize[Acol[i]];
    assert(colsize[Acol[i]] >= 0);
    Aval[pos] = Avalue[i];
    Aindex[pos] = Arow[i];
  }
}

void HPresolve::toCSR(std::vector<double>& ARval, std::vector<int>& ARindex,
                      std::vector<int>& ARstart) {
  // set up the row starts using the row size array
  int numrow = rowsize.size();
  ARstart.resize(numrow + 1);
  int nnz = 0;
  for (int i = 0; i != numrow; ++i) {
    ARstart[i] = nnz;
    nnz += rowsize[i];
  }
  ARstart[numrow] = nnz;

  // now setup the entries of the CSC matrix
  // we reuse the colsize array to count down to zero
  // for determining the position of each nonzero
  ARval.resize(nnz);
  ARindex.resize(nnz);
  for (int i = 0; i != nnz; ++i) {
    if (Avalue[i] == 0.0) continue;
    int pos = ARstart[Arow[i] + 1] - rowsize[Arow[i]];
    --rowsize[Arow[i]];
    assert(rowsize[Arow[i]] >= 0);
    ARval[pos] = Avalue[i];
    ARindex[pos] = Acol[i];
  }
}

HPresolve::Result HPresolve::doubletonEq(HighsPostsolveStack& postSolveStack,
                                         int row) {
  assert(!rowDeleted[row]);
  assert(rowsize[row] == 2);
  assert(model->rowLower_[row] == model->rowUpper_[row]);
  // printf("doubleton equation: ");
  // debugPrintRow(row);
  int nzPos1 = rowroot[row];
  int nzPos2 = ARright[nzPos1] != -1 ? ARright[nzPos1] : ARleft[nzPos1];

  int substcol;
  int staycol;
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
      // both columns continuous, substitutethe one with fewer nonzeros
      if (colsize[Acol[nzPos1]] < colsize[Acol[nzPos2]]) {
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
  if (stayImplLower > oldStayLower + options->primal_feasibility_tolerance)
    changeColLower(staycol, stayImplLower);

  if (stayImplUpper < oldStayUpper - options->primal_feasibility_tolerance)
    changeColUpper(staycol, stayImplUpper);

  postSolveStack.doubletonEquation(
      row, substcol, staycol, substcoef, staycoef, rhs, substLower, substUpper,
      oldStayLower, oldStayUpper, model->colLower_[staycol],
      model->colUpper_[staycol], model->colCost_[substcol],
      getColumnVector(substcol));

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
                                          int row) {
  assert(!rowDeleted[row]);
  assert(rowsize[row] == 1);

  // the tree of nonzeros of this row should just contain the single nonzero
  int nzPos = rowroot[row];
  assert(nzPos != -1);
  // nonzero should have the row in the row array
  assert(Arow[nzPos] == row);
  // tree with one element should not have children
  assert(ARleft[nzPos] == -1);
  assert(ARright[nzPos] == -1);

  int col = Acol[nzPos];
  double val = Avalue[nzPos];

  // printf("singleton row\n");
  // debugPrintRow(row);
  // delete row singleton nonzero directly, we have all information that we need
  // in local variables
  markRowDeleted(row);
  unlink(nzPos);

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

  bool lowerTightened = newColLower > model->colLower_[col];
  bool upperTightened = newColUpper < model->colUpper_[col];
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
    if (ub < lb || (ub > lb && (ub - lb) * col_numerics_threshold[col] *
                                       options->presolve_pivot_threshold <=
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

  // update bounds, or remove as fixed column directly
  if (ub == lb) {
    // does not matter if at upper or at lower
    // try to use bound that is unchanged for fixing
    if (model->colUpper_[col] == ub)
      fixColToUpper(postSolveStack, col);
    else {
      if (model->colLower_[col] != lb) changeColLower(col, lb);

      fixColToLower(postSolveStack, col);
    }
  } else {
    // just update bounds (and row activities)
    if (lowerTightened) changeColLower(col, lb);

    if (upperTightened) changeColUpper(col, ub);
  }

  return checkLimits(postSolveStack);
}

HPresolve::Result HPresolve::singletonCol(HighsPostsolveStack& postSolveStack,
                                          int col) {
  assert(colsize[col] == 1);
  assert(!colDeleted[col]);
  int nzPos = colhead[col];
  int row = Arow[nzPos];
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
    else {
      // todo: forcing column, since this implies colDualLower >= 0 and we
      // already checked that colDualUpper <= 0
      //
    }
    return checkLimits(postSolveStack);
  }
  if (colDualLower >= -options->dual_feasibility_tolerance) {
    if (model->colLower_[col] != -HIGHS_CONST_INF)
      fixColToLower(postSolveStack, col);
    else {
      // todo: forcing column, since this implies colDualUpper <= 0 and we
      // already checked that colDualLower >= 0
      //
    }
    return checkLimits(postSolveStack);
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
    if (rowDualLower[row] >= 0) {
      rhs = model->rowUpper_[row];
      rowType = HighsPostsolveStack::RowType::Leq;
    } else if (rowDualUpper[row] <= 0) {
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

HPresolve::Result HPresolve::rowPresolve(HighsPostsolveStack& postsolveStack,
                                         int row) {
  assert(!rowDeleted[row]);

  // handle special cases directly via a call to the specialized procedure
  switch (rowsize[row]) {
    default:
      break;
    case 1:
      return singletonRow(postsolveStack, row);
    case 2:
      if (model->rowLower_[row] == model->rowUpper_[row])
        return doubletonEq(postsolveStack, row);
  }

  // printf("row presolve: ");
  // debugPrintRow(row);
  double impliedRowUpper = impliedRowBounds.getSumUpper(row);
  double impliedRowLower = impliedRowBounds.getSumLower(row);

#ifndef NDEBUG
  {
    double iRUpper = 0.0;
    double iRLower = 0.0;
    for (const HighsSliceNonzero& nonz : getRowVector(row)) {
      double lb = colLowerSource[nonz.index()] == row
                      ? model->colLower_[nonz.index()]
                      : std::max(implColLower[nonz.index()],
                                 model->colLower_[nonz.index()]);
      double ub = colUpperSource[nonz.index()] == row
                      ? model->colUpper_[nonz.index()]
                      : std::min(implColUpper[nonz.index()],
                                 model->colUpper_[nonz.index()]);
      if (nonz.value() > 0) {
        iRUpper += ub * nonz.value();
        iRLower += lb * nonz.value();
      } else {
        iRUpper += lb * nonz.value();
        iRLower += ub * nonz.value();
      }
    }

    assert(iRUpper == impliedRowUpper ||
           std::abs(iRUpper - impliedRowUpper) <=
               options->primal_feasibility_tolerance);
    assert(iRLower == impliedRowLower ||
           std::abs(iRLower - impliedRowLower) <=
               options->primal_feasibility_tolerance);

    // printf("implied bounds without tightenings: [%g,%g]\n", baseiRLower,
    //        baseiRUpper);
  }
#endif

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
    postsolveStack.redundantRow(row, getRowVector(row));
    removeRow(row);
    return checkLimits(postsolveStack);
  }
  // todo: do additional single row presolve for mip here. It may assume a
  // non-redundant and non-infeasible row when considering variable and implied
  // bounds

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
    ++numForcingRow;
    storeRow(row);
    auto rowVector = getStoredRow();

    postsolveStack.forcingRow(row, rowVector, model->rowLower_[row],
                              HighsPostsolveStack::RowType::Geq);

    // already mark the row as deleted, since otherwise it would be registered
    // as changed/singleton in the process of fixing and removing the contained
    // columns
    int nfixed = 0;
    for (const HighsSliceNonzero& nonzero : rowVector) {
      if (nonzero.value() > 0) {
        if (model->colUpper_[nonzero.index()] <=
                implColUpper[nonzero.index()] ||
            colUpperSource[nonzero.index()] == row) {
          // the upper bound of the column is as tight as the implied upper
          // bound or comes from this row, which means it is not used in the
          // rows implied bounds. Therefore we can fix the variable at its upper
          // bound.
          postsolveStack.fixedColAtUpper(nonzero.index(),
                                         model->colUpper_[nonzero.index()],
                                         model->colCost_[nonzero.index()],
                                         getColumnVector(nonzero.index()));
          if (model->colLower_[nonzero.index()] <
              model->colUpper_[nonzero.index()])
            changeColLower(nonzero.index(), model->colUpper_[nonzero.index()]);
          removeFixedCol(nonzero.index());
          ++nfixed;
        }
      } else {
        if (model->colLower_[nonzero.index()] >=
                implColLower[nonzero.index()] ||
            colLowerSource[nonzero.index()] == row) {
          // the upper bound of the column is as tight as the implied upper
          // bound or comes from this row, which means it is not used in the
          // rows implied bounds. Therefore we can fix the variable at its upper
          // bound.
          postsolveStack.fixedColAtLower(nonzero.index(),
                                         model->colUpper_[nonzero.index()],
                                         model->colCost_[nonzero.index()],
                                         getColumnVector(nonzero.index()));

          if (model->colUpper_[nonzero.index()] >
              model->colLower_[nonzero.index()])
            changeColUpper(nonzero.index(), model->colLower_[nonzero.index()]);
          removeFixedCol(nonzero.index());
          ++nfixed;
        }
      }
    }

    // now the row might be empty, but not necessarily because the implied
    // column bounds might be implied by other rows in which case we cannot fix
    // the column
    if (nfixed != 0) {
      if (rowsize[row] == 0) {
        markRowDeleted(row);
        postsolveStack.redundantRow(row, HighsEmptySlice());
      }
      HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolveStack));

      if (!rowDeleted[row]) return rowPresolve(postsolveStack, row);
    }
    // if there are any new row singletons, also remove them immediately
  } else if (impliedRowLower >=
             model->rowUpper_[row] - options->primal_feasibility_tolerance) {
    // forcing row in the other direction
    ++numForcingRow;
    storeRow(row);
    auto rowVector = getStoredRow();

    postsolveStack.forcingRow(row, rowVector, model->rowUpper_[row],
                              HighsPostsolveStack::RowType::Leq);
    int nfixed = 0;
    for (const HighsSliceNonzero& nonzero : rowVector) {
      if (nonzero.value() < 0) {
        if (model->colUpper_[nonzero.index()] <=
                implColUpper[nonzero.index()] ||
            colUpperSource[nonzero.index()] == row) {
          postsolveStack.fixedColAtUpper(nonzero.index(),
                                         model->colUpper_[nonzero.index()],
                                         model->colCost_[nonzero.index()],
                                         getColumnVector(nonzero.index()));
          if (model->colLower_[nonzero.index()] <
              model->colUpper_[nonzero.index()])
            changeColLower(nonzero.index(), model->colUpper_[nonzero.index()]);

          removeFixedCol(nonzero.index());
          ++nfixed;
        }
      } else {
        if (model->colLower_[nonzero.index()] >=
                implColLower[nonzero.index()] ||
            colLowerSource[nonzero.index()] == row) {
          postsolveStack.fixedColAtLower(nonzero.index(),
                                         model->colUpper_[nonzero.index()],
                                         model->colCost_[nonzero.index()],
                                         getColumnVector(nonzero.index()));
          if (model->colUpper_[nonzero.index()] >
              model->colLower_[nonzero.index()])
            changeColUpper(nonzero.index(), model->colLower_[nonzero.index()]);

          removeFixedCol(nonzero.index());
          ++nfixed;
        }
      }
    }
    if (nfixed != 0) {
      if (rowsize[row] == 0) {
        markRowDeleted(row);
        postsolveStack.redundantRow(row, HighsEmptySlice());
      }
      HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolveStack));

      if (!rowDeleted[row]) return rowPresolve(postsolveStack, row);
    }
  }

  bool hasRowUpper =
      model->rowUpper_[row] != HIGHS_CONST_INF ||
      implRowDualLower[row] > options->dual_feasibility_tolerance;
  bool hasRowLower =
      model->rowLower_[row] != HIGHS_CONST_INF ||
      implRowDualUpper[row] < -options->dual_feasibility_tolerance;

  if ((hasRowUpper && impliedRowBounds.getNumInfSumLower(row) <= 1) ||
      (hasRowLower && impliedRowBounds.getNumInfSumUpper(row) <= 1)) {
    for (const HighsSliceNonzero& nonzero : getRowVector(row))
      updateColImpliedBounds(row, nonzero.index(), nonzero.value());
  }

  return checkLimits(postsolveStack);
}

HPresolve::Result HPresolve::colPresolve(HighsPostsolveStack& postsolveStack,
                                         int col) {
  assert(!colDeleted[col]);

  if (colsize[col] == 1) return singletonCol(postsolveStack, col);

  double colDualUpper =
      impliedDualRowBounds.getSumUpper(col, model->colCost_[col]);
  double colDualLower =
      impliedDualRowBounds.getSumLower(col, model->colCost_[col]);

#if 0
  {
    double cDUpper = model->colCost_[col];
    double cDLower = model->colCost_[col];
    for (const HighsSliceNonzero& nonz : getColumnVector(col)) {
      if (nonz.value() > 0) {
        cDUpper += rowDualUpper[nonz.index()] * nonz.value();
        cDLower += rowDualLower[nonz.index()] * nonz.value();
      } else {
        cDUpper += rowDualLower[nonz.index()] * nonz.value();
        cDLower += rowDualUpper[nonz.index()] * nonz.value();
      }
    }

    assert(cDUpper == colDualUpper ||
           std::abs(cDUpper - colDualUpper) <=
               options->primal_feasibility_tolerance);
    assert(cDLower == colDualLower ||
           std::abs(cDLower - colDualLower) <=
               options->primal_feasibility_tolerance);
  }
#endif

  // check for dominated column
  if (colDualLower > options->dual_feasibility_tolerance) {
    if (model->colLower_[col] == -HIGHS_CONST_INF)
      return Result::DualInfeasible;
    else {
      fixColToLower(postsolveStack, col);
      HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolveStack));
    }
    return checkLimits(postsolveStack);
  }

  if (colDualUpper < -options->dual_feasibility_tolerance) {
    if (model->colUpper_[col] == HIGHS_CONST_INF)
      return Result::DualInfeasible;
    else {
      fixColToUpper(postsolveStack, col);
      HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolveStack));
    }
    return checkLimits(postsolveStack);
  }

  // check for weakly dominated column
  if (colDualUpper <= options->dual_feasibility_tolerance) {
    if (model->colUpper_[col] != HIGHS_CONST_INF) {
      fixColToUpper(postsolveStack, col);
      HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolveStack));
      return checkLimits(postsolveStack);
    } else {
      // todo: forcing column, since this implies colDualLower >= 0 and we
      // already checked that colDualUpper <= 0
      //
      printf("forcing column of size %d\n", colsize[col]);
    }
  } else if (colDualLower >= -options->dual_feasibility_tolerance) {
    // symmetric case for fixing to the lower bound
    if (model->colLower_[col] != HIGHS_CONST_INF) {
      fixColToLower(postsolveStack, col);
      HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolveStack));
      return checkLimits(postsolveStack);
    } else {
      // todo: forcing column, since this implies colDualUpper <= 0 and we
      // already checked that colDualLower >= 0
      printf("forcing column of size %d\n", colsize[col]);
    }
  }

  // column is not (weakly) dominated

  // integer columns cannot be used to tighten bounds on dual multipliers
  if (model->integrality_[col] == HighsVarType::INTEGER) return Result::Ok;

  // the associated dual constraint has an upper bound if there is an infinite
  // or strictly redundant column lower bound as then the reduced cost of the
  // column must not be positive i.e. <= 0
  bool dualConsHasUpper =
      model->colLower_[col] == -HIGHS_CONST_INF ||
      implColLower[col] >
          model->colLower_[col] + options->primal_feasibility_tolerance;

  bool dualConsHasLower =
      model->colUpper_[col] == HIGHS_CONST_INF ||
      implColUpper[col] <
          model->colUpper_[col] - options->primal_feasibility_tolerance;

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
  for (int row = 0; row != model->numRow_; ++row) {
    if (rowDeleted[row]) continue;
    HPRESOLVE_CHECKED_CALL(rowPresolve(postSolveStack, row));
    changedRowFlag[row] = false;
  }

  // same for the columns
  for (int col = 0; col != model->numCol_; ++col) {
    if (colDeleted[col]) continue;
    HPRESOLVE_CHECKED_CALL(colPresolve(postSolveStack, col));
    changedColFlag[col] = false;
  }

  return checkLimits(postSolveStack);
}

HPresolve::Result HPresolve::fastPresolveLoop(
    HighsPostsolveStack& postsolveStack) {
  do {
    storeCurrentProblemSize();

    HPRESOLVE_CHECKED_CALL(presolveChangedRows(postsolveStack));

    HPRESOLVE_CHECKED_CALL(removeDoubletonEquations(postsolveStack));

    HPRESOLVE_CHECKED_CALL(presolveColSingletons(postsolveStack));

    HPRESOLVE_CHECKED_CALL(presolveChangedCols(postsolveStack));

  } while (problemSizeReduction() > 0.05);

  return Result::Ok;
}

HPresolve::Result HPresolve::presolve(HighsPostsolveStack& postsolveStack) {
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

  numForcingRow = 0;

  HPRESOLVE_CHECKED_CALL(initialRowAndColPresolve(postsolveStack));
  printf(
      "after initial presolve scan: %d del rows %d del cols  %d forcing rows\n",
      numDeletedRows, numDeletedCols, numForcingRow);
  while (true) {
    HPRESOLVE_CHECKED_CALL(fastPresolveLoop(postsolveStack));
    printf(
        "after fast presolve loop: %d del rows %d del cols %d forcing rows\n",
        numDeletedRows, numDeletedCols, numForcingRow);

    storeCurrentProblemSize();

    HPRESOLVE_CHECKED_CALL(detectParallelRowsAndCols(postsolveStack));
    printf("after parallel rows/cols: %d del rows %d del cols\n",
           numDeletedRows, numDeletedCols);

    if (problemSizeReduction() > 0) {
      HPRESOLVE_CHECKED_CALL(fastPresolveLoop(postsolveStack));
      printf(
          "after fast presolve loop: %d del rows %d del cols %d forcing rows\n",
          numDeletedRows, numDeletedCols, numForcingRow);
    }

    storeCurrentProblemSize();

    HPRESOLVE_CHECKED_CALL(aggregator(postsolveStack));

    printf("after aggregator: %d del rows %d del cols\n", numDeletedRows,
           numDeletedCols);

    if (problemSizeReduction() > 0.05) continue;

    HPRESOLVE_CHECKED_CALL(fastPresolveLoop(postsolveStack));
    printf(
        "after fast presolve loop: %d del rows %d del cols %d forcing rows\n",
        numDeletedRows, numDeletedCols, numForcingRow);

    // todo sparsify

    break;
  }

  printf("num nonzeros: %lu\n", Avalue.size() - freeslots.size());

  return Result::Ok;
}

HPresolve::Result HPresolve::checkLimits(HighsPostsolveStack& postsolveStack) {
  // todo: check timelimit
#if 0
  for (int row = 0; row != model->numRow_; ++row) {
    if (rowDeleted[row]) continue;

    if (model->rowLower_[row] == model->rowUpper_[row]) {
      assert(eqiters[row] != equations.end());
      assert(eqiters[row]->first == rowsize[row]);
      assert(eqiters[row]->second == row);
    }

    double iRUpper = 0.0;
    double iRLower = 0.0;
    int rowlen = 0;
    for (const HighsSliceNonzero& nonz : getRowVector(row)) {
      ++rowlen;
      if (nonz.value() > 0) {
        iRUpper += model->colUpper_[nonz.index()] * nonz.value();
        iRLower += model->colLower_[nonz.index()] * nonz.value();
      } else {
        iRUpper += model->colLower_[nonz.index()] * nonz.value();
        iRLower += model->colUpper_[nonz.index()] * nonz.value();
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
  }
  for (int col = 0; col != model->numCol_; ++col) {
    if (colDeleted[col]) continue;
    double iDRUpper = 0.0;
    double iDRLower = 0.0;
    int collen = 0;
    for (const HighsSliceNonzero& nonz : getColumnVector(col)) {
      ++collen;
      if (nonz.value() > 0) {
        iDRUpper += rowDualUpper[nonz.index()] * nonz.value();
        iDRLower += rowDualLower[nonz.index()] * nonz.value();
      } else {
        iDRUpper += rowDualLower[nonz.index()] * nonz.value();
        iDRLower += rowDualUpper[nonz.index()] * nonz.value();
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

  return postsolveStack.numReductions() >= reductionLimit ? Result::Stopped
                                                          : Result::Ok;
}

void HPresolve::storeCurrentProblemSize() {
  oldNumCol = model->numCol_ - numDeletedCols;
  oldNumRow = model->numRow_ - numDeletedRows;
}

double HPresolve::problemSizeReduction() {
  double reduction =
      oldNumCol + oldNumRow -
      (model->numCol_ + model->numRow_ - numDeletedRows - numDeletedCols);

  return 100 * reduction / (oldNumRow + oldNumCol);
}

HighsModelStatus HPresolve::run(HighsPostsolveStack& postsolveStack) {
  switch (presolve(postsolveStack)) {
    case Result::Stopped:
    case Result::Ok:
      printf("presolve stopped with status Ok\n");
      break;
    case Result::PrimalInfeasible:
      printf("presolve detected primal infeasible problem\n");
      return HighsModelStatus::PRIMAL_INFEASIBLE;
    case Result::DualInfeasible:
      printf("presolve detected dual infeasible problem\n");
      return HighsModelStatus::DUAL_INFEASIBLE;
  }

  // todo: compress index space

  toCSC(model->Avalue_, model->Aindex_, model->Astart_);

  return HighsModelStatus::NOTSET;
}

HPresolve::Result HPresolve::aggregator(HighsPostsolveStack& postSolveStack) {
  int numcol = colsize.size();
  auto iter = equations.begin();
  int numsubst = 0;
  int numsubstint = 0;
  // substitutionOpportunities.clear();
  std::sort(
      substitutionOpportunities.begin(), substitutionOpportunities.end(),
      [&](const std::pair<int, int>& nz1, const std::pair<int, int>& nz2) {
        return std::min(rowsize[nz1.first], colsize[nz1.second]) <
               std::min(rowsize[nz2.first], colsize[nz2.second]);
      });
  int numSubstOpportunitites = substitutionOpportunities.size();
  for (int i = 0; i != numSubstOpportunitites; ++i) {
    int row = substitutionOpportunities[i].first;
    int col = substitutionOpportunities[i].second;
    if (rowDeleted[row] || colDeleted[col]) continue;

    if (!isImpliedFree(col) || !isDualImpliedFree(row)) continue;

    if (model->integrality_[col] == HighsVarType::INTEGER &&
        !isImpliedIntegral(col))
      continue;

    // in the case where the row has length two or the column has length two
    // we always do the substitution since the fillin can never be problematic
    if (rowsize[row] == 2 || colsize[col] == 2) {
      HighsPostsolveStack::RowType rowType =
          rowDualLower[row] >= 0 ? HighsPostsolveStack::RowType::Leq
                                 : HighsPostsolveStack::RowType::Geq;
      double rhs = rowType == HighsPostsolveStack::RowType::Geq
                       ? model->rowLower_[row]
                       : model->rowUpper_[row];
      ++numsubst;
      if (model->integrality_[col] == HighsVarType::INTEGER) ++numsubstint;
      storeRow(row);

      postSolveStack.freeColSubstitution(row, col, rhs, model->colCost_[col],
                                         rowType, getStoredRow(),
                                         getColumnVector(col));

      substitute(row, col, rhs);
      HPRESOLVE_CHECKED_CALL(removeRowSingletons(postSolveStack));
      HPRESOLVE_CHECKED_CALL(checkLimits(postSolveStack));
      continue;
    }

    int nzPos = findNonzero(row, col);
    if (rowsize[row] < colsize[col]) {
      double maxVal = getMaxAbsRowVal(row);
      if (std::abs(Avalue[nzPos]) <
          maxVal * options->presolve_pivot_threshold) {
        maxVal = getMaxAbsColVal(col);
        if (std::abs(Avalue[nzPos]) <
            maxVal * options->presolve_pivot_threshold)
          continue;
      }
    }

    storeRow(row);
    int fillin = -(rowsize[row] + colsize[col] - 1);
    for (const auto& nz : getColumnVector(col)) {
      if (nz.index() == row) continue;
      fillin += countFillin(nz.index());

      if (fillin > options->presolve_substitution_maxfillin) break;
    }

    if (fillin > options->presolve_substitution_maxfillin) continue;

    ++numsubst;
    if (model->integrality_[col] == HighsVarType::INTEGER) ++numsubstint;
    HighsPostsolveStack::RowType rowType =
        rowDualLower[row] >= 0 ? HighsPostsolveStack::RowType::Leq
                               : HighsPostsolveStack::RowType::Geq;
    double rhs = rowType == HighsPostsolveStack::RowType::Geq
                     ? model->rowLower_[row]
                     : model->rowUpper_[row];

    substitute(row, col, rhs);
    HPRESOLVE_CHECKED_CALL(removeRowSingletons(postSolveStack));
    HPRESOLVE_CHECKED_CALL(checkLimits(postSolveStack));
  }
  substitutionOpportunities.clear();
  return Result::Ok;
  if (numsubst != 0) {
    HPRESOLVE_CHECKED_CALL(fastPresolveLoop(postSolveStack));
  }

  HighsHashTable<int, int> fillinCache;
  std::vector<uint8_t> notimpliedfree(numcol);
  std::vector<std::pair<int, double>> aggr_cands;
  aggr_cands.reserve(colsize.size());

  int nonzerosCheckedSinceLastSubstitution = 0;

  std::pair<int, int> lastCheckedEq;

  while (iter != equations.end()) {
    // if (nonzerosCheckedSinceLastSubstitution > 10000) break;
    // extract sparsest equation
    lastCheckedEq = *iter;
    ++iter;
    int sparsesteq = lastCheckedEq.second;

    // extract aggregation candidates from equation. rule out integers if
    // integrality of coefficients does not work out, then rule out columns that
    // are not implied free
    double minintcoef = HIGHS_CONST_INF;
    int ncont = 0;

    storeRow(sparsesteq);

    nonzerosCheckedSinceLastSubstitution += rowsize[sparsesteq];

    aggr_cands.clear();
    double row_numerics_threshold = 0;
    for (int rowiter : rowpositions) {
      int col = Acol[rowiter];
      double absval = std::abs(Avalue[rowiter]);

      row_numerics_threshold = std::max(row_numerics_threshold, absval);

      if (model->integrality_[col] == HighsVarType::INTEGER) {
        // if there are non-integer variables in the row, no integer variable
        // can be used
        if (ncont != 0) continue;

        // if all variables in a row are integer variables, we still need to
        // check whether their coefficients are all integral
        minintcoef = std::min(absval, minintcoef);
        aggr_cands.emplace_back(col, absval);
      } else {
        // if this is the first continuous variable, we remove all integer
        // candidates that were stored before
        if (ncont == 0) aggr_cands.clear();

        aggr_cands.emplace_back(col, absval);
        ++ncont;
      }
    }

    row_numerics_threshold *= options->presolve_pivot_threshold;
    assert(ncont == 0 || ncont == int(aggr_cands.size()));

    // if all candidates are integral, we check that all coefficients are
    // integral when divided by the minimal absolute coefficient of the row. If
    // that is the case we keep all candidates with a value that is equal to the
    // minimal absolute coefficient value. Otherwise we skip this equation for
    // substitution.
    if (ncont == 0) {
      // all candidates are integer variables so we need to check if
      // all coefficients are integral when divided by the smallest absolute
      // coefficient value
      bool suitable = true;
      for (std::pair<int, double>& cand : aggr_cands) {
        double divval = cand.second / minintcoef;
        double intval = std::floor(divval + 0.5);
        if (std::abs(divval - intval) > options->mip_epsilon) {
          suitable = false;
          break;
        }
      }

      if (!suitable) continue;

      // candidates with the coefficient equal to the minimal absolute
      // coefficient value are suitable for substitution, other candidates are
      // now removed
      double maxintcoef = minintcoef + options->mip_epsilon;
      aggr_cands.erase(std::remove_if(aggr_cands.begin(), aggr_cands.end(),
                                      [&](const std::pair<int, double>& cand) {
                                        return cand.second > maxintcoef;
                                      }),
                       aggr_cands.end());
    }

    // remove candidates that have already been checked to be not implied free,
    // or that do not fulfill the numerics criteria to have their absolute
    // coefficient value in this row above the specified markowitz threshold
    // times the maximal absolute value in the candidates row or column. Note
    // that the "or"-nature of this numerics condition is not accidental.
    aggr_cands.erase(
        std::remove_if(aggr_cands.begin(), aggr_cands.end(),
                       [&](const std::pair<int, double>& cand) {
                         if (notimpliedfree[cand.first]) return true;

                         if (row_numerics_threshold > cand.second &&
                             col_numerics_threshold[cand.first] > cand.second)
                           return true;

                         return false;
                       }),
        aggr_cands.end());

    // check if any candidates are left
    if (aggr_cands.empty()) continue;

    // now sort the candidates to prioritize sparse columns, tiebreak by
    // preferring columns with a larger coefficient in this row which is better
    // for numerics
    std::sort(
        aggr_cands.begin(), aggr_cands.end(),
        [&](const std::pair<int, double>& cand1,
            const std::pair<int, double>& cand2) {
          return std::make_pair(colsize[cand1.first], -std::abs(cand1.second)) <
                 std::make_pair(colsize[cand2.first], -std::abs(cand2.second));
        });
    fillinCache.clear();
    int chosencand = -1;
    for (std::pair<int, double>& cand : aggr_cands) {
      bool isimpliedfree = isImpliedFree(cand.first);

      if (!isimpliedfree) {
        notimpliedfree[cand.first] = true;
        continue;
      }

      if (!checkFillin(fillinCache, sparsesteq, cand.first)) continue;

      // take the first suitable candidate
      chosencand = cand.first;
      break;
    }

    // if we have found no suitable candidate we continue with the next equation
    if (chosencand == -1) continue;

    // finally perform the substitution with the chosen candidate and update the
    // iterator to point to the next sparsest equation
    ++numsubst;
    if (model->integrality_[chosencand] == HighsVarType::INTEGER) ++numsubstint;

    // printf("substituting col %d with row %d\n", chosencand, sparsesteq);
    // debugPrintSubMatrix(sparsesteq, chosencand);

    postSolveStack.freeColSubstitution(
        sparsesteq, chosencand, model->rowUpper_[sparsesteq],
        model->colCost_[chosencand], HighsPostsolveStack::RowType::Eq,
        getStoredRow(), getColumnVector(chosencand));
    substitute(sparsesteq, chosencand, model->rowUpper_[sparsesteq]);

    HPRESOLVE_CHECKED_CALL(removeRowSingletons(postSolveStack));
    // HPRESOLVE_CHECKED_CALL(removeDoubletonEquations(postSolveStack));
    HPRESOLVE_CHECKED_CALL(checkLimits(postSolveStack));

    nonzerosCheckedSinceLastSubstitution = 0;

    iter = equations.upper_bound(lastCheckedEq);
  }

  printf("performed %d(%d int) substitutions\n", numsubst, numsubstint);
  return checkLimits(postSolveStack);
}

void HPresolve::substitute(int substcol, int staycol, double offset,
                           double scale) {
  // substitute the column in each row where it occurs
  for (int coliter = colhead[substcol]; coliter != -1;) {
    int colrow = Arow[coliter];
    double colval = Avalue[coliter];
    // walk to the next position before doing any modifications, because
    // the current position will be deleted in the loop below
    assert(Acol[coliter] == substcol);
    int colpos = coliter;
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

    model->colCost_[staycol] += scale * model->colCost_[substcol];

    if (std::abs(model->colCost_[staycol]) <= options->small_matrix_value)
      model->colCost_[staycol] = 0.0;
    model->colCost_[substcol] = 0.0;
  }
}

void HPresolve::fixColToLower(HighsPostsolveStack& postsolveStack, int col) {
  double fixval = model->colLower_[col];

  // printf("fixing column %d to %.15g\n", col, fixval);

  // mark the column as deleted first so that it is not registered as singleton
  // column upon removing its nonzeros
  postsolveStack.fixedColAtLower(col, fixval, model->colCost_[col],
                                 getColumnVector(col));
  markColDeleted(col);

  for (int coliter = colhead[col]; coliter != -1;) {
    int colrow = Arow[coliter];
    double colval = Avalue[coliter];
    assert(Acol[coliter] == col);

    int colpos = coliter;
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
  model->colCost_[col] = 0;
}

void HPresolve::fixColToUpper(HighsPostsolveStack& postsolveStack, int col) {
  double fixval = model->colUpper_[col];

  // printf("fixing column %d to %.15g\n", col, fixval);

  // mark the column as deleted first so that it is not registered as singleton
  // column upon removing its nonzeros
  postsolveStack.fixedColAtUpper(col, fixval, model->colCost_[col],
                                 getColumnVector(col));
  markColDeleted(col);

  for (int coliter = colhead[col]; coliter != -1;) {
    int colrow = Arow[coliter];
    double colval = Avalue[coliter];
    assert(Acol[coliter] == col);

    int colpos = coliter;
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
  model->colCost_[col] = 0;
}

void HPresolve::removeRow(int row) {
  assert(row < int(rowroot.size()));
  assert(row >= 0);
  // first mark the row as logically deleted, so that it is not register as
  // singleton row upon removing its nonzeros
  markRowDeleted(row);
  storeRow(row);
  for (int rowiter : rowpositions) {
    assert(Arow[rowiter] == row);
    unlink(rowiter);
  }
}

void HPresolve::removeFixedCol(int col) {
  double fixval = model->colLower_[col];

  markColDeleted(col);

  for (int coliter = colhead[col]; coliter != -1;) {
    int colrow = Arow[coliter];
    double colval = Avalue[coliter];
    assert(Acol[coliter] == col);

    int colpos = coliter;
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
  model->colCost_[col] = 0;
}

HPresolve::Result HPresolve::removeRowSingletons(
    HighsPostsolveStack& postSolveStack) {
  for (size_t i = 0; i != singletonRows.size(); ++i) {
    int row = singletonRows[i];
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
    int col = singletonColumns[i];
    if (colDeleted[col]) continue;
    HPRESOLVE_CHECKED_CALL(colPresolve(postSolveStack, col));
  }
  singletonColumns.erase(
      std::remove_if(
          singletonColumns.begin(), singletonColumns.end(),
          [&](int col) { return colDeleted[col] || colsize[col] > 1; }),
      singletonColumns.end());

  return Result::Ok;
}

HPresolve::Result HPresolve::presolveChangedRows(
    HighsPostsolveStack& postSolveStack) {
  std::vector<int> changedRows;
  changedRows.reserve(model->numRow_ - numDeletedRows);
  changedRows.swap(changedRowIndices);
  for (int row : changedRows) {
    if (rowDeleted[row]) continue;
    HPRESOLVE_CHECKED_CALL(rowPresolve(postSolveStack, row));
    changedRowFlag[row] = rowDeleted[row];
  }

  return Result::Ok;
}

HPresolve::Result HPresolve::presolveChangedCols(
    HighsPostsolveStack& postSolveStack) {
  std::vector<int> changedCols;
  changedCols.reserve(model->numCol_ - numDeletedCols);
  changedCols.swap(changedColIndices);
  for (int col : changedCols) {
    if (colDeleted[col]) continue;
    HPRESOLVE_CHECKED_CALL(colPresolve(postSolveStack, col));
    changedColFlag[col] = colDeleted[col];
  }

  return Result::Ok;
}

HPresolve::Result HPresolve::removeDoubletonEquations(
    HighsPostsolveStack& postSolveStack) {
  while (!equations.empty()) {
    auto eq = equations.begin();
    int eqrow = eq->second;
    assert(!rowDeleted[eqrow]);
    assert(eq->first == rowsize[eqrow]);
    assert(model->rowLower_[eqrow] == model->rowUpper_[eqrow]);

    switch (rowsize[eqrow]) {
      case 0:
        if (std::abs(model->rowUpper_[eqrow]) >
            options->primal_feasibility_tolerance)
          return Result::PrimalInfeasible;
        break;
      case 1:
        // handle empty and singleton rows
        HPRESOLVE_CHECKED_CALL(singletonRow(postSolveStack, eqrow));
        break;
      case 2:
        HPRESOLVE_CHECKED_CALL(doubletonEq(postSolveStack, eqrow));
        break;
      default:
        return Result::Ok;
    }
  }

  return Result::Ok;
}

int HPresolve::strengthenInequalities() {
  std::vector<int8_t> complementation;
  std::vector<double> reducedcost;
  std::vector<double> upper;
  std::vector<int> indices;
  std::vector<int> positions;
  std::vector<int> stack;
  std::vector<double> coefs;
  std::vector<int> cover;

  int numstrenghtened = 0;

  for (int row = 0; row != model->numRow_; ++row) {
    if (rowsize[row] <= 1) continue;
    if (model->rowLower_[row] != -HIGHS_CONST_INF &&
        model->rowUpper_[row] != HIGHS_CONST_INF)
      continue;

    // printf("strengthening knapsack of %d vars\n", rowsize[row]);

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
      int pos = stack.back();
      stack.pop_back();

      if (ARright[pos] != -1) stack.push_back(ARright[pos]);
      if (ARleft[pos] != -1) stack.push_back(ARleft[pos]);

      int8_t comp;
      double weight;
      double ub;
      weight = Avalue[pos] * scale;
      int col = Acol[pos];
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

    while (true) {
      if (maxviolation <=
              continuouscontribution + options->mip_feasibility_tolerance ||
          indices.empty())
        break;

      std::sort(indices.begin(), indices.end(), [&](int i1, int i2) {
        return reducedcost[i1] > reducedcost[i2];
      });

      HighsCDouble lambda = maxviolation - continuouscontribution;

      cover.clear();
      cover.reserve(indices.size());

      for (int i = indices.size() - 1; i >= 0; --i) {
        double delta = upper[indices[i]] * reducedcost[indices[i]];

        if (lambda <= delta + options->mip_feasibility_tolerance)
          cover.push_back(indices[i]);
        else
          lambda -= delta;
      }

      if (cover.empty()) break;

      int alpos = *std::min_element(
          cover.begin(), cover.end(),
          [&](int i1, int i2) { return reducedcost[i1] < reducedcost[i2]; });

      int coverend = cover.size();

      double al = reducedcost[alpos];
      coefs.resize(coverend);
      double coverrhs = std::max(
          std::ceil(double(lambda / al - options->mip_feasibility_tolerance)),
          1.0);
      HighsCDouble slackupper = -coverrhs;

      double step = HIGHS_CONST_INF;
      for (int i = 0; i != coverend; ++i) {
        coefs[i] =
            std::ceil(std::min(reducedcost[cover[i]], double(lambda)) / al -
                      options->small_matrix_value);
        slackupper += upper[cover[i]] * coefs[i];
        step = std::min(step, reducedcost[cover[i]] / coefs[i]);
      }
      step = std::min(step, double(maxviolation / coverrhs));
      maxviolation -= step * coverrhs;

      int slackind = reducedcost.size();
      reducedcost.push_back(step);
      upper.push_back(double(slackupper));

      for (int i = 0; i != coverend; ++i)
        reducedcost[cover[i]] -= step * coefs[i];

      indices.erase(std::remove_if(indices.begin(), indices.end(),
                                   [&](int i) {
                                     return reducedcost[i] <=
                                            options->mip_feasibility_tolerance;
                                   }),
                    indices.end());
      indices.push_back(slackind);
    }

    double threshold =
        double(maxviolation + options->mip_feasibility_tolerance);

    indices.erase(std::remove_if(indices.begin(), indices.end(),
                                 [&](int i) {
                                   return i >= (int)positions.size() ||
                                          std::abs(reducedcost[i]) <= threshold;
                                 }),
                  indices.end());
    if (indices.empty()) continue;

    if (scale == -1.0) {
      HighsCDouble lhs = model->rowLower_[row];
      for (int i : indices) {
        double coefdelta = double(reducedcost[i] - maxviolation);
        int pos = positions[i];

        if (complementation[i] == -1) {
          lhs -= coefdelta * model->colLower_[Acol[pos]];
          Avalue[pos] -= coefdelta;
        } else {
          lhs += coefdelta * model->colUpper_[Acol[pos]];
          Avalue[pos] += coefdelta;
        }

        if (std::abs(Avalue[pos]) <= options->small_matrix_value) unlink(pos);
      }

      model->rowLower_[row] = double(lhs);
    } else {
      HighsCDouble rhs = model->rowUpper_[row];
      for (int i : indices) {
        double coefdelta = double(reducedcost[i] - maxviolation);
        int pos = positions[i];

        if (complementation[i] == -1) {
          rhs += coefdelta * model->colLower_[Acol[pos]];
          Avalue[pos] += coefdelta;
        } else {
          rhs -= coefdelta * model->colUpper_[Acol[pos]];
          Avalue[pos] -= coefdelta;
        }

        if (std::abs(Avalue[pos]) <= options->small_matrix_value) unlink(pos);
      }

      model->rowUpper_[row] = double(rhs);
    }

    numstrenghtened += indices.size();
  }

  return numstrenghtened;
}

HPresolve::Result HPresolve::detectParallelRowsAndCols(
    HighsPostsolveStack& postsolveStack) {
  std::vector<std::uint64_t> rowHashes;
  std::vector<std::uint64_t> colHashes;
  std::vector<std::pair<double, int>> rowMax(rowsize.size());
  std::vector<std::pair<double, int>> colMax(colsize.size());

  HighsHashTable<int, int> numRowSingletons;

  int nnz = Avalue.size();
  rowHashes.assign(rowsize.begin(), rowsize.end());
  colHashes.assign(colsize.begin(), colsize.end());

  // Step 1: Determine scales for rows and columns and remove column singletons
  // from the intial row hashes which are initialized with the row sizes
  for (int i = 0; i != nnz; ++i) {
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
  for (int i = 0; i != nnz; ++i) {
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
  int numRowBuckets = 0;
  int numColBuckets = 0;

  std::unordered_multimap<std::uint64_t, int> buckets;

  for (int i = 0; i != model->numCol_; ++i) {
    if (colDeleted[i]) continue;
    if (colsize[i] == 0) {
      colPresolve(postsolveStack, i);
      continue;
    }
    auto it = buckets.find(colHashes[i]);
    decltype(it) last;

    int delCol = -1;
    int parallelColCandidate = -2;

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

      int col = -1;
      int duplicateCol = -1;
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
        return colScale > 0 ? isUpperImplied(col) : isLowerImplied(col);
      };

      auto colLowerInf = [&]() {
        if (!checkColImplBounds) return false;
        return colScale > 0 ? isLowerImplied(col) : isUpperImplied(col);
      };

      auto duplicateColUpperInf = [&]() {
        if (!checkDuplicateColImplBounds) return false;
        return isUpperImplied(duplicateCol);
      };

      auto duplicateColLowerInf = [&]() {
        if (!checkDuplicateColImplBounds) return false;
        return isLowerImplied(duplicateCol);
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
        assert(colScale >= 1.0);
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

      double objDiff =
          model->colCost_[col] * colScale - model->colCost_[duplicateCol];

      constexpr int kMergeParallelCols = 0;
      constexpr int kDominanceColToUpper = 1;
      constexpr int kDominanceColToLower = 2;
      constexpr int kDominanceDuplicateColToUpper = 3;
      constexpr int kDominanceDuplicateColToLower = 4;

      int reductionCase = kMergeParallelCols;
      // now do the case distinctions for dominated columns
      // the cases are a lot simpler due to the helper functions
      // for checking the infinite bounds which automatically
      // incorporate the check for the variable types that allow domination.
      if (objDiff < -options->dual_feasibility_tolerance) {
        // scaled col is better than duplicate col
        if (colUpperInf())
          reductionCase = kDominanceDuplicateColToLower;
        else if (duplicateColLowerInf())
          reductionCase =
              colScale > 0 ? kDominanceColToUpper : kDominanceColToLower;
        else
          continue;
      } else if (objDiff > options->dual_feasibility_tolerance) {
        // duplicate col is better than scaled col
        if (colLowerInf()) reductionCase = kDominanceDuplicateColToUpper;
        if (duplicateColUpperInf())
          reductionCase =
              colScale > 0 ? kDominanceColToLower : kDominanceColToUpper;
        else
          continue;
      } else {
        if (colUpperInf())
          reductionCase = kDominanceDuplicateColToLower;
        else if (colLowerInf())
          reductionCase = kDominanceDuplicateColToUpper;
        else if (duplicateColUpperInf())
          reductionCase =
              colScale > 0 ? kDominanceColToLower : kDominanceColToUpper;
        else if (duplicateColLowerInf())
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
              model->colLower_[col] + colScale * model->colLower_[duplicateCol];
          mergeUpper =
              model->colUpper_[col] + colScale * model->colUpper_[duplicateCol];
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

            int kMax = mergeUpper;
            bool representable = true;
            for (int k = mergeLower; k <= kMax; ++k) {
              // we loop over the domain of the merged variable to check whether
              // there exists a value for col and duplicateCol so that both are
              // within their bounds. since the merged column y is defined as y
              // = col + colScale * duplicateCol, we know that the value of col
              // can be computed as col = y - colScale * duplicateCol. Hence we
              // loop over the domain of col2 until we verify that a suitable
              // value of column 1 exists to yield the desired value for y.
              double mergeVal = mergeLower + k;
              int k2Max = model->colUpper_[duplicateCol];
              assert(k2Max == model->colUpper_[duplicateCol]);
              representable = false;
              for (int k2 = model->colLower_[duplicateCol]; k2 <= k2Max; ++k2) {
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
        int duplicateColRowPos = findNonzero(colNz.index(), duplicateCol);
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
            int row = Arow[colhead[duplicateCol]];
            numRowSingletons[row] -= 1;
          }
          fixColToLower(postsolveStack, duplicateCol);
          break;
        case kDominanceDuplicateColToUpper:
          delCol = duplicateCol;
          if (colsize[duplicateCol] == 1) {
            int row = Arow[colhead[duplicateCol]];
            numRowSingletons[row] -= 1;
          }
          fixColToUpper(postsolveStack, duplicateCol);
          break;
        case kDominanceColToLower:
          delCol = col;
          if (colsize[col] == 1) {
            int row = Arow[colhead[col]];
            numRowSingletons[row] -= 1;
          }
          fixColToLower(postsolveStack, col);
          break;
        case kDominanceColToUpper:
          delCol = col;
          if (colsize[col] == 1) {
            int row = Arow[colhead[col]];
            numRowSingletons[row] -= 1;
          }
          fixColToUpper(postsolveStack, col);
          break;
        case kMergeParallelCols:
          postsolveStack.duplicateColumn(
              colScale, model->colLower_[col], model->colUpper_[col],
              model->colLower_[duplicateCol], model->colUpper_[duplicateCol],
              col, duplicateCol,
              model->integrality_[col] == HighsVarType::INTEGER,
              model->integrality_[duplicateCol] == HighsVarType::INTEGER);

          if (colsize[duplicateCol] == 1) {
            int row = Arow[colhead[duplicateCol]];
            numRowSingletons[row] -= 1;
          }

          // mark duplicate column as deleted
          markColDeleted(duplicateCol);
          // remove all nonzeros of duplicateCol
          for (int coliter = colhead[duplicateCol]; coliter != -1;) {
            assert(Acol[coliter] == duplicateCol);

            int colpos = coliter;
            int colrow = Arow[coliter];
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

          // update bounds of col to bounds of merged column
          if (mergeLower != model->colLower_[col])
            changeColLower(col, mergeLower);

          if (mergeUpper != model->colUpper_[col])
            changeColUpper(col, mergeUpper);

          delCol = duplicateCol;
          break;
      }

      break;
    }

    if (delCol != -1) {
      if (delCol == i)
        buckets.emplace_hint(last, colHashes[i], i);
      else
        buckets.erase(last);

      HPRESOLVE_CHECKED_CALL(checkLimits(postsolveStack));
      // we could have new row singletons since a column was removed. Remove
      // those rows immediately
      HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolveStack));
    } else
      buckets.emplace_hint(last, colHashes[i], i);
  }

  buckets.clear();

  for (int i = 0; i != model->numRow_; ++i) {
    if (rowDeleted[i]) continue;
    if (rowsize[i] <= 1) {
      HPRESOLVE_CHECKED_CALL(rowPresolve(postsolveStack, i));
      ++numRowBuckets;
      continue;
    }
    auto it = buckets.find(rowHashes[i]);
    decltype(it) last;

    const int* numSingletonPtr = numRowSingletons.find(i);
    int numSingleton = numSingletonPtr ? *numSingletonPtr : 0;

    int delRow = -1;
    if (it == buckets.end())
      ++numRowBuckets;
    else
      storeRow(i);
    while (it != buckets.end() && it->first == rowHashes[i]) {
      int parallelRowCand = it->second;
      last = it++;

      numSingletonPtr = numRowSingletons.find(parallelRowCand);
      const int numSingletonCandidate = numSingletonPtr ? *numSingletonPtr : 0;
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
        int nzPos = findNonzero(parallelRowCand, rowNz.index());
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
            if (model->rowUpper_[parallelRowCand] == HIGHS_CONST_INF)
              changeImplRowDualUpper(parallelRowCand, HIGHS_CONST_INF, -1);
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
            rowLowerTightened = true;
            if (model->rowLower_[parallelRowCand] == -HIGHS_CONST_INF)
              changeImplRowDualLower(parallelRowCand, -HIGHS_CONST_INF, -1);
            model->rowLower_[parallelRowCand] = newLower;
          }
        }

        postsolveStack.duplicateRow(parallelRowCand, rowUpperTightened,
                                    rowLowerTightened, i, rowScale);
        delRow = i;
        markRowDeleted(i);
        for (int rowiter : rowpositions) unlink(rowiter);
        break;
      } else if (model->rowLower_[i] == model->rowUpper_[i]) {
        // row i is equation and parallel (except for singletons)
        // add to the row parallelRowCand
        // printf(
        //    "nearly parallel case with %d singletons in eq row and %d "
        //    "singletons in other row(eq=%d)\n",
        //    numSingleton, numSingletonCandidate,
        //    model->rowLower_[parallelRowCand] ==
        //        model->rowUpper_[parallelRowCand]);
        postsolveStack.equalityRowAddition(parallelRowCand, i, -rowScale);
        for (const HighsSliceNonzero& rowNz : getStoredRow()) {
          int pos = findNonzero(parallelRowCand, rowNz.index());
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
        HPRESOLVE_CHECKED_CALL(rowPresolve(postsolveStack, parallelRowCand));

      } else if (model->rowLower_[parallelRowCand] ==
                 model->rowUpper_[parallelRowCand]) {
        // printf(
        //    "nearly parallel case with %d singletons in eq row and %d "
        //    "singletons in other inequality row\n",
        //    numSingletonCandidate, numSingleton);
        // the row parallelRowCand is an equation; add it to the other row
        double scale = -rowMax[i].first / rowMax[parallelRowCand].first;
        postsolveStack.equalityRowAddition(i, parallelRowCand, scale);
        for (const HighsSliceNonzero& rowNz : getRowVector(parallelRowCand)) {
          int pos = findNonzero(i, rowNz.index());
          if (pos != -1)
            unlink(pos);  // all common nonzeros are cancelled, as the rows are
                          // parallel
          else            // might introduce a singleton
            addToMatrix(i, rowNz.index(), scale * rowNz.value());
        }

        if (model->rowUpper_[i] != HIGHS_CONST_INF)
          model->rowUpper_[i] = double(model->rowUpper_[i] -
                                       HighsCDouble(rowScale) *
                                           model->rowUpper_[parallelRowCand]);
        if (model->rowLower_[i] != -HIGHS_CONST_INF)
          model->rowLower_[i] = double(model->rowLower_[i] -
                                       HighsCDouble(rowScale) *
                                           model->rowUpper_[parallelRowCand]);

        HPRESOLVE_CHECKED_CALL(rowPresolve(postsolveStack, i));
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
      }
    }

    if (delRow != -1) {
      if (delRow == i)
        buckets.emplace_hint(last, rowHashes[i], i);
      else
        buckets.erase(last);

      HPRESOLVE_CHECKED_CALL(checkLimits(postsolveStack));
    } else
      buckets.emplace_hint(last, rowHashes[i], i);
  }

  return Result::Ok;
}

}  // namespace presolve