#include "presolve/HighsPostsolveStack.h"

#include <numeric>

#include "HighsCDouble.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsOptions.h"

namespace presolve {

void HighsPostsolveStack::initializeIndexMaps(int numRow, int numCol) {
  origNumRow = numRow;
  origNumCol = numCol;

  origRowIndex.resize(numRow);
  std::iota(origRowIndex.begin(), origRowIndex.end(), 0);

  origColIndex.resize(numCol);
  std::iota(origColIndex.begin(), origColIndex.end(), 0);
}

void HighsPostsolveStack::compressIndexMaps(
    const std::vector<int>& newRowIndex, const std::vector<int>& newColIndex) {
  // loop over rows, decrease row counter for deleted rows (marked with -1),
  // store original index at new index position otherwise
  int numRow = origRowIndex.size();
  for (size_t i = 0; i != newRowIndex.size(); ++i) {
    if (newRowIndex[i] == -1)
      --numRow;
    else
      origRowIndex[newRowIndex[i]] = origRowIndex[i];
  }
  // resize original index array to new size
  origRowIndex.resize(numRow);

  // now compress the column array
  int numCol = origColIndex.size();
  for (size_t i = 0; i != newColIndex.size(); ++i) {
    if (newColIndex[i] == -1)
      --numCol;
    else
      origColIndex[newColIndex[i]] = origColIndex[i];
  }
  origColIndex.resize(numCol);
}

void HighsPostsolveStack::FreeColSubstitution::undo(
    const HighsOptions& options,
    const std::vector<std::pair<int, double>>& rowValues,
    const std::vector<std::pair<int, double>>& colValues,
    HighsSolution& solution, HighsBasis& basis) {
  double colCoef = 0;
  // compute primal values
  HighsCDouble rowValue = 0;
  for (const auto& rowVal : rowValues) {
    if (rowVal.first == col)
      colCoef = rowVal.second;
    else
      rowValue += rowVal.second * solution.col_value[rowVal.first];
  }

  assert(colCoef != 0);
  solution.row_value[row] =
      double(rowValue + colCoef * solution.col_value[col]);
  solution.col_value[col] = double((rhs - rowValue) / colCoef);

  // if no dual values requested, return here
  if (solution.row_dual.empty()) return;

  // compute the row dual value such that reduced cost of basic column is 0
  solution.row_dual[row] = 0;
  HighsCDouble dualval = -colCost;
  for (const auto& colVal : colValues)
    dualval -= colVal.second * solution.row_dual[colVal.first];

  solution.col_dual[col] = 0;
  solution.row_dual[row] = double(dualval / colCoef);

  // set basis status
  basis.col_status[col] = HighsBasisStatus::BASIC;
  if (rowType == RowType::Eq)
    basis.row_status[row] = solution.row_dual[row] > 0
                                ? HighsBasisStatus::UPPER
                                : HighsBasisStatus::LOWER;
  else if (rowType == RowType::Geq)
    basis.row_status[row] = HighsBasisStatus::LOWER;
  else
    basis.row_status[row] = HighsBasisStatus::UPPER;
}

void HighsPostsolveStack::DoubletonEquation::undo(
    const HighsOptions& options,
    const std::vector<std::pair<int, double>>& colValues,
    HighsSolution& solution, HighsBasis& basis) {
  // retrieve the row and column index, the row side and the two
  // coefficients then compute the primal values
  solution.col_value[colSubst] =
      double((rhs - HighsCDouble(coef) * solution.col_value[col]) / coefSubst);

  // can only do primal postsolve, stop here
  if (solution.row_dual.empty()) return;

  if (solution.col_dual[col] > options.dual_feasibility_tolerance)
    basis.col_status[col] = HighsBasisStatus::LOWER;
  else if (solution.col_dual[col] < -options.dual_feasibility_tolerance)
    basis.col_status[col] = HighsBasisStatus::UPPER;

  // compute the current dual values of the row and the substituted column
  // before deciding on which column becomes basic
  // for each entry in a row i of the substituted column we added the doubleton
  // equation row with scale -a_i/substCoef. Therefore the dual multiplier of
  // this row i implicitly increases the dual multiplier of this doubleton
  // equation row with that scale.
  HighsCDouble rowDual = 0.0;
  solution.row_dual[row] = 0;
  for (const auto& colVal : colValues)
    rowDual -= colVal.second * solution.row_dual[colVal.first];
  rowDual /= coefSubst;
  solution.row_dual[row] = double(rowDual);

  // the equation was also added to the objective, so the current values need to
  // be changed
  solution.col_dual[colSubst] = substCost;
  solution.col_dual[col] += substCost * coef / coefSubst;

  if ((upperTightened && basis.col_status[col] == HighsBasisStatus::UPPER) ||
      (lowerTightened && basis.col_status[col] == HighsBasisStatus::LOWER)) {
    // column must get zero reduced cost as the current bound cannot be used
    // so alter the dual multiplier of the row to make the dual multiplier of
    // column zero
    double rowDualDelta = -solution.col_dual[col] / coef;
    solution.row_dual[row] = double(rowDual + rowDualDelta);
    solution.col_dual[col] = 0.0;
    solution.col_dual[colSubst] = double(
        HighsCDouble(solution.col_dual[colSubst]) + rowDualDelta * coefSubst);
    if ((std::signbit(coef) == std::signbit(coefSubst) &&
         basis.col_status[col] == HighsBasisStatus::UPPER) ||
        (std::signbit(coef) != std::signbit(coefSubst) &&
         basis.col_status[col] == HighsBasisStatus::LOWER))
      basis.col_status[colSubst] = HighsBasisStatus::LOWER;
    else
      basis.col_status[colSubst] = HighsBasisStatus::UPPER;
    basis.col_status[col] = HighsBasisStatus::BASIC;
  } else {
    // otherwise make the reduced cost of the subsituted column zero and make
    // that column basic
    double rowDualDelta = -solution.col_dual[colSubst] / coefSubst;
    solution.row_dual[row] = double(rowDual + rowDualDelta);
    solution.col_dual[colSubst] = 0.0;
    solution.col_dual[col] =
        double(HighsCDouble(solution.col_dual[col]) + rowDualDelta * coef);
    basis.col_status[colSubst] = HighsBasisStatus::BASIC;
  }

  if (solution.row_dual[row] < 0)
    basis.row_status[row] = HighsBasisStatus::LOWER;
  else
    basis.row_status[row] = HighsBasisStatus::UPPER;
}

void HighsPostsolveStack::EqualityRowAddition::undo(
    const HighsOptions& options,
    const std::vector<std::pair<int, double>>& eqRowValues,
    int& numMissingBasic, HighsSolution& solution, HighsBasis& basis) {
  // nothing more to do if the row is zero in the dual solution or there is
  // no dual solution
  if (solution.row_dual.empty() || solution.row_dual[row] == 0.0) return;

  // the dual multiplier of the row implicitly increases the dual multiplier
  // of the equation with the scale the equation was added with
  solution.row_dual[addedEqRow] =
      double(HighsCDouble(eqRowScale) * solution.row_dual[row] +
             solution.row_dual[addedEqRow]);

  if (basis.row_status[addedEqRow] == HighsBasisStatus::BASIC &&
      std::abs(solution.row_dual[addedEqRow]) >
          options.dual_feasibility_tolerance) {
    ++numMissingBasic;
    if (solution.row_dual[addedEqRow] > 0)
      basis.row_status[addedEqRow] = HighsBasisStatus::LOWER;
    else
      basis.row_status[addedEqRow] = HighsBasisStatus::UPPER;
#if 0
    // due to redundancy in the linear system it may happen that the equation
    // is basic in the solution. Now it got a nonzero dual multiplier so that
    // we need to make it non-basic. It must however be the case that we have
    // a different column in the equation that we can make basic.
    int bestBasicCol = -1;
    double bestBasicColViol = std::abs(solution.row_value[addedEqRow]);

    for (const auto& entry : eqRowValues) {
      if (basis.col_status[entry.first] != HighsBasisStatus::BASIC &&
          std::abs(solution.col_dual[entry.first]) < bestBasicColViol) {
        bestBasicCol = entry.first;
        bestBasicColViol = std::abs(solution.col_dual[entry.first]);
      }
    }

    if (bestBasicCol != -1) {
      // if we found a column that we can use as basic column instead of the row
      // such that the dual violation gets smaller, we use that column,
      // otherwise we leave the row basic to at least have a consistent basis
      basis.col_status[bestBasicCol] = HighsBasisStatus::BASIC;
      if (solution.row_dual[addedEqRow] > 0)
        basis.row_status[addedEqRow] = HighsBasisStatus::LOWER;
      else
        basis.row_status[addedEqRow] = HighsBasisStatus::UPPER;
    }
#endif
  }
}

void HighsPostsolveStack::EqualityRowAdditions::undo(
    const HighsOptions& options,
    const std::vector<std::pair<int, double>>& eqRowValues,
    const std::vector<std::pair<int, double>>& targetRows, int& numMissingBasic,
    HighsSolution& solution, HighsBasis& basis) {
  // nothing more to do if the row is zero in the dual solution or there is
  // no dual solution
  if (solution.row_dual.empty()) return;

  // the dual multiplier of the rows where the eq row was added implicitly
  // increases the dual multiplier of the equation with the scale that was used
  // for adding the equation
  HighsCDouble eqRowDual = solution.row_dual[addedEqRow];
  for (const auto& targetRow : targetRows)
    eqRowDual +=
        HighsCDouble(targetRow.second) * solution.row_dual[targetRow.first];

  solution.row_dual[addedEqRow] = double(eqRowDual);

  if (basis.row_status[addedEqRow] == HighsBasisStatus::BASIC &&
      std::abs(solution.row_dual[addedEqRow]) >
          options.dual_feasibility_tolerance) {
    ++numMissingBasic;
    if (solution.row_dual[addedEqRow] > 0)
      basis.row_status[addedEqRow] = HighsBasisStatus::LOWER;
    else
      basis.row_status[addedEqRow] = HighsBasisStatus::UPPER;
  }
}

void HighsPostsolveStack::ForcingColumn::undo(
    const HighsOptions& options,
    const std::vector<std::pair<int, double>>& colValues,
    HighsSolution& solution, HighsBasis& basis) {
  int nonbasicRow = -1;
  HighsBasisStatus nonbasicRowStatus = HighsBasisStatus::NONBASIC;
  double colValFromNonbasicRow;

  if (atInfiniteUpper) {
    colValFromNonbasicRow = -HIGHS_CONST_INF;
    // choose largest value as then all rows are feasible
    for (const auto& colVal : colValues) {
      double colValFromRow = solution.row_value[colVal.first] / colVal.second;
      if (colValFromRow > colValFromNonbasicRow) {
        nonbasicRow = colVal.first;
        colValFromNonbasicRow = colValFromRow;
        nonbasicRowStatus = colVal.second > 0 ? HighsBasisStatus::LOWER
                                              : HighsBasisStatus::UPPER;
      }
    }
  } else {
    colValFromNonbasicRow = HIGHS_CONST_INF;
    // choose smallest value, as then all rows are feasible
    for (const auto& colVal : colValues) {
      double colValFromRow = solution.row_value[colVal.first] / colVal.second;
      if (colValFromRow < colValFromNonbasicRow) {
        nonbasicRow = colVal.first;
        colValFromNonbasicRow = colValFromRow;
        nonbasicRowStatus = colVal.second > 0 ? HighsBasisStatus::UPPER
                                              : HighsBasisStatus::LOWER;
      }
    }
  }

  solution.col_value[col] = colValFromNonbasicRow;

  if (solution.col_dual.empty()) return;

  solution.col_dual[col] = 0.0;
  basis.col_status[col] = HighsBasisStatus::BASIC;
  basis.row_status[nonbasicRow] = nonbasicRowStatus;
}

void HighsPostsolveStack::ForcingColumnRemovedRow::undo(
    const HighsOptions& options,
    const std::vector<std::pair<int, double>>& rowValues,
    HighsSolution& solution, HighsBasis& basis) {
  // we use the row value as storage for the scaled value implied on the column
  // dual
  HighsCDouble val = rhs;
  for (const auto& rowVal : rowValues)
    val -= rowVal.second * solution.col_value[rowVal.first];

  solution.row_value[row] = double(val);

  if (!solution.row_dual.empty()) {
    solution.row_dual[row] = 0.0;
    basis.row_status[row] = HighsBasisStatus::BASIC;
  }
}

void HighsPostsolveStack::SingletonRow::undo(const HighsOptions& options,
                                             HighsSolution& solution,
                                             HighsBasis& basis) {
  // nothing to do if the rows dual value is zero in the dual solution or
  // there is no dual solution
  if (solution.row_dual.empty()) return;

  if ((!colLowerTightened ||
       basis.col_status[col] != HighsBasisStatus::LOWER) &&
      (!colUpperTightened ||
       basis.col_status[col] != HighsBasisStatus::UPPER)) {
    // the tightened bound is not used in the basic solution
    // hence we simply make the row basic and give it a dual multiplier of 0
    basis.row_status[row] = HighsBasisStatus::BASIC;
    solution.row_dual[row] = 0;
    return;
  }

  switch (basis.col_status[col]) {
    case HighsBasisStatus::LOWER:
      assert(colLowerTightened);
      if (coef > 0)
        // tightened lower bound comes from row lower bound
        basis.row_status[row] = HighsBasisStatus::LOWER;
      else
        // tightened lower bound comes from row upper bound
        basis.row_status[row] = HighsBasisStatus::UPPER;

      break;
    case HighsBasisStatus::UPPER:
      if (coef > 0)
        // tightened upper bound comes from row lower bound
        basis.row_status[row] = HighsBasisStatus::UPPER;
      else
        // tightened lower bound comes from row upper bound
        basis.row_status[row] = HighsBasisStatus::LOWER;
      break;
    default:
      assert(false);
  }

  // column becomes basic
  basis.col_status[col] = HighsBasisStatus::BASIC;
  // choose the row dual value such that the columns reduced cost becomes
  // zero
  solution.row_dual[row] = -solution.col_dual[col] / coef;
  solution.col_dual[col] = 0;
}

// column fixed to lower or upper bound
void HighsPostsolveStack::FixedCol::undo(
    const HighsOptions& options,
    const std::vector<std::pair<int, double>>& colValues,
    HighsSolution& solution, HighsBasis& basis) {
  // set solution value
  solution.col_value[col] = fixValue;

  if (solution.row_dual.empty()) return;

  // compute reduced cost
  HighsCDouble reducedCost = colCost;
  for (const auto& colVal : colValues)
    reducedCost += colVal.second * solution.row_dual[colVal.first];

  solution.col_dual[col] = double(reducedCost);

  // set basis status
  if (!basis.col_status.empty()) {
    basis.col_status[col] = fixType;
    if (basis.col_status[col] == HighsBasisStatus::NONBASIC)
      basis.col_status[col] = solution.col_dual[col] >= 0
                                  ? HighsBasisStatus::LOWER
                                  : HighsBasisStatus::UPPER;
  }
}

void HighsPostsolveStack::RedundantRow::undo(const HighsOptions& options,
                                             HighsSolution& solution,
                                             HighsBasis& basis) {
  // set row dual to zero if dual solution requested
  if (solution.row_dual.empty()) return;

  solution.row_dual[row] = 0.0;

  if (!basis.row_status.empty())
    basis.row_status[row] = HighsBasisStatus::BASIC;
}

void HighsPostsolveStack::ForcingRow::undo(
    const HighsOptions& options,
    const std::vector<std::pair<int, double>>& rowValues,
    HighsSolution& solution, HighsBasis& basis) {
  if (solution.row_dual.empty()) return;

  // compute the row dual multiplier and determine the new basic column
  int basicCol = -1;
  double dualDelta = 0;
  if (rowType == RowType::Leq) {
    for (const auto& rowVal : rowValues) {
      double colDual =
          solution.col_dual[rowVal.first] + rowVal.second * dualDelta;
      if (colDual * rowVal.second < 0) {
        // column is dual infeasible, increase the row dual such that its
        // reduced cost become zero and remember this column as the new basic
        // column for this row
        dualDelta = -solution.col_dual[rowVal.first] / rowVal.second;
        basicCol = rowVal.first;
      }
    }
  } else {
    for (const auto& rowVal : rowValues) {
      double colDual =
          solution.col_dual[rowVal.first] + rowVal.second * dualDelta;
      if (colDual * rowVal.second > 0) {
        // column is dual infeasible, decrease the row dual such that its
        // reduced cost become zero and remember this column as the new basic
        // column for this row
        dualDelta = -solution.col_dual[rowVal.first] / rowVal.second;
        basicCol = rowVal.first;
      }
    }
  }

  if (basicCol != -1) {
    solution.row_dual[row] = solution.row_dual[row] + dualDelta;
    for (const auto& rowVal : rowValues) {
      solution.col_dual[rowVal.first] =
          double(solution.col_dual[rowVal.first] +
                 HighsCDouble(dualDelta) * rowVal.second);
    }
    solution.col_dual[basicCol] = 0;
    basis.row_status[row] = (rowType == RowType::Geq ? HighsBasisStatus::LOWER
                                                     : HighsBasisStatus::UPPER);

    basis.col_status[basicCol] = HighsBasisStatus::BASIC;
  }
}

void HighsPostsolveStack::DuplicateRow::undo(const HighsOptions& options,
                                             HighsSolution& solution,
                                             HighsBasis& basis) {
  if (!rowUpperTightened && !rowLowerTightened) {
    // simple case of row2 being redundant, in which case it just gets a
    // dual multiplier of 0 and is made basic
    solution.row_dual[duplicateRow] = 0.0;
    basis.row_status[duplicateRow] = HighsBasisStatus::BASIC;
    return;
  }
  if (solution.row_dual[row] > options.dual_feasibility_tolerance)
    basis.row_status[row] = HighsBasisStatus::UPPER;
  else if (solution.row_dual[row] < -options.dual_feasibility_tolerance)
    basis.row_status[row] = HighsBasisStatus::LOWER;

  // at least one bound of the row was tightened by using the bound of the
  // scaled parallel row, hence we might need to make the parallel row
  // nonbasic and the row basic

  switch (basis.row_status[row]) {
    case HighsBasisStatus::BASIC:
      // if row is basic the parallel row is also basic
      solution.row_dual[duplicateRow] = 0.0;
      basis.row_status[duplicateRow] = HighsBasisStatus::BASIC;
      break;
    case HighsBasisStatus::UPPER:
      // if row sits on its upper bound, and the row upper bound was
      // tightened using the parallel row we make the row basic and
      // transfer its dual value to the parallel row with the proper scale
      if (rowUpperTightened) {
        solution.row_dual[duplicateRow] =
            solution.row_dual[row] / duplicateRowScale;
        solution.row_dual[row] = 0.0;
        basis.row_status[row] = HighsBasisStatus::BASIC;
        if (duplicateRowScale > 0)
          basis.row_status[duplicateRow] = HighsBasisStatus::UPPER;
        else
          basis.row_status[duplicateRow] = HighsBasisStatus::LOWER;
      } else {
        solution.row_dual[duplicateRow] = 0.0;
        basis.row_status[duplicateRow] = HighsBasisStatus::BASIC;
      }
      break;
    case HighsBasisStatus::LOWER:
      if (rowLowerTightened) {
        solution.row_dual[duplicateRow] =
            solution.row_dual[row] / duplicateRowScale;
        basis.row_status[row] = HighsBasisStatus::BASIC;
        solution.row_dual[row] = 0.0;
        if (duplicateRowScale > 0)
          basis.row_status[duplicateRow] = HighsBasisStatus::UPPER;
        else
          basis.row_status[duplicateRow] = HighsBasisStatus::LOWER;
      } else {
        solution.row_dual[duplicateRow] = 0.0;
        basis.row_status[duplicateRow] = HighsBasisStatus::BASIC;
      }
      break;
    default:
      assert(false);
  }
}

void HighsPostsolveStack::DuplicateColumn::undo(const HighsOptions& options,
                                                HighsSolution& solution,
                                                HighsBasis& basis) {
  // the column dual of the duplicate column is easily computed by scaling
  // since col * colScale yields the coefficient values and cost of the
  // duplicate column.
  if (!solution.col_dual.empty())
    solution.col_dual[duplicateCol] = solution.col_dual[col] * colScale;

  if (!basis.col_status.empty()) {
    // do postsolve using basis status if a basis is available:
    // if the merged column is nonbasic, we can just set both columns
    // to the corresponding basis status and value
    if (basis.col_status[col] == HighsBasisStatus::LOWER) {
      solution.col_value[col] = colLower;
      if (colScale > 0) {
        basis.col_status[duplicateCol] = HighsBasisStatus::LOWER;
        solution.col_value[duplicateCol] = duplicateColLower;
      } else {
        basis.col_status[duplicateCol] = HighsBasisStatus::UPPER;
        solution.col_value[duplicateCol] = duplicateColUpper;
      }
      // nothing else to do
      return;
    }

    if (basis.col_status[col] == HighsBasisStatus::UPPER) {
      solution.col_value[col] = colUpper;
      if (colScale > 0) {
        basis.col_status[duplicateCol] = HighsBasisStatus::UPPER;
        solution.col_value[duplicateCol] = duplicateColUpper;
      } else {
        basis.col_status[duplicateCol] = HighsBasisStatus::LOWER;
        solution.col_value[duplicateCol] = duplicateColLower;
      }
      // nothing else to do
      return;
    }

    assert(basis.col_status[col] == HighsBasisStatus::BASIC);
  }

  // either no basis for postsolve, or column status is basic. One of
  // the two columns must become nonbasic. In case of integrality it is
  // simpler to choose col, since it has a coefficient of +1 in the equation y
  // = col + colScale * duplicateCol where the merged column is y and is
  // currently using the index of col. The duplicateCol can have a positive or
  // negative coefficient. So for postsolve, we first start out with col
  // sitting at the lower bound and compute the corresponding value for the
  // duplicate column as (y - colLower)/colScale. Then the following things
  // might happen:
  // - case 1: the value computed for duplicateCol is within the bounds
  // - case 1.1: duplicateCol is continuous -> accept value, make col nonbasic
  // at lower and duplicateCol basic
  // - case 1.2: duplicateCol is integer -> accept value if integer feasible,
  // otherwise round down and compute value of col as
  // col = y - colScale * duplicateCol
  // - case 2: the value for duplicateCol violates the column bounds: make it
  // sit at the bound that is violated
  //           and compute the value of col as col = y - colScale *
  //           duplicateCol for basis postsolve col is basic and duplicateCol
  //           nonbasic at lower/upper depending on which bound is violated.

  double mergeVal = solution.col_value[col];
  solution.col_value[duplicateCol] =
      double((HighsCDouble(mergeVal) - colLower) / colScale);

  bool recomputeCol = false;

  if (solution.col_value[duplicateCol] > duplicateColUpper) {
    solution.col_value[duplicateCol] = duplicateColUpper;
    recomputeCol = true;
    if (!basis.col_status.empty())
      basis.col_status[duplicateCol] = HighsBasisStatus::UPPER;
  } else if (solution.col_value[duplicateCol] < duplicateColLower) {
    solution.col_value[duplicateCol] = duplicateColLower;
    recomputeCol = true;
    if (!basis.col_status.empty())
      basis.col_status[duplicateCol] = HighsBasisStatus::LOWER;
  } else if (duplicateColIntegral) {
    double roundVal = std::round(solution.col_value[duplicateCol]);
    if (std::abs(roundVal - solution.col_value[duplicateCol]) >
        options.mip_feasibility_tolerance) {
      solution.col_value[duplicateCol] =
          std::floor(solution.col_value[duplicateCol]);
      recomputeCol = true;
    }
  }

  if (recomputeCol) {
    solution.col_value[col] =
        mergeVal - colScale * solution.col_value[duplicateCol];
  } else if (!basis.col_status.empty()) {
    // setting col to its lower bound yielded a feasible value for
    // duplicateCol
    basis.col_status[duplicateCol] = basis.col_status[col];
    basis.col_status[col] = HighsBasisStatus::LOWER;
    solution.col_value[col] = colLower;
  }
}

}  // namespace presolve