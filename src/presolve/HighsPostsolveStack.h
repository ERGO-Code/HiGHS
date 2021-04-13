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
/**@file HighsPostsolveStack.h
 * @brief Class to hold all information for postsolve and can transform back
 * primal and dual solutions.
 */

#ifndef PRESOLVE_HIGHS_POSTSOLVE_STACK_H_
#define PRESOLVE_HIGHS_POSTSOLVE_STACK_H_

#include <cassert>
#include <cmath>
#include <numeric>
#include <tuple>
#include <vector>

#include "lp_data/HConst.h"
#include "lp_data/HStruct.h"
#include "lp_data/HighsOptions.h"
#include "util/HighsCDouble.h"
#include "util/HighsDataStack.h"
#include "util/HighsMatrixSlice.h"

// class HighsOptions;
namespace presolve {
class HighsPostsolveStack {
  // now a section of individual classes for each type of each transformation
  // step that requires postsolve starts each class gets as first argument the
  // current stack of ReductionValues and custom arguments that contain the
  // necessary information that is required to undo the transformation. The
  // constructor is responsible for storing all necessary information in class
  // members and the reduction value stack. The class members should be as slim
  // as possible and putting values on the reduction value stack should be
  // preferred, because the classes are stored in a discriminated union and the
  // largest size counts. The classes should implement an undo() function which
  // gets the ReductionValues as argument and can be expected to be called such
  // that the stack is in the state as after the constructor has been called.
  // The undo() call must pop all values from the stack that were added during
  // the constructor call, and should restore primal/dual solution values, as
  // well as the basis status as appropriate.
 public:
  enum class RowType {
    kGeq,
    kLeq,
    kEq,
  };
  struct Nonzero {
    HighsInt index;
    double value;

    Nonzero(HighsInt index, double value) : index(index), value(value) {}
    Nonzero() = default;
  };

 private:
  struct FreeColSubstitution {
    double rhs;
    double colCost;
    HighsInt row;
    HighsInt col;
    RowType rowType;

    void undo(const HighsOptions& options,
              const std::vector<Nonzero>& rowValues,
              const std::vector<Nonzero>& colValues, HighsSolution& solution,
              HighsBasis& basis);
  };

  struct DoubletonEquation {
    double coef;
    double coefSubst;
    double rhs;
    double substLower;
    double substUpper;
    double substCost;
    HighsInt row;
    HighsInt colSubst;
    HighsInt col;
    bool lowerTightened;
    bool upperTightened;

    void undo(const HighsOptions& options,
              const std::vector<Nonzero>& colValues, HighsSolution& solution,
              HighsBasis& basis);
  };

  struct EqualityRowAddition {
    HighsInt row;
    HighsInt addedEqRow;
    double eqRowScale;

    void undo(const HighsOptions& options,
              const std::vector<Nonzero>& eqRowValues, HighsSolution& solution,
              HighsBasis& basis);
  };

  struct EqualityRowAdditions {
    HighsInt addedEqRow;

    void undo(const HighsOptions& options,
              const std::vector<Nonzero>& eqRowValues,
              const std::vector<Nonzero>& targetRows, HighsSolution& solution,
              HighsBasis& basis);
  };
  struct SingletonRow {
    double coef;
    HighsInt row;
    HighsInt col;
    bool colLowerTightened;
    bool colUpperTightened;

    void undo(const HighsOptions& options, HighsSolution& solution,
              HighsBasis& basis);
  };

  // column fixed to lower or upper bound
  struct FixedCol {
    double fixValue;
    double colCost;
    HighsInt col;
    HighsBasisStatus fixType;

    void undo(const HighsOptions& options,
              const std::vector<Nonzero>& colValues, HighsSolution& solution,
              HighsBasis& basis);
  };

  struct RedundantRow {
    HighsInt row;

    void undo(const HighsOptions& options, HighsSolution& solution,
              HighsBasis& basis);
  };

  struct ForcingRow {
    double side;
    HighsInt row;
    RowType rowType;

    void undo(const HighsOptions& options,
              const std::vector<Nonzero>& rowValues, HighsSolution& solution,
              HighsBasis& basis);
  };

  struct ForcingColumn {
    double colCost;
    double colBound;
    HighsInt col;
    bool atInfiniteUpper;

    void undo(const HighsOptions& options,
              const std::vector<Nonzero>& colValues, HighsSolution& solution,
              HighsBasis& basis);
  };

  struct ForcingColumnRemovedRow {
    double rhs;
    HighsInt row;
    void undo(const HighsOptions& options,
              const std::vector<Nonzero>& rowValues, HighsSolution& solution,
              HighsBasis& basis);
  };

  struct DuplicateRow {
    double duplicateRowScale;
    HighsInt duplicateRow;
    HighsInt row;
    bool rowLowerTightened;
    bool rowUpperTightened;

    void undo(const HighsOptions& options, HighsSolution& solution,
              HighsBasis& basis);
  };

  struct DuplicateColumn {
    double colScale;
    double colLower;
    double colUpper;
    double duplicateColLower;
    double duplicateColUpper;
    HighsInt col;
    HighsInt duplicateCol;
    bool colIntegral;
    bool duplicateColIntegral;

    void undo(const HighsOptions& options, HighsSolution& solution,
              HighsBasis& basis);
  };

  /// tags for reduction
  enum class ReductionType : uint8_t {
    kFreeColSubstitution,
    kDoubletonEquation,
    kEqualityRowAddition,
    kEqualityRowAdditions,
    kSingletonRow,
    kFixedCol,
    kRedundantRow,
    kForcingRow,
    kForcingColumn,
    kForcingColumnRemovedRow,
    kDuplicateRow,
    kDuplicateColumn,
  };

  HighsDataStack reductionValues;
  std::vector<ReductionType> reductions;
  std::vector<HighsInt> origColIndex;
  std::vector<HighsInt> origRowIndex;

  std::vector<Nonzero> rowValues;
  std::vector<Nonzero> colValues;
  HighsInt origNumCol = -1;
  HighsInt origNumRow = -1;

 public:
  HighsInt getOrigRowIndex(HighsInt row) const {
    assert(row < (HighsInt)origRowIndex.size());
    return origRowIndex[row];
  }

  HighsInt getOrigColIndex(HighsInt col) const {
    assert(col < (HighsInt)origColIndex.size());
    return origColIndex[col];
  }

  void appendCutsToModel(HighsInt numCuts) {
    HighsInt currNumRow = origRowIndex.size();
    HighsInt newNumRow = currNumRow + numCuts;
    origRowIndex.resize(newNumRow);
    for (HighsInt i = currNumRow; i != newNumRow; ++i)
      origRowIndex[i] = origNumRow++;
  }

  void removeCutsFromModel(HighsInt numCuts) {
    origNumRow -= numCuts;

    HighsInt origRowIndexSize = origRowIndex.size();
    for (HighsInt i = origRowIndex.size() - 1; i >= 0; --i) {
      if (origRowIndex[i] < origNumRow) break;
      --origRowIndexSize;
    }

    origRowIndex.resize(origRowIndexSize);
  }

  HighsInt getOrigNumRow() const { return origNumRow; }

  HighsInt getOrigNumCol() const { return origNumCol; }

  void initializeIndexMaps(HighsInt numRow, HighsInt numCol);

  void compressIndexMaps(const std::vector<HighsInt>& newRowIndex,
                         const std::vector<HighsInt>& newColIndex);

  template <typename RowStorageFormat, typename ColStorageFormat>
  void freeColSubstitution(HighsInt row, HighsInt col, double rhs,
                           double colCost, RowType rowType,
                           const HighsMatrixSlice<RowStorageFormat>& rowVec,
                           const HighsMatrixSlice<ColStorageFormat>& colVec) {
    rowValues.clear();
    for (const HighsSliceNonzero& rowVal : rowVec)
      rowValues.emplace_back(origColIndex[rowVal.index()], rowVal.value());

    colValues.clear();
    for (const HighsSliceNonzero& colVal : colVec)
      colValues.emplace_back(origRowIndex[colVal.index()], colVal.value());

    reductionValues.push(FreeColSubstitution{rhs, colCost, origRowIndex[row],
                                             origColIndex[col], rowType});
    reductionValues.push(rowValues);
    reductionValues.push(colValues);
    reductions.push_back(ReductionType::kFreeColSubstitution);
  }

  template <typename ColStorageFormat>
  void doubletonEquation(HighsInt row, HighsInt colSubst, HighsInt col,
                         double coefSubst, double coef, double rhs,
                         double substLower, double substUpper, double substCost,
                         bool lowerTightened, bool upperTightened,
                         const HighsMatrixSlice<ColStorageFormat>& colVec) {
    colValues.clear();
    for (const HighsSliceNonzero& colVal : colVec)
      colValues.emplace_back(origRowIndex[colVal.index()], colVal.value());

    reductionValues.push(DoubletonEquation{
        coef, coefSubst, rhs, substLower, substUpper, substCost,
        row == -1 ? -1 : origRowIndex[row], origColIndex[colSubst],
        origColIndex[col], lowerTightened, upperTightened});
    reductionValues.push(colValues);
    reductions.push_back(ReductionType::kDoubletonEquation);
  }

  template <typename RowStorageFormat>
  void equalityRowAddition(HighsInt row, HighsInt addedEqRow, double eqRowScale,
                           const HighsMatrixSlice<RowStorageFormat>& eqRowVec) {
    rowValues.clear();
    for (const HighsSliceNonzero& rowVal : eqRowVec)
      rowValues.emplace_back(origColIndex[rowVal.index()], rowVal.value());

    reductionValues.push(EqualityRowAddition{
        origRowIndex[row], origRowIndex[addedEqRow], eqRowScale});
    reductionValues.push(rowValues);
    reductions.push_back(ReductionType::kEqualityRowAddition);
  }

  template <typename RowStorageFormat>
  void equalityRowAdditions(HighsInt addedEqRow,
                            const HighsMatrixSlice<RowStorageFormat>& eqRowVec,
                            const std::vector<Nonzero>& targetRows) {
    rowValues.clear();
    for (const HighsSliceNonzero& rowVal : eqRowVec)
      rowValues.emplace_back(origColIndex[rowVal.index()], rowVal.value());

    reductionValues.push(EqualityRowAdditions{origRowIndex[addedEqRow]});
    reductionValues.push(rowValues);
    reductionValues.push(targetRows);
    reductions.push_back(ReductionType::kEqualityRowAdditions);
  }

  void singletonRow(HighsInt row, HighsInt col, double coef,
                    bool tightenedColLower, bool tightenedColUpper) {
    reductionValues.push(SingletonRow{coef, origRowIndex[row],
                                      origColIndex[col], tightenedColLower,
                                      tightenedColUpper});
    reductions.push_back(ReductionType::kSingletonRow);
  }

  template <typename ColStorageFormat>
  void fixedColAtLower(HighsInt col, double fixValue, double colCost,
                       const HighsMatrixSlice<ColStorageFormat>& colVec) {
    assert(std::isfinite(fixValue));
    colValues.clear();
    for (const HighsSliceNonzero& colVal : colVec)
      colValues.emplace_back(origRowIndex[colVal.index()], colVal.value());

    reductionValues.push(FixedCol{fixValue, colCost, origColIndex[col],
                                  HighsBasisStatus::kLower});
    reductionValues.push(colValues);
    reductions.push_back(ReductionType::kFixedCol);
  }

  template <typename ColStorageFormat>
  void fixedColAtUpper(HighsInt col, double fixValue, double colCost,
                       const HighsMatrixSlice<ColStorageFormat>& colVec) {
    assert(std::isfinite(fixValue));
    colValues.clear();
    for (const HighsSliceNonzero& colVal : colVec)
      colValues.emplace_back(origRowIndex[colVal.index()], colVal.value());

    reductionValues.push(FixedCol{fixValue, colCost, origColIndex[col],
                                  HighsBasisStatus::kUpper});
    reductionValues.push(colValues);
    reductions.push_back(ReductionType::kFixedCol);
  }

  template <typename ColStorageFormat>
  void fixedColAtZero(HighsInt col, double colCost,
                      const HighsMatrixSlice<ColStorageFormat>& colVec) {
    colValues.clear();
    for (const HighsSliceNonzero& colVal : colVec)
      colValues.emplace_back(origRowIndex[colVal.index()], colVal.value());

    reductionValues.push(
        FixedCol{0.0, colCost, origColIndex[col], HighsBasisStatus::kZero});
    reductionValues.push(colValues);
    reductions.push_back(ReductionType::kFixedCol);
  }

  template <typename ColStorageFormat>
  void removedFixedCol(HighsInt col, double fixValue, double colCost,
                       const HighsMatrixSlice<ColStorageFormat>& colVec) {
    assert(std::isfinite(fixValue));
    colValues.clear();
    for (const HighsSliceNonzero& colVal : colVec)
      colValues.emplace_back(origRowIndex[colVal.index()], colVal.value());

    reductionValues.push(FixedCol{fixValue, colCost, origColIndex[col],
                                  HighsBasisStatus::kNonbasic});
    reductionValues.push(colValues);
    reductions.push_back(ReductionType::kFixedCol);
  }

  void redundantRow(HighsInt row) {
    reductionValues.push(RedundantRow{origRowIndex[row]});
    reductions.push_back(ReductionType::kRedundantRow);
  }

  template <typename RowStorageFormat>
  void forcingRow(HighsInt row,
                  const HighsMatrixSlice<RowStorageFormat>& rowVec, double side,
                  RowType rowType) {
    rowValues.clear();
    for (const HighsSliceNonzero& rowVal : rowVec)
      rowValues.emplace_back(origColIndex[rowVal.index()], rowVal.value());

    reductionValues.push(ForcingRow{side, origRowIndex[row], rowType});
    reductionValues.push(rowValues);
    reductions.push_back(ReductionType::kForcingRow);
  }

  template <typename ColStorageFormat>
  void forcingColumn(HighsInt col,
                     const HighsMatrixSlice<ColStorageFormat>& colVec,
                     double cost, double boundVal, bool atInfiniteUpper) {
    colValues.clear();
    for (const HighsSliceNonzero& colVal : colVec)
      colValues.emplace_back(origRowIndex[colVal.index()], colVal.value());

    reductionValues.push(
        ForcingColumn{cost, boundVal, origColIndex[col], atInfiniteUpper});
    reductionValues.push(colValues);
    reductions.push_back(ReductionType::kForcingColumn);
  }

  template <typename RowStorageFormat>
  void forcingColumnRemovedRow(
      HighsInt forcingCol, HighsInt row, double rhs,
      const HighsMatrixSlice<RowStorageFormat>& rowVec) {
    rowValues.clear();
    for (const HighsSliceNonzero& rowVal : rowVec)
      if (rowVal.index() != forcingCol)
        rowValues.emplace_back(origColIndex[rowVal.index()], rowVal.value());

    reductionValues.push(ForcingColumnRemovedRow{rhs, origRowIndex[row]});
    reductionValues.push(rowValues);
    reductions.push_back(ReductionType::kForcingColumnRemovedRow);
  }

  void duplicateRow(HighsInt row, bool rowUpperTightened,
                    bool rowLowerTightened, HighsInt duplicateRow,
                    double duplicateRowScale) {
    reductionValues.push(
        DuplicateRow{duplicateRowScale, origRowIndex[duplicateRow],
                     origRowIndex[row], rowLowerTightened, rowUpperTightened});
    reductions.push_back(ReductionType::kDuplicateRow);
  }

  void duplicateColumn(double colScale, double colLower, double colUpper,
                       double duplicateColLower, double duplicateColUpper,
                       HighsInt col, HighsInt duplicateCol, bool colIntegral,
                       bool duplicateColIntegral) {
    reductionValues.push(DuplicateColumn{
        colScale, colLower, colUpper, duplicateColLower, duplicateColUpper,
        origColIndex[col], origColIndex[duplicateCol], colIntegral,
        duplicateColIntegral});
    reductions.push_back(ReductionType::kDuplicateColumn);
  }

  void undo(const HighsOptions& options, HighsSolution& solution,
            HighsBasis& basis) {
    reductionValues.resetPosition();

    if (solution.col_value.size() != origColIndex.size()) return;
    if (solution.row_value.size() != origRowIndex.size()) return;
    bool dualPostSolve = solution.col_dual.size() == solution.col_value.size();

    // expand solution to original index space
    solution.col_value.resize(origNumCol);
    for (HighsInt i = origColIndex.size() - 1; i >= 0; --i) {
      assert(origColIndex[i] >= i);
      solution.col_value[origColIndex[i]] = solution.col_value[i];
    }

    solution.row_value.resize(origNumRow);
    for (HighsInt i = origRowIndex.size() - 1; i >= 0; --i) {
      assert(origRowIndex[i] >= i);
      solution.row_value[origRowIndex[i]] = solution.row_value[i];
    }

    if (dualPostSolve) {
      // if dual solution is given, expand dual solution and basis to original
      // index space
      solution.col_dual.resize(origNumCol);
      basis.col_status.resize(origNumCol);
      for (HighsInt i = origColIndex.size() - 1; i >= 0; --i) {
        basis.col_status[origColIndex[i]] = basis.col_status[i];
        solution.col_dual[origColIndex[i]] = solution.col_dual[i];
      }

      solution.row_dual.resize(origNumRow);
      basis.row_status.resize(origNumRow);
      for (HighsInt i = origRowIndex.size() - 1; i >= 0; --i) {
        basis.row_status[origRowIndex[i]] = basis.row_status[i];
        solution.row_dual[origRowIndex[i]] = solution.row_dual[i];
      }
    }

    // now undo the changes
    for (HighsInt i = reductions.size() - 1; i >= 0; --i) {
      switch (reductions[i]) {
        case ReductionType::kFreeColSubstitution: {
          FreeColSubstitution reduction;
          reductionValues.pop(colValues);
          reductionValues.pop(rowValues);
          reductionValues.pop(reduction);
          reduction.undo(options, rowValues, colValues, solution, basis);
          break;
        }
        case ReductionType::kDoubletonEquation: {
          DoubletonEquation reduction;
          reductionValues.pop(colValues);
          reductionValues.pop(reduction);
          reduction.undo(options, colValues, solution, basis);
          break;
        }
        case ReductionType::kEqualityRowAddition: {
          EqualityRowAddition reduction;
          reductionValues.pop(rowValues);
          reductionValues.pop(reduction);
          reduction.undo(options, rowValues, solution, basis);
          break;
        }
        case ReductionType::kEqualityRowAdditions: {
          EqualityRowAdditions reduction;
          reductionValues.pop(colValues);
          reductionValues.pop(rowValues);
          reductionValues.pop(reduction);
          reduction.undo(options, rowValues, colValues, solution, basis);
          break;
        }
        case ReductionType::kSingletonRow: {
          SingletonRow reduction;
          reductionValues.pop(reduction);
          reduction.undo(options, solution, basis);
          break;
        }
        case ReductionType::kFixedCol: {
          FixedCol reduction;
          reductionValues.pop(colValues);
          reductionValues.pop(reduction);
          reduction.undo(options, colValues, solution, basis);
          break;
        }
        case ReductionType::kRedundantRow: {
          RedundantRow reduction;
          reductionValues.pop(reduction);
          reduction.undo(options, solution, basis);
          break;
        }
        case ReductionType::kForcingRow: {
          ForcingRow reduction;
          reductionValues.pop(rowValues);
          reductionValues.pop(reduction);
          reduction.undo(options, rowValues, solution, basis);
          break;
        }
        case ReductionType::kForcingColumn: {
          ForcingColumn reduction;
          reductionValues.pop(colValues);
          reductionValues.pop(reduction);
          reduction.undo(options, colValues, solution, basis);
          break;
        }
        case ReductionType::kForcingColumnRemovedRow: {
          ForcingColumnRemovedRow reduction;
          reductionValues.pop(rowValues);
          reductionValues.pop(reduction);
          reduction.undo(options, rowValues, solution, basis);
          break;
        }
        case ReductionType::kDuplicateRow: {
          DuplicateRow reduction;
          reductionValues.pop(reduction);
          reduction.undo(options, solution, basis);
          break;
        }
        case ReductionType::kDuplicateColumn: {
          DuplicateColumn reduction;
          reductionValues.pop(reduction);
          reduction.undo(options, solution, basis);
        }
      }
    }
  }

  void undoPrimal(const HighsOptions& options, HighsSolution& solution) {
    reductionValues.resetPosition();

    if (solution.col_value.size() != origColIndex.size()) return;
    if (solution.row_value.size() != origRowIndex.size()) return;

    // expand solution to original index space
    solution.col_value.resize(origNumCol);
    for (HighsInt i = origColIndex.size() - 1; i >= 0; --i) {
      assert(origColIndex[i] >= i);
      solution.col_value[origColIndex[i]] = solution.col_value[i];
    }

    solution.row_value.resize(origNumRow);
    for (HighsInt i = origRowIndex.size() - 1; i >= 0; --i) {
      assert(origRowIndex[i] >= i);
      solution.row_value[origRowIndex[i]] = solution.row_value[i];
    }

    solution.row_dual.clear();
    solution.col_dual.clear();

    HighsBasis basis;
    // now undo the changes
    for (HighsInt i = reductions.size() - 1; i >= 0; --i) {
      switch (reductions[i]) {
        case ReductionType::kFreeColSubstitution: {
          FreeColSubstitution reduction;
          reductionValues.pop(colValues);
          reductionValues.pop(rowValues);
          reductionValues.pop(reduction);
          reduction.undo(options, rowValues, colValues, solution, basis);
          break;
        }
        case ReductionType::kDoubletonEquation: {
          DoubletonEquation reduction;
          reductionValues.pop(colValues);
          reductionValues.pop(reduction);
          reduction.undo(options, colValues, solution, basis);
          break;
        }
        case ReductionType::kEqualityRowAddition: {
          EqualityRowAddition reduction;
          reductionValues.pop(rowValues);
          reductionValues.pop(reduction);
          reduction.undo(options, rowValues, solution, basis);
          break;
        }
        case ReductionType::kEqualityRowAdditions: {
          EqualityRowAdditions reduction;
          reductionValues.pop(colValues);
          reductionValues.pop(rowValues);
          reductionValues.pop(reduction);
          reduction.undo(options, rowValues, colValues, solution, basis);
          break;
        }
        case ReductionType::kSingletonRow: {
          SingletonRow reduction;
          reductionValues.pop(reduction);
          reduction.undo(options, solution, basis);
          break;
        }
        case ReductionType::kFixedCol: {
          FixedCol reduction;
          reductionValues.pop(colValues);
          reductionValues.pop(reduction);
          reduction.undo(options, colValues, solution, basis);
          break;
        }
        case ReductionType::kRedundantRow: {
          RedundantRow reduction;
          reductionValues.pop(reduction);
          reduction.undo(options, solution, basis);
          break;
        }
        case ReductionType::kForcingRow: {
          ForcingRow reduction;
          reductionValues.pop(rowValues);
          reductionValues.pop(reduction);
          reduction.undo(options, rowValues, solution, basis);
          break;
        }
        case ReductionType::kForcingColumn: {
          ForcingColumn reduction;
          reductionValues.pop(colValues);
          reductionValues.pop(reduction);
          reduction.undo(options, colValues, solution, basis);
          break;
        }
        case ReductionType::kForcingColumnRemovedRow: {
          ForcingColumnRemovedRow reduction;
          reductionValues.pop(rowValues);
          reductionValues.pop(reduction);
          reduction.undo(options, rowValues, solution, basis);
          break;
        }
        case ReductionType::kDuplicateRow: {
          DuplicateRow reduction;
          reductionValues.pop(reduction);
          reduction.undo(options, solution, basis);
          break;
        }
        case ReductionType::kDuplicateColumn: {
          DuplicateColumn reduction;
          reductionValues.pop(reduction);
          reduction.undo(options, solution, basis);
        }
      }
    }
  }

  void undoUntil(const HighsOptions& options,
                 const std::vector<HighsInt>& flagRow,
                 const std::vector<HighsInt>& flagCol, HighsSolution& solution,
                 HighsBasis& basis, HighsInt numReductions) {
    reductionValues.resetPosition();

    if (solution.col_value.size() != origColIndex.size()) return;
    if (solution.row_value.size() != origRowIndex.size()) return;

    bool dualPostSolve = solution.col_dual.size() == solution.col_value.size();

    // expand solution to original index space
    solution.col_value.resize(origNumCol);
    for (HighsInt i = origColIndex.size() - 1; i >= 0; --i) {
      assert(origColIndex[i] >= i);
      solution.col_value[origColIndex[i]] = solution.col_value[i];
    }

    solution.row_value.resize(origNumRow);
    for (HighsInt i = origRowIndex.size() - 1; i >= 0; --i) {
      assert(origRowIndex[i] >= i);
      solution.row_value[origRowIndex[i]] = solution.row_value[i];
    }

    if (dualPostSolve) {
      // if dual solution is given, expand dual solution and basis to original
      // index space
      solution.col_dual.resize(origNumCol);
      basis.col_status.resize(origNumCol);
      for (HighsInt i = origColIndex.size() - 1; i >= 0; --i) {
        basis.col_status[origColIndex[i]] = basis.col_status[i];
        solution.col_dual[origColIndex[i]] = solution.col_dual[i];
      }

      solution.row_dual.resize(origNumRow);
      basis.row_status.resize(origNumRow);
      for (HighsInt i = origRowIndex.size() - 1; i >= 0; --i) {
        basis.row_status[origRowIndex[i]] = basis.row_status[i];
        solution.row_dual[origRowIndex[i]] = solution.row_dual[i];
      }
    }

    // now undo the changes
    for (HighsInt i = reductions.size() - 1; i >= numReductions; --i) {
      switch (reductions[i]) {
        case ReductionType::kFreeColSubstitution: {
          FreeColSubstitution reduction;
          reductionValues.pop(colValues);
          reductionValues.pop(rowValues);
          reductionValues.pop(reduction);
          reduction.undo(options, rowValues, colValues, solution, basis);
          break;
        }
        case ReductionType::kDoubletonEquation: {
          DoubletonEquation reduction;
          reductionValues.pop(colValues);
          reductionValues.pop(reduction);
          reduction.undo(options, colValues, solution, basis);
          break;
        }
        case ReductionType::kEqualityRowAddition: {
          EqualityRowAddition reduction;
          reductionValues.pop(rowValues);
          reductionValues.pop(reduction);
          reduction.undo(options, rowValues, solution, basis);
          break;
        }
        case ReductionType::kEqualityRowAdditions: {
          EqualityRowAdditions reduction;
          reductionValues.pop(colValues);
          reductionValues.pop(rowValues);
          reductionValues.pop(reduction);
          reduction.undo(options, rowValues, colValues, solution, basis);
          break;
        }
        case ReductionType::kSingletonRow: {
          SingletonRow reduction;
          reductionValues.pop(reduction);
          reduction.undo(options, solution, basis);
          break;
        }
        case ReductionType::kFixedCol: {
          FixedCol reduction;
          reductionValues.pop(colValues);
          reductionValues.pop(reduction);
          reduction.undo(options, colValues, solution, basis);
          break;
        }
        case ReductionType::kRedundantRow: {
          RedundantRow reduction;
          reductionValues.pop(reduction);
          reduction.undo(options, solution, basis);
          break;
        }
        case ReductionType::kForcingRow: {
          ForcingRow reduction;
          reductionValues.pop(rowValues);
          reductionValues.pop(reduction);
          reduction.undo(options, rowValues, solution, basis);
          break;
        }
        case ReductionType::kForcingColumn: {
          ForcingColumn reduction;
          reductionValues.pop(colValues);
          reductionValues.pop(reduction);
          reduction.undo(options, colValues, solution, basis);
          break;
        }
        case ReductionType::kForcingColumnRemovedRow: {
          ForcingColumnRemovedRow reduction;
          reductionValues.pop(rowValues);
          reductionValues.pop(reduction);
          reduction.undo(options, rowValues, solution, basis);
          break;
        }
        case ReductionType::kDuplicateRow: {
          DuplicateRow reduction;
          reductionValues.pop(reduction);
          reduction.undo(options, solution, basis);
          break;
        }
        case ReductionType::kDuplicateColumn: {
          DuplicateColumn reduction;
          reductionValues.pop(reduction);
          reduction.undo(options, solution, basis);
        }
      }
    }
  }

  size_t numReductions() const { return reductions.size(); }
};

}  // namespace presolve

#endif
