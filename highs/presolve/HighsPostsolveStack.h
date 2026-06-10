/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file HighsPostsolveStack.h
 * @brief Class to hold all information for postsolve and can transform back
 * primal and dual solutions.
 */

#ifndef PRESOLVE_HIGHS_POSTSOLVE_STACK_H_
#define PRESOLVE_HIGHS_POSTSOLVE_STACK_H_

#include <algorithm>
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
  enum class OrigRowType : uint8_t { kOriginal, kCut, kAppended };

  enum class RowType {
    kGeq,
    kLeq,
    kEq,
  };
  struct Nonzero {
    HighsInt index;
    double value;

    Nonzero(HighsInt index_, double value_) : index(index_), value(value_) {}
    Nonzero() = default;
  };

  template <typename RowStorageFormat>
  struct FmeRowData {
    HighsInt row;
    double rowLower;
    double rowUpper;
    HighsMatrixSlice<RowStorageFormat> rowVec;
  };

  struct FmeRowHeader {
    HighsInt row;
    double rowLower;
    double rowUpper;
  };

  struct FmeStepHeader {
    double colLower;
    double colUpper;
    double colCost;
    HighsInt col;
    HighsInt numPlus;
    HighsInt numMinus;
    HighsInt numNewRows;
  };

  struct FmeDescendant {
    HighsInt row;
    double scaleFactor;
  };

  struct FmeNewRow {
    HighsInt row;
    HighsInt plusParentIdx;   // index into plus parents (-1 if bound row)
    HighsInt minusParentIdx;  // index into minus parents (-1 if bound row)
  };

  size_t debug_prev_numreductions = 0;
  double debug_prev_col_lower = 0;
  double debug_prev_col_upper = 0;
  double debug_prev_row_lower = 0;
  double debug_prev_row_upper = 0;

 private:
  /// transform a column x by a linear mapping with a new column x'.
  /// I.e. substitute x = a * x' + b
  struct LinearTransform {
    double scale;
    double constant;
    HighsInt col;

    void undo(const HighsOptions& options, HighsSolution& solution) const;

    void transformToPresolvedSpace(std::vector<double>& primalSol) const;
  };

  struct FourierMotzkinObjCol {
    double offset;
    HighsInt col;

    void transformToPresolvedSpace(const std::vector<Nonzero>& costEntries,
                                   std::vector<double>& primalSol) const;

    void undo(const std::vector<Nonzero>& costEntries,
              HighsSolution& solution) const;
  };

  struct FreeColSubstitution {
    double rhs;
    double colCost;
    HighsInt row;
    HighsInt col;
    RowType rowType;

    void undo(const HighsPostsolveStack& postsolveStack,
              const HighsOptions& options,
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
    RowType rowType;

    void undo(const HighsPostsolveStack& postsolveStack,
              const HighsOptions& options,
              const std::vector<Nonzero>& colValues, HighsSolution& solution,
              HighsBasis& basis) const;
  };

  struct EqualityRowAddition {
    HighsInt row;
    HighsInt addedEqRow;
    double eqRowScale;

    void undo(const HighsPostsolveStack& postsolveStack,
              const HighsOptions& options,
              const std::vector<Nonzero>& eqRowValues, HighsSolution& solution,
              HighsBasis& basis) const;
  };

  struct EqualityRowAdditions {
    HighsInt addedEqRow;

    void undo(const HighsPostsolveStack& postsolveStack,
              const HighsOptions& options,
              const std::vector<Nonzero>& eqRowValues,
              const std::vector<Nonzero>& targetRows, HighsSolution& solution,
              HighsBasis& basis) const;
  };
  struct SingletonRow {
    double coef;
    HighsInt row;
    HighsInt col;
    bool colLowerTightened;
    bool colUpperTightened;

    void undo(const HighsPostsolveStack& postsolveStack,
              const HighsOptions& options, HighsSolution& solution,
              HighsBasis& basis) const;
  };

  // column fixed to lower or upper bound
  struct FixedCol {
    double fixValue;
    double colCost;
    HighsInt col;
    HighsBasisStatus fixType;

    void undo(const HighsPostsolveStack& postsolveStack,
              const HighsOptions& options,
              const std::vector<Nonzero>& colValues, HighsSolution& solution,
              HighsBasis& basis) const;
  };

  struct RedundantRow {
    HighsInt row;

    void undo(const HighsPostsolveStack& postsolveStack,
              const HighsOptions& options, HighsSolution& solution,
              HighsBasis& basis) const;
  };

  struct ForcingRow {
    double side;
    HighsInt row;
    RowType rowType;

    void undo(const HighsPostsolveStack& postsolveStack,
              const HighsOptions& options,
              const std::vector<Nonzero>& rowValues, HighsSolution& solution,
              HighsBasis& basis) const;
  };

  struct ForcingColumn {
    double colCost;
    double colBound;
    HighsInt col;
    bool atInfiniteUpper;
    bool colIntegral;

    void undo(const HighsPostsolveStack& postsolveStack,
              const HighsOptions& options,
              const std::vector<Nonzero>& colValues, HighsSolution& solution,
              HighsBasis& basis) const;
  };

  struct ForcingColumnRemovedRow {
    double rhs;
    HighsInt row;
    void undo(const HighsPostsolveStack& postsolveStack,
              const HighsOptions& options,
              const std::vector<Nonzero>& rowValues, HighsSolution& solution,
              HighsBasis& basis) const;
  };

  struct DuplicateRow {
    double duplicateRowScale;
    HighsInt duplicateRow;
    HighsInt row;
    bool rowLowerTightened;
    bool rowUpperTightened;

    void undo(const HighsPostsolveStack& postsolveStack,
              const HighsOptions& options, HighsSolution& solution,
              HighsBasis& basis) const;
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
              HighsBasis& basis) const;
    bool okMerge(const double tolerance) const;
    void undoFix(const HighsOptions& options, HighsSolution& solution,
                 const double mergeValue) const;
    void transformToPresolvedSpace(std::vector<double>& primalSol) const;
  };

  struct SlackColSubstitution {
    double rhs;
    HighsInt row;
    HighsInt col;

    void undo(const HighsPostsolveStack& postsolveStack,
              const HighsOptions& options,
              const std::vector<Nonzero>& rowValues, HighsSolution& solution,
              HighsBasis& basis);
  };

  static void undoFourierMotzkinBlock(HighsDataStack& stack,
                                      const HighsOptions& options,
                                      HighsSolution& solution,
                                      HighsBasis& basis);

  /// tags for reduction
  enum class ReductionType : uint8_t {
    kLinearTransform,
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
    kSlackColSubstitution,
    kFourierMotzkinBlock,
    kFourierMotzkinObjCol,
  };

  HighsDataStack reductionValues;
  std::vector<std::pair<ReductionType, size_t>> reductions;
  std::vector<HighsInt> origColIndex;
  std::vector<HighsInt> origRowIndex;
  std::vector<OrigRowType> origRowType;
  std::vector<uint8_t> linearlyTransformable;

  std::vector<Nonzero> rowValues;
  std::vector<Nonzero> colValues;
  HighsInt origNumCol = -1;
  HighsInt origNumRow = -1;
  HighsInt numAppendedRows = 0;
  HighsInt nextRowIndex = -1;
  HighsInt nextColIndex = -1;

  void reductionAdded(ReductionType type) {
    size_t position = reductionValues.getCurrentDataSize();
    reductions.emplace_back(type, position);
  }

  bool isModelRow(HighsInt row) const { return row < nextRowIndex; }

 public:
  const std::vector<HighsInt>& getOrigColIndex() const { return origColIndex; }

  const std::vector<HighsInt>& getOrigRowIndex() const { return origRowIndex; }

  bool isOrigCol(HighsInt col) const { return origColIndex[col] < origNumCol; }

  bool isOrigRow(HighsInt row) const {
    return origRowType[row] == OrigRowType::kOriginal;
  }

  bool isAppendedRow(HighsInt row) const {
    return origRowType[row] == OrigRowType::kAppended;
  }

  bool isCutRow(HighsInt row) const {
    return origRowType[row] == OrigRowType::kCut;
  }

  bool hasAppendedRows() const { return numAppendedRows > 0; }

  void appendToModel(HighsInt& numRows, HighsInt numRowsToAppend,
                     OrigRowType rowType) {
    if (numRowsToAppend <= 0) return;
    size_t currNumRow = origRowIndex.size();
    size_t newNumRow = currNumRow + numRowsToAppend;
    origRowIndex.resize(newNumRow);
    origRowType.resize(newNumRow, rowType);
    for (size_t i = currNumRow; i != newNumRow; ++i)
      origRowIndex[i] = nextRowIndex++;
    numRows += numRowsToAppend;
  }

  void appendCutsToModel(HighsInt numCuts) {
    appendToModel(origNumRow, numCuts, OrigRowType::kCut);
  }

  void appendRowsToModel(HighsInt numRows) {
    appendToModel(numAppendedRows, numRows, OrigRowType::kAppended);
  }

  void removeCutsFromModel(HighsInt numCuts) {
    if (numCuts <= 0) return;
    origNumRow -= numCuts;
    size_t newSize = 0;
    for (size_t i = 0; i < origRowIndex.size(); ++i) {
      if (origRowType[i] != OrigRowType::kCut) {
        if (i != newSize) {
          origRowIndex[newSize] = origRowIndex[i];
          origRowType[newSize] = origRowType[i];
        }
        ++newSize;
      }
    }
    origRowIndex.resize(newSize);
    origRowType.resize(newSize);
  }

  HighsInt getOrigNumRow() const { return origNumRow; }

  HighsInt getOrigNumCol() const { return origNumCol; }

  HighsInt getNextRowIndex() const { return nextRowIndex; }

  HighsInt getNextColIndex() const { return nextColIndex; }

  void appendColToModel() {
    origColIndex.push_back(nextColIndex++);
    linearlyTransformable.push_back(false);
  }

  void initializeIndexMaps(HighsInt numRow, HighsInt numCol);

  void compressIndexMaps(const std::vector<HighsInt>& newRowIndex,
                         const std::vector<HighsInt>& newColIndex);

  /// transform a column x by a linear mapping with a new column x'.
  /// I.e. substitute x = scale * x' + constant
  void linearTransform(HighsInt col, double scale, double constant) {
    reductionValues.push(LinearTransform{scale, constant, origColIndex[col]});
    reductionAdded(ReductionType::kLinearTransform);
  }

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
    reductionAdded(ReductionType::kFreeColSubstitution);
  }

  template <typename RowStorageFormat>
  void slackColSubstitution(HighsInt row, HighsInt col, double rhs,
                            const HighsMatrixSlice<RowStorageFormat>& rowVec) {
    rowValues.clear();
    for (const HighsSliceNonzero& rowVal : rowVec)
      rowValues.emplace_back(origColIndex[rowVal.index()], rowVal.value());

    reductionValues.push(
        SlackColSubstitution{rhs, origRowIndex[row], origColIndex[col]});
    reductionValues.push(rowValues);
    reductionAdded(ReductionType::kSlackColSubstitution);
  }

  template <typename ColStorageFormat>
  void doubletonEquation(HighsInt row, HighsInt colSubst, HighsInt col,
                         double coefSubst, double coef, double rhs,
                         double substLower, double substUpper, double substCost,
                         bool lowerTightened, bool upperTightened,
                         RowType rowType,
                         const HighsMatrixSlice<ColStorageFormat>& colVec) {
    colValues.clear();
    for (const HighsSliceNonzero& colVal : colVec)
      colValues.emplace_back(origRowIndex[colVal.index()], colVal.value());

    reductionValues.push(DoubletonEquation{
        coef, coefSubst, rhs, substLower, substUpper, substCost,
        row == -1 ? -1 : origRowIndex[row], origColIndex[colSubst],
        origColIndex[col], lowerTightened, upperTightened, rowType});
    reductionValues.push(colValues);
    reductionAdded(ReductionType::kDoubletonEquation);
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
    reductionAdded(ReductionType::kEqualityRowAddition);
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
    reductionAdded(ReductionType::kEqualityRowAdditions);
  }

  void singletonRow(HighsInt row, HighsInt col, double coef,
                    bool tightenedColLower, bool tightenedColUpper) {
    reductionValues.push(SingletonRow{coef, origRowIndex[row],
                                      origColIndex[col], tightenedColLower,
                                      tightenedColUpper});
    reductionAdded(ReductionType::kSingletonRow);
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
    reductionAdded(ReductionType::kFixedCol);
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
    reductionAdded(ReductionType::kFixedCol);
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
    reductionAdded(ReductionType::kFixedCol);
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
    reductionAdded(ReductionType::kFixedCol);
  }

  void redundantRow(HighsInt row) {
    reductionValues.push(RedundantRow{origRowIndex[row]});
    reductionAdded(ReductionType::kRedundantRow);
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
    reductionAdded(ReductionType::kForcingRow);
  }

  template <typename ColStorageFormat>
  void forcingColumn(HighsInt col,
                     const HighsMatrixSlice<ColStorageFormat>& colVec,
                     double cost, double boundVal, bool atInfiniteUpper,
                     bool colIntegral) {
    colValues.clear();
    for (const HighsSliceNonzero& colVal : colVec)
      colValues.emplace_back(origRowIndex[colVal.index()], colVal.value());

    reductionValues.push(ForcingColumn{cost, boundVal, origColIndex[col],
                                       atInfiniteUpper, colIntegral});
    reductionValues.push(colValues);
    reductionAdded(ReductionType::kForcingColumn);
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
    reductionAdded(ReductionType::kForcingColumnRemovedRow);
  }

  // Serialization layout for FM block (push order, so pop is reversed):
  // For each step (first eliminated to last):
  //   For each plus row: vector<Nonzero> entries
  //   vector<double> plusCoefs
  //   vector<FmeRowHeader> plusHeaders
  //   For each minus row: vector<Nonzero> entries
  //   vector<double> minusCoefs
  //   vector<FmeRowHeader> minusHeaders
  // Then (after all steps):
  //   For each step, for each parent: vector<FmeDescendant>
  //   For each step: FmeStepHeader
  //   numSteps (HighsInt)

  // Push one step's row data. Must be called before addToMatrix invalidates
  // the row slices. Returns (numPlus, numMinus) for later use.
  template <typename RowStorageFormat>
  std::pair<HighsInt, HighsInt> fourierMotzkinBlockPushStep(
      HighsInt col, const std::vector<FmeRowData<RowStorageFormat>>& plusRows,
      const std::vector<FmeRowData<RowStorageFormat>>& minusRows) {
    // push plus row entries
    std::vector<FmeRowHeader> plusHeaders;
    std::vector<double> plusCoefs;
    plusHeaders.reserve(plusRows.size());
    plusCoefs.reserve(plusRows.size());
    for (const auto& rd : plusRows) {
      std::vector<Nonzero> translated;
      double coef = 0.0;
      for (const HighsSliceNonzero& nz : rd.rowVec) {
        if (nz.index() == col)
          coef = nz.value();
        else
          translated.push_back({origColIndex[nz.index()], nz.value()});
      }
      reductionValues.push(translated);
      plusCoefs.push_back(coef);
      plusHeaders.push_back({origRowIndex[rd.row], rd.rowLower, rd.rowUpper});
    }
    reductionValues.push(plusCoefs);
    reductionValues.push(plusHeaders);

    // push minus row entries
    std::vector<FmeRowHeader> minusHeaders;
    std::vector<double> minusCoefs;
    minusHeaders.reserve(minusRows.size());
    minusCoefs.reserve(minusRows.size());
    for (const auto& rd : minusRows) {
      std::vector<Nonzero> translated;
      double coef = 0.0;
      for (const HighsSliceNonzero& nz : rd.rowVec) {
        if (nz.index() == col)
          coef = nz.value();
        else
          translated.push_back({origColIndex[nz.index()], nz.value()});
      }
      reductionValues.push(translated);
      minusCoefs.push_back(coef);
      minusHeaders.push_back({origRowIndex[rd.row], rd.rowLower, rd.rowUpper});
    }
    reductionValues.push(minusCoefs);
    reductionValues.push(minusHeaders);

    return {static_cast<HighsInt>(plusRows.size()),
            static_cast<HighsInt>(minusRows.size())};
  }

  // Finalize the FM block: push descendants mapping, new row origins,
  // and step headers. Called once after all elimination steps are complete.
  void fourierMotzkinBlockFinalize(
      const std::vector<HighsInt>& eliminatedCols,
      const std::vector<double>& colLowers,
      const std::vector<double>& colUppers, const std::vector<double>& colCosts,
      const std::vector<HighsInt>& numPlusPerStep,
      const std::vector<HighsInt>& numMinusPerStep,
      const std::vector<std::vector<std::vector<FmeDescendant>>>&
          descendantsAll,
      const std::vector<std::vector<FmeNewRow>>& newRowsAll) {
    HighsInt numSteps = static_cast<HighsInt>(eliminatedCols.size());

    // push descendants for each step's parents
    for (HighsInt s = 0; s < numSteps; ++s) {
      HighsInt numParents = numPlusPerStep[s] + numMinusPerStep[s];
      assert(static_cast<HighsInt>(descendantsAll[s].size()) == numParents);
      for (HighsInt p = 0; p < numParents; ++p)
        reductionValues.push(descendantsAll[s][p]);
    }

    // push new row origins for each step (translate row to orig space)
    for (HighsInt s = 0; s < numSteps; ++s) {
      std::vector<FmeNewRow> translated;
      translated.reserve(newRowsAll[s].size());
      for (const auto& nr : newRowsAll[s])
        translated.push_back(
            {origRowIndex[nr.row], nr.plusParentIdx, nr.minusParentIdx});
      reductionValues.push(translated);
    }

    // push step headers
    for (HighsInt s = 0; s < numSteps; ++s) {
      FmeStepHeader header{colLowers[s],
                           colUppers[s],
                           colCosts[s],
                           origColIndex[eliminatedCols[s]],
                           numPlusPerStep[s],
                           numMinusPerStep[s],
                           static_cast<HighsInt>(newRowsAll[s].size())};
      reductionValues.push(header);
    }

    reductionValues.push(numSteps);
    reductionAdded(ReductionType::kFourierMotzkinBlock);
  }

  void fourierMotzkinObjCol(HighsInt col, double offset,
                            const std::vector<Nonzero>& costEntries) {
    reductionValues.push(FourierMotzkinObjCol{offset, origColIndex[col]});
    std::vector<Nonzero> translatedEntries;
    translatedEntries.reserve(costEntries.size());
    for (const Nonzero& entry : costEntries)
      if (entry.index != col)
        translatedEntries.emplace_back(origColIndex[entry.index], entry.value);
    reductionValues.push(translatedEntries);
    reductionAdded(ReductionType::kFourierMotzkinObjCol);
  }

  void duplicateRow(HighsInt row, bool rowUpperTightened,
                    bool rowLowerTightened, HighsInt duplicateRow,
                    double duplicateRowScale) {
    reductionValues.push(
        DuplicateRow{duplicateRowScale, origRowIndex[duplicateRow],
                     origRowIndex[row], rowLowerTightened, rowUpperTightened});
    reductionAdded(ReductionType::kDuplicateRow);
  }

  bool duplicateColumn(double colScale, double colLower, double colUpper,
                       double duplicateColLower, double duplicateColUpper,
                       HighsInt col, HighsInt duplicateCol, bool colIntegral,
                       bool duplicateColIntegral,
                       const double ok_merge_tolerance) {
    const HighsInt origCol = origColIndex[col];
    const HighsInt origDuplicateCol = origColIndex[duplicateCol];
    DuplicateColumn debug_values = {
        colScale,          colLower,          colUpper,
        duplicateColLower, duplicateColUpper, origCol,
        origDuplicateCol,  colIntegral,       duplicateColIntegral};
    const bool ok_merge = debug_values.okMerge(ok_merge_tolerance);
    const bool prevent_illegal_merge = true;
    if (!ok_merge && prevent_illegal_merge) return false;
    reductionValues.push(debug_values);
    //    reductionValues.push(DuplicateColumn{
    //        colScale, colLower, colUpper, duplicateColLower,
    //        duplicateColUpper, origCol, origDuplicateCol, colIntegral,
    //        duplicateColIntegral});

    reductionAdded(ReductionType::kDuplicateColumn);

    // mark columns as not linearly transformable
    linearlyTransformable[origCol] = false;
    linearlyTransformable[origDuplicateCol] = false;
    return true;
  }

  std::vector<double> getReducedPrimalSolution(
      const std::vector<double>& origPrimalSolution) {
    std::vector<double> reducedSolution = origPrimalSolution;
    reducedSolution.resize(nextColIndex, 0.0);

    for (const std::pair<ReductionType, size_t>& primalColTransformation :
         reductions) {
      switch (primalColTransformation.first) {
        case ReductionType::kDuplicateColumn: {
          DuplicateColumn duplicateColReduction;
          reductionValues.setPosition(primalColTransformation.second);
          reductionValues.pop(duplicateColReduction);
          duplicateColReduction.transformToPresolvedSpace(reducedSolution);
          break;
        }
        case ReductionType::kLinearTransform: {
          reductionValues.setPosition(primalColTransformation.second);
          LinearTransform linearTransform;
          reductionValues.pop(linearTransform);
          linearTransform.transformToPresolvedSpace(reducedSolution);
          break;
        }
        case ReductionType::kFourierMotzkinObjCol: {
          reductionValues.setPosition(primalColTransformation.second);
          std::vector<Nonzero> costEntries;
          reductionValues.pop(costEntries);
          FourierMotzkinObjCol fmObjCol;
          reductionValues.pop(fmObjCol);
          fmObjCol.transformToPresolvedSpace(costEntries, reducedSolution);
          break;
        }
        default:
          continue;
      }
    }

    size_t reducedNumCol = origColIndex.size();
    for (size_t i = 0; i < reducedNumCol; ++i)
      reducedSolution[i] = reducedSolution[origColIndex[i]];

    reducedSolution.resize(reducedNumCol);
    return reducedSolution;
  }

  bool isColLinearlyTransformable(HighsInt col) const {
    assert(col >= 0);
    assert(static_cast<size_t>(col) < origColIndex.size());
    return (linearlyTransformable[origColIndex[col]] != 0);
  }

  template <typename T>
  void undoIterateBackwards(std::vector<T>& values,
                            const std::vector<HighsInt>& index,
                            HighsInt origSize) {
    values.resize(origSize);
#ifdef DEBUG_EXTRA
    // Fill vector with NaN for debugging purposes
    std::vector<T> valuesNew;
    valuesNew.resize(origSize, std::numeric_limits<T>::signaling_NaN());
    for (size_t i = index.size(); i > 0; --i) {
      assert(static_cast<size_t>(index[i - 1]) >= i - 1);
      valuesNew[index[i - 1]] = values[i - 1];
    }
    std::copy(valuesNew.cbegin(), valuesNew.cend(), values.begin());
#else
    for (size_t i = index.size(); i > 0; --i) {
      assert(static_cast<size_t>(index[i - 1]) >= i - 1);
      values[index[i - 1]] = values[i - 1];
    }
#endif
  }

  /// check if vector contains NaN or Inf
  bool containsNanOrInf(const std::vector<double>& v) const {
    return std::find_if(v.cbegin(), v.cend(), [](const double& d) {
             return (std::isnan(d) || std::isinf(d));
           }) != v.cend();
  }

  /// undo presolve steps for primal dual solution and basis
  void undo(const HighsOptions& options, HighsSolution& solution,
            HighsBasis& basis, size_t numReductions = 0,
            const HighsInt report_col = -1) {
    reductionValues.resetPosition();

    // Verify that undo can be performed
    assert(solution.value_valid);
    bool perform_dual_postsolve = solution.dual_valid;
    bool perform_basis_postsolve = basis.valid;

    // expand solution to original index space
    assert(nextColIndex > 0);
    undoIterateBackwards(solution.col_value, origColIndex, nextColIndex);

    assert(nextRowIndex >= 0);
    undoIterateBackwards(solution.row_value, origRowIndex, nextRowIndex);

    if (perform_dual_postsolve) {
      // if dual solution is given, expand dual solution and basis to original
      // index space
      undoIterateBackwards(solution.col_dual, origColIndex, nextColIndex);

      undoIterateBackwards(solution.row_dual, origRowIndex, nextRowIndex);
    }

    if (perform_basis_postsolve) {
      // if basis is given, expand basis status values to original index space
      undoIterateBackwards(basis.col_status, origColIndex, nextColIndex);

      undoIterateBackwards(basis.row_status, origRowIndex, nextRowIndex);
    }

    // now undo the changes
    for (size_t i = reductions.size(); i > numReductions; --i) {
      if (report_col >= 0)
        printf("Before  reduction %2d (type %2d): col_value[%2d] = %g\n",
               int(i - 1), int(reductions[i - 1].first), int(report_col),
               solution.col_value[report_col]);
      switch (reductions[i - 1].first) {
        case ReductionType::kLinearTransform: {
          LinearTransform reduction;
          reductionValues.pop(reduction);
          reduction.undo(options, solution);
          break;
        }
        case ReductionType::kFreeColSubstitution: {
          FreeColSubstitution reduction;
          reductionValues.pop(colValues);
          reductionValues.pop(rowValues);
          reductionValues.pop(reduction);
          reduction.undo(*this, options, rowValues, colValues, solution, basis);
          break;
        }
        case ReductionType::kDoubletonEquation: {
          DoubletonEquation reduction;
          reductionValues.pop(colValues);
          reductionValues.pop(reduction);
          reduction.undo(*this, options, colValues, solution, basis);
          break;
        }
        case ReductionType::kEqualityRowAddition: {
          EqualityRowAddition reduction;
          reductionValues.pop(rowValues);
          reductionValues.pop(reduction);
          reduction.undo(*this, options, rowValues, solution, basis);
          break;
        }
        case ReductionType::kEqualityRowAdditions: {
          EqualityRowAdditions reduction;
          reductionValues.pop(colValues);
          reductionValues.pop(rowValues);
          reductionValues.pop(reduction);
          reduction.undo(*this, options, rowValues, colValues, solution, basis);
          break;
        }
        case ReductionType::kSingletonRow: {
          SingletonRow reduction;
          reductionValues.pop(reduction);
          reduction.undo(*this, options, solution, basis);
          break;
        }
        case ReductionType::kFixedCol: {
          FixedCol reduction;
          reductionValues.pop(colValues);
          reductionValues.pop(reduction);
          reduction.undo(*this, options, colValues, solution, basis);
          break;
        }
        case ReductionType::kRedundantRow: {
          RedundantRow reduction;
          reductionValues.pop(reduction);
          reduction.undo(*this, options, solution, basis);
          break;
        }
        case ReductionType::kForcingRow: {
          ForcingRow reduction;
          reductionValues.pop(rowValues);
          reductionValues.pop(reduction);
          reduction.undo(*this, options, rowValues, solution, basis);
          break;
        }
        case ReductionType::kForcingColumn: {
          ForcingColumn reduction;
          reductionValues.pop(colValues);
          reductionValues.pop(reduction);
          reduction.undo(*this, options, colValues, solution, basis);
          break;
        }
        case ReductionType::kForcingColumnRemovedRow: {
          ForcingColumnRemovedRow reduction;
          reductionValues.pop(rowValues);
          reductionValues.pop(reduction);
          reduction.undo(*this, options, rowValues, solution, basis);
          break;
        }
        case ReductionType::kDuplicateRow: {
          DuplicateRow reduction;
          reductionValues.pop(reduction);
          reduction.undo(*this, options, solution, basis);
          break;
        }
        case ReductionType::kDuplicateColumn: {
          DuplicateColumn reduction;
          reductionValues.pop(reduction);
          reduction.undo(options, solution, basis);
          break;
        }
        case ReductionType::kSlackColSubstitution: {
          SlackColSubstitution reduction;
          reductionValues.pop(rowValues);
          reductionValues.pop(reduction);
          reduction.undo(*this, options, rowValues, solution, basis);
          break;
        }
        case ReductionType::kFourierMotzkinBlock: {
          undoFourierMotzkinBlock(reductionValues, options, solution, basis);
          break;
        }
        case ReductionType::kFourierMotzkinObjCol: {
          std::vector<Nonzero> costEntries;
          reductionValues.pop(costEntries);
          FourierMotzkinObjCol reduction;
          reductionValues.pop(reduction);
          reduction.undo(costEntries, solution);
          break;
        }
        default:
          printf("Reduction case %d not handled\n",
                 int(reductions[i - 1].first));
          if (kAllowDeveloperAssert) assert(1 == 0);
      }
    }
    if (report_col >= 0)
      printf("After last reduction: col_value[%2d] = %g\n", int(report_col),
             solution.col_value[report_col]);

    solution.col_value.resize(origNumCol);
    if (perform_dual_postsolve) solution.col_dual.resize(origNumCol);
    if (perform_basis_postsolve) basis.col_status.resize(origNumCol);

    solution.row_value.resize(origNumRow);
    if (perform_dual_postsolve) solution.row_dual.resize(origNumRow);
    if (perform_basis_postsolve) basis.row_status.resize(origNumRow);

#ifdef DEBUG_EXTRA
    // solution should not contain NaN or Inf
    assert(!containsNanOrInf(solution.col_value));
    // row values are not determined by postsolve
    // assert(!containsNanOrInf(solution.row_value));
    assert(!containsNanOrInf(solution.col_dual));
    assert(!containsNanOrInf(solution.row_dual));
#endif
  }

  /// undo presolve steps for primal solution
  void undoPrimal(const HighsOptions& options, HighsSolution& solution,
                  const HighsInt report_col = -1) {
    // Call to reductionValues.resetPosition(); seems unnecessary as
    // it's the first thing done in undo
    reductionValues.resetPosition();
    HighsBasis basis;
    basis.valid = false;
    solution.dual_valid = false;
    undo(options, solution, basis, 0, report_col);
  }

  /*
    // Not used
  /// undo presolve steps for primal and dual solution
  void undoPrimalDual(const HighsOptions& options, HighsSolution& solution) {
    reductionValues.resetPosition();
    HighsBasis basis;
    basis.valid = false;
    assert(solution.value_valid);
    assert(solution.dual_valid);
    undo(options, solution, basis);
  }
  */

  size_t numReductions() const { return reductions.size(); }
};

}  // namespace presolve

#endif
