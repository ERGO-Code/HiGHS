/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file HighsSymmetry.cpp
 * @brief Facilities for symmetry detection
 * @author Leona Gottwald
 */

#include "presolve/HighsSymmetry.h"

#include <algorithm>
#include <numeric>

void HighsSymmetryDetection::removeFixPoints() {
  colPartition.erase(std::remove_if(colPartition.begin(), colPartition.end(),
                                    [&](HighsInt col) {
                                      return colCellSize(colToCell[col]) == 1;
                                    }),
                     colPartition.end());

  if ((HighsInt)colPartition.size() < numCol) {
    numCol = colPartition.size();
    colPartitionLinks.resize(numCol);

    HighsInt cellStart = 0;
    HighsInt cellNumber = 0;
    for (HighsInt i = 0; i < numCol; ++i) {
      HighsInt col = colPartition[i];
      // if the cell number is different to the current cell number this is the
      // start of a new cell
      if (cellNumber != colToCell[col]) {
        // remember the number of this cell to indetify its end
        cellNumber = colToCell[col];
        // set the link of the cell start to point to its end
        colPartitionLinks[cellStart] = i;
        // remember start of this cell
        cellStart = i;
      }

      // correct the colToCell array to not store the start index of the
      // cell, not its number
      colToCell[col] = cellStart;
      // set the link of the column to the cellStart
      colPartitionLinks[i] = cellStart;
    }

    // set the column partition link of the last started cell to point past the
    // end
    colPartitionLinks[cellStart] = numCol;
  }

  rowPartition.erase(std::remove_if(rowPartition.begin(), rowPartition.end(),
                                    [&](HighsInt col) {
                                      return rowCellSize(rowToCell[col]) == 1;
                                    }),
                     rowPartition.end());

  if ((HighsInt)rowPartition.size() < numRow) {
    numRow = rowPartition.size();
    rowPartitionLinks.resize(numRow);
    HighsInt cellStart = 0;
    HighsInt cellNumber = 0;
    for (HighsInt i = 0; i < numRow; ++i) {
      HighsInt row = rowPartition[i];
      // if the cell number is different to the current cell number this is the
      // start of a new cell
      if (cellNumber != rowToCell[row]) {
        // remember the number of this cell to indetify its end
        cellNumber = rowToCell[row];
        // set the link of the cell start to point to its end
        rowPartitionLinks[cellStart] = i;
        // remember start of this cell
        cellStart = i;
      }

      // correct the rowToCell array to not store the start index of the cell,
      // not its number
      rowToCell[row] = cellStart;
      // set the link of the column to the cellStart
      rowPartitionLinks[i] = cellStart;
    }

    // set the column partition link of the last started cell to point past the
    // end
    rowPartitionLinks[cellStart] = numRow;
  }
}

void HighsSymmetryDetection::initializeGroundSet() {
  groundSet = colPartition;
  std::sort(groundSet.begin(), groundSet.end());
  colPosition.resize(colToCell.size(), -1);
  for (HighsInt i = 0; i < numCol; ++i) colPosition[groundSet[i]] = i;

  orbitPartition.resize(numCol);
  std::iota(orbitPartition.begin(), orbitPartition.end(), 0);
  orbitSize.resize(numCol, 1);

  permutation.resize(numCol);
  currNodeCertificate.resize(numCol + numRow);
}

bool HighsSymmetryDetection::mergeOrbits(HighsInt col1, HighsInt col2) {
  if (col1 == col2) return false;

  HighsInt orbit1 = getOrbit(col1);
  HighsInt orbit2 = getOrbit(col2);

  if (orbit1 == orbit2) return false;

  if (orbitSize[orbit2] <= orbitSize[orbit1]) {
    orbitPartition[orbit2] = orbitPartition[orbit1];
    orbitSize[orbit1] += orbitSize[orbit2];
  } else {
    orbitPartition[orbit1] = orbitPartition[orbit2];
    orbitSize[orbit2] += orbitSize[orbit1];
  }

  return true;
}

HighsInt HighsSymmetryDetection::getOrbit(HighsInt col) {
  HighsInt i = colPosition[col];
  HighsInt orbit = orbitPartition[i];
  if (orbit != orbitPartition[orbit]) {
    do {
      linkCompressionStack.push_back(i);
      i = orbit;
      orbit = orbitPartition[orbit];
    } while (orbit != orbitPartition[orbit]);

    do {
      i = linkCompressionStack.back();
      linkCompressionStack.pop_back();
      orbitPartition[i] = orbit;
    } while (!linkCompressionStack.empty());
  }

  return orbit;
}

void HighsSymmetryDetection::initializeHashValues() {
  colHashes.resize(numCol);
  rowHashes.resize(numRow);
  for (HighsInt i = 0; i != numCol; ++i) {
    HighsInt cell = colToCell[i];
    HighsInt cellSize = colCellSize(cell);
    if (cellSize == 1) continue;

    if (cellsOnRefinementStack.insert(cell, kColCell))
      refinementStack.emplace_back(cell, kColCell);
    for (HighsInt j = Astart[i]; j != Astart[i + 1]; ++j) {
      HighsInt rowCell = rowToCell[Aindex[j]];
      HighsHashHelpers::sparse_combine(colHashes[i], rowCell, Acolor[j]);
      HighsHashHelpers::sparse_combine(rowHashes[Aindex[j]], cell, Acolor[j]);
    }
  }

  for (HighsInt i = 0; i != numRow; ++i) {
    HighsInt cell = rowToCell[i];
    HighsInt cellSize = rowPartitionLinks[cell] - cell;
    if (cellSize == 1) continue;
    if (cellsOnRefinementStack.insert(cell, kRowCell))
      refinementStack.emplace_back(cell, kRowCell);
  }
}

void HighsSymmetryDetection::partitionRefinement() {
  while (!refinementStack.empty()) {
    auto cell = refinementStack.back();
    refinementStack.pop_back();

    HighsInt cellStart = cell.first;

    if (cell.second) {
      // col cell
      cellsOnRefinementStack.erase(cell);
      if (colCellSize(cellStart) == 1) continue;
      HighsInt cellEnd = colPartitionLinks[cellStart];
      assert(cellEnd >= cellStart);

      std::sort(&colPartition[cellStart], &colPartition[cellEnd],
                [&](HighsInt col1, HighsInt col2) {
                  return colHashes[col1] < colHashes[col2];
                });

      for (HighsInt i = cellStart + 1; i < cellEnd; ++i) {
        HighsInt col = colPartition[i];
        if (colHashes[col] != colHashes[colPartition[i - 1]]) {
          colPartitionLinks[cellStart] = i;
          cellStart = i;
          // remember that a new column cell was created
          cellCreationStack.emplace_back(i, kColCell);
        }

        if (cellStart != colPartitionLinks[i]) {
          // set new cell id
          assert(colToCell[col] == cell.first);
          colToCell[col] = cellStart;
          colPartitionLinks[i] = cellStart;

          // update hashes of affected rows
          for (HighsInt j = Astart[col]; j != Astart[col + 1]; ++j) {
            HighsInt row = Aindex[j];
            // remove hash contribution with old cell index and add with new
            // cell index
            HighsHashHelpers::sparse_inverse_combine(rowHashes[row], cell.first,
                                                     Acolor[j]);
            HighsHashHelpers::sparse_combine(rowHashes[row], cellStart,
                                             Acolor[j]);

            HighsInt rowCell = rowToCell[row];
            if (rowCellSize(rowCell) > 1 &&
                cellsOnRefinementStack.insert(rowCell, kRowCell))
              refinementStack.emplace_back(rowCell, kRowCell);
          }
        }
      }

      // set the link of the last started cell to point to the cell end
      colPartitionLinks[cellStart] = cellEnd;
    } else {
      // row partition
      cellsOnRefinementStack.erase(cell);
      if (rowCellSize(cellStart) == 1) continue;

      HighsInt cellEnd = rowPartitionLinks[cellStart];
      assert(cellEnd >= cellStart);

      std::sort(&rowPartition[cellStart], &rowPartition[cellEnd],
                [&](HighsInt row1, HighsInt row2) {
                  return rowHashes[row1] < rowHashes[row2];
                });

      for (HighsInt i = cellStart + 1; i < cellEnd; ++i) {
        HighsInt row = rowPartition[i];
        if (rowHashes[row] != rowHashes[rowPartition[i - 1]]) {
          rowPartitionLinks[cellStart] = i;
          cellStart = i;
          // remember that a new row cell was created
          cellCreationStack.emplace_back(i, kRowCell);
        }

        if (cellStart != rowPartitionLinks[i]) {
          // set new cell id
          assert(rowToCell[row] == cell.first);
          rowToCell[row] = cellStart;
          rowPartitionLinks[i] = cellStart;

          // update hashes of affected columns
          for (HighsInt j = ARstart[row]; j != ARstart[row + 1]; ++j) {
            HighsInt col = ARindex[j];
            // remove hash contribution with old cell index and add with new
            // cell index
            HighsHashHelpers::sparse_inverse_combine(colHashes[col], cell.first,
                                                     ARcolor[j]);
            HighsHashHelpers::sparse_combine(colHashes[col], cellStart,
                                             ARcolor[j]);

            HighsInt colCell = colToCell[col];
            if (colCellSize(colCell) > 1 &&
                cellsOnRefinementStack.insert(colCell, kColCell))
              refinementStack.emplace_back(colCell, kColCell);
          }
        }
      }

      // set the link of the last started cell to point to the cell end
      rowPartitionLinks[cellStart] = cellEnd;
    }
  }
}

std::pair<HighsInt, HighsSymmetryDetection::CellType>
HighsSymmetryDetection::selectTargetCell() {
  HighsInt i = 0;
  if (nodeStack.size() > 1) {
    i = nodeStack[nodeStack.size() - 2].targetCell + 1;
    if (nodeStack[nodeStack.size() - 2].targetCellType == kRowCell)
      goto rowTargetCell;
  }
  while (i < numCol) {
    currNodeCertificate[i] = colHashes[colPartition[i]];

    // employ prefix pruning scheme as in bliss
    if (!firstLeaveCertificate.empty() &&
        currNodeCertificate[i] != firstLeaveCertificate[i] &&
        currNodeCertificate[i] > smallestLeaveCertificate[i])
      return std::make_pair(-2, CellType());

    HighsInt cellEnd = colPartitionLinks[i];
    if (cellEnd - i > 1) return std::make_pair(i, kColCell);

    ++i;
  }

  i = 0;
rowTargetCell:
  while (i < numRow) {
    HighsInt k = i + numCol;
    currNodeCertificate[k] = rowHashes[rowPartition[i]];

    if (!firstLeaveCertificate.empty() &&
        currNodeCertificate[k] != firstLeaveCertificate[k] &&
        currNodeCertificate[k] > smallestLeaveCertificate[k])
      return std::make_pair(-2, CellType());

    HighsInt cellEnd = rowPartitionLinks[i];
    if (cellEnd - i > 1) return std::make_pair(i, kRowCell);

    ++i;
  }

  return std::make_pair(-1, CellType());
}

bool HighsSymmetryDetection::determineNextToDistinguish() {
  Node& currNode = nodeStack.back();
  distinguishCands.clear();
  std::vector<HighsInt>::iterator cellStart;
  std::vector<HighsInt>::iterator cellEnd;
  if (currNode.targetCellType == kColCell) {
    cellStart = colPartition.begin() + currNode.targetCell;
    cellEnd = colPartition.begin() + colPartitionLinks[currNode.targetCell];
  } else {
    cellStart = rowPartition.begin() + currNode.targetCell;
    cellEnd = rowPartition.begin() + rowPartitionLinks[currNode.targetCell];
  }

  if (currNode.lastDistiguished == -1) {
    auto nextDistinguishPos = std::min_element(cellStart, cellEnd);
    distinguishCands.push_back(&*nextDistinguishPos);
  } else if (!currNode.onFirstPath) {
    for (auto i = cellStart; i != cellEnd; ++i) {
      if (*i > currNode.lastDistiguished) distinguishCands.push_back(&*i);
    }
    if (distinguishCands.empty()) return false;
    auto nextDistinguishPos =
        std::min_element(distinguishCands.begin(), distinguishCands.end(),
                         [](HighsInt* a, HighsInt* b) { return *a < *b; });
    std::swap(*distinguishCands.begin(), *nextDistinguishPos);
    distinguishCands.resize(1);
  } else {
    HighsHashTable<HighsInt> prunedOrbits;
    for (auto i = cellStart; i != cellEnd; ++i) {
      if (*i > currNode.lastDistiguished)
        distinguishCands.push_back(&*i);
      else
        prunedOrbits.insert(getOrbit(*i));
    }
    if (distinguishCands.empty()) return false;
    std::make_heap(distinguishCands.begin(), distinguishCands.end(),
                   [](HighsInt* a, HighsInt* b) { return *a > *b; });
    while (true) {
      // check if an already distinguished vertex is in the same orbit as
      // this vertex
      HighsInt orbit = getOrbit(*distinguishCands.front());
      if (prunedOrbits.find(orbit)) {
        // if such a vertex was found we prune this sub tree und look for
        // the next node to distinguish
        std::pop_heap(distinguishCands.begin(), distinguishCands.end(),
                      [](HighsInt* a, HighsInt* b) { return *a > *b; });
        distinguishCands.pop_back();
        if (distinguishCands.empty()) return false;
        continue;
      }

      break;
    }

    distinguishCands.resize(1);
  }

  return true;
}

void HighsSymmetryDetection::distinguishVertex(HighsInt targetCell,
                                               CellType cellType) {
  if (cellType == kColCell) {
    HighsInt targetCellEnd = colPartitionLinks[targetCell];
    std::swap(*distinguishCands[0], colPartition[targetCell]);
    nodeStack.back().lastDistiguished = colPartition[targetCell];
    currNodeCertificate[targetCell] = colHashes[colPartition[targetCell]];

    HighsInt newCell = targetCell + 1;
    colPartitionLinks[targetCell] = newCell;
    cellCreationStack.emplace_back(newCell, kColCell);

    for (HighsInt i = newCell; i < targetCellEnd; ++i) {
      HighsInt col = colPartition[i];
      assert(colToCell[col] == targetCell);
      colToCell[col] = newCell;
      colPartitionLinks[i] = newCell;

      // update hashes of affected rows
      for (HighsInt j = Astart[col]; j != Astart[col + 1]; ++j) {
        HighsInt row = Aindex[j];
        // remove hash contribution with old cell index and add with new
        // cell index
        HighsHashHelpers::sparse_inverse_combine(rowHashes[row], targetCell,
                                                 Acolor[j]);
        HighsHashHelpers::sparse_combine(rowHashes[row], newCell, Acolor[j]);

        HighsInt rowCell = rowToCell[row];
        if (rowCellSize(rowCell) > 1 &&
            cellsOnRefinementStack.insert(rowCell, kRowCell))
          refinementStack.emplace_back(rowCell, kRowCell);
      }
    }

    colPartitionLinks[newCell] = targetCellEnd;
  } else {
    HighsInt targetCellEnd = rowPartitionLinks[targetCell];
    std::swap(*distinguishCands[0], rowPartition[targetCell]);
    nodeStack.back().lastDistiguished = rowPartition[targetCell];
    currNodeCertificate[numCol + targetCell] =
        rowHashes[rowPartition[targetCell]];

    HighsInt newCell = targetCell + 1;
    rowPartitionLinks[targetCell] = newCell;
    cellCreationStack.emplace_back(newCell, kRowCell);

    for (HighsInt i = newCell; i < targetCellEnd; ++i) {
      HighsInt row = rowPartition[i];
      assert(rowToCell[row] == targetCell);
      rowToCell[row] = newCell;
      rowPartitionLinks[i] = newCell;

      // update hashes of affected columns
      for (HighsInt j = ARstart[row]; j != ARstart[row + 1]; ++j) {
        HighsInt col = ARindex[j];
        // remove hash contribution with old cell index and add with new
        // cell index
        HighsHashHelpers::sparse_inverse_combine(colHashes[col], targetCell,
                                                 ARcolor[j]);
        HighsHashHelpers::sparse_combine(colHashes[col], newCell, ARcolor[j]);

        HighsInt colCell = colToCell[col];
        if (colCellSize(colCell) > 1 &&
            cellsOnRefinementStack.insert(colCell, kColCell))
          refinementStack.emplace_back(colCell, kColCell);
      }
    }
    rowPartitionLinks[newCell] = targetCellEnd;
  }
}

void HighsSymmetryDetection::backtrack(HighsInt backtrackStackNewEnd,
                                       HighsInt backtrackStackEnd) {
  // we assume that we always backtrack from a leave node, i.e. a discrete
  // partition therefore we do not need to remember the values of the hash
  // contributions as it is the indentity for each position and all new cells
  // are on the cell creation stack.
  for (HighsInt stackPos = backtrackStackEnd - 1;
       stackPos >= backtrackStackNewEnd; --stackPos) {
    HighsInt cell = cellCreationStack[stackPos].first;
    CellType cellType = cellCreationStack[stackPos].second;
    if (cellType == kColCell) {
      // look up the cell start of the preceding cell with link compression
      HighsInt newStart = getColCellStart(cell - 1);
      // remember the current end
      HighsInt currEnd = colPartitionLinks[cell];
      // change the link to point to the start of the preceding cell
      colPartitionLinks[cell] = newStart;
      // change the link of the start pointer of the preceding cell to point to
      // the end of this cell
      colPartitionLinks[newStart] = currEnd;
    } else {
      HighsInt newStart = getRowCellStart(cell - 1);
      HighsInt currEnd = rowPartitionLinks[cell];
      rowPartitionLinks[cell] = newStart;
      rowPartitionLinks[newStart] = currEnd;
    }
  }
}

void HighsSymmetryDetection::cleanupBacktrack(HighsInt cellCreationStackPos) {
  // the links have been updated. Even though they might still not be fully
  // compressed the cell starts will all point to the correct cell end and the
  // lookup with path compression will give the correct start
  for (HighsInt stackPos = cellCreationStack.size() - 1;
       stackPos >= cellCreationStackPos; --stackPos) {
    HighsInt cell = cellCreationStack[stackPos].first;
    CellType cellType = cellCreationStack[stackPos].second;
    if (cellType == kColCell) {
      HighsInt cellStart = getColCellStart(cell);
      HighsInt cellEnd = colPartitionLinks[cellStart];

      for (HighsInt c = cell; c < cellEnd && colToCell[colPartition[c]] == cell;
           ++c) {
        colPartitionLinks[c] = cellStart;

        HighsInt col = colPartition[c];
        colToCell[col] = cellStart;

        // update hashes of affected columns
        for (HighsInt j = Astart[col]; j != Astart[col + 1]; ++j) {
          HighsInt row = Aindex[j];
          // remove hash contribution with old cell index and add with new
          // cell index
          HighsHashHelpers::sparse_inverse_combine(rowHashes[row], cell,
                                                   Acolor[j]);
          HighsHashHelpers::sparse_combine(rowHashes[row], cellStart,
                                           Acolor[j]);
        }
      }
    } else {
      HighsInt cellStart = getRowCellStart(cell);
      HighsInt cellEnd = rowPartitionLinks[cellStart];
      for (HighsInt c = cell; c < cellEnd && rowToCell[rowPartition[c]] == cell;
           ++c) {
        rowPartitionLinks[c] = cellStart;

        HighsInt row = rowPartition[c];
        rowToCell[row] = cellStart;

        // update hashes of affected columns
        for (HighsInt j = ARstart[row]; j != ARstart[row + 1]; ++j) {
          HighsInt col = ARindex[j];
          // remove hash contribution with old cell index and add with new
          // cell index
          HighsHashHelpers::sparse_inverse_combine(colHashes[col], cell,
                                                   ARcolor[j]);
          HighsHashHelpers::sparse_combine(colHashes[col], cellStart,
                                           ARcolor[j]);
        }
      }
    }
  }

  cellCreationStack.resize(cellCreationStackPos);
}

HighsInt HighsSymmetryDetection::getColCellStart(HighsInt pos) {
  HighsInt startPos = colPartitionLinks[pos];
  if (startPos > pos) return pos;
  if (colPartitionLinks[startPos] < startPos) {
    do {
      linkCompressionStack.push_back(pos);
      pos = startPos;
      startPos = colPartitionLinks[startPos];
    } while (colPartitionLinks[startPos] < startPos);

    do {
      colPartitionLinks[linkCompressionStack.back()] = startPos;
      linkCompressionStack.pop_back();
    } while (!linkCompressionStack.empty());
  }

  return startPos;
}

HighsInt HighsSymmetryDetection::getRowCellStart(HighsInt pos) {
  HighsInt startPos = rowPartitionLinks[pos];
  if (startPos > pos) return pos;
  if (rowPartitionLinks[startPos] < startPos) {
    do {
      linkCompressionStack.push_back(pos);
      pos = startPos;
      startPos = rowPartitionLinks[startPos];
    } while (rowPartitionLinks[startPos] < startPos);

    do {
      rowPartitionLinks[linkCompressionStack.back()] = startPos;
      linkCompressionStack.pop_back();
    } while (!linkCompressionStack.empty());
  }

  return startPos;
}

void HighsSymmetryDetection::createNode() {
  nodeStack.emplace_back();
  nodeStack.back().stackStart = cellCreationStack.size();
  nodeStack.back().targetCell = -1;
  nodeStack.back().lastDistiguished = -1;
  nodeStack.back().onFirstPath = firstLeaveColPartition.empty();
}

struct MatrixColumn {
  uint32_t cost;
  uint32_t lb;
  uint32_t ub;
  uint32_t integral;
  uint32_t len;

  bool operator==(const MatrixColumn& other) const {
    return std::memcmp(this, &other, sizeof(MatrixColumn)) == 0;
  }
};

struct MatrixRow {
  uint32_t lb;
  uint32_t ub;
  uint32_t len;

  bool operator==(const MatrixRow& other) const {
    return std::memcmp(this, &other, sizeof(MatrixRow)) == 0;
  }
};

void HighsSymmetryDetection::loadModelAsGraph(const HighsLp& lp,
                                              double epsilon) {
  numCol = lp.numCol_;
  numRow = lp.numRow_;
  HighsHashTable<MatrixColumn, HighsInt> columnSet;
  HighsHashTable<MatrixRow, HighsInt> rowSet;
  HighsMatrixColoring coloring(epsilon);

  // set up row and column based incidence matrix
  Aindex = lp.Aindex_;
  Astart = lp.Astart_;

  HighsInt numNz = lp.Aindex_.size();
  Acolor.resize(numNz);

  // set up the column colors and count row sizes
  std::vector<HighsInt> rowSizes(numRow);
  for (HighsInt i = 0; i < numCol; ++i) {
    for (HighsInt j = Astart[i]; j < Astart[i + 1]; ++j) {
      Acolor[j] = coloring.color(lp.Avalue_[j]);
      rowSizes[Aindex[j]] += 1;
    }
  }

  // next set up the row starts using the computed row sizes
  ARstart.resize(numRow + 1);

  HighsInt offset = 0;
  for (HighsInt i = 0; i < numRow; ++i) {
    ARstart[i] = offset;
    offset += rowSizes[i];
  }
  ARstart[numRow] = offset;

  // finally add the nonzeros to the row major matrix
  ARindex.resize(numNz);
  ARcolor.resize(numNz);
  for (HighsInt i = 0; i < numCol; ++i) {
    for (HighsInt j = Astart[i]; j < Astart[i + 1]; ++j) {
      HighsInt row = Aindex[j];
      HighsInt ARpos = ARstart[row + 1] - rowSizes[row];
      rowSizes[row] -= 1;
      ARcolor[ARpos] = Acolor[j];
      ARindex[ARpos] = i;
    }
  }

  colToCell.resize(numCol);

  // loop over the columns and assign them a number that is distinct based on
  // their upper/lower bounds, cost, and integrality status. Use the columnSet
  // hash table to look up whether a column with similar properties exists and
  // use the previous number in that case. The number is stored in the
  // colToCell array which is subsequently used to sort an initial column
  // permutation.
  for (HighsInt i = 0; i < numCol; ++i) {
    MatrixColumn matrixCol;

    matrixCol.cost = coloring.color(lp.colCost_[i]);
    matrixCol.lb = coloring.color(lp.colLower_[i]);
    matrixCol.ub = coloring.color(lp.colUpper_[i]);
    matrixCol.integral = (u32)lp.integrality_[i];
    matrixCol.len = Astart[i + 1] - Astart[i];

    HighsInt* columnCell = &columnSet[matrixCol];

    if (*columnCell == 0) *columnCell = columnSet.size();

    colToCell[i] = *columnCell;
  }

  // set up the initial partition array, sort by the colToCell value
  // assigned above
  colPartition.resize(numCol);
  std::iota(colPartition.begin(), colPartition.end(), 0);
  std::sort(colPartition.begin(), colPartition.end(),
            [&](HighsInt col1, HighsInt col2) {
              return colToCell[col1] < colToCell[col2];
            });

  // now set up partition links and correct the colToCell array to the
  // correct cell index
  colPartitionLinks.resize(numCol);
  HighsInt cellStart = 0;
  HighsInt cellNumber = 0;
  for (HighsInt i = 0; i < numCol; ++i) {
    HighsInt col = colPartition[i];
    // if the cell number is different to the current cell number this is the
    // start of a new cell
    if (cellNumber != colToCell[col]) {
      // remember the number of this cell to indetify its end
      cellNumber = colToCell[col];
      // set the link of the cell start to point to its end
      colPartitionLinks[cellStart] = i;
      // remember start of this cell
      cellStart = i;
    }

    // correct the colToCell array to not store the start index of the
    // cell, not its number
    colToCell[col] = cellStart;
    // set the link of the column to the cellStart
    colPartitionLinks[i] = cellStart;
  }

  // set the column partition link of the last started cell to point past the
  // end
  colPartitionLinks[cellStart] = numCol;

  // loop over the rows and assign them a number that is distinct based on
  // their upper and lower bounds. Use the rowSet hash table to look up
  // whether a row with similar bounds exists and use the previous number in
  // that case. The number is stored in the colToCell
  // array which is subsequently used to sort an initial row permutation.
  rowToCell.resize(numRow);

  for (HighsInt i = 0; i < numRow; ++i) {
    MatrixRow matrixRow;

    matrixRow.lb = coloring.color(lp.rowLower_[i]);
    matrixRow.ub = coloring.color(lp.rowUpper_[i]);
    matrixRow.len = ARstart[i + 1] - ARstart[i];

    HighsInt* rowCell = &rowSet[matrixRow];

    if (*rowCell == 0) *rowCell = columnSet.size();

    rowToCell[i] = *rowCell;
  }

  // set up the initial partition array, sort by the rowToCell value
  // assigned above
  rowPartition.resize(numRow);
  std::iota(rowPartition.begin(), rowPartition.end(), 0);
  std::sort(rowPartition.begin(), rowPartition.end(),
            [&](HighsInt row1, HighsInt row2) {
              return rowToCell[row1] < rowToCell[row2];
            });

  // now set up partition links and correct the rowToCell array to the correct
  // cell index
  rowPartitionLinks.resize(numRow);
  cellStart = 0;
  cellNumber = 0;
  for (HighsInt i = 0; i < numRow; ++i) {
    HighsInt row = rowPartition[i];
    // if the cell number is different to the current cell number this is the
    // start of a new cell
    if (cellNumber != rowToCell[row]) {
      // remember the number of this cell to indetify its end
      cellNumber = rowToCell[row];
      // set the link of the cell start to point to its end
      rowPartitionLinks[cellStart] = i;
      // remember start of this cell
      cellStart = i;
    }

    // correct the rowToCell array to not store the start index of the cell,
    // not its number
    rowToCell[row] = cellStart;
    // set the link of the column to the cellStart
    rowPartitionLinks[i] = cellStart;
  }

  // set the column partition link of the last started cell to point past the
  // end
  rowPartitionLinks[cellStart] = numRow;
#if 0
  numCol = lp.numCol_;
  numRow = lp.numRow_;
  HighsHashTable<MatrixColumn, HighsInt> columnSet;
  HighsHashTable<MatrixRow, HighsInt> rowSet;
  HighsMatrixColoring coloring(epsilon);

  // set up row and column based incidence matrix
  HighsInt numNz = lp.Aindex_.size();
  Gindex.resize(2 * numNz);
  std::transform(lp.Aindex_.begin(), lp.Aindex_.end(), Gindex.begin(),
                 [&](HighsInt rowIndex) { return rowIndex + numCol; });

  Gstart.resize(numRow + numCol + 1);
  std::copy(lp.Astart_.begin(), lp.Astart_.end(), Gstart.begin());

  Gcolor.resize(2*numNz);

  // set up the column colors and count row sizes
  std::vector<HighsInt> rowSizes(numRow);
  for (HighsInt i = 0; i < numCol; ++i) {
    for (HighsInt j = Gstart[i]; j < Gstart[i + 1]; ++j) {
      Gcolor[j] = coloring.color(lp.Avalue_[j]);
      rowSizes[lp.Aindex_[j]] += 1;
    }
  }

  // next set up the row starts using the computed row sizes
  HighsInt offset = numNz;
  for (HighsInt i = 0; i < numRow; ++i) {
    Gstart[numCol+i] = offset;
    offset += rowSizes[i];
  }
  Gstart[numCol+numRow] = offset;

  // finally add the nonzeros to the row major matrix
  for (HighsInt i = 0; i < numCol; ++i) {
    for (HighsInt j = Gstart[i]; j < Gstart[i + 1]; ++j) {
      HighsInt row = lp.Aindex_[j];
      HighsInt ARpos = Gstart[numCol + row + 1] - rowSizes[row];
      rowSizes[row] -= 1;
      Gcolor[ARpos] = Gcolor[j];
      Gindex[ARpos] = i;
    }
  }

  vertexToCell.resize(numCol + numRow);

  // loop over the columns and assign them a number that is distinct based on
  // their upper/lower bounds, cost, and integrality status. Use the columnSet
  // hash table to look up whether a column with similar properties exists and
  // use the previous number in that case. The number is stored in the
  // colToCell array which is subsequently used to sort an initial column
  // permutation.
  for (HighsInt i = 0; i < numCol; ++i) {
    MatrixColumn matrixCol;

    matrixCol.cost = coloring.color(lp.colCost_[i]);
    matrixCol.lb = coloring.color(lp.colLower_[i]);
    matrixCol.ub = coloring.color(lp.colUpper_[i]);
    matrixCol.integral = (u32)lp.integrality_[i];
    matrixCol.len = Astart[i + 1] - Astart[i];

    HighsInt* columnCell = &columnSet[matrixCol];

    if (*columnCell == 0) *columnCell = columnSet.size();

    vertexToCell[i] = *columnCell;
  }

  for (HighsInt i = 0; i < numRow; ++i) {
    MatrixRow matrixRow;

    matrixRow.lb = coloring.color(lp.rowLower_[i]);
    matrixRow.ub = coloring.color(lp.rowUpper_[i]);
    matrixRow.len = Gstart[numCol+ i + 1] - Gstart[numCol + i];

    HighsInt* rowCell = &rowSet[matrixRow];

    if (*rowCell == 0) *rowCell = columnSet.size();

    vertexToCell[numCol + i] = *rowCell;
  }

  // set up the initial partition array, sort by the colToCell value
  // assigned above
  currentPartition.resize(numCol + numRow);
  std::iota(currentPartition.begin(), currentPartition.end(), 0);
  std::sort(currentPartition.begin(), currentPartition.end(),
            [&](HighsInt v1, HighsInt v2) {
              return vertexToCell[v1] < vertexToCell[v2];
            });

  // now set up partition links and correct the colToCell array to the
  // correct cell index
  HighsInt cellStart = 0;
  HighsInt cellNumber = 0;
  for (HighsInt i = 0; i < numCol; ++i) {
    HighsInt vertex = currentPartition[i];
    // if the cell number is different to the current cell number this is the
    // start of a new cell
    if (cellNumber != vertexToCell[vertex]) {
      // remember the number of this cell to indetify its end
      cellNumber = vertexToCell[vertex];
      // set the link of the cell start to point to its end
      currentPartitionLinks[cellStart] = i;
      // remember start of this cell
      cellStart = i;
    }

    // correct the colToCell array to not store the start index of the
    // cell, not its number
    vertexToCell[vertex] = cellStart;
    // set the link of the column to the cellStart
    currentPartitionLinks[i] = cellStart;
  }

  // set the column partition link of the last started cell to point past the
  // end
  currentPartitionLinks[cellStart] = numCol;
#endif
}

std::vector<std::tuple<HighsInt, HighsInt, HighsUInt>>
HighsSymmetryDetection::dumpCurrentGraph() {
  std::vector<std::tuple<HighsInt, HighsInt, HighsUInt>> graphTriplets;

  HighsInt numColTotal = Astart.size() - 1;
  graphTriplets.reserve(Aindex.size());

  for (HighsInt i = 0; i < numColTotal; ++i) {
    HighsInt colCell = colToCell[i];
    for (HighsInt j = Astart[i]; j != Astart[i + 1]; ++j)
      graphTriplets.emplace_back(rowToCell[Aindex[j]], colCell, Acolor[j]);
  }

  std::sort(graphTriplets.begin(), graphTriplets.end());
  return graphTriplets;
}

bool HighsSymmetryDetection::isomorphicToFirstLeave() {
  // verify result on hashes, even though a collision on all row and column
  // hashes is incredibly unlikely

  HighsInt certificateLen = numRow + numCol;

  if (std::memcmp(firstLeaveCertificate.data(), currNodeCertificate.data(),
                  certificateLen * sizeof(HighsInt)) != 0) {
#ifndef NDEBUG
    assert(std::memcmp(smallestLeaveCertificate.data(),
                       currNodeCertificate.data(),
                       certificateLen * sizeof(HighsInt)) != 0);
    HighsInt i = 0;
    for (HighsInt i = 0; i < certificateLen; ++i)
      assert(currNodeCertificate[i] <= smallestLeaveCertificate[i]);
#endif
    smallestLeaveCertificate = currNodeCertificate;

    return false;
  }

  // verify result on hashes, even though a collision on all row and column
  // hashes is incredibly unlikely
  auto graph = dumpCurrentGraph();
  return std::equal(firstLeaveGraph.begin(), firstLeaveGraph.end(),
                    graph.begin());
}

void HighsSymmetryDetection::run(HighsSymmetries& symmetries) {
  initializeHashValues();
  partitionRefinement();
  removeFixPoints();
  initializeGroundSet();
  if (numCol == 0) return;
  printf("numCol: %d\n", numCol);
  createNode();
  bool backtrackToFirstPath = false;
  while (!nodeStack.empty()) {
    HighsInt stackEnd = cellCreationStack.size();
    if (stackEnd > nodeStack.back().stackStart) {
      // we need to backtrack the datastructures
      do {
        Node& currNode = nodeStack.back();
        backtrack(currNode.stackStart, stackEnd);
        stackEnd = currNode.stackStart;
        if ((backtrackToFirstPath && !currNode.onFirstPath) ||
            !determineNextToDistinguish()) {
          nodeStack.pop_back();
          continue;
        }
        break;
      } while (!nodeStack.empty());

      if (nodeStack.empty()) break;

      Node& currNode = nodeStack.back();
      // call cleanup backtrack with the final stackEnd
      // so that all hashes are up to date and the link arrays do not contain
      // chains anymore
      cleanupBacktrack(stackEnd);
      HighsInt targetCell = currNode.targetCell;
      CellType targetCellType = currNode.targetCellType;

      distinguishVertex(targetCell, targetCellType);
      partitionRefinement();
      createNode();
    }
    backtrackToFirstPath = false;
    std::pair<HighsInt, CellType> targetCell = selectTargetCell();
    if (targetCell.first == -1) {
      if (firstLeaveColPartition.empty()) {
        firstLeaveColPartition = colPartition;
        firstLeaveCertificate = currNodeCertificate;
        smallestLeaveCertificate = currNodeCertificate;
        firstLeaveGraph = dumpCurrentGraph();
        backtrackToFirstPath = true;
        printf("dicovered first leave\n");
      } else if (isomorphicToFirstLeave()) {
        std::vector<HighsInt> permutation(numCol);

        for (HighsInt i = 0; i < numCol; ++i) {
          HighsInt firstLeaveCol = firstLeaveColPartition[i];
          permutation[colPosition[colPartition[i]]] = firstLeaveCol;
        }

        bool report = false;
        for (HighsInt i = 0; i < numCol; ++i) {
          if (mergeOrbits(permutation[i], groundSet[i])) report = true;
        }

        if (report) {
          symmetries.permutations.insert(symmetries.permutations.end(),
                                         permutation.begin(),
                                         permutation.end());
          ++symmetries.numPerms;
          printf("discovered useful generator nr. %d\n", symmetries.numPerms);
        } else
          printf("discovered redundant generator\n");

        backtrackToFirstPath = true;
      } else
        printf("dicovered bad leave\n");
      nodeStack.pop_back();
      continue;
    } else if (targetCell.first == -2) {
      // printf("pruned node due to bad invariant\n");
      nodeStack.pop_back();
      continue;
    }
    Node& currNode = nodeStack.back();
    currNode.targetCell = targetCell.first;
    currNode.targetCellType = targetCell.second;
    if (currNode.onFirstPath)
      currNode.onFirstPath = currNode.targetCellType == kColCell;
    determineNextToDistinguish();
    distinguishVertex(targetCell.first, targetCell.second);
    partitionRefinement();
    createNode();
  }

  symmetries.permutationColumns = std::move(groundSet);
}