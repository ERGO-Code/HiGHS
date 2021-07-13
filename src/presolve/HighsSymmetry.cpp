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

#include "mip/HighsDomain.h"

void HighsSymmetryDetection::removeFixPoints() {
  Gend.resize(numVertices);
  for (HighsInt i = 0; i < numVertices; ++i) {
    Gend[i] =
        std::partition(Gedge.begin() + Gstart[i], Gedge.begin() + Gstart[i + 1],
                       [&](const std::pair<HighsInt, HighsUInt>& edge) {
                         return cellSize(vertexToCell[edge.first]) > 1;
                       }) -
        Gedge.begin();
    assert(Gend[i] >= Gstart[i] && Gend[i] <= Gstart[i + 1]);
  }

  HighsInt unitCellIndex = numVertices;
  currentPartition.erase(
      std::remove_if(currentPartition.begin(), currentPartition.end(),
                     [&](HighsInt vertex) {
                       if (cellSize(vertexToCell[vertex]) == 1) {
                         --unitCellIndex;
                         HighsInt oldCellStart = vertexToCell[vertex];
                         vertexToCell[vertex] = unitCellIndex;
                         return true;
                       }
                       return false;
                     }),
      currentPartition.end());

  hashValid.assign(numVertices, false);
  for (HighsInt i = 0; i < numVertices; ++i) {
    if (Gend[i] == Gstart[i + 1]) continue;

    for (HighsInt j = Gend[i]; j < Gstart[i + 1]; ++j)
      Gedge[j].first = vertexToCell[Gedge[j].first];

    std::sort(Gedge.data() + Gend[i], Gedge.data() + Gstart[i + 1]);

    u64 vecHash = HighsHashHelpers::vector_hash(Gedge.data() + Gend[i],
                                                Gstart[i + 1] - Gend[i]);
    vertexHashes[i] = HighsHashHelpers::hash(vecHash);
  }

  if ((HighsInt)currentPartition.size() < numVertices) {
    numVertices = currentPartition.size();
    currentPartitionLinks.resize(numVertices);
    cellInRefinementQueue.assign(numVertices, false);
    assert(refinementQueue.empty());
    refinementQueue.clear();
    HighsInt cellStart = 0;
    HighsInt cellNumber = 0;
    for (HighsInt i = 0; i < numVertices; ++i) {
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

      if (!cellInRefinementQueue[cellStart]) {
        cellInRefinementQueue[cellStart] = true;
        refinementQueue.push_back(cellStart);
        std::push_heap(refinementQueue.begin(), refinementQueue.end(),
                       std::greater<HighsInt>());
      }

      // correct the vertexToCell array to not store the start index of the
      // cell, not its number
      updateCellMembership(i, cellStart, false);
    }

    // set the column partition link of the last started cell to point past the
    // end
    currentPartitionLinks[cellStart] = numVertices;

    numActiveCols =
        std::partition_point(currentPartition.begin(), currentPartition.end(),
                             [&](HighsInt v) { return v < numCol; }) -
        currentPartition.begin();

    partitionRefinement();
  } else
    numActiveCols = numCol;
}

void HighsSymmetries::mergeOrbits(HighsInt v1, HighsInt v2) {
  if (v1 == v2) return;

  HighsInt orbit1 = getOrbit(v1);
  HighsInt orbit2 = getOrbit(v2);

  if (orbit1 == orbit2) return;

  if (orbitSize[orbit2] < orbitSize[orbit1]) {
    orbitPartition[orbit2] = orbit1;
    orbitSize[orbit1] += orbitSize[orbit2];
  } else {
    orbitPartition[orbit1] = orbit2;
    orbitSize[orbit2] += orbitSize[orbit1];
  }

  return;
}

HighsInt HighsSymmetries::getOrbit(HighsInt col) {
  HighsInt i = columnPosition[col];
  if (i == -1) return -1;
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

std::shared_ptr<const StabilizerOrbits>
HighsSymmetries::computeStabilizerOrbits(const HighsDomain& localdom) {
  const auto& domchgStack = localdom.getDomainChangeStack();
  const auto& branchingPos = localdom.getBranchingPositions();
  const auto& prevBounds = localdom.getPreviousBounds();

  StabilizerOrbits stabilizerOrbits;
  stabilizerOrbits.stabilizedCols.reserve(permutationColumns.size());
  for (HighsInt i : branchingPos) {
    HighsInt col = domchgStack[i].column;
    if (columnPosition[col] == -1) continue;

    assert(localdom.variableType(col) != HighsVarType::kContinuous);

    // if we branch a variable upwards it is either binary and branched to one
    // and needs to be stabilized or it is a general integer and needs to be
    // stabilized regardless of branching direction.
    // if we branch downwards we only need to stabilize on this
    // branching column if it is a general integer
    if (domchgStack[i].boundtype == HighsBoundType::kLower ||
        !localdom.isGlobalBinary(col))
      stabilizerOrbits.stabilizedCols.push_back(columnPosition[col]);
  }

  HighsInt permLength = permutationColumns.size();
  orbitPartition.resize(permLength);
  std::iota(orbitPartition.begin(), orbitPartition.end(), 0);
  orbitSize.assign(permLength, 1);

  for (HighsInt i = 0; i < numPerms; ++i) {
    const HighsInt* perm = permutations.data() + i * permutationColumns.size();

    bool permRespectsBranchings = true;
    for (HighsInt i : stabilizerOrbits.stabilizedCols) {
      if (permutationColumns[i] != perm[i]) {
        permRespectsBranchings = false;
        break;
      }
    }

    if (!permRespectsBranchings) continue;

    for (HighsInt j = 0; j < permLength; ++j) {
      mergeOrbits(permutationColumns[j], perm[j]);
    }
  }

  stabilizerOrbits.stabilizedCols.clear();

  stabilizerOrbits.orbitCols.reserve(permLength);
  for (HighsInt i = 0; i < permLength; ++i) {
    if (localdom.variableType(permutationColumns[i]) ==
        HighsVarType::kContinuous)
      continue;
    HighsInt orbit = getOrbit(permutationColumns[i]);
    if (orbitSize[orbit] == 1)
      stabilizerOrbits.stabilizedCols.push_back(permutationColumns[i]);
    else
      stabilizerOrbits.orbitCols.push_back(permutationColumns[i]);
  }
  stabilizerOrbits.columnPosition = columnPosition.data();
  std::sort(stabilizerOrbits.stabilizedCols.begin(),
            stabilizerOrbits.stabilizedCols.end());
  if (!stabilizerOrbits.orbitCols.empty()) {
    std::sort(stabilizerOrbits.orbitCols.begin(),
              stabilizerOrbits.orbitCols.end(),
              [&](HighsInt col1, HighsInt col2) {
                return getOrbit(col1) < getOrbit(col2);
              });
    HighsInt numOrbitCols = stabilizerOrbits.orbitCols.size();
    stabilizerOrbits.orbitStarts.reserve(numOrbitCols + 1);
    stabilizerOrbits.orbitStarts.push_back(0);

    for (HighsInt i = 1; i < numOrbitCols; ++i) {
      if (getOrbit(stabilizerOrbits.orbitCols[i]) !=
          getOrbit(stabilizerOrbits.orbitCols[i - 1]))
        stabilizerOrbits.orbitStarts.push_back(i);
    }
    stabilizerOrbits.orbitStarts.push_back(numOrbitCols);
  }

  return std::make_shared<const StabilizerOrbits>(std::move(stabilizerOrbits));
}

HighsInt StabilizerOrbits::orbitalFixing(HighsDomain& domain) const {
  if (orbitCols.empty()) return 0;

  HighsInt numFixed = 0;
  HighsInt numOrbits = orbitStarts.size() - 1;
  for (HighsInt i = 0; i < numOrbits; ++i) {
    bool containsNonBin = false;
    HighsInt fixcol = -1;
    for (HighsInt j = orbitStarts[i]; j < orbitStarts[i + 1]; ++j) {
      if (domain.isFixed(orbitCols[j])) {
        fixcol = orbitCols[j];
        break;
      }
    }

    if (fixcol != -1) {
      HighsInt oldNumFixed = numFixed;
      double fixVal = domain.colLower_[fixcol];
      auto oldSize = domain.getDomainChangeStack().size();
      if (domain.colLower_[fixcol] == 1.0) {
        for (HighsInt j = orbitStarts[i]; j < orbitStarts[i + 1]; ++j) {
          if (domain.colLower_[orbitCols[j]] == 1.0) continue;
          ++numFixed;
          domain.changeBound(HighsBoundType::kLower, orbitCols[j], 1.0,
                             HighsDomain::Reason::unspecified());
          if (domain.infeasible()) return numFixed;
        }
      } else {
        for (HighsInt j = orbitStarts[i]; j < orbitStarts[i + 1]; ++j) {
          if (domain.colUpper_[orbitCols[j]] == 0.0) continue;
          ++numFixed;
          domain.changeBound(HighsBoundType::kUpper, orbitCols[j], 0.0,
                             HighsDomain::Reason::unspecified());
          if (domain.infeasible()) return numFixed;
        }
      }

      HighsInt newFixed = numFixed - oldNumFixed;

      if (newFixed != 0) {
        domain.propagate();
        if (domain.infeasible()) return numFixed;
        if (domain.getDomainChangeStack().size() - oldSize > newFixed) i = -1;
      }
    }
  }

  return numFixed;
}

void HighsSymmetryDetection::initializeGroundSet() {
  vertexGroundSet = currentPartition;
  std::sort(vertexGroundSet.begin(), vertexGroundSet.end());
  vertexPosition.resize(vertexToCell.size(), -1);
  for (HighsInt i = 0; i < numVertices; ++i)
    vertexPosition[vertexGroundSet[i]] = i;

  orbitPartition.resize(numVertices);
  std::iota(orbitPartition.begin(), orbitPartition.end(), 0);
  orbitSize.assign(numVertices, 1);

  automorphisms.resize(numVertices * 64);
  numAutomorphisms = 0;
  currNodeCertificate.reserve(numVertices);
}

bool HighsSymmetryDetection::mergeOrbits(HighsInt v1, HighsInt v2) {
  if (v1 == v2) return false;

  HighsInt orbit1 = getOrbit(v1);
  HighsInt orbit2 = getOrbit(v2);

  if (orbit1 == orbit2) return false;

  if (orbit1 < orbit2) {
    orbitPartition[orbit2] = orbit1;
    orbitSize[orbit1] += orbitSize[orbit2];
  } else {
    orbitPartition[orbit1] = orbit2;
    orbitSize[orbit2] += orbitSize[orbit1];
  }

  return true;
}

HighsInt HighsSymmetryDetection::getOrbit(HighsInt vertex) {
  HighsInt i = vertexPosition[vertex];
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
  vertexHashes.resize(numVertices);
  hashValid.resize(numVertices);

  for (HighsInt i = 0; i != numVertices; ++i) {
    HighsInt cell = vertexToCell[i];
    markCellForRefinement(cell);
  }
}

bool HighsSymmetryDetection::updateCellMembership(HighsInt i, HighsInt cell,
                                                  bool markForRefinement) {
  HighsInt vertex = currentPartition[i];
  if (vertexToCell[vertex] != cell) {
    // set new cell id
    HighsInt oldCellStart = vertexToCell[vertex];
    vertexToCell[vertex] = cell;
    if (i != cell) currentPartitionLinks[i] = cell;

    // update hashes of affected rows
    for (HighsInt j = Gstart[vertex]; j != Gend[vertex]; ++j) {
      hashValid[Gedge[j].first] = false;
      if (markForRefinement)
        markCellForRefinement(vertexToCell[Gedge[j].first]);
    }

    return true;
  }

  return false;
}

bool HighsSymmetryDetection::splitCell(HighsInt cell, HighsInt splitPoint) {
  u64 hSplit = getVertexHash(currentPartition[splitPoint]);
  u64 hCell = getVertexHash(currentPartition[cell]);

  u32 certificateVal =
      (HighsHashHelpers::pair_hash<0>(hSplit, hSplit >> 32) +
       HighsHashHelpers::pair_hash<1>(hCell, hCell >> 32) +
       HighsHashHelpers::pair_hash<2>(
           cell, currentPartitionLinks[cell] - splitPoint) +
       HighsHashHelpers::pair_hash<3>(splitPoint, splitPoint - cell)) >>
      32;

  // employ prefix pruning scheme as in bliss
  if (!firstLeaveCertificate.empty()) {
    firstLeavePrefixLen +=
        (firstLeavePrefixLen == currNodeCertificate.size()) *
        (certificateVal == firstLeaveCertificate[currNodeCertificate.size()]);
    bestLeavePrefixLen +=
        (bestLeavePrefixLen == currNodeCertificate.size()) *
        (certificateVal == bestLeaveCertificate[currNodeCertificate.size()]);

    // if the node certificate is not a prefix of the first leave's certificate
    // and it comes lexicographically after the certificate value of the
    // lexicographically smallest leave certificate we prune the node
    if (firstLeavePrefixLen <= currNodeCertificate.size() &&
        bestLeavePrefixLen == currNodeCertificate.size() &&
        certificateVal > bestLeaveCertificate[currNodeCertificate.size()])
      return false;
  }

  currentPartitionLinks[splitPoint] = currentPartitionLinks[cell];
  currentPartitionLinks[cell] = splitPoint;
  cellCreationStack.push_back(splitPoint);
  currNodeCertificate.push_back(certificateVal);

  return true;
}

void HighsSymmetryDetection::markCellForRefinement(HighsInt cell) {
  if (cellSize(cell) == 1 || cellInRefinementQueue[cell]) return;

  cellInRefinementQueue[cell] = true;
  refinementQueue.push_back(cell);
  std::push_heap(refinementQueue.begin(), refinementQueue.end(),
                 std::greater<HighsInt>());
}

HighsSymmetryDetection::u64 HighsSymmetryDetection::getVertexHash(HighsInt v) {
  if (!hashValid[v]) {
    HighsInt dynamicLen = Gend[v] - Gstart[v];
    if (dynamicLen != 0) {
      std::transform(
          Gedge.begin() + Gstart[v], Gedge.begin() + Gend[v],
          edgeBuffer.begin(), [&](const std::pair<HighsInt, HighsUInt>& edge) {
            return std::make_pair(vertexToCell[edge.first], edge.second);
          });

      std::sort(edgeBuffer.begin(), edgeBuffer.begin() + dynamicLen);

      u64 vecHash =
          HighsHashHelpers::vector_hash(edgeBuffer.data(), dynamicLen);

      // the vertex hash already contains the hash of the stable part in the
      // lower bits so make sure to keep it and only overwrite the upper part
      vertexHashes[v] =
          HighsHashHelpers::hash(vecHash) << 32 | u32(vertexHashes[v]);
    }
    hashValid[v] = true;
  }

  return vertexHashes[v];
}

bool HighsSymmetryDetection::partitionRefinement() {
  while (!refinementQueue.empty()) {
    std::pop_heap(refinementQueue.begin(), refinementQueue.end(),
                  std::greater<HighsInt>());

    HighsInt cellStart = refinementQueue.back();
    HighsInt firstCellStart = cellStart;
    refinementQueue.pop_back();
    cellInRefinementQueue[cellStart] = false;

    if (cellSize(cellStart) == 1) continue;
    HighsInt cellEnd = currentPartitionLinks[cellStart];
    assert(cellEnd >= cellStart);

    for (HighsInt i = cellStart; i < cellEnd; ++i) {
      HighsInt v = currentPartition[i];
      getVertexHash(v);
    }

    std::sort(currentPartition.begin() + cellStart,
              currentPartition.begin() + cellEnd,
              [&](HighsInt v1, HighsInt v2) {
                assert(hashValid[v1]);
                assert(hashValid[v2]);
                return vertexHashes[v1] < vertexHashes[v2];
              });

    bool prune = false;
    HighsInt i;
    for (i = cellStart + 1; i < cellEnd; ++i) {
      HighsInt vertex = currentPartition[i];
      if (vertexHashes[vertex] != vertexHashes[currentPartition[i - 1]]) {
        if (!splitCell(cellStart, i)) {
          prune = true;
          break;
        }
        cellStart = i;
      }

      updateCellMembership(i, cellStart);
    }

    if (prune) {
      for (HighsInt c : refinementQueue) cellInRefinementQueue[c] = false;
      refinementQueue.clear();

      for (--i; i > firstCellStart; --i) {
        HighsInt vertex = currentPartition[i];
        if (!updateCellMembership(i, firstCellStart, false)) break;
      }
      return false;
    }

    // set the link of the last started cell to point to the cell end
    assert(currentPartitionLinks[cellStart] == cellEnd);
  }

  return true;
}

HighsInt HighsSymmetryDetection::selectTargetCell() {
  HighsInt i = 0;
  if (nodeStack.size() > 1) i = nodeStack[nodeStack.size() - 2].targetCell + 1;

  while (i < numVertices) {
    if (cellSize(i) > 1) return i;

    ++i;
  }

  return -1;
}

bool HighsSymmetryDetection::checkStoredAutomorphism(HighsInt vertex) {
  HighsInt numCheck = std::min(numAutomorphisms, 64);

  for (HighsInt i = 0; i < numCheck; ++i) {
    const HighsInt* automorphism = automorphisms.data() + i * numVertices;
    bool automorphismUseful = true;
    for (HighsInt j = nodeStack.size() - 2; j >= firstPathDepth; --j) {
      HighsInt fixPos = vertexPosition[nodeStack[j].lastDistiguished];

      if (automorphism[fixPos] != vertexGroundSet[fixPos]) {
        automorphismUseful = false;
        break;
      }
    }

    if (!automorphismUseful) continue;

    if (automorphism[vertexPosition[vertex]] < vertex) return false;
  }

  return true;
}

bool HighsSymmetryDetection::determineNextToDistinguish() {
  Node& currNode = nodeStack.back();
  distinguishCands.clear();
  std::vector<HighsInt>::iterator cellStart;
  std::vector<HighsInt>::iterator cellEnd;
  cellStart = currentPartition.begin() + currNode.targetCell;
  cellEnd =
      currentPartition.begin() + currentPartitionLinks[currNode.targetCell];

  if (currNode.lastDistiguished == -1) {
    auto nextDistinguishPos = std::min_element(cellStart, cellEnd);
    distinguishCands.push_back(&*nextDistinguishPos);
  } else if ((HighsInt)nodeStack.size() > firstPathDepth) {
    for (auto i = cellStart; i != cellEnd; ++i) {
      if (*i > currNode.lastDistiguished && checkStoredAutomorphism(*i))
        distinguishCands.push_back(&*i);
    }
    if (distinguishCands.empty()) return false;
    auto nextDistinguishPos =
        std::min_element(distinguishCands.begin(), distinguishCands.end(),
                         [](HighsInt* a, HighsInt* b) { return *a < *b; });
    std::swap(*distinguishCands.begin(), *nextDistinguishPos);
    distinguishCands.resize(1);
  } else {
    for (auto i = cellStart; i != cellEnd; ++i) {
      if (*i > currNode.lastDistiguished && vertexGroundSet[getOrbit(*i)] == *i)
        distinguishCands.push_back(&*i);
    }
    if (distinguishCands.empty()) return false;
    auto nextDistinguishPos =
        std::min_element(distinguishCands.begin(), distinguishCands.end(),
                         [](HighsInt* a, HighsInt* b) { return *a < *b; });
    std::swap(*distinguishCands.begin(), *nextDistinguishPos);
    distinguishCands.resize(1);
  }

  return true;
}

bool HighsSymmetryDetection::distinguishVertex(HighsInt targetCell) {
  assert(distinguishCands.size() == 1u);
  HighsInt targetCellEnd = currentPartitionLinks[targetCell];
  HighsInt newCell = targetCellEnd - 1;
  std::swap(*distinguishCands[0], currentPartition[newCell]);
  nodeStack.back().lastDistiguished = currentPartition[newCell];

  if (!splitCell(targetCell, newCell)) return false;

  updateCellMembership(newCell, newCell);

  return true;
}

void HighsSymmetryDetection::backtrack(HighsInt backtrackStackNewEnd,
                                       HighsInt backtrackStackEnd) {
  // we assume that we always backtrack from a leave node, i.e. a discrete
  // partition therefore we do not need to remember the values of the hash
  // contributions as it is the indentity for each position and all new cells
  // are on the cell creation stack.
  for (HighsInt stackPos = backtrackStackEnd - 1;
       stackPos >= backtrackStackNewEnd; --stackPos) {
    HighsInt cell = cellCreationStack[stackPos];
    // look up the cell start of the preceding cell with link compression
    HighsInt newStart = getCellStart(cell - 1);
    // remember the current end
    HighsInt currEnd = currentPartitionLinks[cell];
    // change the link to point to the start of the preceding cell
    currentPartitionLinks[cell] = newStart;
    // change the link of the start pointer of the preceding cell to point to
    // the end of this cell
    currentPartitionLinks[newStart] = currEnd;
  }
}

void HighsSymmetryDetection::cleanupBacktrack(HighsInt cellCreationStackPos) {
  // the links have been updated. Even though they might still not be fully
  // compressed the cell starts will all point to the correct cell end and the
  // lookup with path compression will give the correct start
  for (HighsInt stackPos = cellCreationStack.size() - 1;
       stackPos >= cellCreationStackPos; --stackPos) {
    HighsInt cell = cellCreationStack[stackPos];

    HighsInt cellStart = getCellStart(cell);
    HighsInt cellEnd = currentPartitionLinks[cellStart];

    for (HighsInt v = cell;
         v < cellEnd && vertexToCell[currentPartition[v]] == cell; ++v)
      updateCellMembership(v, cellStart, false);
  }

  cellCreationStack.resize(cellCreationStackPos);
}

HighsInt HighsSymmetryDetection::getCellStart(HighsInt pos) {
  HighsInt startPos = currentPartitionLinks[pos];
  if (startPos > pos) return pos;
  if (currentPartitionLinks[startPos] < startPos) {
    do {
      linkCompressionStack.push_back(pos);
      pos = startPos;
      startPos = currentPartitionLinks[startPos];
    } while (currentPartitionLinks[startPos] < startPos);

    do {
      currentPartitionLinks[linkCompressionStack.back()] = startPos;
      linkCompressionStack.pop_back();
    } while (!linkCompressionStack.empty());
  }

  return startPos;
}

void HighsSymmetryDetection::createNode() {
  nodeStack.emplace_back();
  nodeStack.back().stackStart = cellCreationStack.size();
  nodeStack.back().certificateEnd = currNodeCertificate.size();
  nodeStack.back().targetCell = -1;
  nodeStack.back().lastDistiguished = -1;
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

void HighsSymmetryDetection::loadModelAsGraph(const HighsLp& model,
                                              double epsilon) {
  this->model = &model;
  numCol = model.numCol_;
  numRow = model.numRow_;
  numVertices = numRow + numCol;

  cellInRefinementQueue.resize(numVertices);
  vertexToCell.resize(numVertices);
  refinementQueue.reserve(numVertices);
  currNodeCertificate.reserve(numVertices);

  HighsHashTable<MatrixColumn, HighsInt> columnSet;
  HighsHashTable<MatrixRow, HighsInt> rowSet;
  HighsMatrixColoring coloring(epsilon);
  edgeBuffer.resize(numVertices);
  // set up row and column based incidence matrix
  HighsInt numNz = model.Aindex_.size();
  Gedge.resize(2 * numNz);
  std::transform(model.Aindex_.begin(), model.Aindex_.end(), Gedge.begin(),
                 [&](HighsInt rowIndex) {
                   return std::make_pair(rowIndex + numCol, HighsUInt{0});
                 });

  Gstart.resize(numVertices + 1);
  std::copy(model.Astart_.begin(), model.Astart_.end(), Gstart.begin());

  // set up the column colors and count row sizes
  std::vector<HighsInt> rowSizes(numRow);
  for (HighsInt i = 0; i < numCol; ++i) {
    for (HighsInt j = Gstart[i]; j < Gstart[i + 1]; ++j) {
      Gedge[j].second = coloring.color(model.Avalue_[j]);
      rowSizes[model.Aindex_[j]] += 1;
    }
  }

  // next set up the row starts using the computed row sizes
  HighsInt offset = numNz;
  for (HighsInt i = 0; i < numRow; ++i) {
    Gstart[numCol + i] = offset;
    offset += rowSizes[i];
  }
  Gstart[numCol + numRow] = offset;

  Gend.assign(Gstart.begin() + 1, Gstart.end());

  // finally add the nonzeros to the row major matrix
  for (HighsInt i = 0; i < numCol; ++i) {
    for (HighsInt j = Gstart[i]; j < Gstart[i + 1]; ++j) {
      HighsInt row = model.Aindex_[j];
      HighsInt ARpos = Gstart[numCol + row + 1] - rowSizes[row];
      rowSizes[row] -= 1;
      Gedge[ARpos].first = i;
      Gedge[ARpos].second = Gedge[j].second;
    }
  }

  // loop over the columns and assign them a number that is distinct based on
  // their upper/lower bounds, cost, and integrality status. Use the columnSet
  // hash table to look up whether a column with similar properties exists and
  // use the previous number in that case. The number is stored in the
  // colToCell array which is subsequently used to sort an initial column
  // permutation.
  HighsInt indexOffset = numCol + 1;
  for (HighsInt i = 0; i < numCol; ++i) {
    MatrixColumn matrixCol;

    matrixCol.cost = coloring.color(model.colCost_[i]);
    matrixCol.lb = coloring.color(model.colLower_[i]);
    matrixCol.ub = coloring.color(model.colUpper_[i]);
    matrixCol.integral = (u32)model.integrality_[i];
    matrixCol.len = Gstart[i + 1] - Gstart[i];

    HighsInt* columnCell = &columnSet[matrixCol];

    if (*columnCell == 0) {
      *columnCell = columnSet.size();
      if (model.colLower_[i] != 0.0 || model.colUpper_[i] != 1.0 ||
          model.integrality_[i] == HighsVarType::kContinuous)
        *columnCell += indexOffset;
    }

    vertexToCell[i] = *columnCell;
  }

  indexOffset = 2 * numCol + 1;
  for (HighsInt i = 0; i < numRow; ++i) {
    MatrixRow matrixRow;

    matrixRow.lb = coloring.color(model.rowLower_[i]);
    matrixRow.ub = coloring.color(model.rowUpper_[i]);
    matrixRow.len = Gstart[numCol + i + 1] - Gstart[numCol + i];

    HighsInt* rowCell = &rowSet[matrixRow];

    if (*rowCell == 0) *rowCell = rowSet.size();

    vertexToCell[numCol + i] = indexOffset + *rowCell;
  }

  // set up the initial partition array, sort by the colToCell value
  // assigned above
  currentPartition.resize(numVertices);
  std::iota(currentPartition.begin(), currentPartition.end(), 0);
  std::sort(currentPartition.begin(), currentPartition.end(),
            [&](HighsInt v1, HighsInt v2) {
              return vertexToCell[v1] < vertexToCell[v2];
            });

  // now set up partition links and correct the colToCell array to the
  // correct cell index
  currentPartitionLinks.resize(numVertices);
  HighsInt cellStart = 0;
  HighsInt cellNumber = 0;
  for (HighsInt i = 0; i < numVertices; ++i) {
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
  currentPartitionLinks[cellStart] = numVertices;
}

HighsHashTable<std::tuple<HighsInt, HighsInt, HighsUInt>>
HighsSymmetryDetection::dumpCurrentGraph() {
  HighsHashTable<std::tuple<HighsInt, HighsInt, HighsUInt>> graphTriplets;

  for (HighsInt i = 0; i < numCol; ++i) {
    HighsInt colCell = vertexToCell[i];
    for (HighsInt j = Gstart[i]; j != Gend[i]; ++j)
      graphTriplets.insert(vertexToCell[Gedge[j].first], colCell,
                           Gedge[j].second);
    for (HighsInt j = Gend[i]; j != Gstart[i + 1]; ++j)
      graphTriplets.insert(Gedge[j].first, colCell, Gedge[j].second);
  }

  return graphTriplets;
}

void HighsSymmetryDetection::switchToNextNode(HighsInt backtrackDepth) {
  HighsInt stackEnd = cellCreationStack.size();
  // we need to backtrack the datastructures

  nodeStack.resize(backtrackDepth);
  if (backtrackDepth == 0) return;
  do {
    Node& currNode = nodeStack.back();
    backtrack(currNode.stackStart, stackEnd);
    stackEnd = currNode.stackStart;
    firstPathDepth = std::min((HighsInt)nodeStack.size(), firstPathDepth);
    bestPathDepth = std::min((HighsInt)nodeStack.size(), bestPathDepth);
    firstLeavePrefixLen =
        std::min(currNode.certificateEnd, firstLeavePrefixLen);
    bestLeavePrefixLen = std::min(currNode.certificateEnd, bestLeavePrefixLen);
    currNodeCertificate.resize(currNode.certificateEnd);
    if (!determineNextToDistinguish()) {
      nodeStack.pop_back();
      continue;
    }

    // call cleanup backtrack with the final stackEnd
    // so that all hashes are up to date and the link arrays do not contain
    // chains anymore
    cleanupBacktrack(stackEnd);
    HighsInt targetCell = currNode.targetCell;

    if (!distinguishVertex(targetCell)) {
      // if distinguishing the next vertex fails, it means that its certificate
      // value is lexicographically larger than that of the best leave
      nodeStack.pop_back();
      continue;
    }

    if (!partitionRefinement()) continue;

    createNode();
    break;
  } while (!nodeStack.empty());
}

bool HighsSymmetryDetection::compareCurrentGraph(
    const HighsHashTable<std::tuple<HighsInt, HighsInt, HighsUInt>>&
        otherGraph) {
  for (HighsInt i = 0; i < numCol; ++i) {
    HighsInt colCell = vertexToCell[i];

    for (HighsInt j = Gstart[i]; j != Gend[i]; ++j)
      if (!otherGraph.find(std::make_tuple(vertexToCell[Gedge[j].first],
                                           colCell, Gedge[j].second)))
        return false;
    for (HighsInt j = Gend[i]; j != Gstart[i + 1]; ++j)
      if (!otherGraph.find(
              std::make_tuple(Gedge[j].first, colCell, Gedge[j].second)))
        return false;
  }

  return true;
}

bool HighsSymmetryDetection::isFromBinaryColumn(HighsInt pos) const {
  if (pos >= numActiveCols) return false;

  HighsInt col = currentPartition[pos];

  if (model->colLower_[col] != 0.0 || model->colUpper_[col] != 1.0 ||
      model->integrality_[col] == HighsVarType::kContinuous)
    return false;

  return true;
}

void HighsSymmetryDetection::run(HighsSymmetries& symmetries) {
  initializeHashValues();
  partitionRefinement();
  removeFixPoints();
  initializeGroundSet();
  if (numActiveCols == 0) return;
  currNodeCertificate.clear();
  cellCreationStack.clear();
  createNode();
  while (!nodeStack.empty()) {
    HighsInt targetCell = selectTargetCell();
    if (targetCell == -1) {
      if (firstLeavePartition.empty()) {
        firstLeavePartition = currentPartition;
        firstLeaveCertificate = currNodeCertificate;
        bestLeaveCertificate = currNodeCertificate;
        firstLeaveGraph = dumpCurrentGraph();
        firstPathDepth = nodeStack.size();
        bestPathDepth = nodeStack.size();
        firstLeavePrefixLen = currNodeCertificate.size();
        bestLeavePrefixLen = currNodeCertificate.size();

        HighsInt backtrackDepth = firstPathDepth - 1;
        while (backtrackDepth > 0 &&
               !isFromBinaryColumn(nodeStack[backtrackDepth - 1].targetCell))
          --backtrackDepth;
        switchToNextNode(backtrackDepth);
      } else {
        HighsInt backtrackDepth = nodeStack.size() - 1;
        assert(currNodeCertificate.size() == firstLeaveCertificate.size());
        if (firstLeavePrefixLen == currNodeCertificate.size() ||
            bestLeavePrefixLen == currNodeCertificate.size()) {
          if (firstLeavePrefixLen == currNodeCertificate.size() &&
              compareCurrentGraph(firstLeaveGraph)) {
            HighsInt k = (numAutomorphisms++) & 63;
            HighsInt* permutation = automorphisms.data() + k * numVertices;
            for (HighsInt i = 0; i < numVertices; ++i) {
              HighsInt firstLeaveCol = firstLeavePartition[i];
              permutation[vertexPosition[currentPartition[i]]] = firstLeaveCol;
            }

            bool report = false;
            for (HighsInt i = 0; i < numVertices; ++i) {
              if (mergeOrbits(permutation[i], vertexGroundSet[i]) &&
                  i < numActiveCols) {
                assert(permutation[i] < numCol);
                report = true;
              }
            }

            if (report) {
              symmetries.permutations.insert(symmetries.permutations.end(),
                                             permutation,
                                             permutation + numActiveCols);
              ++symmetries.numPerms;
            }
            backtrackDepth = std::min(backtrackDepth, firstPathDepth);
          } else if (!bestLeavePartition.empty() &&
                     bestLeavePrefixLen == currNodeCertificate.size() &&
                     compareCurrentGraph(bestLeaveGraph)) {
            HighsInt k = (numAutomorphisms++) & 63;
            HighsInt* permutation = automorphisms.data() + k * numVertices;
            for (HighsInt i = 0; i < numVertices; ++i) {
              HighsInt bestLeaveCol = bestLeavePartition[i];
              permutation[vertexPosition[currentPartition[i]]] = bestLeaveCol;
            }

            bool report = false;
            for (HighsInt i = 0; i < numVertices; ++i) {
              if (mergeOrbits(permutation[i], vertexGroundSet[i]) &&
                  i < numActiveCols) {
                assert(permutation[i] < numCol);
                report = true;
              }
            }

            if (report) {
              symmetries.permutations.insert(symmetries.permutations.end(),
                                             permutation,
                                             permutation + numActiveCols);
              ++symmetries.numPerms;
            }

            backtrackDepth = std::min(backtrackDepth, bestPathDepth);
          } else if (bestLeavePrefixLen < currNodeCertificate.size() &&
                     currNodeCertificate[bestLeavePrefixLen] >
                         bestLeaveCertificate[bestLeavePrefixLen]) {
            // certificate value is lexicographically above the smallest one
            // seen so far, so we might be able to backtrack to a higher level
            HighsInt possibleBacktrackDepth = firstPathDepth - 1;
            while (nodeStack[possibleBacktrackDepth].certificateEnd <=
                   bestLeavePrefixLen)
              ++possibleBacktrackDepth;

            backtrackDepth = std::min(possibleBacktrackDepth, backtrackDepth);
          }
        } else {
          // leave must have a lexicographically smaller certificate value
          // than the current best leave, because its prefix length is smaller
          // than the best leaves and it would have been already pruned if
          // it's certificate value was larger unless it is equal to the first
          // leave nodes certificate value which is caught by the first case
          // of the if confition. Hence, having a lexicographically smaller
          // certificate value than the best leave is the only way to get
          // here.
          assert(bestLeaveCertificate[bestLeavePrefixLen] >
                     currNodeCertificate[bestLeavePrefixLen] &&
                 std::memcmp(bestLeaveCertificate.data(),
                             currNodeCertificate.data(),
                             bestLeavePrefixLen * sizeof(u32)) == 0);
          bestLeaveCertificate = currNodeCertificate;
          bestLeaveGraph = dumpCurrentGraph();
          bestLeavePartition = currentPartition;
          bestPathDepth = nodeStack.size();
          bestLeavePrefixLen = currNodeCertificate.size();
        }

        switchToNextNode(backtrackDepth);
      }
    } else {
      Node& currNode = nodeStack.back();
      currNode.targetCell = targetCell;
      bool success = determineNextToDistinguish();
      assert(success);
      if (!distinguishVertex(targetCell)) {
        switchToNextNode(nodeStack.size() - 1);
        continue;
      }
      if (!partitionRefinement()) {
        switchToNextNode(nodeStack.size());
        continue;
      }

      createNode();
    }
  }

  if (symmetries.numPerms > 0) {
    // check which columns have non-trivial orbits
    HighsInt numFixed = 0;
    for (HighsInt i = 0; i < numActiveCols; ++i) {
      if (orbitSize[getOrbit(vertexGroundSet[i])] == 1) {
        vertexPosition[vertexGroundSet[i]] = -1;
        vertexGroundSet[i] = -1;
        numFixed += 1;
      }
    }

    vertexPosition.resize(numCol);

    if (numFixed != 0) {
      // now compress symmetries and the groundset to only contain the unfixed
      // columns
      HighsInt* perms = symmetries.permutations.data();
      HighsInt* permEnd =
          symmetries.permutations.data() + symmetries.numPerms * numActiveCols;
      HighsInt* permOutput = symmetries.permutations.data();
      while (perms != permEnd) {
        for (HighsInt i = 0; i < numActiveCols; ++i) {
          if (vertexGroundSet[i] == -1) continue;

          *permOutput = perms[i];
          ++permOutput;
        }

        perms += numActiveCols;
      }

      HighsInt outPos = 0;
      for (HighsInt i = 0; i < numActiveCols; ++i) {
        if (vertexGroundSet[i] == -1) continue;

        vertexGroundSet[outPos] = vertexGroundSet[i];
        vertexPosition[vertexGroundSet[outPos]] = outPos;
        outPos += 1;
      }

      numActiveCols -= numFixed;
    }

    vertexGroundSet.resize(numActiveCols);
    symmetries.permutationColumns = std::move(vertexGroundSet);
    symmetries.columnPosition = std::move(vertexPosition);
    symmetries.permutations.resize(symmetries.numPerms * numActiveCols);
  }
}