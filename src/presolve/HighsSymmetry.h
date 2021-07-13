/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file HighsSymmetry.h
 * @brief Facilities for symmetry detection
 * @author Leona Gottwald
 */

#ifndef PRESOLVE_HIGHS_SYMMETRY_H_
#define PRESOLVE_HIGHS_SYMMETRY_H_

#include <algorithm>
#include <map>
#include <vector>

#include "lp_data/HighsLp.h"
#include "util/HighsHash.h"
#include "util/HighsInt.h"

/// class that is responsible for assiging distinct colors for each distinct
/// double value
class HighsMatrixColoring {
  using u32 = std::uint32_t;

  std::map<double, u32> colorMap;
  double tolerance;

 public:
  // initialize with exact 0.0 and 1.0, to not have differing results due tiny
  // numerical differences on those values
  HighsMatrixColoring(double tolerance)
      : colorMap({{0.0, 1}, {1.0, 2}, {-kHighsInf, 3}, {kHighsInf, 4}}),
        tolerance(tolerance) {}

  u32 color(double value) {
    // iterator points to smallest element in map which fulfills key >= value -
    // tolerance
    auto it = colorMap.lower_bound(value - tolerance);
    u32 color;
    // check if there is no such element, or if this element has a key value +
    // tolerance in which case we create a new color and store it with the key
    // value
    if (it == colorMap.end() || it->first > value + tolerance)
      it = colorMap.emplace_hint(it, value, colorMap.size() + 1);
    return it->second;
  }
};

class HighsDomain;

struct StabilizerOrbits {
  std::vector<HighsInt> orbitCols;
  std::vector<HighsInt> orbitStarts;
  std::vector<HighsInt> stabilizedCols;
  const HighsInt* columnPosition;

  HighsInt orbitalFixing(HighsDomain& domain) const;

  bool isStabilized(HighsInt col) const {
    return columnPosition[col] == -1 ||
           std::binary_search(stabilizedCols.begin(), stabilizedCols.end(),
                              col);
  }
};

struct HighsSymmetries {
  std::vector<HighsInt> permutationColumns;
  std::vector<HighsInt> permutations;
  std::vector<HighsInt> orbitPartition;
  std::vector<HighsInt> orbitSize;
  std::vector<HighsInt> columnPosition;
  std::vector<HighsInt> linkCompressionStack;
  HighsInt numPerms = 0;

  void mergeOrbits(HighsInt col1, HighsInt col2);
  HighsInt getOrbit(HighsInt col);
  std::shared_ptr<const StabilizerOrbits> computeStabilizerOrbits(
      const HighsDomain& localdom);
};

class HighsSymmetryDetection {
  using u64 = std::uint64_t;
  using u32 = std::uint32_t;

  const HighsLp* model;
  // compressed graph storage
  std::vector<HighsInt> Gstart;
  std::vector<HighsInt> Gend;
  std::vector<std::pair<HighsInt, HighsUInt>> Gedge;

  std::vector<std::pair<HighsInt, HighsUInt>> edgeBuffer;

  std::vector<HighsInt> currentPartition;
  std::vector<HighsInt> currentPartitionLinks;
  std::vector<HighsInt> vertexToCell;
  std::vector<HighsInt> vertexPosition;
  std::vector<HighsInt> vertexGroundSet;
  std::vector<u64> vertexHashes;
  std::vector<bool> hashValid;
  std::vector<HighsInt> orbitPartition;
  std::vector<HighsInt> orbitSize;

  std::vector<HighsInt> cellCreationStack;
  std::vector<std::uint8_t> cellInRefinementQueue;
  std::vector<HighsInt> refinementQueue;
  std::vector<HighsInt*> distinguishCands;
  std::vector<HighsInt> automorphisms;

  std::vector<HighsInt> linkCompressionStack;

  std::vector<u32> currNodeCertificate;
  std::vector<u32> firstLeaveCertificate;
  std::vector<u32> bestLeaveCertificate;
  std::vector<HighsInt> firstLeavePartition;
  std::vector<HighsInt> bestLeavePartition;

  HighsHashTable<std::tuple<HighsInt, HighsInt, HighsUInt>> firstLeaveGraph;
  HighsHashTable<std::tuple<HighsInt, HighsInt, HighsUInt>> bestLeaveGraph;

  HighsInt firstLeavePrefixLen;
  HighsInt bestLeavePrefixLen;
  HighsInt firstPathDepth;
  HighsInt bestPathDepth;

  HighsInt numAutomorphisms;
  HighsInt numCol;
  HighsInt numRow;
  HighsInt numVertices;
  HighsInt numActiveCols;

  // node in the search tree for finding automorphisms
  struct Node {
    HighsInt stackStart;
    HighsInt certificateEnd;
    HighsInt targetCell;
    HighsInt lastDistiguished;
  };

  std::vector<Node> nodeStack;

  HighsInt getCellStart(HighsInt pos);

  void backtrack(HighsInt backtrackStackNewEnd, HighsInt backtrackStackEnd);
  void cleanupBacktrack(HighsInt cellCreationStackPos);

  void switchToNextNode(HighsInt backtrackDepth);

  bool compareCurrentGraph(
      const HighsHashTable<std::tuple<HighsInt, HighsInt, HighsUInt>>&
          otherGraph);

  void removeFixPoints();
  void initializeGroundSet();
  HighsHashTable<std::tuple<HighsInt, HighsInt, HighsUInt>> dumpCurrentGraph();
  bool mergeOrbits(HighsInt v1, HighsInt v2);
  HighsInt getOrbit(HighsInt vertex);

  void initializeHashValues();
  bool isomorphicToFirstLeave();
  bool partitionRefinement();
  bool checkStoredAutomorphism(HighsInt vertex);
  u64 getVertexHash(HighsInt vertex);
  HighsInt selectTargetCell();

  bool updateCellMembership(HighsInt vertex, HighsInt cell,
                            bool markForRefinement = true);
  bool splitCell(HighsInt cell, HighsInt splitPoint);
  void markCellForRefinement(HighsInt cell);

  bool distinguishVertex(HighsInt targetCell);
  bool determineNextToDistinguish();
  void createNode();

  HighsInt cellSize(HighsInt cell) const {
    return currentPartitionLinks[cell] - cell;
  }

  bool isFromBinaryColumn(HighsInt vertex) const;

 public:
  void loadModelAsGraph(const HighsLp& model, double epsilon);

  void run(HighsSymmetries& symmetries);
};

#endif