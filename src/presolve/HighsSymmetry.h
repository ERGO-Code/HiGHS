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
      : colorMap({{0.0, 0}, {1.0, 1}, {-kHighsInf, 2}, {kHighsInf, 3}}),
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
      it = colorMap.emplace_hint(it, value, colorMap.size());
    return it->second;
  }
};

struct HighsSymmetries
{
  std::vector<HighsInt> permutationColumns;
  std::vector<HighsInt> permutations;
  HighsInt numPerms = 0;
};

class HighsSymmetryDetection {
  using u64 = std::uint64_t;
  using u32 = std::uint32_t;

  std::vector<HighsInt> currentPartition;
  std::vector<HighsInt> currentPartitionLinks;
  std::vector<HighsInt> vertexToCell;
  std::vector<HighsInt> vertexPosition;
  std::vector<HighsInt> vertexGroundSet;
  std::vector<HighsInt> vertexHashes;

  // compressed graph storage
  std::vector<HighsInt> Gstart;
  std::vector<HighsInt> Gindex;
  std::vector<u32> Gcolor;


  // vector that stores the column indices in the order of the partition
  std::vector<HighsInt> colPartition;
  // vector that stores links for each column, where the link i corresponds to
  // column at position i in the partition the link points to the start of its
  // cell in the partition, except for the representant at the start of the
  // partition which points to the end of the cell to quickly compute cell sizes
  std::vector<HighsInt> colPartitionLinks;
  // stores at position i the cell index of column i
  std::vector<HighsInt> colToCell;

  // row partition data structure as above for the columns
  std::vector<HighsInt> rowPartition;
  std::vector<HighsInt> rowPartitionLinks;
  std::vector<HighsInt> rowToCell;

  // row and column major storage for fast updating of hashes
  std::vector<HighsInt> Astart;
  std::vector<HighsInt> Aindex;
  std::vector<u32> Acolor;

  std::vector<HighsInt> ARstart;
  std::vector<HighsInt> ARindex;
  std::vector<u32> ARcolor;

  // hashes for rows and columns
  std::vector<u64> colHashes;
  std::vector<u64> rowHashes;

  std::vector<HighsInt> linkCompressionStack;

  std::vector<HighsInt> colPosition;
  std::vector<HighsInt> groundSet;
  std::vector<HighsInt> permutation;
  std::vector<HighsInt> orbitPartition;
  std::vector<HighsInt> orbitSize;
  std::vector<HighsInt*> distinguishCands;

  std::vector<u64> currNodeCertificate;
  std::vector<u64> smallestLeaveCertificate;

  std::vector<HighsInt> firstLeaveColPartition;
  std::vector<u64> firstLeaveCertificate;
  std::vector<std::tuple<HighsInt,HighsInt,HighsUInt>> firstLeaveGraph;

  HighsInt numCol;
  HighsInt numRow;

  enum CellType : HighsUInt {
    kRowCell = 0,
    kColCell = 1,
  };

  std::vector<std::pair<HighsInt, CellType>> cellCreationStack;

  HighsHashTable<std::pair<HighsInt, CellType>> cellsOnRefinementStack;
  std::vector<std::pair<HighsInt, CellType>> refinementStack;

  // node in the search tree for finding automorphisms
  struct Node {
    HighsInt stackStart;
    HighsInt certificateEnd;
    CellType targetCellType;
    HighsInt targetCell;
    HighsInt lastDistiguished;
    bool onFirstPath;
  };

  std::vector<Node> nodeStack;

  HighsInt getColCellStart(HighsInt pos);
  HighsInt getRowCellStart(HighsInt pos);

  void backtrack(HighsInt backtrackStackNewEnd, HighsInt backtrackStackEnd);
  void cleanupBacktrack(HighsInt cellCreationStackPos);
  void removeFixPoints();
  void initializeGroundSet();
  std::vector<std::tuple<HighsInt,HighsInt,HighsUInt>> dumpCurrentGraph();
  bool mergeOrbits(HighsInt col1, HighsInt col2);
  HighsInt getOrbit(HighsInt col1);
  void initializeHashValues();
  bool isomorphicToFirstLeave();
  void partitionRefinement();
  std::pair<HighsInt,CellType> selectTargetCell();
  void distinguishVertex(HighsInt targetCell, CellType cellType);
  bool determineNextToDistinguish();
  void createNode();

  HighsInt rowCellSize(HighsInt rowCell) const {
    return rowPartitionLinks[rowCell] - rowCell;
  }
  HighsInt colCellSize(HighsInt colCell) const {
    return colPartitionLinks[colCell] - colCell;
  }

 public:
  void loadModelAsGraph(const HighsLp& lp, double epsilon);

  void run(HighsSymmetries& symmetries);
};

#endif