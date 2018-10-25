#ifndef LP_DATA_H
#define LP_DATA_H

#include <string>
#include <vector>

#include "HConst.h"

enum class LpError {
  none,
  matrix_dimensions,
  matrix_indices,
  matrix_start,
  matrix_value,
  col_bounds,
  row_bounds,
  objective
};

class LpData {
 public:
  // Model data
  int numCol;
  int numRow;

  std::vector<int> Astart;
  std::vector<int> Aindex;
  std::vector<double> Avalue;
  std::vector<double> colCost;
  std::vector<double> colLower;
  std::vector<double> colUpper;
  std::vector<double> rowLower;
  std::vector<double> rowUpper;

  static LpError checkLp(const LpData& lp) const;
};


struct Options {
  string filename = "";
  int presolve = 0, int crash = 0;
  int edgeWeight = 0;
  int price = 0;
  int pami = 0;
  int sip = 0;
  int scip = 0;
  int timeLimit = 0;
  double cut = 0;
  const char *fileName = "";
  const char *presolveMode = "";
  const char *edWtMode = "";
  const char *priceMode = "";
  const char *crashMode = "";
}

enum Status {
  OK,
  InputError,
  FileNotFound,
  ParseError,
  ProblemReduced,
  ProblemReducedToEmpty,
  Presolved,
  ReducedSolution,
  Postsolved,
  SimplexCleanUpFinished,
  Unknown,
  Infeasible,
  Unbounded,
  Optimal,
  NotImplemented
};

void printStatus(Status status) {
  switch (status) {
    case Status::OK:
      std::cout << "OK";
      break;
    case Status::FileNotFound:
      std::cout << "Error: File not found.";
      break;
    case Status::ParseError:
      std::cout << "Parse error.";
      break;
    case Status::Presolved:
      std::cout << "Presolve finished." break;
    case Status::ReducedSolution:
      std::cout << "Reduced problem solved";
      break;
    case Status::Postsolved:
      std::cout << "Postsolved";
      break;
  }
}

checkStatus(Status status) {
  if (status != Status::OK) {
    printStatus(status);
    if (status == Status::InputError) printHelp(argv[0]);
    exit(0);
  }
}

LpError checkLp(const LpData& lp) const {
  // Check dimensions.
  if (lp.numCol <= 0 || lp.numRow <= 0) return LpError::matrix_dimensions;

  // Check vectors.
  if (lp.colCost.size() != lp.numCol) return LpError::objective;

  if (lp.colLower.size() != lp.numCol || lp.colUpper.size() != lp.numCol)
    return LpError::col_bounds;
  if (lp.rowLower.size() != lp.numRow || lp.rowUpper.size() != lp.numRow)
    return LpError::row_bounds;

  for (int i = 0; i < numRow; i++)
    if (lp.rowLower[i] < -HSOL_CONST_INF || lp.rowUpper[i] > HSOL_CONST_INF)
      return LpError::row_bounds;

  for (int j = 0; j < numCol; j++) {
    if (lp.colCost[j] < -HSOL_CONST_INF || lp.colCost[j] > HSOL_CONST_INF)
      return LpError::objective;

    if (lp.colLower[j] < -HSOL_CONST_INF || lp.colUpper[j] > HSOL_CONST_INF)
      return LpError::col_bounds;
    if (lp.colLower[j] > lp.colUpper[j] + kBoundTolerance)
      return LpError::col_bounds;
  }

  // Check matrix.
  const int nnz = lp.Avalue.size();
  if (nnz <= 0) return LpError::matrix_value;
  if (lp.Aindex.size() != nnz) return LpError::matrix_indices;

  if (lp.Astart.size() != numCol + 1) return LpError::matrix_start;
  for (int i = 0; i < numCol; i++) {
    if (lp.Astart[i] > lp.Astart[i + 1] || lp.Astart[i] >= nnz ||
        lp.Astart[i] < 0)
      return LpError::matrix_start;
  }

  for (int k = 0; k < nnz; k++) {
    if (lp.Aindex[k] < 0 || lp.Aindex[k] >= lp.numRow)
      return LpError::matrix_indices;
    if (lp.Avalue[k] < -HSOL_CONST_INF || lp.Avalue[k] > HSOL_CONST_INF)
      return LpError::matrix_value;
  }

  return LpError::none;
}

#endif