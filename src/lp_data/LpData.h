#ifndef LP_DATA_H
#define LP_DATA_H

#include <string>
#include <vector>

#include "HConst.h"
// The free parser also reads fixed format MPS files but the fixed
// parser does not read free mps files.
enum class MpsParserType { free, fixed };

// For now, but later change so HiGHS properties are string based so that new
// options (for debug and testing too) can be added easily.
struct Options {
  string filename = "";
  int presolve = 0;
  int crash = 0;
  int edgeWeight = 0;
  int price = 0;
  int pami = 0;
  int sip = 0;
  int scip = 0;

  double timeLimit = 0;
  double cut = 0;

  MpsParserType parser_type = MpsParserType::free;

  const char* fileName = "";
  const char* presolveMode = "";
  const char* edWtMode = "";
  const char* priceMode = "";
  const char* crashMode = "";
  const char* partitionFile = "";
};

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
};

// HiGHS status
enum Status {
  OK,
  InputError,
  FileNotFound,
  ParseError,
  ProblemReduced,
  ProblemReducedToEmpty,
  ReducedSolution,
  Postsolved,
  Infeasible,
  Unbounded,
  Optimal,
  NotImplemented
};

struct Solution {
  std::vector<double> colValue;
  std::vector<double> colDual;
  std::vector<double> rowValue;
};

// Return a string representation of status.
string toString(Status status) {
  switch (status) {
    case Status::OK:
      return "OK.";
      break;
    case Status::FileNotFound:
      return "Error: File not found.";
      break;
    case Status::ParseError:
      return "Parse error.";
      break;
    case Status::InputError:
      return "Input error.";
      break;
    case Status::ProblemReduced:
      return "Problem reduced.";
      break;
    case Status::ProblemReducedToEmpty:
      return "Problem reduced to empty.";
      break;
    case Status::ReducedSolution:
      return "Reduced problem solved.";
      break;
    case Status::Postsolved:
      return "Postsolved.";
      break;
    case Status::Infeasible:
      return "Infeasible.";
      break;
    case Status::Unbounded:
      return "Unbounded.";
      break;
    case Status::Optimal:
      return "Optimal.";
      break;
    case Status::NotImplemented:
      return "Not implemented.";
      break;
  }
  return "";
}

// If debug this method terminates the program when the status is not OK. If
// standard build it only prints a message.
void checkStatus(Status status) {
  assert(status == Status::OK);
  if (status != Status::OK)
    std::cout << "Unexpected status: " << toString(status);
}

LpError checkLp(const LpData& lp) {
  // Check dimensions.
  if (lp.numCol <= 0 || lp.numRow <= 0) return LpError::matrix_dimensions;

  // Check vectors.
  if ((int)lp.colCost.size() != lp.numCol) return LpError::objective;

  if ((int)lp.colLower.size() != lp.numCol ||
      (int)lp.colUpper.size() != lp.numCol)
    return LpError::col_bounds;
  if ((int)lp.rowLower.size() != lp.numRow ||
      (int)lp.rowUpper.size() != lp.numRow)
    return LpError::row_bounds;

  for (int i = 0; i < lp.numRow; i++)
    if (lp.rowLower[i] < -HSOL_CONST_INF || lp.rowUpper[i] > HSOL_CONST_INF)
      return LpError::row_bounds;

  for (int j = 0; j < lp.numCol; j++) {
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
  if ((int)lp.Aindex.size() != nnz) return LpError::matrix_indices;

  if ((int)lp.Astart.size() != lp.numCol + 1) return LpError::matrix_start;
  for (int i = 0; i < lp.numCol; i++) {
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