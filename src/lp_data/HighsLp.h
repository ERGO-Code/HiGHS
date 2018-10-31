#ifndef LP_DATA_LP_DATA_H
#define LP_DATA_LP_DATA_H

#include <cassert>
#include <string>
#include <vector>
#include <iostream>

#include "HConst.h"
// The free parser also reads fixed format MPS files but the fixed
// parser does not read free mps files.
enum class HighsMpsParserType { free, fixed };

// For now, but later change so HiGHS properties are string based so that new
// options (for debug and testing too) can be added easily. The options below
// are just what has been used to parse options from argv.
// todo: when creating the new options don't forget underscores for class
// variables but no underscores for struct
struct HighsOptions {
  int filename = 0;
  int presolve = 0;
  int crash = 0;
  int edgeWeight = 0;
  int price = 0;
  int pami = 0;
  int sip = 0;
  int scip = 0;

  double timeLimit = 0;
  double cut = 0;

  HighsMpsParserType parser_type = HighsMpsParserType::free;

  const char* fileName = "";
  const char* presolveMode = "";
  const char* edWtMode = "";
  const char* priceMode = "";
  const char* crashMode = "";
  const char* partitionFile = "";
};

class HighsLp {
 public:
  // Model data
  int numCol_;
  int numRow_;

  std::vector<int> Astart_;
  std::vector<int> Aindex_;
  std::vector<double> Avalue_;
  std::vector<double> colCost_;
  std::vector<double> colLower_;
  std::vector<double> colUpper_;
  std::vector<double> rowLower_;
  std::vector<double> rowUpper_;
};

// HiGHS status
enum class HighsStatus {
  OK,
  LpError,
  OptionsError,
  PresolveError,
  SolutionError,
  PostsolveError,
  NotImplemented
};

enum class HighsSolutionStatus {
  Unset,
  Unbounded,
  Infeasible,
  Feasible,
  Optimal,
};

enum class HighsInputStatus {
  OK,
  FileNotFound,
  ErrorMatrixDimensions,
  ErrorMatrixIndices,
  ErrorMatrixStart,
  ErrorMatrixValue,
  ErrorColBounds,
  ErrorRowBounds,
  ErrorObjective
};

struct HighsSolution {
  std::vector<double> colValue;
  std::vector<double> colDual;
  std::vector<double> rowValue;
};

// Return a string representation of SolutionStatus.
// Capitalized because it is ClassNameToString for the following three methods.
std::string HighsSolutionStatusToString(HighsSolutionStatus status) {
  switch (status) {
    case HighsSolutionStatus::Unset:
      return "Unset.";
      break;
    case HighsSolutionStatus::Unbounded:
      return "Unbounded.";
      break;
    case HighsSolutionStatus::Infeasible:
      return "Infeasible.";
      break;
    case HighsSolutionStatus::Feasible:
      return "Feasible.";
      break;
    case HighsSolutionStatus::Optimal:
      return "Optimal.";
      break;
  }
  return "";
}

// Return a string representation of HighsStatus.
std::string HighsStatusToString(HighsStatus status) {
  switch (status) {
    case HighsStatus::OK:
      return "OK.";
      break;
    case HighsStatus::LpError:
      return "Lp Error.";
      break;
    case HighsStatus::OptionsError:
      return "Options Error.";
      break;
    case HighsStatus::PresolveError:
      return "Presolve Error.";
      break;
    case HighsStatus::SolutionError:
      return "Solution Error.";
      break;
    case HighsStatus::PostsolveError:
      return "PostsolveError.";
      break;
    case HighsStatus::NotImplemented:
      return "Not implemented.";
      break;
  }
  return "";
}

// Return a string representation of ParseStatus.
std::string HighsInputStatusToString(HighsInputStatus status) {
  switch (status) {
    case HighsInputStatus::OK:
      return "OK.";
      break;
    case HighsInputStatus::FileNotFound:
      return "Error: File not found.";
      break;
    case HighsInputStatus::ErrorMatrixDimensions:
      return "Error Matrix Dimensions";
      break;
    case HighsInputStatus::ErrorMatrixIndices:
      return "Error Matrix Indices";
      break;
    case HighsInputStatus::ErrorMatrixStart:
      return "Error Matrix Start";
      break;
    case HighsInputStatus::ErrorMatrixValue:
      return "Error Matrix Value";
      break;
    case HighsInputStatus::ErrorColBounds:
      return "Error Col Bound";
      break;
    case HighsInputStatus::ErrorRowBounds:
      return "Error Row Bounds";
      break;
    case HighsInputStatus::ErrorObjective:
      return "Error Objective";
      break;
  }
  return "";
}

// If debug this method terminates the program when the status is not OK. If
// standard build it only prints a message.
void checkStatus(HighsStatus status) {
  assert(status == HighsStatus::OK);
  if (status != HighsStatus::OK)
    std::cout << "Unexpected status: " << HighsStatusToString(status);
}

HighsInputStatus checkLp(const HighsLp& lp) {
  // Check dimensions.
  if (lp.numCol_ <= 0 || lp.numRow_ <= 0)
    return HighsInputStatus::ErrorMatrixDimensions;

  // Check vectors.
  if ((int)lp.colCost_.size() != lp.numCol_)
    return HighsInputStatus::ErrorObjective;

  if ((int)lp.colLower_.size() != lp.numCol_ ||
      (int)lp.colUpper_.size() != lp.numCol_)
    return HighsInputStatus::ErrorColBounds;
  if ((int)lp.rowLower_.size() != lp.numRow_ ||
      (int)lp.rowUpper_.size() != lp.numRow_)
    return HighsInputStatus::ErrorRowBounds;

  for (int i = 0; i < lp.numRow_; i++)
    if (lp.rowLower_[i] < -HSOL_CONST_INF || lp.rowUpper_[i] > HSOL_CONST_INF)
      return HighsInputStatus::ErrorRowBounds;

  for (int j = 0; j < lp.numCol_; j++) {
    if (lp.colCost_[j] < -HSOL_CONST_INF || lp.colCost_[j] > HSOL_CONST_INF)
      return HighsInputStatus::ErrorObjective;

    if (lp.colLower_[j] < -HSOL_CONST_INF || lp.colUpper_[j] > HSOL_CONST_INF)
      return HighsInputStatus::ErrorColBounds;
    if (lp.colLower_[j] > lp.colUpper_[j] + kBoundTolerance)
      return HighsInputStatus::ErrorColBounds;
  }

  // Check matrix.
  const int nnz = lp.Avalue_.size();
  if (nnz <= 0) return HighsInputStatus::ErrorMatrixValue;
  if ((int)lp.Aindex_.size() != nnz)
    return HighsInputStatus::ErrorMatrixIndices;

  if ((int)lp.Astart_.size() != lp.numCol_ + 1)
    return HighsInputStatus::ErrorMatrixStart;
  for (int i = 0; i < lp.numCol_; i++) {
    if (lp.Astart_[i] > lp.Astart_[i + 1] || lp.Astart_[i] >= nnz ||
        lp.Astart_[i] < 0)
      return HighsInputStatus::ErrorMatrixStart;
  }

  for (int k = 0; k < nnz; k++) {
    if (lp.Aindex_[k] < 0 || lp.Aindex_[k] >= lp.numRow_)
      return HighsInputStatus::ErrorMatrixIndices;
    if (lp.Avalue_[k] < -HSOL_CONST_INF || lp.Avalue_[k] > HSOL_CONST_INF)
      return HighsInputStatus::ErrorRowBounds;
  }

  return HighsInputStatus::OK;
}

#endif