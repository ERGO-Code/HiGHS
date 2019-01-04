/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsLp.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "HighsLp.h"
#include "HConst.h"

// If debug this method terminates the program when the status is not OK. If
// standard build it only prints a message.
void checkStatus(HighsStatus status) {
  assert(status == HighsStatus::OK);
  if (status != HighsStatus::OK)
    std::cout << "Unexpected status: " << HighsStatusToString(status);
}

bool isSolutionConsistent(const HighsLp& lp, const HighsSolution& solution) {
  if (solution.colDual_.size() == (size_t)lp.numCol_ ||
      solution.colValue_.size() == (size_t)lp.numCol_ ||
      solution.rowDual_.size() == (size_t)lp.numRow_ ||
      solution.rowValue_.size() == (size_t)lp.numRow_)
    return true;
  return false;
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
    if (lp.rowLower_[i] < -HIGHS_CONST_INF || lp.rowUpper_[i] > HIGHS_CONST_INF)
      return HighsInputStatus::ErrorRowBounds;

  for (int j = 0; j < lp.numCol_; j++) {
    if (lp.colCost_[j] < -HIGHS_CONST_INF || lp.colCost_[j] > HIGHS_CONST_INF)
      return HighsInputStatus::ErrorObjective;

    if (lp.colLower_[j] < -HIGHS_CONST_INF || lp.colUpper_[j] > HIGHS_CONST_INF)
      return HighsInputStatus::ErrorColBounds;
    if (lp.colLower_[j] > lp.colUpper_[j] + kBoundTolerance)
      return HighsInputStatus::ErrorColBounds;
  }

  // Check matrix.
  if ((size_t)lp.nnz_ != lp.Avalue_.size())
    return HighsInputStatus::ErrorMatrixValue;
  if (lp.nnz_ <= 0) return HighsInputStatus::ErrorMatrixValue;
  if ((int)lp.Aindex_.size() != lp.nnz_)
    return HighsInputStatus::ErrorMatrixIndices;

  if ((int)lp.Astart_.size() != lp.numCol_ + 1)
    return HighsInputStatus::ErrorMatrixStart;
  for (int i = 0; i < lp.numCol_; i++) {
    if (lp.Astart_[i] > lp.Astart_[i + 1] || lp.Astart_[i] >= lp.nnz_ ||
        lp.Astart_[i] < 0)
      return HighsInputStatus::ErrorMatrixStart;
  }

  for (int k = 0; k < lp.nnz_; k++) {
    if (lp.Aindex_[k] < 0 || lp.Aindex_[k] >= lp.numRow_)
      return HighsInputStatus::ErrorMatrixIndices;
    if (lp.Avalue_[k] < -HIGHS_CONST_INF || lp.Avalue_[k] > HIGHS_CONST_INF)
      return HighsInputStatus::ErrorRowBounds;
  }

  return HighsInputStatus::OK;
}

// Return a string representation of HighsStatus.
std::string HighsStatusToString(HighsStatus status) {
  switch (status) {
    case HighsStatus::OK:
      return "OK";
      break;
    case HighsStatus::Init:
      return "Init";
      break;
    case HighsStatus::LpError:
      return "Lp Error";
      break;
    case HighsStatus::OptionsError:
      return "Options Error";
      break;
    case HighsStatus::PresolveError:
      return "Presolve Error";
      break;
    case HighsStatus::SolutionError:
      return "Solution Error";
      break;
    case HighsStatus::PostsolveError:
      return "Postsolve Error";
      break;
    case HighsStatus::NotImplemented:
      return "Not implemented";
      break;
    case HighsStatus::Unbounded:
      return "Unbounded";
      break;
    case HighsStatus::Infeasible:
      return "Infeasible";
      break;
    case HighsStatus::Feasible:
      return "Feasible";
      break;
    case HighsStatus::Optimal:
      return "Optimal";
      break;
    case HighsStatus::Timeout:
      return "Timeout";
      break;
  }
  return "";
}

// Return a string representation of ParseStatus.
std::string HighsInputStatusToString(HighsInputStatus status) {
  switch (status) {
    case HighsInputStatus::OK:
      return "OK";
      break;
    case HighsInputStatus::FileNotFound:
      return "Error: File not found";
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
