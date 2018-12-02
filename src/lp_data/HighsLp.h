/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsLp.h
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef LP_DATA_HIGHS_LP_H_
#define LP_DATA_HIGHS_LP_H_

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

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
  int nnz_;

  std::vector<int> Astart_;
  std::vector<int> Aindex_;
  std::vector<double> Avalue_;
  std::vector<double> colCost_;
  std::vector<double> colLower_;
  std::vector<double> colUpper_;
  std::vector<double> rowLower_;
  std::vector<double> rowUpper_;

  // 1 is default, -1 is maximize
  int sense_ = 1;
  double offset_ = 0;
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
  std::vector<double> rowDual;
};

bool isSolutionConsistent(const HighsLp& lp, const HighsSolution& solution);

// Return a string representation of SolutionStatus.
// Capitalized because it is ClassNameToString for the following three methods.
std::string HighsSolutionStatusToString(HighsSolutionStatus status);

// Return a string representation of HighsStatus.
std::string HighsStatusToString(HighsStatus status);

// Return a string representation of ParseStatus.
std::string HighsInputStatusToString(HighsInputStatus status);

// If debug this method terminates the program when the status is not OK. If
// standard build it only prints a message.
void checkStatus(HighsStatus status);

HighsInputStatus checkLp(const HighsLp& lp);

#endif
