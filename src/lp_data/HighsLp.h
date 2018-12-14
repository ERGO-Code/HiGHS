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

/** SCIP/HiGHS Objective sense */
enum objSense
{
  OBJSENSE_MINIMIZE = 1,
  OBJSENSE_MAXIMIZE = -1
};


// For now, but later change so HiGHS properties are string based so that new
// options (for debug and testing too) can be added easily. The options below
// are just what has been used to parse options from argv.
// todo: when creating the new options don't forget underscores for class
// variables but no underscores for struct
struct HighsOptions {
  std::string filenames = "";

  bool pami = 0;
  bool sip = 0;
  bool scip = 0;

  double timeLimit = 0;

  HighsMpsParserType parser_type = HighsMpsParserType::free;

  std::string fileName = "";
  std::string presolveMode = "";
  std::string edWtMode = "";
  std::string priceMode = "";
  std::string crashMode = "";
  std::string partitionFile = "";

  bool clean_up = false;
};

class HighsLp {
 public:
  // Model data
  int numCol_ = 0;
  int numRow_ = 0;
  int nnz_ = 0;

  std::vector<int> Astart_;
  std::vector<int> Aindex_;
  std::vector<double> Avalue_;
  std::vector<double> colCost_;
  std::vector<double> colLower_;
  std::vector<double> colUpper_;
  std::vector<double> rowLower_;
  std::vector<double> rowUpper_;

  // sense 1 = minimize, -1 = maximize
  int sense_ = 1;
  double offset_ = 0;
  std::string model_name_ = "";

void reportLp();
void reportLpBrief();
void reportLpDimensions();
void reportLpObjSense();
void reportLpColVec();
void reportLpRowVec();
void reportLpColMtx();


};

// HiGHS status
enum class HighsStatus {
  OK,
  Init,
  LpError,
  OptionsError,
  PresolveError,
  SolutionError,
  PostsolveError,
  NotImplemented,
  Unbounded,
  Infeasible,
  Feasible,
  Optimal,
  Timeout
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
  std::vector<double> colValue_;
  std::vector<double> colDual_;
  std::vector<double> rowValue_;
  std::vector<double> rowDual_;
};

struct HighsRanging {
  std::vector<double> colCostRangeUpValue_;
  std::vector<double> colCostRangeUpObjective_;
  std::vector<int>    colCostRangeUpInCol_;
  std::vector<int>    colCostRangeUpOutCol_;
  std::vector<double> colCostRangeDnValue_;
  std::vector<double> colCostRangeDnObjective_;
  std::vector<int>    colCostRangeDnInCol_;
  std::vector<int>    colCostRangeDnOutCol_;
  std::vector<double> rowBoundRangeUpValue_;
  std::vector<double> rowBoundRangeUpObjective_;
  std::vector<int>    rowBoundRangeUpInCol_;
  std::vector<int>    rowBoundRangeUpOutCol_;
  std::vector<double> rowBoundRangeDnValue_;
  std::vector<double> rowBoundRangeDnObjective_;
  std::vector<int>    rowBoundRangeDnInCol_;
  std::vector<int>    rowBoundRangeDnOutCol_;
};

// Make sure the dimensions of solution are the same as numRow_ and numCol_.
bool isSolutionConsistent(const HighsLp& lp, const HighsSolution& solution);

// Return a string representation of HighsStatus.
std::string HighsStatusToString(HighsStatus status);

// Return a string representation of ParseStatus.
std::string HighsInputStatusToString(HighsInputStatus status);

// If debug this method terminates the program when the status is not OK. If
// standard build it only prints a message.
void checkStatus(HighsStatus status);

HighsInputStatus checkLp(const HighsLp& lp);

#endif
