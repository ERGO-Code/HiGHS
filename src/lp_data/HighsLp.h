/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
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

#include "HConfig.h"
#include "lp_data/HConst.h"        // For HiGHS strategy options
#include "simplex/SimplexConst.h"  // For simplex strategy options
#include "simplex/SimplexStruct.h"  // For SimplexBasis

enum class LpAction {
  DUALISE = 0,
  PERMUTE,
  SCALE,
  NEW_COSTS,
  NEW_BOUNDS,
  NEW_BASIS,
  NEW_COLS,
  NEW_ROWS,
  DEL_COLS,
  DEL_ROWS,
  DEL_ROWS_BASIS_OK,
  SCALED_COL,
  SCALED_ROW,
  BACKTRACKING
};

class HighsLp;

class HighsLp {
 public:
  // Model data
  int numCol_ = 0;
  int numRow_ = 0;

  std::vector<int> Astart_;
  std::vector<int> Aindex_;
  std::vector<double> Avalue_;
  std::vector<double> colCost_;
  std::vector<double> colLower_;
  std::vector<double> colUpper_;
  std::vector<double> rowLower_;
  std::vector<double> rowUpper_;

  ObjSense sense_ = ObjSense::MINIMIZE;
  double offset_ = 0;

  std::string model_name_ = "";
  std::string lp_name_ = "";

  std::vector<std::string> row_names_;
  std::vector<std::string> col_names_;

  std::vector<int> integrality_;

  bool equalButForNames(const HighsLp& lp);
  bool operator==(const HighsLp& lp);
};

// Cost, column and row scaling factors
struct HighsScale {
  bool is_scaled_ = false;
  double cost_;
  std::vector<double> col_;
  std::vector<double> row_;
};

struct HighsSolutionParams {
  // Input to solution analysis method
  double primal_feasibility_tolerance;
  double dual_feasibility_tolerance;
  int primal_status = PrimalDualStatus::STATUS_NOTSET;
  int dual_status = PrimalDualStatus::STATUS_NOTSET;
  // Output from solution analysis method
  double objective_function_value;
  int num_primal_infeasibilities;
  double sum_primal_infeasibilities;
  double max_primal_infeasibility;
  int num_dual_infeasibilities;
  double sum_dual_infeasibilities;
  double max_dual_infeasibility;
};

struct HighsIterationCounts {
  int simplex = 0;
  int ipm = 0;
  int crossover = 0;
};

struct HighsSolution {
  std::vector<double> col_value;
  std::vector<double> col_dual;
  std::vector<double> row_value;
  std::vector<double> row_dual;
};

// To be the basis representation given back to the user. Values of
// HighsBasisStatus are defined in HConst.h
struct HighsBasis {
  bool valid_ = false;
  std::vector<HighsBasisStatus> col_status;
  std::vector<HighsBasisStatus> row_status;
};

// Set a basis to be logical for the LP
void setLogicalBasis(const HighsLp& lp, HighsBasis& basis);

// Make sure the sizes of solution and basis vectors are consistent
// with numRow_ and numCol_
bool isBasisConsistent(const HighsLp& lp, const HighsBasis& basis);
bool isSolutionRightSize(const HighsLp& lp, const HighsSolution& solution);
bool isBasisRightSize(const HighsLp& lp, const HighsBasis& basis);

// If debug this method terminates the program when the status is not OK. If
// standard build it only prints a message.
// void checkStatus(HighsStatus status);

void clearSolutionUtil(HighsSolution& solution);
void clearBasisUtil(HighsBasis& solution);
void clearLp(HighsLp& lp);

#endif
