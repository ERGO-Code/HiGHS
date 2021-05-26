/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HStruct.h
 * @brief Structs for HiGHS
 */
#ifndef LP_DATA_HSTRUCT_H_
#define LP_DATA_HSTRUCT_H_

#include <vector>

#include "lp_data/HConst.h"

struct HighsIterationCounts {
  HighsInt simplex = 0;
  HighsInt ipm = 0;
  HighsInt crossover = 0;
};

struct HighsScale {
  bool is_scaled = false;
  double cost;
  std::vector<double> col;
  std::vector<double> row;
};

struct HighsSolution {
  bool value_valid = false;
  bool dual_valid = false;
  std::vector<double> col_value;
  std::vector<double> col_dual;
  std::vector<double> row_value;
  std::vector<double> row_dual;
};

struct HighsBasis {
  bool valid = false;
  std::vector<HighsBasisStatus> col_status;
  std::vector<HighsBasisStatus> row_status;
};

struct HighsSolutionParams {
  // Input to solution analysis method
  double primal_feasibility_tolerance;
  double dual_feasibility_tolerance;
  HighsInt primal_solution_status;
  HighsInt dual_solution_status;
  // Output from solution analysis method
  double objective_function_value;
  HighsInt num_primal_infeasibility;
  double sum_primal_infeasibility;
  double max_primal_infeasibility;
  HighsInt num_dual_infeasibility;
  double sum_dual_infeasibility;
  double max_dual_infeasibility;
};

#endif /* LP_DATA_HSTRUCT_H_ */
