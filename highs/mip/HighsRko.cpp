/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsRko.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

// This method should return true if an integer feasible solution
// (returned as solution) has been found
bool rkoHeuristic(const HighsLp* lp, std::vector<double>& solution) {
  if (lp->mip_type_ != kMipTypeKnapsack) return false;
  printf("Calling the RKO heuristic for a knapsack problem with %d items\n",
         int(lp->num_col_));

  const HighsInt num_col = lp->num_col_;
  if (num_col <= 0 || lp->num_row_ != 1) return false;

  solution.assign(num_col, 0.0);

  std::vector<double> weight(num_col, 0.0);
  const HighsSparseMatrix& matrix = lp->a_matrix_;
  if (matrix.isRowwise()) {
    if (matrix.start_.size() < 2) return false;
    for (HighsInt iEl = matrix.start_[0]; iEl < matrix.start_[1]; ++iEl) {
      const HighsInt iCol = matrix.index_[iEl];
      if (iCol >= 0 && iCol < num_col) weight[iCol] = matrix.value_[iEl];
    }
  } else if (matrix.isColwise()) {
    for (HighsInt iCol = 0; iCol < num_col; ++iCol) {
      const HighsInt end = matrix.p_end_.empty() ? matrix.start_[iCol + 1]
                                                 : matrix.p_end_[iCol];
      for (HighsInt iEl = matrix.start_[iCol]; iEl < end; ++iEl) {
        if (matrix.index_[iEl] == 0) weight[iCol] = matrix.value_[iEl];
      }
    }
  } else {
    return false;
  }

  double capacity = lp->row_upper_[0];
  if (capacity >= kHighsInf && lp->row_lower_[0] > -kHighsInf)
    capacity = -lp->row_lower_[0];
  if (capacity < 0 || capacity >= kHighsInf) return false;

  std::vector<double> profit(num_col, 0.0);
  HighsInt num_negative_cost = 0;
  for (HighsInt iCol = 0; iCol < num_col; ++iCol) {
    if (lp->col_cost_[iCol] < 0) ++num_negative_cost;
  }
  const bool costs_are_negated = num_negative_cost > num_col / 2;
  for (HighsInt iCol = 0; iCol < num_col; ++iCol) {
    profit[iCol] = costs_are_negated ? -lp->col_cost_[iCol] : lp->col_cost_[iCol];
  }

  std::vector<double> best_solution(num_col, 0.0);
  double best_profit = 0.0;

  auto tryKeys = [&](const std::vector<double>& key) {
    std::vector<HighsInt> order(num_col);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](const HighsInt a, const HighsInt b) {
                if (key[a] == key[b]) return a < b;
                return key[a] > key[b];
              });

    std::vector<double> candidate_solution(num_col, 0.0);
    double used_capacity = 0.0;
    double candidate_profit = 0.0;
    for (const HighsInt iCol : order) {
      if (profit[iCol] <= 0 || weight[iCol] <= 0) continue;
      if (used_capacity + weight[iCol] <= capacity + 1e-9) {
        candidate_solution[iCol] = 1.0;
        used_capacity += weight[iCol];
        candidate_profit += profit[iCol];
      }
    }
    if (candidate_profit > best_profit + 1e-9) {
      best_profit = candidate_profit;
      best_solution = candidate_solution;
    }
  };

  std::vector<double> key(num_col, 0.0);

  for (HighsInt iCol = 0; iCol < num_col; ++iCol)
    key[iCol] = weight[iCol] > 0 ? profit[iCol] / weight[iCol] : 0.0;
  tryKeys(key);

  for (HighsInt iCol = 0; iCol < num_col; ++iCol) key[iCol] = profit[iCol];
  tryKeys(key);

  for (HighsInt pass = 0; pass < 64; ++pass) {
    for (HighsInt iCol = 0; iCol < num_col; ++iCol) {
      const HighsInt hash =
          (1103515245 * (iCol + 1) + 12345 * (pass + 1)) & 0x7fffffff;
      const double random_key = double(hash % 10000) / 10000.0;
      const double density = weight[iCol] > 0 ? profit[iCol] / weight[iCol] : 0.0;
      key[iCol] = density * (0.75 + 0.5 * random_key);
    }
    tryKeys(key);
  }

  solution = best_solution;
  return true;
}
