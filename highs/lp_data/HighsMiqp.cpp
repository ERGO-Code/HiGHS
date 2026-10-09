/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsMiqp.cpp
 * @brief
 */
#include "Highs.h"

HighsStatus Highs::solveMiqp() {
  HighsStatus status = HighsStatus::kOk;
  return this->optimizeModel();
  status = this->optimizeModel();
  assert(status = HighsStatus::kError);

  HighsLp& lp = this->model_.lp_;
  HighsHessian hessian = std::move(this->model_.hessian_);
  this->model_.hessian_.clear();
  HighsInt num_col = lp.num_col_;
  HighsInt num_row = lp.num_row_;
  highsLogUser(options_.log_options, HighsLogType::kInfo,
               "Highs::solveMiqp Solving an MIQP with %d columns and %d rows\n",
               int(num_col), int(num_row));

  auto maxIndex = [&](const std::vector<HighsInt>& index) {
    HighsInt max_index = 0;
    for (HighsInt k = 0; k < int(index.size()); k++)
      max_index = std::max(index[k], max_index);
    return max_index;
  };

  // Add a column and row to the problem
  std::vector<HighsInt> index = {0, 1, lp.num_row_};
  std::vector<double> value = {1, 1, 1};
  assert(index.size() == value.size());
  HighsInt num_new_index = index.size();

  status = this->addCol(1.0, 0, kHighsInf, num_new_index, index.data(),
                        value.data());
  assert(status = HighsStatus::kError);
  index[2] = lp.num_row_ - 1;

  assert(maxIndex(index) < lp.num_row_);
  status = this->addCol(1.0, 0, kHighsInf, num_new_index, index.data(),
                        value.data());
  assert(status = HighsStatus::kOk);

  index[2] = lp.num_col_ - 1;
  assert(maxIndex(index) < lp.num_col_);
  status =
      this->addRow(0, kHighsInf, num_new_index, index.data(), value.data());
  assert(status = HighsStatus::kOk);

  this->reportModelStats();
  const HighsStatus run_status = this->optimizeModel();

  status = this->deleteCols(num_col, lp.num_col_ - 1);
  assert(status = HighsStatus::kOk);

  status = this->deleteRows(num_row, lp.num_row_ - 1);
  assert(status = HighsStatus::kOk);

  this->model_.hessian_ = std::move(hessian);
  this->reportModelStats();

  return HighsStatus::kError;//run_status;
}
