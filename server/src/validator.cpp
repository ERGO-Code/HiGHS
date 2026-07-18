// validator.cpp
#include "validator.h"
#include <sstream>

grpc::Status ValidateRequest(const highsserver::v1::SolveRequest& req) {
  const auto n_col = req.col_cost_size();
  if (n_col == 0)
    return {grpc::INVALID_ARGUMENT, "col_cost is empty"};
  if (req.col_lower_size() != n_col || req.col_upper_size() != n_col) {
    std::ostringstream os;
    os << "col arrays size mismatch: cost=" << n_col
       << " lower=" << req.col_lower_size()
       << " upper=" << req.col_upper_size();
    return {grpc::INVALID_ARGUMENT, os.str()};
  }
  if (req.col_integrality_size() != 0 && req.col_integrality_size() != n_col)
    return {grpc::INVALID_ARGUMENT, "col_integrality size mismatch"};

  const auto n_row = req.row_lower_size();
  if (req.row_upper_size() != n_row) {
    std::ostringstream os;
    os << "row arrays size mismatch: lower=" << n_row
       << " upper=" << req.row_upper_size();
    return {grpc::INVALID_ARGUMENT, os.str()};
  }

  // CSC: start length must be num_col + 1
  if (req.a_format_start_size() != n_col + 1) {
    std::ostringstream os;
    os << "a_format_start size must be num_col+1=" << (n_col + 1)
       << " got=" << req.a_format_start_size();
    return {grpc::INVALID_ARGUMENT, os.str()};
  }
  const auto nnz = req.a_format_index_size();
  if (req.a_format_value_size() != nnz)
    return {grpc::INVALID_ARGUMENT, "a_format_index/value size mismatch"};
  for (int i = 0; i <= n_col; ++i) {
    const auto s = req.a_format_start(i);
    if (s < 0 || s > nnz)
      return {grpc::INVALID_ARGUMENT, "a_format_start out of range"};
  }
  for (int k = 0; k < nnz; ++k) {
    const auto idx = req.a_format_index(k);
    if (idx < 0 || idx >= n_row) {
      std::ostringstream os;
      os << "a_format_index[" << k << "]=" << idx << " out of [0," << n_row << ")";
      return {grpc::INVALID_ARGUMENT, os.str()};
    }
  }
  return grpc::Status::OK;
}
