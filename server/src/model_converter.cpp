// model_converter.cpp
#include "model_converter.h"

grpc::Status FillHighsModel(const highsserver::v1::SolveRequest& req,
                            HighsModel& model) {
  auto& lp = model.lp_;
  lp.num_col_ = req.col_cost_size();
  lp.num_row_ = req.row_lower_size();

  lp.col_cost_.assign(req.col_cost().begin(), req.col_cost().end());
  lp.col_lower_.assign(req.col_lower().begin(), req.col_lower().end());
  lp.col_upper_.assign(req.col_upper().begin(), req.col_upper().end());
  lp.row_lower_.assign(req.row_lower().begin(), req.row_lower().end());
  lp.row_upper_.assign(req.row_upper().begin(), req.row_upper().end());

  // Objective sense: use directly, no "negate cost" hack
  lp.sense_ = (req.sense() == highsserver::v1::OBJ_SENSE_MAX)
                  ? ObjSense::kMaximize : ObjSense::kMinimize;

  // CSC matrix (must set a_matrix_ dimensions, else HiGHS asserts)
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = lp.num_col_;
  lp.a_matrix_.num_row_ = lp.num_row_;
  // proto int64 → HighsInt, assign converts automatically
  lp.a_matrix_.start_.assign(req.a_format_start().begin(),
                             req.a_format_start().end());
  lp.a_matrix_.index_.assign(req.a_format_index().begin(),
                             req.a_format_index().end());
  lp.a_matrix_.value_.assign(req.a_format_value().begin(),
                             req.a_format_value().end());

  // Integer variables (MIP support)
  if (req.col_integrality_size() > 0) {
    lp.integrality_.resize(lp.num_col_);
    for (int i = 0; i < lp.num_col_; ++i) {
      using VT = highsserver::v1::VarType;
      switch (req.col_integrality(i)) {
        case VT::VAR_INTEGER:
        case VT::VAR_BINARY:
          lp.integrality_[i] = HighsVarType::kInteger; break;
        default:
          lp.integrality_[i] = HighsVarType::kContinuous;
      }
    }
  }
  return grpc::Status::OK;
}
