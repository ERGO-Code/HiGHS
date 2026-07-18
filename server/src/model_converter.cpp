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

  // 目标方向：直接用 sense，不再让客户端取负 cost
  lp.sense_ = (req.sense() == highsserver::v1::OBJ_SENSE_MAX)
                  ? ObjSense::kMaximize : ObjSense::kMinimize;

  // 约束矩阵 CSC（关键：补齐 a_matrix_ 自身维度，否则 HiGHS 内部断言失败）
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = lp.num_col_;
  lp.a_matrix_.num_row_ = lp.num_row_;
  // proto int64 → HighsInt，assign 自动转换
  lp.a_matrix_.start_.assign(req.a_format_start().begin(),
                             req.a_format_start().end());
  lp.a_matrix_.index_.assign(req.a_format_index().begin(),
                             req.a_format_index().end());
  lp.a_matrix_.value_.assign(req.a_format_value().begin(),
                             req.a_format_value().end());

  // 整数变量（MIP 支持）
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
