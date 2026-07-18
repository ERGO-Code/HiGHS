// model_converter.h
// 将 SolveRequest 转换为 HighsModel。修复 v1 漏设 a_matrix_.num_col_/num_row_、
// 目标方向 hack、integrality 缺失。
#pragma once
#include <grpcpp/grpcpp.h>
#include "Highs.h"
#include "solver.pb.h"

grpc::Status FillHighsModel(const highsserver::v1::SolveRequest& req,
                            HighsModel& model);
