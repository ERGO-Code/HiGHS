// validator.h
// 输入校验：所有非法输入统一返回 grpc::Status(INVALID_ARGUMENT, ...)
#pragma once
#include <grpcpp/grpcpp.h>
#include "solver.pb.h"

grpc::Status ValidateRequest(const highsserver::v1::SolveRequest& req);
