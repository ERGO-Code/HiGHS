// validator.h
// Input validation: all invalid inputs return grpc::Status(INVALID_ARGUMENT, ...)
#pragma once
#include <grpcpp/grpcpp.h>
#include "solver.pb.h"

grpc::Status ValidateRequest(const highsserver::v1::SolveRequest& req);
