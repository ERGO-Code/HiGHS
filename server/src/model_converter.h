// model_converter.h
// Converts SolveRequest to HighsModel. Fixes v1 issues:
//   - sets a_matrix_.num_col_/num_row_ (v1 omitted, caused HiGHS assertion)
//   - uses sense field directly (no "negate cost" hack)
//   - maps col_integrality for MIP support
#pragma once
#include <grpcpp/grpcpp.h>
#include "Highs.h"
#include "solver.pb.h"

grpc::Status FillHighsModel(const highsserver::v1::SolveRequest& req,
                            HighsModel& model);
