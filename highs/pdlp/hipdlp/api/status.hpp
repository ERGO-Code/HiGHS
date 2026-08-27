#ifndef PDLP_HIPDLP_API_STATUS_HPP_
#define PDLP_HIPDLP_API_STATUS_HPP_

#include <string>

#include "lp_data/HighsSolution.h"

namespace hipdlp {
namespace api {

enum class ErrorCode {
  kNone,
  kInvalidOption,
  kInvalidArgument,
  kNotReady,
  kBackendUnavailable,
  kUnsupportedCapability,
  kAllocationFailure,
  kBackendFailure,
  kInvalidModel,
  kInternal
};

struct Error {
  ErrorCode code;
  std::string message;

  Error();
  Error(ErrorCode code, const std::string& message);

  bool ok() const;
};

enum class SolverStatus {
  kNotSet,
  kOptimal,
  kInfeasible,
  kUnbounded,
  kIterationLimit,
  kTimeLimit,
  kError
};

struct SolveResult {
  SolverStatus status;
  Error error;
  HighsInt iteration_count;
  bool has_solution;
  HighsSolution solution;

  SolveResult();

  static SolveResult failure(ErrorCode code, const std::string& message);
  static SolveResult success(SolverStatus status, const HighsSolution& solution,
                             HighsInt iteration_count);

  bool ok() const;
};

HighsStatus applySolveResult(const SolveResult& result, HighsBasis& highs_basis,
                             HighsSolution& highs_solution,
                             HighsModelStatus& model_status,
                             HighsInfo& highs_info);

}  // namespace api
}  // namespace hipdlp

#endif