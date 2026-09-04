#include "pdlp/hipdlp/api/status.hpp"

namespace hipdlp {
namespace api {

Error::Error() : code(ErrorCode::kNone) {}

Error::Error(const ErrorCode code, const std::string& message)
    : code(code), message(message) {}

bool Error::ok() const { return code == ErrorCode::kNone; }

CertificateMetrics::CertificateMetrics()
    : direction_norm(0.0),
      max_normalized_violation(0.0),
      max_normalized_residual(0.0),
      margin(0.0),
      worst_row(-1),
      worst_col(-1) {}

Certificate::Certificate() : kind(CertificateKind::kNone), valid(false) {}

SolveResult::SolveResult()
    : status(SolverStatus::kNotSet), iteration_count(0), has_solution(false) {}

SolveResult SolveResult::failure(const ErrorCode code,
                                 const std::string& message) {
  SolveResult result;
  result.status = SolverStatus::kError;
  result.error = Error(code, message);
  return result;
}

SolveResult SolveResult::success(const SolverStatus status,
                                 const HighsSolution& solution,
                                 const HighsInt iteration_count) {
  SolveResult result;
  result.status = status;
  result.iteration_count = iteration_count;
  result.has_solution = true;
  result.solution = solution;
  return result;
}

bool SolveResult::ok() const {
  return error.ok() && status != SolverStatus::kNotSet &&
         status != SolverStatus::kError;
}

HighsStatus applySolveResult(const SolveResult& result, HighsBasis& highs_basis,
                             HighsSolution& highs_solution,
                             HighsModelStatus& model_status,
                             HighsInfo& highs_info) {
  resetModelStatusAndHighsInfo(model_status, highs_info);
  highs_basis.valid = false;
  highs_solution.clear();
  highs_info.pdlp_iteration_count = result.iteration_count;

  if (!result.ok() || !result.has_solution || !result.solution.value_valid) {
    model_status = HighsModelStatus::kSolveError;
    return HighsStatus::kError;
  }

  switch (result.status) {
    case SolverStatus::kOptimal:
      model_status = HighsModelStatus::kOptimal;
      break;
    case SolverStatus::kInfeasible:
      model_status = HighsModelStatus::kInfeasible;
      break;
    case SolverStatus::kUnbounded:
      model_status = HighsModelStatus::kUnbounded;
      break;
    case SolverStatus::kUnboundedOrInfeasible:
      model_status = HighsModelStatus::kUnboundedOrInfeasible;
      break;
    case SolverStatus::kIterationLimit:
      model_status = HighsModelStatus::kIterationLimit;
      break;
    case SolverStatus::kTimeLimit:
      model_status = HighsModelStatus::kTimeLimit;
      break;
    default:
      model_status = HighsModelStatus::kSolveError;
      return HighsStatus::kError;
  }

  highs_solution = result.solution;
  return HighsStatus::kOk;
}

}  // namespace api
}  // namespace hipdlp