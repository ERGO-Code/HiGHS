#ifndef PDLP_HIPDLP_API_STATUS_HPP_
#define PDLP_HIPDLP_API_STATUS_HPP_

#include <string>
#include <vector>

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
  kUnboundedOrInfeasible,
  kIterationLimit,
  kTimeLimit,
  kError
};

// What a verified direction proves about the original LP. This is solver
// neutral: it never reaches HighsSolution and never reaches simplex ray state.
enum class CertificateKind {
  kNone,
  kPrimalInfeasibility,
  kPrimalUnboundedness,
  kAmbiguous
};

// Summary of what the independent verifier measured on the recovered
// direction. Residual is the stationarity residual of a Farkas certificate and
// stays zero for a primal ray.
struct CertificateMetrics {
  double direction_norm;
  double max_normalized_violation;
  double max_normalized_residual;
  double margin;
  HighsInt worst_row;
  HighsInt worst_col;

  CertificateMetrics();
};

struct Certificate {
  CertificateKind kind;
  bool valid;
  std::vector<double> column_direction;
  std::vector<double> row_direction;
  CertificateMetrics metrics;

  Certificate();
};

struct SolveResult {
  SolverStatus status;
  Error error;
  HighsInt iteration_count;
  bool has_solution;
  HighsSolution solution;
  Certificate certificate;

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