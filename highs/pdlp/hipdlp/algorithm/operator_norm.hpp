#ifndef PDLP_HIPDLP_ALGORITHM_OPERATOR_NORM_HPP_
#define PDLP_HIPDLP_ALGORITHM_OPERATOR_NORM_HPP_

#include <cstddef>

#include "pdlp/hipdlp/api/status.hpp"
#include "pdlp/hipdlp/backend/types.hpp"

namespace hipdlp {
namespace backend {
class Backend;
}
namespace model {
struct SparseMatrix;
}
namespace algorithm {

struct OperatorNormParameters {
  HighsInt iteration_budget;
  double relative_tolerance;

  OperatorNormParameters();
};

struct OperatorNormDiagnostics {
  std::size_t resource_creation_calls_before_iterations;
  std::size_t resource_creation_calls_during_iterations;
  std::size_t device_to_host_scalar_readbacks_per_iteration;
  std::size_t device_to_host_scalar_readbacks_during_iterations;

  OperatorNormDiagnostics();
};

struct OperatorNormResult {
  api::Error error;
  double operator_norm_squared;
  HighsInt performed_iterations;
  double final_relative_change;
  bool tolerance_met;
  bool empty_matrix;
  OperatorNormDiagnostics diagnostics;

  OperatorNormResult();

  bool ok() const;
};

api::Error computeOperatorNormSquaredUpperBound(
    const model::SparseMatrix& matrix, double& upper_bound_squared);

OperatorNormResult estimateOperatorNormSquared(
    backend::Backend& backend, const backend::ProblemHandle& problem,
    HighsInt num_rows, HighsInt num_cols,
    const OperatorNormParameters& parameters = OperatorNormParameters());

}  // namespace algorithm
}  // namespace hipdlp

#endif