#ifndef PDLP_HIPDLP_ALGORITHM_FIXED_POINT_ERROR_HPP_
#define PDLP_HIPDLP_ALGORITHM_FIXED_POINT_ERROR_HPP_

#include <cstddef>

#include "pdlp/hipdlp/api/status.hpp"
#include "pdlp/hipdlp/backend/types.hpp"

namespace hipdlp {
namespace backend {
class Backend;
}
namespace algorithm {

struct FixedPointErrorInput {
  const backend::VectorHandle* next_primal;
  const backend::VectorHandle* next_dual;
  const backend::VectorHandle* reflected_primal;
  const backend::VectorHandle* reflected_dual;
  double omega;
  double eta;
  backend::VectorHandle* delta_primal;
  backend::VectorHandle* delta_dual;
  backend::VectorHandle* aty_delta;

  FixedPointErrorInput();
};

struct FixedPointErrorDiagnostics {
  std::size_t device_to_host_scalar_readbacks;

  FixedPointErrorDiagnostics();
};

struct FixedPointErrorResult {
  api::Error error;
  bool computed;
  double fixed_point_error;
  FixedPointErrorDiagnostics diagnostics;

  FixedPointErrorResult();

  bool ok() const;
};

  api::Error prepareFixedPointError(
    backend::Backend& backend, const backend::ProblemHandle& problem,
    const FixedPointErrorInput& input);

  FixedPointErrorResult computeFixedPointError(double primal_norm,
                        double dual_norm,
                        double cross_term, double omega,
                        double eta);

FixedPointErrorResult evaluateFixedPointError(
    backend::Backend& backend, const backend::ProblemHandle& problem,
    const FixedPointErrorInput& input);

}  // namespace algorithm
}  // namespace hipdlp

#endif