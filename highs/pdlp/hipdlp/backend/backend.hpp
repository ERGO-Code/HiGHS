#ifndef PDLP_HIPDLP_BACKEND_BACKEND_HPP_
#define PDLP_HIPDLP_BACKEND_BACKEND_HPP_

#include <vector>

#include "pdlp/hipdlp/api/status.hpp"
#include "pdlp/hipdlp/backend/types.hpp"
#include "pdlp/hipdlp/model/problem.hpp"

namespace hipdlp {
namespace backend {

struct Descriptor {
  Kind kind;
  const char* name;
  bool build_enabled;
  bool implemented;
  PrecisionDescriptor precision;
  Capabilities capabilities;
};

class Backend {
 public:
  virtual ~Backend() {}

  virtual Descriptor descriptor() const = 0;
  virtual api::Error queryDevice(DeviceInfo& device) const = 0;
  virtual api::Error requireCapability(Capability capability) const = 0;
  virtual api::Error bindProblem(const model::Problem& problem,
                                 ProblemHandle& handle) = 0;
  virtual api::Error createWorkspace(const ProblemHandle& problem,
                                     WorkspaceHandle& workspace) = 0;
  virtual api::Error createVector(HighsInt size, VectorHandle& vector) = 0;
  virtual api::Error createBatchedVector(HighsInt size, HighsInt batch_size,
                                         VectorHandle& vector) = 0;
  virtual api::Error uploadVector(const std::vector<double>& values,
                                  VectorHandle& vector) = 0;
  virtual api::Error initialize(VectorHandle& primal, VectorHandle& dual) = 0;
  virtual api::Error fill(VectorHandle& vector, double value) = 0;
  virtual api::Error copy(const VectorHandle& source,
                          VectorHandle& destination) = 0;
  virtual api::Error scale(double factor, VectorHandle& vector) = 0;
  virtual api::Error axpy(double factor, const VectorHandle& input,
                          VectorHandle& output) = 0;
  virtual api::Error multiply(const ProblemHandle& problem,
                              const VectorHandle& input,
                              VectorHandle& output) = 0;
  virtual api::Error multiplyTranspose(const ProblemHandle& problem,
                                       const VectorHandle& input,
                                       VectorHandle& output) = 0;
  virtual api::Error multiplyBatched(const ProblemHandle& problem,
                                     const VectorHandle& input,
                                     VectorHandle& output) = 0;
  virtual api::Error fusedProjectedPrimalUpdate(
      const ProblemHandle& problem, WorkspaceHandle& workspace,
      const VectorHandle& current_primal, const VectorHandle& transpose_product,
      double step_size, VectorHandle& next_primal) = 0;
  virtual api::Error fusedProjectedDualUpdate(
      const ProblemHandle& problem, WorkspaceHandle& workspace,
      const VectorHandle& current_dual,
      const VectorHandle& extrapolated_product, double step_size,
      VectorHandle& next_dual) = 0;
  virtual api::Error fusedVectorBlend(
      WorkspaceHandle& workspace, const VectorHandle& primal_left,
      const VectorHandle& dual_left, const VectorHandle& primal_right,
      const VectorHandle& dual_right, double left_weight,
      VectorHandle& primal_output, VectorHandle& dual_output) = 0;
  virtual api::Error fusedAverageUpdate(WorkspaceHandle& workspace,
                                        const VectorHandle& primal,
                                        const VectorHandle& dual, double weight,
                                        VectorHandle& primal_sum,
                                        VectorHandle& dual_sum) = 0;
  virtual api::Error dot(const VectorHandle& left, const VectorHandle& right,
                         double& result) = 0;
  virtual api::Error norm2(const VectorHandle& vector, double& result) = 0;
  virtual api::Error differenceNorm2(const VectorHandle& left,
                                     const VectorHandle& right,
                                     double& result) = 0;
  virtual api::Error reduceConvergenceMetrics(
      const VectorHandle& primal_residual, const VectorHandle& dual_residual,
      double primal_objective, double dual_objective,
      ConvergenceMetrics& metrics) = 0;
  virtual api::Error synchronize() = 0;
  virtual api::Error downloadVector(const VectorHandle& vector,
                                    std::vector<double>& values) = 0;
  virtual api::Error downloadSolution(const VectorHandle& primal,
                                      const VectorHandle& dual,
                                      std::vector<double>& primal_values,
                                      std::vector<double>& dual_values) = 0;
};

}  // namespace backend
}  // namespace hipdlp

#endif