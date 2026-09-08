#ifndef PDLP_HIPDLP_BACKEND_CUDA_CUDA_BACKEND_HPP_
#define PDLP_HIPDLP_BACKEND_CUDA_CUDA_BACKEND_HPP_

#include <memory>

#include "pdlp/hipdlp/backend/backend.hpp"
#include "pdlp/hipdlp/backend/command_buffer.hpp"

namespace hipdlp {
namespace backend {
namespace cuda {

struct CudaProblemResource;
struct CudaVectorResource;
struct CudaWorkspaceResource;
struct CudaCommandBufferResource;

class CudaBackend : public Backend, public CommandBufferRecorder {
 public:
  ~CudaBackend() override;

  bool timingUsesEvents() const override;
  api::Error prepareTiming() override;
  api::Error beginTiming() override;
  api::Error endTiming(double& seconds) override;

  Descriptor descriptor() const override;
  api::Error queryDevice(DeviceInfo& device) const override;
  api::Error requireCapability(Capability capability) const override;
  api::Error bindProblem(const model::Problem& problem,
                         ProblemHandle& handle) override;
  api::Error createWorkspace(const ProblemHandle& problem,
                             WorkspaceHandle& workspace) override;
  api::Error prepareCheckWorkspace(const model::Problem& problem,
                                   const HighsLp& original_lp,
                                   WorkspaceHandle& workspace) override;
  api::Error prepareCheckFixedPoint(
      const ProblemHandle& problem, WorkspaceHandle& workspace,
      const VectorHandle& primal, const VectorHandle& dual,
      const VectorHandle& reflected_primal, const VectorHandle& reflected_dual,
      VectorHandle& delta_primal, VectorHandle& delta_dual,
      VectorHandle& aty_delta) override;
  api::Error evaluateCheck(const model::Problem& problem,
                           WorkspaceHandle& workspace, const CheckInput& input,
                           CheckResult& result) override;
  api::Error createVector(HighsInt size, VectorHandle& vector) override;
  api::Error createBatchedVector(HighsInt size, HighsInt batch_size,
                                 VectorHandle& vector) override;
  api::Error uploadVector(const std::vector<double>& values,
                          VectorHandle& vector) override;
  api::Error initialize(VectorHandle& primal, VectorHandle& dual) override;
  api::Error fill(VectorHandle& vector, double value) override;
  api::Error copy(const VectorHandle& source,
                  VectorHandle& destination) override;
  api::Error scale(double factor, VectorHandle& vector) override;
  api::Error axpy(double factor, const VectorHandle& input,
                  VectorHandle& output) override;
  api::Error multiply(const ProblemHandle& problem, const VectorHandle& input,
                      VectorHandle& output) override;
  api::Error multiplyTranspose(const ProblemHandle& problem,
                               const VectorHandle& input,
                               VectorHandle& output) override;
  api::Error multiplyBatched(const ProblemHandle& problem,
                             const VectorHandle& input,
                             VectorHandle& output) override;
  api::Error fusedProjectedPrimalUpdate(const ProblemHandle& problem,
                                        WorkspaceHandle& workspace,
                                        const VectorHandle& current_primal,
                                        const VectorHandle& transpose_product,
                                        double step_size,
                                        VectorHandle& next_primal) override;
  api::Error fusedProjectedDualUpdate(const ProblemHandle& problem,
                                      WorkspaceHandle& workspace,
                                      const VectorHandle& current_dual,
                                      const VectorHandle& extrapolated_product,
                                      double step_size,
                                      VectorHandle& next_dual) override;
  api::Error fusedVectorBlend(WorkspaceHandle& workspace,
                              const VectorHandle& primal_left,
                              const VectorHandle& dual_left,
                              const VectorHandle& primal_right,
                              const VectorHandle& dual_right,
                              double left_weight, VectorHandle& primal_output,
                              VectorHandle& dual_output) override;
  api::Error fusedAverageUpdate(WorkspaceHandle& workspace,
                                const VectorHandle& primal,
                                const VectorHandle& dual, double weight,
                                VectorHandle& primal_sum,
                                VectorHandle& dual_sum) override;
  api::Error dot(const VectorHandle& left, const VectorHandle& right,
                 double& result) override;
  api::Error norm2(const VectorHandle& vector, double& result) override;
  api::Error differenceNorm2(const VectorHandle& left,
                             const VectorHandle& right,
                             double& result) override;
  api::Error reduceConvergenceMetrics(const VectorHandle& primal_residual,
                                      const VectorHandle& dual_residual,
                                      double primal_objective,
                                      double dual_objective,
                                      ConvergenceMetrics& metrics) override;
  api::Error beginRecording() override;
  api::Error endRecording(CommandBufferHandle& command_buffer) override;
  api::Error cancelRecording() override;
  api::Error replay(const CommandBufferHandle& command_buffer,
                    CommandBufferStatistics& statistics) override;
  api::Error synchronize() override;
    api::Error prepareVectorDownloads(const VectorDownload* downloads,
                                                                        std::size_t count) override;
    api::Error downloadVectors(const VectorDownload* downloads,
                                                         std::size_t count) override;
  api::Error downloadVector(const VectorHandle& vector,
                            std::vector<double>& values) override;
  api::Error downloadSolution(const VectorHandle& primal,
                              const VectorHandle& dual,
                              std::vector<double>& primal_values,
                              std::vector<double>& dual_values) override;

 private:
  struct Impl;

  explicit CudaBackend(std::unique_ptr<Impl> impl);

  CudaProblemResource* problemResource(ProblemHandle& handle) const;
  const CudaProblemResource* problemResource(const ProblemHandle& handle) const;
  CudaVectorResource* vectorResource(VectorHandle& handle) const;
  const CudaVectorResource* vectorResource(const VectorHandle& handle) const;
  CudaWorkspaceResource* workspaceResource(WorkspaceHandle& handle) const;
  const CudaCommandBufferResource* commandBufferResource(
      const CommandBufferHandle& handle) const;
    api::Error prepareDownloadBuffers(const VectorDownload* downloads,
                                                                        std::size_t count);

  std::unique_ptr<Impl> impl_;

  friend api::Error createCudaBackend(std::unique_ptr<Backend>& backend);
};

Descriptor cudaDescriptor();
api::Error createCudaBackend(std::unique_ptr<Backend>& backend);

}  // namespace cuda
}  // namespace backend
}  // namespace hipdlp

#endif