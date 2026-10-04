#ifndef PDLP_HIPDLP_BACKEND_BACKEND_HPP_
#define PDLP_HIPDLP_BACKEND_BACKEND_HPP_

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <new>
#include <stdexcept>
#include <vector>

#include "pdlp/hipdlp/api/status.hpp"
#include "pdlp/hipdlp/backend/check.hpp"
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

struct VectorDownload {
    const VectorHandle* vector;
    std::vector<double>* values;
};

class Backend {
 public:
  virtual ~Backend() {}

    virtual bool timingUsesEvents() const { return false; }
    virtual api::Error prepareTiming() {
        if (timing_active_)
            return api::Error(api::ErrorCode::kInvalidArgument,
                                                "cannot prepare an active backend timer");
        timing_prepared_ = true;
        return api::Error();
    }
    virtual api::Error beginTiming() {
        if (!timing_prepared_ || timing_active_)
            return api::Error(api::ErrorCode::kInvalidArgument,
                                                "backend timer must be prepared and inactive");
        api::Error error = synchronize();
        if (!error.ok()) return error;
        timing_start_ = std::chrono::steady_clock::now();
        timing_active_ = true;
        return api::Error();
    }
    virtual api::Error endTiming(double& seconds) {
        if (!timing_active_)
            return api::Error(api::ErrorCode::kInvalidArgument,
                                                "backend timer has not begun");
        timing_active_ = false;
        api::Error error = synchronize();
        if (!error.ok()) return error;
        const double elapsed = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - timing_start_).count();
        if (!std::isfinite(elapsed) || elapsed < 0.0)
            return api::Error(api::ErrorCode::kBackendFailure,
                                                "backend timer returned invalid elapsed time");
        seconds = elapsed;
        return api::Error();
    }

  virtual Descriptor descriptor() const = 0;
  virtual api::Error queryDevice(DeviceInfo& device) const = 0;
  virtual api::Error requireCapability(Capability capability) const = 0;
  virtual api::Error bindProblem(const model::Problem& problem,
                                 ProblemHandle& handle) = 0;
  virtual api::Error createWorkspace(const ProblemHandle& problem,
                                     WorkspaceHandle& workspace) = 0;
    virtual api::Error prepareCheckWorkspace(const model::Problem&,
                                                                                     const HighsLp&,
                                                                                     WorkspaceHandle&) {
        return api::Error(api::ErrorCode::kUnsupportedCapability,
                                            "backend does not support resident checks");
    }
    virtual api::Error prepareCheckFixedPoint(
            const ProblemHandle& problem, WorkspaceHandle& workspace,
            const VectorHandle& primal, const VectorHandle& dual,
            const VectorHandle& reflected_primal, const VectorHandle& reflected_dual,
            VectorHandle& delta_primal, VectorHandle& delta_dual,
            VectorHandle& aty_delta) {
        if (!problem.valid() || !workspace.valid())
            return api::Error(api::ErrorCode::kInvalidArgument,
                                                "HiPDLP check preparation problem or workspace is invalid");
        const VectorHandle* handles[] = {
                &primal, &dual, &reflected_primal, &reflected_dual,
                &delta_primal, &delta_dual, &aty_delta};
        for (std::size_t index = 0; index < 7; ++index) {
            if (!handles[index]->valid())
                return api::Error(api::ErrorCode::kInvalidArgument,
                                                    "HiPDLP check preparation vector is invalid");
            for (std::size_t previous = 0; previous < index; ++previous) {
                if (handles[index] == handles[previous])
                    return api::Error(api::ErrorCode::kInvalidArgument,
                                                        "HiPDLP check preparation vectors must be distinct");
            }
        }
        api::Error error = copy(primal, delta_primal);
        if (!error.ok()) return error;
        error = axpy(-1.0, reflected_primal, delta_primal);
        if (!error.ok()) return error;
        error = copy(dual, delta_dual);
        if (!error.ok()) return error;
        error = axpy(-1.0, reflected_dual, delta_dual);
        if (!error.ok()) return error;
        return multiplyTranspose(problem, delta_dual, aty_delta);
    }
    virtual api::Error evaluateCheck(const model::Problem&, WorkspaceHandle&,
                                                                     const CheckInput&, CheckResult& result) {
        result = CheckResult();
        return api::Error(api::ErrorCode::kUnsupportedCapability,
                                            "backend does not support resident checks");
    }
  // Bounds the original-model KKT metrics of a scaled candidate on the
  // backend. The default has no screen and leaves the result invalid.
  virtual api::Error screenOriginalKkt(const ProblemHandle&, WorkspaceHandle&,
                                       const VectorHandle&, const VectorHandle&,
                                       double, KktScreenResult& result) {
    result = KktScreenResult();
    return api::Error();
  }
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
    virtual api::Error difference(const VectorHandle& left, const VectorHandle& right,
                                                                VectorHandle& output) {
        if (!left.valid() || !right.valid() || !output.valid() ||
                &output == &left || &output == &right)
            return api::Error(api::ErrorCode::kInvalidArgument,
                                                "difference inputs and output must be valid and distinct");
        api::Error error = copy(left, output);
        if (!error.ok()) return error;
        return axpy(-1.0, right, output);
    }
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
    virtual api::Error fusedProjectedPrimalUpdateAndReflect(
            const ProblemHandle& problem, WorkspaceHandle& workspace,
            const VectorHandle& current, const VectorHandle& product,
            double step_size, VectorHandle& candidate, VectorHandle& reflected) {
        api::Error error = validateReflectionHandles(current, product, candidate,
                                                                                                 reflected);
        if (!error.ok()) return error;
        error = fusedProjectedPrimalUpdate(problem, workspace, current, product,
                                                                             step_size, candidate);
        if (!error.ok()) return error;
        error = copy(candidate, reflected);
        if (!error.ok()) return error;
        error = scale(2.0, reflected);
        if (!error.ok()) return error;
        return axpy(-1.0, current, reflected);
    }
    virtual api::Error fusedProjectedDualUpdateAndReflect(
            const ProblemHandle& problem, WorkspaceHandle& workspace,
            const VectorHandle& current, const VectorHandle& product,
            double step_size, VectorHandle& candidate, VectorHandle& reflected) {
        api::Error error = validateReflectionHandles(current, product, candidate,
                                                                                                 reflected);
        if (!error.ok()) return error;
        error = fusedProjectedDualUpdate(problem, workspace, current, product,
                                                                         step_size, candidate);
        if (!error.ok()) return error;
        error = copy(candidate, reflected);
        if (!error.ok()) return error;
        error = scale(2.0, reflected);
        if (!error.ok()) return error;
        return axpy(-1.0, current, reflected);
    }
  virtual api::Error fusedVectorBlend(
      WorkspaceHandle& workspace, const VectorHandle& primal_left,
      const VectorHandle& dual_left, const VectorHandle& primal_right,
      const VectorHandle& dual_right, double left_weight,
      VectorHandle& primal_output, VectorHandle& dual_output) = 0;
    virtual api::Error fusedProjectedPrimalUpdateAndReflect(
            const ProblemHandle& problem, WorkspaceHandle& workspace,
            const VectorHandle& current, const VectorHandle& product,
            const ScalarParameter& step, VectorHandle& candidate, VectorHandle& reflected) {
        if (step.vector != nullptr)
            return api::Error(api::ErrorCode::kUnsupportedCapability,
                                                "backend does not support vector-backed scalars");
        return fusedProjectedPrimalUpdateAndReflect(
                problem, workspace, current, product, step.value, candidate, reflected);
    }
    virtual api::Error fusedProjectedDualUpdateAndReflect(
            const ProblemHandle& problem, WorkspaceHandle& workspace,
            const VectorHandle& current, const VectorHandle& product,
            const ScalarParameter& step, VectorHandle& candidate, VectorHandle& reflected) {
        if (step.vector != nullptr)
            return api::Error(api::ErrorCode::kUnsupportedCapability,
                                                "backend does not support vector-backed scalars");
        return fusedProjectedDualUpdateAndReflect(
                problem, workspace, current, product, step.value, candidate, reflected);
    }
    virtual api::Error fusedVectorBlend(
            WorkspaceHandle& workspace, const VectorHandle& primal_left,
            const VectorHandle& dual_left, const VectorHandle& primal_right,
            const VectorHandle& dual_right, const ScalarParameter& weight,
            VectorHandle& primal_output, VectorHandle& dual_output) {
        if (weight.vector != nullptr)
            return api::Error(api::ErrorCode::kUnsupportedCapability,
                                                "backend does not support vector-backed scalars");
        return fusedVectorBlend(workspace, primal_left, dual_left, primal_right,
                                                        dual_right, weight.value, primal_output, dual_output);
    }
  virtual api::Error multiplyAndProjectedDualUpdateAndReflect(
      const ProblemHandle& problem, WorkspaceHandle& workspace,
      const VectorHandle& input, const VectorHandle& current,
      const ScalarParameter& step, VectorHandle& product,
      VectorHandle& candidate, VectorHandle& reflected) {
    api::Error error = multiply(problem, input, product);
    if (!error.ok()) return error;
    return step.vector != nullptr
        ? fusedProjectedDualUpdateAndReflect(problem, workspace, current,
            product, step, candidate, reflected)
        : fusedProjectedDualUpdateAndReflect(problem, workspace, current,
            product, step.value, candidate, reflected);
  }
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
    virtual api::Error reduceFixedPointMetrics(
            const VectorHandle& delta_primal, const VectorHandle& delta_dual,
            const VectorHandle& aty_delta, double& primal_norm, double& dual_norm,
            double& cross_term) {
        double primal = 0.0, dual = 0.0, cross = 0.0;
        api::Error error = norm2(delta_primal, primal);
        if (!error.ok()) return error;
        error = norm2(delta_dual, dual);
        if (!error.ok()) return error;
        error = dot(delta_primal, aty_delta, cross);
        if (!error.ok()) return error;
        primal_norm = primal;
        dual_norm = dual;
        cross_term = cross;
        return api::Error();
    }
  virtual api::Error reduceConvergenceMetrics(
      const VectorHandle& primal_residual, const VectorHandle& dual_residual,
      double primal_objective, double dual_objective,
      ConvergenceMetrics& metrics) = 0;
  virtual api::Error synchronize() = 0;
    virtual api::Error prepareVectorDownloads(const VectorDownload* downloads,
                                                                                        std::size_t count) {
        return validateVectorDownloads(downloads, count);
    }
    virtual api::Error downloadVectors(const VectorDownload* downloads,
                                                                         std::size_t count) {
        api::Error error = validateVectorDownloads(downloads, count);
        if (!error.ok() || count == 0) return error;
        try {
            std::vector<std::vector<double>> downloaded(count);
            for (std::size_t index = 0; index < count; ++index) {
                error = downloadVector(*downloads[index].vector, downloaded[index]);
                if (!error.ok()) return error;
            }
            DownloadOutputs outputs(downloads, count);
            for (std::size_t index = 0; index < count; ++index) {
                error = outputs.reserve(index, downloaded[index].size());
                if (!error.ok()) return error;
            }
            for (std::size_t index = 0; index < count; ++index)
                outputs.publish(index, downloaded[index].data(), downloaded[index].size());
        } catch (const std::bad_alloc&) {
            return api::Error(api::ErrorCode::kAllocationFailure,
                                                "allocating vector download outputs failed");
        } catch (const std::length_error&) {
            return api::Error(api::ErrorCode::kAllocationFailure,
                                                "vector download outputs exceed host capacity");
        }
        return api::Error();
    }
  virtual api::Error downloadVector(const VectorHandle& vector,
                                    std::vector<double>& values) = 0;
  virtual api::Error downloadSolution(const VectorHandle& primal,
                                      const VectorHandle& dual,
                                      std::vector<double>& primal_values,
                                      std::vector<double>& dual_values) = 0;

 protected:
    static api::Error validateReflectionHandles(
            const VectorHandle& current, const VectorHandle& product,
            const VectorHandle& candidate, const VectorHandle& reflected) {
        if (!current.valid() || !product.valid() || !candidate.valid() ||
                !reflected.valid() || &candidate == &reflected ||
                &candidate == &current || &candidate == &product ||
                &reflected == &current || &reflected == &product)
            return api::Error(api::ErrorCode::kInvalidArgument,
                                                "projection/reflection outputs must be valid and "
                                                "distinct from inputs");
        return api::Error();
    }
    static api::Error validateVectorDownloads(const VectorDownload* downloads,
                                                                                        const std::size_t count) {
        if (count != 0 && downloads == nullptr)
            return api::Error(api::ErrorCode::kInvalidArgument,
                                                "vector download array is null");
        for (std::size_t index = 0; index < count; ++index) {
            if (downloads[index].vector == nullptr ||
                    !downloads[index].vector->valid() || downloads[index].values == nullptr)
                return api::Error(api::ErrorCode::kInvalidArgument,
                                                    "vector download handle or output is invalid");
            for (std::size_t previous = 0; previous < index; ++previous) {
                if (downloads[index].values == downloads[previous].values)
                    return api::Error(api::ErrorCode::kInvalidArgument,
                                                        "vector download outputs must be distinct");
            }
        }
        return api::Error();
    }

    class DownloadOutputs {
     public:
        DownloadOutputs(const VectorDownload* downloads, const std::size_t count)
                : downloads_(downloads), count_(count) {}

        api::Error reserve(const std::size_t index, const std::size_t size) {
            const std::vector<double>& destination = *downloads_[index].values;
            if (destination.capacity() >= size) return api::Error();
            try {
                if (replacements_.empty()) replacements_.resize(count_);
                replacements_[index].reserve(size);
                replacements_[index].assign(destination.begin(), destination.end());
            } catch (const std::bad_alloc&) {
                return api::Error(api::ErrorCode::kAllocationFailure,
                                                    "allocating vector download output failed");
            } catch (const std::length_error&) {
                return api::Error(api::ErrorCode::kAllocationFailure,
                                                    "vector download output exceeds host capacity");
            }
            return api::Error();
        }

        void publishReservations() {
            for (std::size_t index = 0; index < replacements_.size(); ++index) {
                if (replacements_[index].capacity() != 0)
                    downloads_[index].values->swap(replacements_[index]);
            }
        }

        void publish(const std::size_t index, const double* values,
                                 const std::size_t size) {
            const bool replaced =
                    !replacements_.empty() && replacements_[index].capacity() != 0;
            std::vector<double>& destination =
                    replaced ? replacements_[index] : *downloads_[index].values;
            destination.resize(size);
            if (size != 0) std::copy_n(values, size, destination.begin());
            if (replaced) downloads_[index].values->swap(destination);
        }

     private:
        const VectorDownload* downloads_;
        std::size_t count_;
        std::vector<std::vector<double>> replacements_;
    };

 private:
  std::chrono::steady_clock::time_point timing_start_;
  bool timing_prepared_ = false;
  bool timing_active_ = false;
};

}  // namespace backend
}  // namespace hipdlp

#endif