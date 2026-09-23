#include "pdlp/hipdlp/backend/cpu/cpu_backend.hpp"

#include <algorithm>
#include <cmath>
#include <new>
#include <string>

#include "pdlp/hipdlp/backend/cpu/sparse_kernels.hpp"
#include "pdlp/hipdlp/backend/cpu/vector_kernels.hpp"

namespace hipdlp {
namespace backend {
namespace cpu {

namespace {

const unsigned long long kCpuCapabilityMask =
    static_cast<unsigned long long>(Capability::kSparseMatrixVector) |
    static_cast<unsigned long long>(Capability::kFusedPrimalDualUpdate) |
    static_cast<unsigned long long>(Capability::kFusedAverageUpdate) |
    static_cast<unsigned long long>(Capability::kConvergenceReduction) |
    static_cast<unsigned long long>(Capability::kVectorOperations) |
    static_cast<unsigned long long>(Capability::kReductions) |
    static_cast<unsigned long long>(Capability::kHostDownload);

api::Error invalidArgument(const std::string& message) {
  return api::Error(api::ErrorCode::kInvalidArgument, message);
}

api::Error backendFailure(const std::string& message) {
  return api::Error(api::ErrorCode::kBackendFailure, message);
}

api::Error allocationFailure(const std::string& operation) {
  return api::Error(api::ErrorCode::kAllocationFailure,
                    "serial CPU allocation failed while " + operation);
}

api::Error unsupported(const char* capability) {
  return api::Error(
      api::ErrorCode::kUnsupportedCapability,
      std::string("serial CPU backend does not support ") + capability);
}

bool finiteVector(const std::vector<double>& values) {
  for (std::size_t index = 0; index < values.size(); ++index) {
    if (!std::isfinite(values[index])) return false;
  }
  return true;
}

bool finiteScaleResult(const std::vector<double>& values, const double factor) {
  for (std::size_t index = 0; index < values.size(); ++index) {
    if (!std::isfinite(values[index] * factor)) return false;
  }
  return true;
}

bool finiteAxpyResult(const std::vector<double>& input,
                      const std::vector<double>& output, const double factor) {
  for (std::size_t index = 0; index < input.size(); ++index) {
    if (!std::isfinite(output[index] + factor * input[index])) return false;
  }
  return true;
}

}  // namespace

struct CpuProblemResource : Resource {
  const CpuBackend* owner;
  HighsInt num_rows;
  HighsInt num_cols;
  model::IndexWidth index_width;
  std::vector<double> costs;
  std::vector<double> col_lower;
  std::vector<double> col_upper;
  std::vector<double> row_lower;
  std::vector<HighsBool> equality_rows;
  model::SparseMatrix matrix;
  model::SparseMatrix transpose;

  CpuProblemResource(const CpuBackend* backend, const model::Problem& problem)
      : owner(backend),
        num_rows(problem.numRows()),
        num_cols(problem.numCols()),
        index_width(problem.indexWidth()),
        costs(problem.costs()),
        col_lower(problem.colLower()),
        col_upper(problem.colUpper()),
        row_lower(problem.rowLower()),
        equality_rows(problem.equalityRows()),
        matrix(problem.matrix()),
        transpose(problem.transpose()) {}
};

struct CpuVectorResource : Resource {
  const CpuBackend* owner;
  std::vector<double> values;

  CpuVectorResource(const CpuBackend* backend, const HighsInt size)
      : owner(backend), values(size, 0.0) {}
};

struct CpuWorkspaceResource : Resource {
  const CpuBackend* owner;
  const model::Problem* check_problem;
  HighsInt num_rows;
  HighsInt num_cols;
  std::vector<double> row_scratch;
  std::vector<double> col_scratch;

  CpuWorkspaceResource(const CpuBackend* backend, const HighsInt rows,
                       const HighsInt cols)
      : owner(backend),
        check_problem(nullptr),
        num_rows(rows),
        num_cols(cols),
        row_scratch(rows, 0.0),
        col_scratch(cols, 0.0) {}
};

CpuProblemResource* CpuBackend::problemResource(ProblemHandle& handle) const {
  CpuProblemResource* resource =
      dynamic_cast<CpuProblemResource*>(handle.resource_.get());
  return resource != nullptr && resource->owner == this ? resource : nullptr;
}

const CpuProblemResource* CpuBackend::problemResource(
    const ProblemHandle& handle) const {
  const CpuProblemResource* resource =
      dynamic_cast<const CpuProblemResource*>(handle.resource_.get());
  return resource != nullptr && resource->owner == this ? resource : nullptr;
}

CpuVectorResource* CpuBackend::vectorResource(VectorHandle& handle) const {
  CpuVectorResource* resource =
      dynamic_cast<CpuVectorResource*>(handle.resource_.get());
  return resource != nullptr && resource->owner == this ? resource : nullptr;
}

const CpuVectorResource* CpuBackend::vectorResource(
    const VectorHandle& handle) const {
  const CpuVectorResource* resource =
      dynamic_cast<const CpuVectorResource*>(handle.resource_.get());
  return resource != nullptr && resource->owner == this ? resource : nullptr;
}

CpuWorkspaceResource* CpuBackend::workspaceResource(
    WorkspaceHandle& handle) const {
  CpuWorkspaceResource* resource =
      dynamic_cast<CpuWorkspaceResource*>(handle.resource_.get());
  return resource != nullptr && resource->owner == this ? resource : nullptr;
}

Descriptor CpuBackend::descriptor() const {
  return {Kind::kCpu,
          "cpu",
          true,
          true,
          PrecisionDescriptor(),
          Capabilities(kCpuCapabilityMask)};
}

api::Error CpuBackend::queryDevice(DeviceInfo& device) const {
  device.available = true;
  device.name = "deterministic serial host CPU";
  return api::Error();
}

api::Error CpuBackend::requireCapability(const Capability capability) const {
  if (!descriptor().capabilities.supports(capability))
    return unsupported(capability == Capability::kBatchedVectors
                           ? "batched vectors"
                       : capability == Capability::kSparseMatrixMatrix
                           ? "SpMM"
                           : "the requested capability");
  return api::Error();
}

api::Error CpuBackend::bindProblem(const model::Problem& problem,
                                   ProblemHandle& handle) {
  if (!problem.initialized())
    return invalidArgument("cannot bind an uninitialized HiPDLP problem");
  try {
    handle = ProblemHandle(
      std::unique_ptr<Resource>(new CpuProblemResource(this, problem)));
  } catch (const std::bad_alloc&) {
    return allocationFailure("binding a problem");
  }
  return api::Error();
}

api::Error CpuBackend::createWorkspace(const ProblemHandle& problem,
                                       WorkspaceHandle& workspace) {
  const CpuProblemResource* cpu_problem = problemResource(problem);
  if (cpu_problem == nullptr)
    return invalidArgument("serial CPU problem handle is invalid");
  try {
    workspace =
        WorkspaceHandle(std::unique_ptr<Resource>(new CpuWorkspaceResource(
          this, cpu_problem->num_rows, cpu_problem->num_cols)));
  } catch (const std::bad_alloc&) {
    return allocationFailure("creating a workspace");
  }
  return api::Error();
}

api::Error CpuBackend::prepareCheckWorkspace(const model::Problem& problem,
                                            const HighsLp& original_lp,
                                            WorkspaceHandle& workspace) {
  (void)original_lp;
  CpuWorkspaceResource* cpu_workspace = workspaceResource(workspace);
  if (cpu_workspace == nullptr)
    return invalidArgument("serial CPU check workspace is invalid");
  cpu_workspace->check_problem = nullptr;
  if (!problem.initialized() || cpu_workspace->num_rows != problem.numRows() ||
      cpu_workspace->num_cols != problem.numCols())
    return invalidArgument("serial CPU check problem dimensions are invalid");
  cpu_workspace->check_problem = &problem;
  return api::Error();
}

api::Error CpuBackend::evaluateCheck(const model::Problem& problem,
                                    WorkspaceHandle& workspace,
                                    const CheckInput& input,
                                    CheckResult& result) {
  result = CheckResult();
  if (!input.validTolerances())
    return invalidArgument("HiPDLP check tolerances must be finite and nonnegative");
  const CpuWorkspaceResource* cpu_workspace = workspaceResource(workspace);
  if (cpu_workspace == nullptr || cpu_workspace->check_problem != &problem ||
      !problem.initialized() || cpu_workspace->num_rows != problem.numRows() ||
      cpu_workspace->num_cols != problem.numCols())
    return invalidArgument("serial CPU check workspace is not prepared");

  const VectorHandle* handles[] = {
      input.primal, input.dual, input.ax, input.aty, input.delta_primal,
      input.delta_dual, input.aty_delta, input.primal_anchor, input.dual_anchor};
  const HighsInt sizes[] = {
      problem.numCols(), problem.numRows(), problem.numRows(),
      problem.numCols(), problem.numCols(), problem.numRows(),
      problem.numCols(), problem.numCols(), problem.numRows()};
  const CpuVectorResource* vectors[9] = {};
  for (std::size_t index = 0; index < 9; ++index) {
    const bool required = index < 4 ||
                          (index < 7 ? input.fixed_point : input.anchor_distances);
    if (!required) continue;
    vectors[index] = handles[index] == nullptr ? nullptr : vectorResource(*handles[index]);
    if (vectors[index] == nullptr || vectors[index]->values.size() !=
                                         static_cast<std::size_t>(sizes[index]))
      return invalidArgument("serial CPU check vector dimensions or owner are invalid");
    if (!finiteVector(vectors[index]->values))
      return backendFailure("serial CPU check input contains a non-finite value");
  }

  CheckResult computed;
  api::Error error = algorithm::computeRawConvergenceMetrics(
      problem, vectors[0]->values, vectors[1]->values, vectors[2]->values,
      vectors[3]->values, computed.raw);
  if (!error.ok()) return error;
  if (!std::isfinite(computed.raw.unscaled_rhs_norm) ||
      !std::isfinite(computed.raw.unscaled_cost_norm))
    return backendFailure("serial CPU check model norms are non-finite");
  if (input.fixed_point) {
    computed.primal_delta_norm = norm2Vector(vectors[4]->values);
    computed.dual_delta_norm = norm2Vector(vectors[5]->values);
    computed.delta_dot = dotVector(vectors[4]->values, vectors[6]->values);
  }
  if (input.anchor_distances) {
    computed.primal_anchor_distance =
        differenceNorm2Vector(vectors[0]->values, vectors[7]->values);
    computed.dual_anchor_distance =
        differenceNorm2Vector(vectors[1]->values, vectors[8]->values);
  }
  if (!std::isfinite(computed.primal_delta_norm) ||
      !std::isfinite(computed.dual_delta_norm) ||
      !std::isfinite(computed.delta_dot) ||
      !std::isfinite(computed.primal_anchor_distance) ||
      !std::isfinite(computed.dual_anchor_distance))
    return backendFailure("serial CPU check reduction produced a non-finite value");
  computed.scalar_count = 7 + (input.fixed_point ? 3 : 0) +
                         (input.anchor_distances ? 2 : 0);
  computed.valid = true;
  computed.fixed_point_valid = input.fixed_point;
  computed.fixed_point_prepared = input.fixed_point;
  computed.anchor_distances_valid = input.anchor_distances;
  result = computed;
  return api::Error();
}

api::Error CpuBackend::createVector(const HighsInt size, VectorHandle& vector) {
  if (size < 0) return invalidArgument("vector size must be nonnegative");
  try {
    vector =
      VectorHandle(std::unique_ptr<Resource>(new CpuVectorResource(this, size)));
  } catch (const std::bad_alloc&) {
    return allocationFailure("creating a vector");
  }
  return api::Error();
}

api::Error CpuBackend::createBatchedVector(const HighsInt size,
                                           const HighsInt batch_size,
                                           VectorHandle& vector) {
  (void)size;
  (void)batch_size;
  (void)vector;
  return unsupported("batched vectors");
}

api::Error CpuBackend::uploadVector(const std::vector<double>& values,
                                    VectorHandle& vector) {
  CpuVectorResource* cpu_vector = vectorResource(vector);
  if (cpu_vector == nullptr || cpu_vector->values.size() != values.size())
    return invalidArgument("serial CPU upload vector dimension mismatch");
  if (!finiteVector(values))
    return invalidArgument("serial CPU upload contains a non-finite value");
  copyVector(values, cpu_vector->values);
  return api::Error();
}

api::Error CpuBackend::initialize(VectorHandle& primal, VectorHandle& dual) {
  api::Error error = fill(primal, 0.0);
  if (!error.ok()) return error;
  return fill(dual, 0.0);
}

api::Error CpuBackend::fill(VectorHandle& vector, const double value) {
  CpuVectorResource* cpu_vector = vectorResource(vector);
  if (cpu_vector == nullptr)
    return invalidArgument("serial CPU vector handle is invalid");
  if (!std::isfinite(value))
    return invalidArgument("serial CPU fill value must be finite");
  fillVector(cpu_vector->values, value);
  return api::Error();
}

api::Error CpuBackend::copy(const VectorHandle& source,
                            VectorHandle& destination) {
  const CpuVectorResource* cpu_source = vectorResource(source);
  CpuVectorResource* cpu_destination = vectorResource(destination);
  if (cpu_source == nullptr || cpu_destination == nullptr ||
      cpu_source->values.size() != cpu_destination->values.size())
    return invalidArgument("serial CPU copy vector dimension mismatch");
  copyVector(cpu_source->values, cpu_destination->values);
  return api::Error();
}

api::Error CpuBackend::scale(const double factor, VectorHandle& vector) {
  CpuVectorResource* cpu_vector = vectorResource(vector);
  if (cpu_vector == nullptr || !std::isfinite(factor))
    return invalidArgument("serial CPU scale arguments are invalid");
  if (!finiteScaleResult(cpu_vector->values, factor))
    return backendFailure("serial CPU scale produced a non-finite value");
  scaleVector(factor, cpu_vector->values);
  return api::Error();
}

api::Error CpuBackend::axpy(const double factor, const VectorHandle& input,
                            VectorHandle& output) {
  const CpuVectorResource* cpu_input = vectorResource(input);
  CpuVectorResource* cpu_output = vectorResource(output);
  if (cpu_input == nullptr || cpu_output == nullptr ||
      cpu_input->values.size() != cpu_output->values.size() ||
      !std::isfinite(factor))
    return invalidArgument("serial CPU axpy arguments are invalid");
  if (!finiteAxpyResult(cpu_input->values, cpu_output->values, factor))
    return backendFailure("serial CPU axpy produced a non-finite value");
  axpyVector(factor, cpu_input->values, cpu_output->values);
  return api::Error();
}

api::Error CpuBackend::multiply(const ProblemHandle& problem,
                                const VectorHandle& input,
                                VectorHandle& output) {
  const CpuProblemResource* cpu_problem = problemResource(problem);
  const CpuVectorResource* cpu_input = vectorResource(input);
  CpuVectorResource* cpu_output = vectorResource(output);
  if (cpu_problem == nullptr || cpu_input == nullptr || cpu_output == nullptr ||
      cpu_input->values.size() !=
          static_cast<std::size_t>(cpu_problem->num_cols) ||
      cpu_output->values.size() !=
          static_cast<std::size_t>(cpu_problem->num_rows))
    return invalidArgument("serial CPU Ax dimensions are inconsistent");
  if (cpu_input == cpu_output)
    return invalidArgument("serial CPU Ax does not support aliased vectors");
  multiplySparse(cpu_problem->matrix, cpu_input->values, cpu_output->values);
  if (!finiteVector(cpu_output->values))
    return backendFailure("serial CPU Ax produced a non-finite value");
  return api::Error();
}

api::Error CpuBackend::multiplyTranspose(const ProblemHandle& problem,
                                         const VectorHandle& input,
                                         VectorHandle& output) {
  const CpuProblemResource* cpu_problem = problemResource(problem);
  const CpuVectorResource* cpu_input = vectorResource(input);
  CpuVectorResource* cpu_output = vectorResource(output);
  if (cpu_problem == nullptr || cpu_input == nullptr || cpu_output == nullptr ||
      cpu_input->values.size() !=
          static_cast<std::size_t>(cpu_problem->num_rows) ||
      cpu_output->values.size() !=
          static_cast<std::size_t>(cpu_problem->num_cols))
    return invalidArgument("serial CPU A^T y dimensions are inconsistent");
  if (cpu_input == cpu_output)
    return invalidArgument("serial CPU A^T y does not support aliased vectors");
  multiplySparse(cpu_problem->transpose, cpu_input->values, cpu_output->values);
  if (!finiteVector(cpu_output->values))
    return backendFailure("serial CPU A^T y produced a non-finite value");
  return api::Error();
}

api::Error CpuBackend::multiplyBatched(const ProblemHandle& problem,
                                       const VectorHandle& input,
                                       VectorHandle& output) {
  (void)problem;
  (void)input;
  (void)output;
  return unsupported("SpMM");
}

api::Error CpuBackend::fusedProjectedPrimalUpdate(
    const ProblemHandle& problem, WorkspaceHandle& workspace,
    const VectorHandle& current_primal, const VectorHandle& transpose_product,
    const double step_size, VectorHandle& next_primal) {
  const CpuProblemResource* cpu_problem = problemResource(problem);
  CpuWorkspaceResource* cpu_workspace = workspaceResource(workspace);
  const CpuVectorResource* current = vectorResource(current_primal);
  const CpuVectorResource* product = vectorResource(transpose_product);
  CpuVectorResource* next = vectorResource(next_primal);
  if (cpu_problem == nullptr || cpu_workspace == nullptr ||
      current == nullptr || product == nullptr || next == nullptr ||
      !std::isfinite(step_size) || step_size < 0.0 ||
      cpu_workspace->num_cols != cpu_problem->num_cols ||
      current->values.size() !=
          static_cast<std::size_t>(cpu_problem->num_cols) ||
      product->values.size() != current->values.size() ||
      next->values.size() != current->values.size())
    return invalidArgument(
        "serial CPU fused primal update arguments are "
        "invalid");
  projectedPrimalUpdate(current->values, product->values, cpu_problem->costs,
                        cpu_problem->col_lower, cpu_problem->col_upper,
                        step_size, next->values);
  if (!finiteVector(next->values))
    return backendFailure(
        "serial CPU fused primal update produced a non-finite value");
  return api::Error();
}

api::Error CpuBackend::fusedProjectedDualUpdate(
    const ProblemHandle& problem, WorkspaceHandle& workspace,
    const VectorHandle& current_dual, const VectorHandle& extrapolated_product,
    const double step_size, VectorHandle& next_dual) {
  const CpuProblemResource* cpu_problem = problemResource(problem);
  CpuWorkspaceResource* cpu_workspace = workspaceResource(workspace);
  const CpuVectorResource* current = vectorResource(current_dual);
  const CpuVectorResource* product = vectorResource(extrapolated_product);
  CpuVectorResource* next = vectorResource(next_dual);
  if (cpu_problem == nullptr || cpu_workspace == nullptr ||
      current == nullptr || product == nullptr || next == nullptr ||
      !std::isfinite(step_size) || step_size < 0.0 ||
      cpu_workspace->num_rows != cpu_problem->num_rows ||
      current->values.size() !=
          static_cast<std::size_t>(cpu_problem->num_rows) ||
      product->values.size() != current->values.size() ||
      next->values.size() != current->values.size())
    return invalidArgument(
        "serial CPU fused dual update arguments are "
        "invalid");
  projectedDualUpdate(current->values, product->values, cpu_problem->row_lower,
                      cpu_problem->equality_rows, step_size, next->values);
  if (!finiteVector(next->values))
    return backendFailure(
        "serial CPU fused dual update produced a non-finite value");
  return api::Error();
}

api::Error CpuBackend::fusedProjectedPrimalUpdateAndReflect(
    const ProblemHandle& problem, WorkspaceHandle& workspace,
    const VectorHandle& current, const VectorHandle& product,
    const double step_size, VectorHandle& candidate, VectorHandle& reflected) {
  const api::Error handles = validateReflectionHandles(current, product,
                                                       candidate, reflected);
  if (!handles.ok()) return handles;
  const CpuProblemResource* cpu_problem = problemResource(problem);
  const CpuWorkspaceResource* cpu_workspace = workspaceResource(workspace);
  const CpuVectorResource* cpu_current = vectorResource(current);
  const CpuVectorResource* cpu_product = vectorResource(product);
  CpuVectorResource* cpu_candidate = vectorResource(candidate);
  CpuVectorResource* cpu_reflected = vectorResource(reflected);
  if (!cpu_problem || !cpu_workspace || !cpu_current || !cpu_product ||
      !cpu_candidate || !cpu_reflected || !std::isfinite(step_size) ||
      step_size < 0.0 || cpu_workspace->num_cols != cpu_problem->num_cols ||
      cpu_current->values.size() != static_cast<std::size_t>(cpu_problem->num_cols) ||
      cpu_product->values.size() != cpu_current->values.size() ||
      cpu_candidate->values.size() != cpu_current->values.size() ||
      cpu_reflected->values.size() != cpu_current->values.size())
    return invalidArgument("serial CPU primal projection/reflection arguments are invalid");
  for (std::size_t col = 0; col < cpu_current->values.size(); ++col) {
    const double projected = std::max(cpu_problem->col_lower[col],
        std::min(cpu_current->values[col] - step_size *
            (cpu_problem->costs[col] - cpu_product->values[col]),
            cpu_problem->col_upper[col]));
    const double doubled = 2.0 * projected;
    const double reflection = doubled - cpu_current->values[col];
    if (!std::isfinite(projected) || !std::isfinite(doubled) ||
        !std::isfinite(reflection))
      return backendFailure("serial CPU primal projection/reflection produced a non-finite value");
    cpu_candidate->values[col] = projected;
    cpu_reflected->values[col] = reflection;
  }
  return api::Error();
}

api::Error CpuBackend::fusedProjectedDualUpdateAndReflect(
    const ProblemHandle& problem, WorkspaceHandle& workspace,
    const VectorHandle& current, const VectorHandle& product,
    const double step_size, VectorHandle& candidate, VectorHandle& reflected) {
  const api::Error handles = validateReflectionHandles(current, product,
                                                       candidate, reflected);
  if (!handles.ok()) return handles;
  const CpuProblemResource* cpu_problem = problemResource(problem);
  const CpuWorkspaceResource* cpu_workspace = workspaceResource(workspace);
  const CpuVectorResource* cpu_current = vectorResource(current);
  const CpuVectorResource* cpu_product = vectorResource(product);
  CpuVectorResource* cpu_candidate = vectorResource(candidate);
  CpuVectorResource* cpu_reflected = vectorResource(reflected);
  if (!cpu_problem || !cpu_workspace || !cpu_current || !cpu_product ||
      !cpu_candidate || !cpu_reflected || !std::isfinite(step_size) ||
      step_size < 0.0 || cpu_workspace->num_rows != cpu_problem->num_rows ||
      cpu_current->values.size() != static_cast<std::size_t>(cpu_problem->num_rows) ||
      cpu_product->values.size() != cpu_current->values.size() ||
      cpu_candidate->values.size() != cpu_current->values.size() ||
      cpu_reflected->values.size() != cpu_current->values.size())
    return invalidArgument("serial CPU dual projection/reflection arguments are invalid");
  for (std::size_t row = 0; row < cpu_current->values.size(); ++row) {
    const double value = cpu_current->values[row] + step_size *
        (cpu_problem->row_lower[row] - cpu_product->values[row]);
    const double projected = cpu_problem->equality_rows[row]
        ? value : std::max(0.0, value);
    const double doubled = 2.0 * projected;
    const double reflection = doubled - cpu_current->values[row];
    if (!std::isfinite(projected) || !std::isfinite(doubled) ||
        !std::isfinite(reflection))
      return backendFailure("serial CPU dual projection/reflection produced a non-finite value");
    cpu_candidate->values[row] = projected;
    cpu_reflected->values[row] = reflection;
  }
  return api::Error();
}

api::Error CpuBackend::fusedVectorBlend(
    WorkspaceHandle& workspace, const VectorHandle& primal_left,
    const VectorHandle& dual_left, const VectorHandle& primal_right,
    const VectorHandle& dual_right, const double left_weight,
    VectorHandle& primal_output, VectorHandle& dual_output) {
  CpuWorkspaceResource* cpu_workspace = workspaceResource(workspace);
  const CpuVectorResource* p_left = vectorResource(primal_left);
  const CpuVectorResource* d_left = vectorResource(dual_left);
  const CpuVectorResource* p_right = vectorResource(primal_right);
  const CpuVectorResource* d_right = vectorResource(dual_right);
  CpuVectorResource* p_output = vectorResource(primal_output);
  CpuVectorResource* d_output = vectorResource(dual_output);
  if (cpu_workspace == nullptr || p_left == nullptr || d_left == nullptr ||
      p_right == nullptr || d_right == nullptr || p_output == nullptr ||
      d_output == nullptr || !std::isfinite(left_weight) || left_weight < 0.0 ||
      left_weight > 1.0 || p_left->values.size() != p_right->values.size() ||
      p_left->values.size() != p_output->values.size() ||
      d_left->values.size() != d_right->values.size() ||
      d_left->values.size() != d_output->values.size() ||
      p_left->values.size() !=
          static_cast<std::size_t>(cpu_workspace->num_cols) ||
      d_left->values.size() !=
          static_cast<std::size_t>(cpu_workspace->num_rows))
    return invalidArgument(
        "serial CPU fused vector blend arguments are "
        "invalid");
  blendVector(p_left->values, p_right->values, left_weight, p_output->values);
  blendVector(d_left->values, d_right->values, left_weight, d_output->values);
  if (!finiteVector(p_output->values) || !finiteVector(d_output->values))
    return backendFailure(
        "serial CPU fused vector blend produced a non-finite value");
  return api::Error();
}

api::Error CpuBackend::fusedAverageUpdate(WorkspaceHandle& workspace,
                                          const VectorHandle& primal,
                                          const VectorHandle& dual,
                                          const double weight,
                                          VectorHandle& primal_sum,
                                          VectorHandle& dual_sum) {
  CpuWorkspaceResource* cpu_workspace = workspaceResource(workspace);
  const CpuVectorResource* p = vectorResource(primal);
  const CpuVectorResource* d = vectorResource(dual);
  CpuVectorResource* p_sum = vectorResource(primal_sum);
  CpuVectorResource* d_sum = vectorResource(dual_sum);
  if (cpu_workspace == nullptr || p == nullptr || d == nullptr ||
      p_sum == nullptr || d_sum == nullptr || !std::isfinite(weight) ||
      p->values.size() != p_sum->values.size() ||
      d->values.size() != d_sum->values.size() ||
      p->values.size() != static_cast<std::size_t>(cpu_workspace->num_cols) ||
      d->values.size() != static_cast<std::size_t>(cpu_workspace->num_rows) ||
      !finiteAxpyResult(p->values, p_sum->values, weight) ||
      !finiteAxpyResult(d->values, d_sum->values, weight))
    return invalidArgument(
        "serial CPU fused average update arguments are "
        "invalid");
  averageUpdate(p->values, weight, p_sum->values);
  averageUpdate(d->values, weight, d_sum->values);
  return api::Error();
}

api::Error CpuBackend::dot(const VectorHandle& left, const VectorHandle& right,
                           double& result) {
  const CpuVectorResource* cpu_left = vectorResource(left);
  const CpuVectorResource* cpu_right = vectorResource(right);
  if (cpu_left == nullptr || cpu_right == nullptr ||
      cpu_left->values.size() != cpu_right->values.size())
    return invalidArgument("serial CPU dot dimensions are inconsistent");
  const double reduced = dotVector(cpu_left->values, cpu_right->values);
  if (!std::isfinite(reduced))
    return backendFailure("serial CPU dot produced a non-finite value");
  result = reduced;
  return api::Error();
}

api::Error CpuBackend::norm2(const VectorHandle& vector, double& result) {
  const CpuVectorResource* cpu_vector = vectorResource(vector);
  if (cpu_vector == nullptr)
    return invalidArgument("serial CPU norm vector is invalid");
  const double reduced = norm2Vector(cpu_vector->values);
  if (!std::isfinite(reduced))
    return backendFailure("serial CPU norm produced a non-finite value");
  result = reduced;
  return api::Error();
}

api::Error CpuBackend::differenceNorm2(const VectorHandle& left,
                                       const VectorHandle& right,
                                       double& result) {
  const CpuVectorResource* cpu_left = vectorResource(left);
  const CpuVectorResource* cpu_right = vectorResource(right);
  if (cpu_left == nullptr || cpu_right == nullptr ||
      cpu_left->values.size() != cpu_right->values.size())
    return invalidArgument(
        "serial CPU difference norm dimensions are inconsistent");
  const double reduced =
      differenceNorm2Vector(cpu_left->values, cpu_right->values);
  if (!std::isfinite(reduced))
    return backendFailure(
        "serial CPU difference norm produced a non-finite value");
  result = reduced;
  return api::Error();
}

api::Error CpuBackend::reduceConvergenceMetrics(
    const VectorHandle& primal_residual, const VectorHandle& dual_residual,
    const double primal_objective, const double dual_objective,
    ConvergenceMetrics& metrics) {
  if (!std::isfinite(primal_objective) || !std::isfinite(dual_objective))
    return invalidArgument("convergence objectives must be finite");
  double primal_norm = 0.0;
  api::Error error = norm2(primal_residual, primal_norm);
  if (!error.ok()) return error;
  double dual_norm = 0.0;
  error = norm2(dual_residual, dual_norm);
  if (!error.ok()) return error;
  const double gap = std::abs(primal_objective - dual_objective);
  if (!std::isfinite(gap))
    return backendFailure("convergence gap is non-finite");
  metrics.primal_residual = primal_norm;
  metrics.dual_residual = dual_norm;
  metrics.objective_gap = gap;
  return api::Error();
}

api::Error CpuBackend::synchronize() { return api::Error(); }

api::Error CpuBackend::prepareVectorDownloads(const VectorDownload* downloads,
                                              const std::size_t count) {
  api::Error error = validateVectorDownloads(downloads, count);
  if (!error.ok()) return error;
  for (std::size_t index = 0; index < count; ++index) {
    if (vectorResource(*downloads[index].vector) == nullptr)
      return invalidArgument("serial CPU download vector is invalid");
  }
  DownloadOutputs outputs(downloads, count);
  for (std::size_t index = 0; index < count; ++index) {
    error = outputs.reserve(
        index, vectorResource(*downloads[index].vector)->values.size());
    if (!error.ok()) return error;
  }
  outputs.publishReservations();
  return api::Error();
}

api::Error CpuBackend::downloadVectors(const VectorDownload* downloads,
                                       const std::size_t count) {
  const api::Error error = prepareVectorDownloads(downloads, count);
  if (!error.ok()) return error;
  for (std::size_t index = 0; index < count; ++index) {
    const std::vector<double>& source =
        vectorResource(*downloads[index].vector)->values;
    downloads[index].values->assign(source.begin(), source.end());
  }
  return api::Error();
}

api::Error CpuBackend::downloadVector(const VectorHandle& vector,
                                      std::vector<double>& values) {
  const VectorDownload download = {&vector, &values};
  return downloadVectors(&download, 1);
}

api::Error CpuBackend::downloadSolution(const VectorHandle& primal,
                                        const VectorHandle& dual,
                                        std::vector<double>& primal_values,
                                        std::vector<double>& dual_values) {
  const VectorDownload downloads[] = {{&primal, &primal_values},
                                       {&dual, &dual_values}};
  return downloadVectors(downloads, 2);
}

std::unique_ptr<Backend> createCpuBackend() {
  return std::unique_ptr<Backend>(new CpuBackend());
}

}  // namespace cpu
}  // namespace backend
}  // namespace hipdlp