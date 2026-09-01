#include <sys/stat.h>

#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "HCheckConfig.h"
#include "HConfig.h"
#include "catch.hpp"
#include "lp_data/HConst.h"
#include "pdlp/hipdlp/backend/backend_factory.hpp"
#include "pdlp/hipdlp/backend/command_buffer.hpp"
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
#include "pdlp/hipdlp/backend/cuda/cuda_backend_testing.hpp"
#endif
#include "pdlp/hipdlp/model/problem.hpp"

namespace {

HighsLp makeEqualityLp(const HighsInt rows, const HighsInt cols,
                       const std::vector<HighsInt>& starts,
                       const std::vector<HighsInt>& indices,
                       const std::vector<double>& values) {
  HighsLp lp;
  lp.num_row_ = rows;
  lp.num_col_ = cols;
  lp.col_cost_.assign(cols, 0.0);
  lp.col_lower_.assign(cols, -kHighsInf);
  lp.col_upper_.assign(cols, kHighsInf);
  lp.row_lower_.assign(rows, 0.0);
  lp.row_upper_.assign(rows, 0.0);
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_row_ = rows;
  lp.a_matrix_.num_col_ = cols;
  lp.a_matrix_.start_ = starts;
  lp.a_matrix_.index_ = indices;
  lp.a_matrix_.value_ = values;
  lp.sense_ = ObjSense::kMinimize;
  lp.offset_ = 0.0;
  return lp;
}

std::unique_ptr<hipdlp::backend::Backend> makeCpuBackend() {
  hipdlp::backend::FactoryResult result =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kCpu);
  INFO(result.error.message);
  REQUIRE(result.ok());
  return std::move(result.backend);
}

std::vector<double> download(hipdlp::backend::Backend& backend,
                             const hipdlp::backend::VectorHandle& vector) {
  std::vector<double> values;
  REQUIRE(backend.downloadVector(vector, values).ok());
  return values;
}

std::string readSourceFile(const char* relative_path) {
  const std::string path = std::string(HIGHS_DIR) + "/" + relative_path;
  std::ifstream stream(path.c_str(), std::ios::in | std::ios::binary);
  INFO(path);
  REQUIRE(stream.good());
  std::ostringstream contents;
  contents << stream.rdbuf();
  return contents.str();
}

bool sourcePathExists(const char* relative_path) {
  const std::string source_path = std::string(HIGHS_DIR) + "/" + relative_path;
#ifdef _WIN32
  struct _stat path_status;
  return _stat(source_path.c_str(), &path_status) == 0;
#else
  struct stat path_status;
  return stat(source_path.c_str(), &path_status) == 0;
#endif
}

std::size_t countText(const std::string& contents, const std::string& text) {
  std::size_t count = 0;
  std::size_t position = 0;
  while ((position = contents.find(text, position)) != std::string::npos) {
    ++count;
    position += text.size();
  }
  return count;
}

void checkApproxVector(const std::vector<double>& actual,
                       const std::vector<double>& expected,
                       const double epsilon = 1e-14) {
  REQUIRE(actual.size() == expected.size());
  for (std::size_t index = 0; index < actual.size(); ++index)
    CHECK(actual[index] == Approx(expected[index]).epsilon(epsilon));
}

void checkProducts(const HighsLp& lp, const std::vector<double>& x,
                   const std::vector<double>& expected_ax,
                   const std::vector<double>& y,
                   const std::vector<double>& expected_aty) {
  hipdlp::model::ProblemBuildResult built = hipdlp::model::buildProblem(lp);
  INFO(built.error.message);
  REQUIRE(built.ok());
  std::unique_ptr<hipdlp::backend::Backend> backend = makeCpuBackend();
  hipdlp::backend::ProblemHandle problem;
  REQUIRE(backend->bindProblem(built.problem, problem).ok());
  hipdlp::backend::VectorHandle x_handle;
  hipdlp::backend::VectorHandle ax_handle;
  hipdlp::backend::VectorHandle y_handle;
  hipdlp::backend::VectorHandle aty_handle;
  REQUIRE(backend->createVector(lp.num_col_, x_handle).ok());
  REQUIRE(backend->createVector(lp.num_row_, ax_handle).ok());
  REQUIRE(backend->createVector(lp.num_row_, y_handle).ok());
  REQUIRE(backend->createVector(lp.num_col_, aty_handle).ok());
  REQUIRE(backend->uploadVector(x, x_handle).ok());
  REQUIRE(backend->uploadVector(y, y_handle).ok());
  REQUIRE(backend->multiply(problem, x_handle, ax_handle).ok());
  REQUIRE(backend->multiplyTranspose(problem, y_handle, aty_handle).ok());
  checkApproxVector(download(*backend, ax_handle), expected_ax, 1e-13);
  checkApproxVector(download(*backend, aty_handle), expected_aty, 1e-13);
}

struct TrackingResource : hipdlp::backend::Resource {
  explicit TrackingResource(int& live) : live_(live) { ++live_; }
  ~TrackingResource() override { --live_; }
  int& live_;
};

#ifdef HIGHS_HIPDLP_CUDA_DETECTED
std::unique_ptr<hipdlp::backend::Backend> makeCudaBackend() {
  hipdlp::backend::FactoryResult result =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda);
  INFO(result.error.message);
  REQUIRE(result.ok());
  return std::move(result.backend);
}

struct ProductResult {
  std::vector<double> ax;
  std::vector<double> aty;
};

ProductResult runProducts(const hipdlp::backend::Kind kind, const HighsLp& lp,
                          const std::vector<double>& x,
                          const std::vector<double>& y) {
  hipdlp::model::ProblemBuildResult built = hipdlp::model::buildProblem(lp);
  INFO(built.error.message);
  REQUIRE(built.ok());
  hipdlp::backend::FactoryResult factory = hipdlp::backend::createBackend(kind);
  INFO(factory.error.message);
  REQUIRE(factory.ok());
  hipdlp::backend::ProblemHandle problem;
  REQUIRE(factory.backend->bindProblem(built.problem, problem).ok());
  hipdlp::backend::VectorHandle x_handle;
  hipdlp::backend::VectorHandle ax_handle;
  hipdlp::backend::VectorHandle y_handle;
  hipdlp::backend::VectorHandle aty_handle;
  REQUIRE(factory.backend->createVector(lp.num_col_, x_handle).ok());
  REQUIRE(factory.backend->createVector(lp.num_row_, ax_handle).ok());
  REQUIRE(factory.backend->createVector(lp.num_row_, y_handle).ok());
  REQUIRE(factory.backend->createVector(lp.num_col_, aty_handle).ok());
  REQUIRE(factory.backend->uploadVector(x, x_handle).ok());
  REQUIRE(factory.backend->uploadVector(y, y_handle).ok());
  REQUIRE(factory.backend->multiply(problem, x_handle, ax_handle).ok());
  REQUIRE(
      factory.backend->multiplyTranspose(problem, y_handle, aty_handle).ok());
  ProductResult result;
  REQUIRE(factory.backend->downloadVector(ax_handle, result.ax).ok());
  REQUIRE(factory.backend->downloadVector(aty_handle, result.aty).ok());
  return result;
}

struct PrimitiveResult {
  std::vector<double> vector_update;
  double dot;
  double norm;
  double difference_norm;
  hipdlp::backend::ConvergenceMetrics metrics;
  std::vector<double> primal_update;
  std::vector<double> dual_update;
  std::vector<double> halpern_primal;
  std::vector<double> halpern_dual;
  std::vector<double> average_primal;
  std::vector<double> average_dual;
};

PrimitiveResult runPrimitiveSuite(const hipdlp::backend::Kind kind) {
  hipdlp::backend::FactoryResult factory = hipdlp::backend::createBackend(kind);
  INFO(factory.error.message);
  REQUIRE(factory.ok());
  hipdlp::backend::Backend& backend = *factory.backend;
  PrimitiveResult result;

  hipdlp::backend::VectorHandle left;
  hipdlp::backend::VectorHandle right;
  hipdlp::backend::VectorHandle output;
  REQUIRE(backend.createVector(3, left).ok());
  REQUIRE(backend.createVector(3, right).ok());
  REQUIRE(backend.createVector(3, output).ok());
  REQUIRE(backend.uploadVector({1.0, -2.0, 3.0}, left).ok());
  REQUIRE(backend.uploadVector({4.0, 5.0, -6.0}, right).ok());
  REQUIRE(backend.copy(left, output).ok());
  REQUIRE(backend.scale(2.0, output).ok());
  REQUIRE(backend.axpy(-0.5, right, output).ok());
  REQUIRE(backend.downloadVector(output, result.vector_update).ok());
  REQUIRE(backend.dot(left, right, result.dot).ok());
  REQUIRE(backend.norm2(left, result.norm).ok());
  REQUIRE(backend.differenceNorm2(left, right, result.difference_norm).ok());
  REQUIRE(
      backend.reduceConvergenceMetrics(left, right, 10.0, 7.0, result.metrics)
          .ok());

  HighsLp lp = makeEqualityLp(2, 2, {0, 1, 2}, {0, 1}, {1.0, 1.0});
  lp.col_cost_ = {1.0, -2.0};
  lp.col_lower_ = {0.0, -kHighsInf};
  lp.col_upper_ = {2.0, kHighsInf};
  lp.row_lower_ = {1.0, 2.0};
  lp.row_upper_ = {1.0, kHighsInf};
  hipdlp::model::ProblemBuildResult built = hipdlp::model::buildProblem(lp);
  INFO(built.error.message);
  REQUIRE(built.ok());
  hipdlp::backend::ProblemHandle problem;
  hipdlp::backend::WorkspaceHandle workspace;
  REQUIRE(backend.bindProblem(built.problem, problem).ok());
  REQUIRE(backend.createWorkspace(problem, workspace).ok());

  hipdlp::backend::VectorHandle current_primal;
  hipdlp::backend::VectorHandle transpose_product;
  hipdlp::backend::VectorHandle next_primal;
  hipdlp::backend::VectorHandle current_dual;
  hipdlp::backend::VectorHandle extrapolated_product;
  hipdlp::backend::VectorHandle next_dual;
  hipdlp::backend::VectorHandle primal_anchor;
  hipdlp::backend::VectorHandle dual_anchor;
  hipdlp::backend::VectorHandle primal_sum;
  hipdlp::backend::VectorHandle dual_sum;
  for (hipdlp::backend::VectorHandle* handle :
       {&current_primal, &transpose_product, &next_primal, &current_dual,
        &extrapolated_product, &next_dual, &primal_anchor, &dual_anchor,
        &primal_sum, &dual_sum})
    REQUIRE(backend.createVector(2, *handle).ok());

  REQUIRE(backend.uploadVector({1.0, 3.0}, current_primal).ok());
  REQUIRE(backend.uploadVector({4.0, -1.0}, transpose_product).ok());
  REQUIRE(backend
              .fusedProjectedPrimalUpdate(problem, workspace, current_primal,
                                          transpose_product, 0.5, next_primal)
              .ok());
  REQUIRE(backend.downloadVector(next_primal, result.primal_update).ok());

  REQUIRE(backend.uploadVector({0.5, -1.0}, current_dual).ok());
  REQUIRE(backend.uploadVector({0.0, 3.0}, extrapolated_product).ok());
  REQUIRE(backend
              .fusedProjectedDualUpdate(problem, workspace, current_dual,
                                        extrapolated_product, 2.0, next_dual)
              .ok());
  REQUIRE(backend.downloadVector(next_dual, result.dual_update).ok());

  REQUIRE(backend.uploadVector({0.0, 10.0}, primal_anchor).ok());
  REQUIRE(backend.uploadVector({4.0, 8.0}, dual_anchor).ok());
  REQUIRE(backend
              .fusedVectorBlend(workspace, primal_anchor, dual_anchor,
                                next_primal, next_dual, 0.25, primal_sum,
                                dual_sum)
              .ok());
  REQUIRE(backend.downloadVector(primal_sum, result.halpern_primal).ok());
  REQUIRE(backend.downloadVector(dual_sum, result.halpern_dual).ok());

  REQUIRE(backend.fill(primal_sum, 1.0).ok());
  REQUIRE(backend.fill(dual_sum, 1.0).ok());
  REQUIRE(backend
              .fusedAverageUpdate(workspace, next_primal, next_dual, 0.5,
                                  primal_sum, dual_sum)
              .ok());
  REQUIRE(backend.downloadVector(primal_sum, result.average_primal).ok());
  REQUIRE(backend.downloadVector(dual_sum, result.average_dual).ok());
  return result;
}

struct CudaCommandBufferFixture {
  std::unique_ptr<hipdlp::backend::Backend> backend;
  hipdlp::backend::CommandBufferRecorder* recorder;
  hipdlp::backend::VectorHandle source;
  hipdlp::backend::VectorHandle destination;
  hipdlp::backend::CommandBufferHandle command_buffer;

  CudaCommandBufferFixture() : recorder(nullptr) {}
};

CudaCommandBufferFixture makeCudaCommandBufferFixture() {
  CudaCommandBufferFixture fixture;
  fixture.backend = makeCudaBackend();
  fixture.recorder = hipdlp::backend::commandBufferRecorder(*fixture.backend);
  REQUIRE(fixture.recorder != nullptr);
  REQUIRE(fixture.backend->createVector(1, fixture.source).ok());
  REQUIRE(fixture.backend->createVector(1, fixture.destination).ok());
  REQUIRE(fixture.backend->uploadVector({2.0}, fixture.source).ok());
  REQUIRE(fixture.backend->fill(fixture.destination, 0.0).ok());
  REQUIRE(fixture.backend->synchronize().ok());
  return fixture;
}

hipdlp::api::Error recordCudaCommandBuffer(CudaCommandBufferFixture& fixture) {
  hipdlp::api::Error error = fixture.recorder->beginRecording();
  if (!error.ok()) return error;
  error = fixture.backend->copy(fixture.source, fixture.destination);
  if (!error.ok()) {
    fixture.recorder->cancelRecording();
    return error;
  }
  error = fixture.backend->scale(3.0, fixture.destination);
  if (!error.ok()) {
    fixture.recorder->cancelRecording();
    return error;
  }
  return fixture.recorder->endRecording(fixture.command_buffer);
}

bool bitwiseVectorEqual(const std::vector<double>& left,
                        const std::vector<double>& right) {
  return left.size() == right.size() &&
         (left.empty() || std::memcmp(left.data(), right.data(),
                                      left.size() * sizeof(double)) == 0);
}
#endif

}  // namespace

TEST_CASE("hi-pdlp-cpu-backend-factory-and-capabilities",
          "[pdlp][hipdlp][phase3][phase4][phase5][backend]") {
  hipdlp::backend::FactoryResult explicit_cpu =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kCpu);
  REQUIRE(explicit_cpu.ok());
  CHECK(explicit_cpu.backend->descriptor().kind == hipdlp::backend::Kind::kCpu);
  CHECK(explicit_cpu.backend->descriptor().precision.isFp64());
  CHECK(explicit_cpu.backend->descriptor().capabilities.supports(
      hipdlp::backend::Capability::kSparseMatrixVector));
  CHECK(explicit_cpu.backend->descriptor().capabilities.supports(
      hipdlp::backend::Capability::kVectorOperations));
  CHECK_FALSE(explicit_cpu.backend->descriptor().capabilities.supports(
      hipdlp::backend::Capability::kBatchedVectors));
  CHECK_FALSE(explicit_cpu.backend->descriptor().capabilities.supports(
      hipdlp::backend::Capability::kSparseMatrixMatrix));
  CHECK_FALSE(explicit_cpu.backend->descriptor().capabilities.supports(
      hipdlp::backend::Capability::kCommandBuffer));
  CHECK(explicit_cpu.backend
            ->requireCapability(hipdlp::backend::Capability::kCommandBuffer)
            .code == hipdlp::api::ErrorCode::kUnsupportedCapability);
  CHECK(hipdlp::backend::commandBufferRecorder(*explicit_cpu.backend) ==
        nullptr);

  hipdlp::backend::DeviceInfo device;
  REQUIRE(explicit_cpu.backend->queryDevice(device).ok());
  CHECK(device.available);
  CHECK_FALSE(device.name.empty());

  CHECK(hipdlp::backend::createBackend(hipdlp::backend::Kind::kOpenmp)
            .error.code == hipdlp::api::ErrorCode::kBackendUnavailable);
  CHECK(
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kRocm).error.code ==
      hipdlp::api::ErrorCode::kBackendUnavailable);
  const hipdlp::backend::Descriptor* cuda =
      hipdlp::backend::findRegistration(hipdlp::backend::Kind::kCuda);
  REQUIRE(cuda != nullptr);
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
  CHECK(cuda->implemented);
#else
  CHECK_FALSE(cuda->implemented);
#endif
  const hipdlp::backend::FactoryResult cuda_result =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda);
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
  CHECK((cuda_result.ok() || cuda_result.error.code ==
                                 hipdlp::api::ErrorCode::kBackendUnavailable));
  if (!cuda_result.ok()) {
    CHECK(cuda_result.error.code ==
          hipdlp::api::ErrorCode::kBackendUnavailable);
    CHECK(cuda_result.backend.get() == nullptr);
  }
#else
  CHECK_FALSE(cuda_result.ok());
  CHECK(cuda_result.error.code ==
        (cuda->build_enabled ? hipdlp::api::ErrorCode::kNotReady
                             : hipdlp::api::ErrorCode::kBackendUnavailable));
  CHECK(cuda_result.backend.get() == nullptr);
  const hipdlp::backend::FactoryResult repeated_cuda_result =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda);
  CHECK_FALSE(repeated_cuda_result.ok());
  CHECK(repeated_cuda_result.error.code == cuda_result.error.code);
  CHECK(repeated_cuda_result.error.message == cuda_result.error.message);
  CHECK(repeated_cuda_result.backend.get() == nullptr);
#endif

  hipdlp::backend::FactoryResult automatic =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kAuto);
  REQUIRE(automatic.ok());
  CHECK(automatic.backend->descriptor().kind ==
        (cuda_result.ok() ? hipdlp::backend::Kind::kCuda
                          : hipdlp::backend::Kind::kCpu));

  hipdlp::backend::VectorHandle vector;
  CHECK(explicit_cpu.backend
            ->requireCapability(hipdlp::backend::Capability::kBatchedVectors)
            .code == hipdlp::api::ErrorCode::kUnsupportedCapability);
  CHECK(explicit_cpu.backend->createBatchedVector(3, 2, vector).code ==
        hipdlp::api::ErrorCode::kUnsupportedCapability);
  hipdlp::backend::ProblemHandle problem;
  CHECK(explicit_cpu.backend->multiplyBatched(problem, vector, vector).code ==
        hipdlp::api::ErrorCode::kUnsupportedCapability);
}

TEST_CASE("hi-pdlp-rocm-backend-remains-unavailable",
          "[pdlp][hipdlp][phase6][rocm][factory]") {
  CHECK(hipdlp::backend::isValidKind(hipdlp::backend::Kind::kRocm));
  const hipdlp::backend::Descriptor* registration =
      hipdlp::backend::findRegistration(hipdlp::backend::Kind::kRocm);
  REQUIRE(registration != nullptr);
  CHECK(registration->kind == hipdlp::backend::Kind::kRocm);
  CHECK(std::string(registration->name) == "rocm");
  CHECK_FALSE(registration->build_enabled);
  CHECK_FALSE(registration->implemented);
  CHECK(registration->precision.isFp64());
  CHECK(registration->capabilities.mask == 0);

  hipdlp::backend::FactoryResult first =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kRocm);
  CHECK_FALSE(first.ok());
  CHECK(first.error.code == hipdlp::api::ErrorCode::kBackendUnavailable);
  CHECK(first.error.message == "requested HiPDLP backend is unavailable");
  CHECK(first.backend.get() == nullptr);

  hipdlp::backend::FactoryResult second =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kRocm);
  CHECK_FALSE(second.ok());
  CHECK(second.error.code == first.error.code);
  CHECK(second.error.message == first.error.message);
  CHECK(second.backend.get() == nullptr);

#ifdef HIGHS_HIPDLP_CUDA_DETECTED
  const hipdlp::backend::FactoryResult cuda =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda);
  const bool cuda_usable = cuda.ok();
#else
  const bool cuda_usable = false;
#endif
  for (int request = 0; request < 3; ++request) {
    hipdlp::backend::FactoryResult automatic =
        hipdlp::backend::createBackend(hipdlp::backend::Kind::kAuto);
    INFO(automatic.error.message);
    REQUIRE(automatic.ok());
    CHECK(automatic.backend->descriptor().kind ==
          (cuda_usable ? hipdlp::backend::Kind::kCuda
                       : hipdlp::backend::Kind::kCpu));
    CHECK(automatic.backend->descriptor().kind != hipdlp::backend::Kind::kRocm);
  }
}

TEST_CASE("hi-pdlp-rocm-reservation-has-no-implementation-surface",
          "[pdlp][hipdlp][phase6][rocm][structure]") {
  CHECK_FALSE(sourcePathExists("highs/pdlp/hipdlp/backend/rocm"));

  const std::string backend_types =
      readSourceFile("highs/pdlp/hipdlp/backend/types.hpp");
  const std::string backend_type_definitions =
      readSourceFile("highs/pdlp/hipdlp/backend/types.cc");
  const std::string backend_factory =
      readSourceFile("highs/pdlp/hipdlp/backend/backend_factory.cc");
  CHECK(countText(backend_types, "kRocm") == 1);
  CHECK(countText(backend_type_definitions, "kRocm") == 1);
  CHECK(countText(backend_factory, "kRocm") == 1);
  CHECK(countText(backend_factory, "\"rocm\"") == 1);
  CHECK(backend_factory.find("RocmBackend") == std::string::npos);
  CHECK(backend_factory.find("createRocmBackend") == std::string::npos);
  CHECK(backend_factory.find("namespace rocm") == std::string::npos);

  const char* rocm_neutral_files[] = {
      "highs/pdlp/hipdlp/algorithm/fixed_trace.cc",
      "highs/pdlp/hipdlp/algorithm/fixed_trace.hpp",
      "highs/pdlp/hipdlp/algorithm/halpern_pdhg_step.cc",
      "highs/pdlp/hipdlp/algorithm/halpern_pdhg_step.hpp",
      "highs/pdlp/hipdlp/algorithm/solver.cc",
      "highs/pdlp/hipdlp/algorithm/solver.hpp",
      "highs/pdlp/hipdlp/api/solver.cc",
      "highs/pdlp/hipdlp/api/solver.hpp",
      "highs/pdlp/hipdlp/api/status.cc",
      "highs/pdlp/hipdlp/api/status.hpp",
      "highs/pdlp/hipdlp/backend/backend.hpp",
      "highs/pdlp/hipdlp/backend/backend_factory.hpp",
      "highs/pdlp/hipdlp/backend/command_buffer.cc",
      "highs/pdlp/hipdlp/backend/command_buffer.hpp",
      "highs/pdlp/hipdlp/backend/cpu/cpu_backend.cc",
      "highs/pdlp/hipdlp/backend/cpu/cpu_backend.hpp",
      "highs/pdlp/hipdlp/backend/cpu/sparse_kernels.cc",
      "highs/pdlp/hipdlp/backend/cpu/sparse_kernels.hpp",
      "highs/pdlp/hipdlp/backend/cpu/vector_kernels.cc",
      "highs/pdlp/hipdlp/backend/cpu/vector_kernels.hpp",
      "highs/pdlp/hipdlp/backend/cuda/cuda_backend.cu",
      "highs/pdlp/hipdlp/backend/cuda/cuda_backend.hpp",
      "highs/pdlp/hipdlp/backend/cuda/cuda_backend_testing.hpp",
      "highs/pdlp/hipdlp/backend/cuda/cuda_kernels.cu",
      "highs/pdlp/hipdlp/backend/cuda/cuda_kernels.cuh",
      "highs/pdlp/hipdlp/legacy/legacy_solver_adapter.cc",
      "highs/pdlp/hipdlp/legacy/legacy_solver_adapter.hpp",
      "highs/pdlp/hipdlp/model/problem.cc",
      "highs/pdlp/hipdlp/model/problem.hpp",
      "highs/pdlp/hipdlp/model/scaling.cc",
      "highs/pdlp/hipdlp/model/scaling.hpp",
      "highs/pdlp/hipdlp/verification/kkt_verifier.cc",
      "highs/pdlp/hipdlp/verification/kkt_verifier.hpp",
      "highs/pdlp/HiPdlpWrapper.cpp",
      "highs/pdlp/HiPdlpWrapper.h",
      "highs/pdlp/hipdlp/pdhg.cc",
      "highs/pdlp/hipdlp/pdhg.hpp",
      "highs/lp_data/HighsOptions.h"};
  const char* rocm_symbols[] = {
      "kRocm",      "RocmBackend",   "rocm_backend", "namespace rocm",
      "ROCm",       "rocblas",       "rocsparse",    "hip_runtime",
      "hipblas",    "hipsparse",     "hipStream_t",  "hipEvent_t",
      "hipError_t", "hipDevice",     "hipMalloc",    "hipFree",
      "hipMemcpy",  "#include <hip/"};
  for (const char* path : rocm_neutral_files) {
    CAPTURE(path);
    const std::string contents = readSourceFile(path);
    for (const char* symbol : rocm_symbols) {
      CAPTURE(symbol);
      CHECK(contents.find(symbol) == std::string::npos);
    }
  }

  // Option translation maps the reserved public name and nothing else.
  const char* rocm_name_mapping_files[] = {"highs/pdlp/hipdlp/api/options.cc",
                                           "highs/pdlp/hipdlp/api/options.hpp"};
  const char* rocm_implementation_symbols[] = {
      "RocmBackend",   "rocm_backend", "namespace rocm", "ROCm",
      "rocblas",       "rocsparse",    "hip_runtime",    "hipblas",
      "hipsparse",     "hipStream_t",  "hipEvent_t",     "hipError_t",
      "hipDevice",     "hipMalloc",    "hipFree",        "hipMemcpy",
      "#include <hip/"};
  for (const char* path : rocm_name_mapping_files) {
    CAPTURE(path);
    const std::string contents = readSourceFile(path);
    for (const char* symbol : rocm_implementation_symbols) {
      CAPTURE(symbol);
      CHECK(contents.find(symbol) == std::string::npos);
    }
  }
  CHECK(countText(readSourceFile("highs/pdlp/hipdlp/api/options.cc"),
                  "kRocm") == 1);

  const char* build_files[] = {"CMakeLists.txt",
                               "cmake/sources.cmake",
                               "highs/CMakeLists.txt",
                               "highs/HConfig.h.in",
                               "highs/HConfig.h.bazel.in",
                               "highs/HConfig.h.meson.in",
                               "BUILD.bazel",
                               "MODULE.bazel",
                               "WORKSPACE",
                               "meson.build",
                               "meson_options.txt",
                               "highs/meson.build",
                               "check/meson.build"};
  const char* rocm_build_terms[] = {"HIGHS_HIPDLP_ROCM",
                                    "HIGHS_HIPDLP_HIP",
                                    "enable_language(HIP)",
                                    "check_language(HIP)",
                                    "CMAKE_HIP",
                                    "HIP::",
                                    "find_package(HIP",
                                    "find_package(hip",
                                    "hip_runtime",
                                    "hipblas",
                                    "hipsparse",
                                    "backend/rocm",
                                    "RocmBackend",
                                    "hipcc",
                                    "/opt/rocm",
                                    "rocblas",
                                    "rocsparse",
                                    "rocBLAS",
                                    "rocSPARSE",
                                    "ROCm"};
  for (const char* path : build_files) {
    CAPTURE(path);
    const std::string contents = readSourceFile(path);
    for (const char* term : rocm_build_terms) {
      CAPTURE(term);
      CHECK(contents.find(term) == std::string::npos);
    }
  }
}

TEST_CASE("hi-pdlp-phase5h-shared-step-ownership-boundary",
          "[pdlp][hipdlp][phase5h][phase6][structure]") {
  const std::string backend_types =
      readSourceFile("highs/pdlp/hipdlp/backend/types.hpp");
  const std::string backend_contract =
      readSourceFile("highs/pdlp/hipdlp/backend/backend.hpp");
  const char* removed_contract_symbols[] = {
      "FixedIterationPlanVectors", "createFixedIterationPlan",
      "executeFixedIterationPlan", "downloadFixedIterationTrace"};
  for (const char* symbol : removed_contract_symbols) {
    CAPTURE(symbol);
    CHECK(backend_types.find(symbol) == std::string::npos);
    CHECK(backend_contract.find(symbol) == std::string::npos);
  }
  const char* algorithm_names[] = {"primal_anchor", "reflected_dual",
                                   "current_aty", "Halpern"};
  for (const char* name : algorithm_names) {
    CAPTURE(name);
    CHECK(backend_contract.find(name) == std::string::npos);
  }

  const std::string cpu_header =
      readSourceFile("highs/pdlp/hipdlp/backend/cpu/cpu_backend.hpp");
  const std::string cpu_source =
      readSourceFile("highs/pdlp/hipdlp/backend/cpu/cpu_backend.cc");
  const char* recorder_methods[] = {"beginRecording", "endRecording",
                                    "cancelRecording", "replay("};
  for (const char* method : recorder_methods) {
    CAPTURE(method);
    CHECK(cpu_header.find(method) == std::string::npos);
    CHECK(cpu_source.find(method) == std::string::npos);
  }

  const std::string cuda_source =
      readSourceFile("highs/pdlp/hipdlp/backend/cuda/cuda_backend.cu");
  const std::string cuda_kernels =
      readSourceFile("highs/pdlp/hipdlp/backend/cuda/cuda_kernels.cu");
  const char* cuda_algorithm_symbols[] = {"enqueuePlannedIteration",
                                          "CudaExecutionPlanResource",
                                          "HalpernPdhgState",
                                          "enqueueHalpernPdhgStep",
                                          "primal_anchor",
                                          "reflected_dual",
                                          "current_aty",
                                          "Halpern iteration"};
  for (const char* symbol : cuda_algorithm_symbols) {
    CAPTURE(symbol);
    CHECK(cuda_source.find(symbol) == std::string::npos);
    CHECK(cuda_kernels.find(symbol) == std::string::npos);
  }

  const std::string fixed_trace =
      readSourceFile("highs/pdlp/hipdlp/algorithm/fixed_trace.cc");
  const std::string shared_step =
      readSourceFile("highs/pdlp/hipdlp/algorithm/halpern_pdhg_step.cc");
  CHECK(countText(shared_step, "api::Error enqueueHalpernPdhgStep(") == 1);
  CHECK(countText(fixed_trace, "enqueueHalpernPdhgStep(*backend, state)") == 2);
  const std::size_t direct_path = fixed_trace.find(
      "options.execution_mode == FixedTraceExecutionMode::kDirect");
  const std::size_t direct_call =
      fixed_trace.find("enqueueHalpernPdhgStep(*backend, state)", direct_path);
  const std::size_t recording = fixed_trace.find("recorder->beginRecording()");
  const std::size_t recorded_call =
      fixed_trace.find("enqueueHalpernPdhgStep(*backend, state)", recording);
  CHECK(direct_path != std::string::npos);
  CHECK(direct_call != std::string::npos);
  CHECK(recording != std::string::npos);
  CHECK(recorded_call != std::string::npos);

  const char* backend_independent_files[] = {
      "highs/pdlp/hipdlp/algorithm/fixed_trace.cc",
      "highs/pdlp/hipdlp/algorithm/fixed_trace.hpp",
      "highs/pdlp/hipdlp/algorithm/halpern_pdhg_step.cc",
      "highs/pdlp/hipdlp/algorithm/halpern_pdhg_step.hpp",
      "highs/pdlp/hipdlp/algorithm/solver.cc",
      "highs/pdlp/hipdlp/algorithm/solver.hpp",
      "highs/pdlp/hipdlp/api/solver.cc",
      "highs/pdlp/hipdlp/api/solver.hpp",
      "highs/pdlp/hipdlp/api/status.cc",
      "highs/pdlp/hipdlp/api/status.hpp",
      "highs/pdlp/hipdlp/model/problem.cc",
      "highs/pdlp/hipdlp/model/problem.hpp",
      "highs/pdlp/hipdlp/model/scaling.cc",
      "highs/pdlp/hipdlp/model/scaling.hpp",
      "highs/pdlp/hipdlp/legacy/legacy_solver_adapter.cc",
      "highs/pdlp/hipdlp/legacy/legacy_solver_adapter.hpp"};
  for (const char* path : backend_independent_files) {
    CAPTURE(path);
    const std::string contents = readSourceFile(path);
    CHECK(contents.find("cuda") == std::string::npos);
    CHECK(contents.find("CUDA") == std::string::npos);
    CHECK(contents.find("graph") == std::string::npos);
    CHECK(contents.find("Graph") == std::string::npos);
  }
  const std::string public_options =
      readSourceFile("highs/lp_data/HighsOptions.h");
  // The public backend option names cuda without importing any CUDA surface.
  const char* option_name_mapping_files[] = {
      "highs/pdlp/hipdlp/api/options.cc", "highs/pdlp/hipdlp/api/options.hpp",
      "highs/lp_data/HighsOptions.h"};
  const char* cuda_surface_symbols[] = {
      "cuda_runtime",  "cublas",      "cusparse",       "cudaStream",
      "cudaError",     "cudaMalloc",  "cudaFree",       "cudaMemcpy",
      "cudaGraph",     "CudaBackend", "namespace cuda", "#include <cuda",
      "command_buffer"};
  for (const char* path : option_name_mapping_files) {
    CAPTURE(path);
    const std::string contents = readSourceFile(path);
    for (const char* symbol : cuda_surface_symbols) {
      CAPTURE(symbol);
      CHECK(contents.find(symbol) == std::string::npos);
    }
  }
  CHECK(public_options.find("CUDA") == std::string::npos);
  CHECK(public_options.find("cudaGraph") == std::string::npos);
  CHECK(public_options.find("command_buffer") == std::string::npos);
  std::string options_without_existing_graph_words = public_options;
  const char* existing_graph_words[] = {"conflict graph", "lexicographically"};
  for (const char* word : existing_graph_words) {
    const std::size_t position =
        options_without_existing_graph_words.find(word);
    REQUIRE(position != std::string::npos);
    options_without_existing_graph_words.erase(position, std::strlen(word));
  }
  CHECK(options_without_existing_graph_words.find("graph") ==
        std::string::npos);
  CHECK(options_without_existing_graph_words.find("Graph") ==
        std::string::npos);
}

#ifdef HIGHS_HIPDLP_CUDA_DETECTED
TEST_CASE("hi-pdlp-cuda-backend-factory-and-device",
          "[pdlp][hipdlp][phase5][cuda][factory]") {
  const hipdlp::backend::Descriptor* registration =
      hipdlp::backend::findRegistration(hipdlp::backend::Kind::kCuda);
  REQUIRE(registration != nullptr);
  CHECK(registration->build_enabled);
  CHECK(registration->implemented);
  CHECK(registration->precision.isFp64());

  hipdlp::backend::FactoryResult first =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda);
  if (!first.ok()) {
    CHECK_FALSE(first.ok());
    CHECK(first.backend.get() == nullptr);
    CHECK(first.error.code == hipdlp::api::ErrorCode::kBackendUnavailable);
    CHECK_FALSE(first.error.message.empty());
    const hipdlp::api::Error first_error = first.error;
    for (int repetition = 1; repetition < 4; ++repetition) {
      hipdlp::backend::FactoryResult explicit_cuda =
          hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda);
      CHECK_FALSE(explicit_cuda.ok());
      CHECK(explicit_cuda.backend.get() == nullptr);
      CHECK(explicit_cuda.error.code ==
            hipdlp::api::ErrorCode::kBackendUnavailable);
      CHECK_FALSE(explicit_cuda.error.message.empty());
      CHECK(explicit_cuda.error.code == first_error.code);
      CHECK(explicit_cuda.error.message == first_error.message);
    }
    hipdlp::backend::FactoryResult automatic =
        hipdlp::backend::createBackend(hipdlp::backend::Kind::kAuto);
    INFO(automatic.error.message);
    REQUIRE(automatic.ok());
    CHECK(automatic.backend->descriptor().kind == hipdlp::backend::Kind::kCpu);
    return;
  }

  for (int repetition = 0; repetition < 4; ++repetition) {
    hipdlp::backend::FactoryResult explicit_cuda =
        repetition == 0
            ? std::move(first)
            : hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda);
    INFO(explicit_cuda.error.message);
    REQUIRE(explicit_cuda.ok());
    REQUIRE(explicit_cuda.backend.get() != nullptr);
    CHECK(explicit_cuda.backend->descriptor().kind ==
          hipdlp::backend::Kind::kCuda);
    hipdlp::backend::DeviceInfo device;
    REQUIRE(explicit_cuda.backend->queryDevice(device).ok());
    CHECK(device.available);
    CHECK_FALSE(device.name.empty());
    REQUIRE(explicit_cuda.backend->synchronize().ok());
  }

  hipdlp::backend::FactoryResult automatic =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kAuto);
  INFO(automatic.error.message);
  REQUIRE(automatic.ok());
  CHECK(automatic.backend->descriptor().kind == hipdlp::backend::Kind::kCuda);
}

TEST_CASE("hi-pdlp-cuda-command-buffer-factory-capability",
          "[pdlp][hipdlp][phase5g][cuda][graph][factory]") {
  hipdlp::backend::FactoryResult factory =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda);
  INFO(factory.error.message);
  REQUIRE(factory.ok());
  REQUIRE(factory.backend.get() != nullptr);
  CHECK(factory.backend->descriptor().capabilities.supports(
      hipdlp::backend::Capability::kCommandBuffer));
  CHECK(factory.backend
            ->requireCapability(hipdlp::backend::Capability::kCommandBuffer)
            .ok());
  CHECK(hipdlp::backend::commandBufferRecorder(*factory.backend) != nullptr);
}

TEST_CASE("hi-pdlp-cuda-command-buffer-reuse-is-deterministic",
          "[pdlp][hipdlp][phase5g][cuda][graph][reuse][determinism]") {
  CudaCommandBufferFixture fixture = makeCudaCommandBufferFixture();
  REQUIRE(recordCudaCommandBuffer(fixture).ok());
  REQUIRE(fixture.command_buffer.valid());

  hipdlp::backend::CommandBufferStatistics first_statistics;
  REQUIRE(
      fixture.recorder->replay(fixture.command_buffer, first_statistics).ok());
  CHECK(first_statistics.replay_calls == 1);
  CHECK(first_statistics.direct_operation_enqueue_calls == 0);
  CHECK(first_statistics.allocation_calls == 0);
  CHECK(first_statistics.descriptor_creation_calls == 0);
  CHECK(first_statistics.host_to_device_transfer_calls == 0);
  CHECK(first_statistics.device_to_host_transfer_calls == 0);
  CHECK(first_statistics.synchronization_calls == 0);
  REQUIRE(fixture.backend->synchronize().ok());
  const std::vector<double> first =
      download(*fixture.backend, fixture.destination);
  CHECK(first == std::vector<double>({6.0}));

  REQUIRE(fixture.backend->fill(fixture.destination, 0.0).ok());
  REQUIRE(fixture.backend->synchronize().ok());
  hipdlp::backend::CommandBufferStatistics second_statistics;
  REQUIRE(
      fixture.recorder->replay(fixture.command_buffer, second_statistics).ok());
  REQUIRE(fixture.backend->synchronize().ok());
  const std::vector<double> second =
      download(*fixture.backend, fixture.destination);
  CHECK(second_statistics.replay_calls == 1);
  CHECK(bitwiseVectorEqual(first, second));
}

TEST_CASE("hi-pdlp-cuda-command-buffer-construction-is-transactional",
          "[pdlp][hipdlp][phase5g][cuda][graph][resources][failure]") {
  const hipdlp::backend::cuda::testing::FailurePoint failure_points[] = {
      hipdlp::backend::cuda::testing::FailurePoint::kGraphCapture,
      hipdlp::backend::cuda::testing::FailurePoint::kGraphInstantiation};
  for (const hipdlp::backend::cuda::testing::FailurePoint failure_point :
       failure_points) {
    CAPTURE(static_cast<int>(failure_point));
    CudaCommandBufferFixture fixture = makeCudaCommandBufferFixture();
    hipdlp::backend::cuda::testing::resetResourceAudit();
    hipdlp::backend::cuda::testing::setFailurePoint(failure_point);
    const hipdlp::api::Error error = recordCudaCommandBuffer(fixture);
    hipdlp::backend::cuda::testing::setFailurePoint(
        hipdlp::backend::cuda::testing::FailurePoint::kNone);
    CHECK(error.code == hipdlp::api::ErrorCode::kBackendFailure);
    CHECK_FALSE(error.message.empty());
    CHECK_FALSE(fixture.command_buffer.valid());
    const hipdlp::backend::cuda::testing::ResourceAudit failed_audit =
        hipdlp::backend::cuda::testing::resourceAudit();
    CHECK(failed_audit.allocation_calls == failed_audit.free_calls);
    CHECK(failed_audit.graph_creations == failed_audit.graph_destructions);
    CHECK(failed_audit.graph_exec_creations ==
          failed_audit.graph_exec_destructions);

    const hipdlp::api::Error recovered = recordCudaCommandBuffer(fixture);
    INFO(recovered.message);
    REQUIRE(recovered.ok());
    REQUIRE(fixture.command_buffer.valid());
  }
}

TEST_CASE("hi-pdlp-cuda-command-buffer-validates-resource-lifetimes",
          "[pdlp][hipdlp][phase5g][cuda][graph][resources][lifetime]") {
  CudaCommandBufferFixture fixture = makeCudaCommandBufferFixture();
  REQUIRE(recordCudaCommandBuffer(fixture).ok());
  fixture.source.reset();
  hipdlp::backend::cuda::testing::resetResourceAudit();
  hipdlp::backend::CommandBufferStatistics statistics;
  const hipdlp::api::Error error =
      fixture.recorder->replay(fixture.command_buffer, statistics);
  CHECK(error.code == hipdlp::api::ErrorCode::kInvalidArgument);
  CHECK_FALSE(error.message.empty());
  CHECK(statistics.replay_calls == 0);
  const hipdlp::backend::cuda::testing::ResourceAudit audit =
      hipdlp::backend::cuda::testing::resourceAudit();
  CHECK(audit.graph_launch_calls == 0);
  CHECK(audit.direct_primitive_enqueue_calls == 0);
}
#endif

TEST_CASE("hi-pdlp-openmp-backend-remains-unavailable",
          "[pdlp][hipdlp][phase4][backend][openmp]") {
  const hipdlp::backend::Descriptor* registration =
      hipdlp::backend::findRegistration(hipdlp::backend::Kind::kOpenmp);
  REQUIRE(registration != nullptr);
  CHECK(registration->kind == hipdlp::backend::Kind::kOpenmp);
  CHECK(std::string(registration->name) == "openmp");
  CHECK_FALSE(registration->build_enabled);
  CHECK_FALSE(registration->implemented);

  hipdlp::backend::FactoryResult first =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kOpenmp);
  CHECK_FALSE(first.ok());
  CHECK(first.error.code == hipdlp::api::ErrorCode::kBackendUnavailable);
  CHECK(first.error.message == "requested HiPDLP backend is unavailable");
  CHECK(first.backend.get() == nullptr);

  hipdlp::backend::FactoryResult second =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kOpenmp);
  CHECK_FALSE(second.ok());
  CHECK(second.error.code == first.error.code);
  CHECK(second.error.message == first.error.message);
  CHECK(second.backend.get() == nullptr);

#ifdef HIGHS_HIPDLP_CUDA_DETECTED
  const hipdlp::backend::FactoryResult cuda =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda);
  const bool cuda_usable = cuda.ok();
#else
  const bool cuda_usable = false;
#endif
  for (int request = 0; request < 3; ++request) {
    hipdlp::backend::FactoryResult automatic =
        hipdlp::backend::createBackend(hipdlp::backend::Kind::kAuto);
    INFO(automatic.error.message);
    REQUIRE(automatic.ok());
    CHECK(automatic.backend->descriptor().kind ==
          (cuda_usable ? hipdlp::backend::Kind::kCuda
                       : hipdlp::backend::Kind::kCpu));
    CHECK(automatic.backend->descriptor().precision.isFp64());
  }
}

TEST_CASE("hi-pdlp-cpu-handle-ownership",
          "[pdlp][hipdlp][phase3][phase4][backend]") {
  CHECK_FALSE(
      std::is_copy_constructible<hipdlp::backend::ProblemHandle>::value);
  CHECK(std::is_nothrow_move_constructible<
        hipdlp::backend::ProblemHandle>::value);
  CHECK_FALSE(std::is_copy_constructible<hipdlp::backend::VectorHandle>::value);
  CHECK(
      std::is_nothrow_move_constructible<hipdlp::backend::VectorHandle>::value);
  CHECK_FALSE(
      std::is_copy_constructible<hipdlp::backend::WorkspaceHandle>::value);
  CHECK(std::is_nothrow_move_constructible<
        hipdlp::backend::WorkspaceHandle>::value);
  CHECK_FALSE(
      std::is_copy_constructible<hipdlp::backend::CommandBufferHandle>::value);
  CHECK(std::is_nothrow_move_constructible<
        hipdlp::backend::CommandBufferHandle>::value);

  int live = 0;
  {
    hipdlp::backend::ProblemHandle first(
        std::unique_ptr<hipdlp::backend::Resource>(new TrackingResource(live)));
    REQUIRE(first.valid());
    CHECK(live == 1);
    hipdlp::backend::ProblemHandle second(std::move(first));
    CHECK_FALSE(first.valid());
    CHECK(second.valid());
    second.reset();
    CHECK(live == 0);
  }
  {
    hipdlp::backend::VectorHandle vector(
        std::unique_ptr<hipdlp::backend::Resource>(new TrackingResource(live)));
    hipdlp::backend::WorkspaceHandle workspace(
        std::unique_ptr<hipdlp::backend::Resource>(new TrackingResource(live)));
    CHECK(live == 2);
    hipdlp::backend::VectorHandle moved_vector(std::move(vector));
    hipdlp::backend::WorkspaceHandle moved_workspace(std::move(workspace));
    CHECK_FALSE(vector.valid());
    CHECK_FALSE(workspace.valid());
    CHECK(moved_vector.valid());
    CHECK(moved_workspace.valid());
  }
  {
    hipdlp::backend::CommandBufferHandle command_buffer(
        std::unique_ptr<hipdlp::backend::Resource>(new TrackingResource(live)));
    REQUIRE(command_buffer.valid());
    hipdlp::backend::CommandBufferHandle moved_command_buffer(
        std::move(command_buffer));
    CHECK_FALSE(command_buffer.valid());
    CHECK(moved_command_buffer.valid());
  }
  CHECK(live == 0);
}

TEST_CASE("hi-pdlp-cpu-sparse-products",
          "[pdlp][hipdlp][phase3][phase4][backend]") {
  SECTION("diagonal") {
    checkProducts(
        makeEqualityLp(3, 3, {0, 1, 2, 3}, {0, 1, 2}, {2.0, -3.0, 4.0}),
        {1.0, 2.0, -1.0}, {2.0, -6.0, -4.0}, {1.0, 2.0, -1.0},
        {2.0, -6.0, -4.0});
  }

  SECTION("rectangular") {
    checkProducts(
        makeEqualityLp(2, 3, {0, 2, 3, 4}, {0, 1, 1, 0}, {1.0, -1.0, 3.0, 2.0}),
        {2.0, -1.0, 4.0}, {10.0, -5.0}, {3.0, -2.0}, {5.0, -6.0, 6.0});
  }

  SECTION("empty") {
    checkProducts(makeEqualityLp(2, 3, {0, 0, 0, 0}, {}, {}), {2.0, -1.0, 4.0},
                  {0.0, 0.0}, {3.0, -2.0}, {0.0, 0.0, 0.0});
  }

  SECTION("skewed") {
    checkProducts(makeEqualityLp(2, 2, {0, 1, 2}, {0, 1}, {1e-200, 1e200}),
                  {1e200, 1e-200}, {1.0, 1.0}, {1e200, 1e-200}, {1.0, 1.0});
  }
}

TEST_CASE("hi-pdlp-cpu-vector-and-reduction-primitives",
          "[pdlp][hipdlp][phase3][phase4][backend]") {
  std::unique_ptr<hipdlp::backend::Backend> backend = makeCpuBackend();
  hipdlp::backend::VectorHandle left;
  hipdlp::backend::VectorHandle right;
  hipdlp::backend::VectorHandle output;
  REQUIRE(backend->createVector(3, left).ok());
  REQUIRE(backend->createVector(3, right).ok());
  REQUIRE(backend->createVector(3, output).ok());
  REQUIRE(backend->uploadVector({1.0, -2.0, 3.0}, left).ok());
  REQUIRE(backend->uploadVector({4.0, 5.0, -6.0}, right).ok());
  REQUIRE(backend->copy(left, output).ok());
  REQUIRE(backend->scale(2.0, output).ok());
  REQUIRE(backend->axpy(-0.5, right, output).ok());
  checkApproxVector(download(*backend, output), {0.0, -6.5, 9.0});

  double reduced = 0.0;
  REQUIRE(backend->dot(left, right, reduced).ok());
  CHECK(reduced == -24.0);
  REQUIRE(backend->norm2(left, reduced).ok());
  CHECK(reduced == Approx(std::sqrt(14.0)));
  REQUIRE(backend->differenceNorm2(left, right, reduced).ok());
  CHECK(reduced == Approx(std::sqrt(139.0)));

  REQUIRE(backend->uploadVector({1e200, 1e-200, 0.0}, left).ok());
  REQUIRE(backend->norm2(left, reduced).ok());
  CHECK(reduced == Approx(1e200));

  hipdlp::backend::ConvergenceMetrics metrics;
  REQUIRE(backend->uploadVector({1.0, -2.0, 3.0}, left).ok());
  REQUIRE(
      backend->reduceConvergenceMetrics(left, right, 10.0, 7.0, metrics).ok());
  CHECK(metrics.primal_residual == Approx(std::sqrt(14.0)));
  CHECK(metrics.dual_residual == Approx(std::sqrt(77.0)));
  CHECK(metrics.objective_gap == 3.0);
  REQUIRE(backend->synchronize().ok());

  hipdlp::backend::VectorHandle empty_left;
  hipdlp::backend::VectorHandle empty_right;
  REQUIRE(backend->createVector(0, empty_left).ok());
  REQUIRE(backend->createVector(0, empty_right).ok());
  REQUIRE(backend->dot(empty_left, empty_right, reduced).ok());
  CHECK(reduced == 0.0);
  REQUIRE(backend->norm2(empty_left, reduced).ok());
  CHECK(reduced == 0.0);
}

TEST_CASE("hi-pdlp-cpu-fused-primitives",
          "[pdlp][hipdlp][phase3][phase4][backend]") {
  HighsLp lp = makeEqualityLp(2, 2, {0, 1, 2}, {0, 1}, {1.0, 1.0});
  lp.col_cost_ = {1.0, -2.0};
  lp.col_lower_ = {0.0, -kHighsInf};
  lp.col_upper_ = {2.0, kHighsInf};
  lp.row_lower_ = {1.0, 2.0};
  lp.row_upper_ = {1.0, kHighsInf};
  hipdlp::model::ProblemBuildResult built = hipdlp::model::buildProblem(lp);
  REQUIRE(built.ok());
  std::unique_ptr<hipdlp::backend::Backend> backend = makeCpuBackend();
  hipdlp::backend::ProblemHandle problem;
  hipdlp::backend::WorkspaceHandle workspace;
  REQUIRE(backend->bindProblem(built.problem, problem).ok());
  REQUIRE(backend->createWorkspace(problem, workspace).ok());

  hipdlp::backend::VectorHandle current_primal;
  hipdlp::backend::VectorHandle transpose_product;
  hipdlp::backend::VectorHandle next_primal;
  hipdlp::backend::VectorHandle current_dual;
  hipdlp::backend::VectorHandle extrapolated_product;
  hipdlp::backend::VectorHandle next_dual;
  for (hipdlp::backend::VectorHandle* handle :
       {&current_primal, &transpose_product, &next_primal, &current_dual,
        &extrapolated_product, &next_dual})
    REQUIRE(backend->createVector(2, *handle).ok());
  REQUIRE(backend->uploadVector({1.0, 3.0}, current_primal).ok());
  REQUIRE(backend->uploadVector({4.0, -1.0}, transpose_product).ok());
  REQUIRE(backend
              ->fusedProjectedPrimalUpdate(problem, workspace, current_primal,
                                           transpose_product, 0.5, next_primal)
              .ok());
  checkApproxVector(download(*backend, next_primal), {2.0, 3.5});

  REQUIRE(backend->uploadVector({0.5, -1.0}, current_dual).ok());
  REQUIRE(backend->uploadVector({0.0, 3.0}, extrapolated_product).ok());
  REQUIRE(backend
              ->fusedProjectedDualUpdate(problem, workspace, current_dual,
                                         extrapolated_product, 2.0, next_dual)
              .ok());
  checkApproxVector(download(*backend, next_dual), {2.5, 0.0});

  hipdlp::backend::VectorHandle primal_anchor;
  hipdlp::backend::VectorHandle dual_anchor;
  hipdlp::backend::VectorHandle primal_sum;
  hipdlp::backend::VectorHandle dual_sum;
  for (hipdlp::backend::VectorHandle* handle :
       {&primal_anchor, &dual_anchor, &primal_sum, &dual_sum})
    REQUIRE(backend->createVector(2, *handle).ok());
  REQUIRE(backend->uploadVector({0.0, 10.0}, primal_anchor).ok());
  REQUIRE(backend->uploadVector({4.0, 8.0}, dual_anchor).ok());
  REQUIRE(backend
              ->fusedVectorBlend(workspace, primal_anchor, dual_anchor,
                                 next_primal, next_dual, 0.25, primal_sum,
                                 dual_sum)
              .ok());
  checkApproxVector(download(*backend, primal_sum), {1.5, 5.125});
  checkApproxVector(download(*backend, dual_sum), {2.875, 2.0});

  REQUIRE(backend->fill(primal_sum, 1.0).ok());
  REQUIRE(backend->fill(dual_sum, 1.0).ok());
  REQUIRE(backend
              ->fusedAverageUpdate(workspace, next_primal, next_dual, 0.5,
                                   primal_sum, dual_sum)
              .ok());
  checkApproxVector(download(*backend, primal_sum), {2.0, 2.75});
  checkApproxVector(download(*backend, dual_sum), {2.25, 1.0});
}

#ifdef HIGHS_HIPDLP_CUDA_DETECTED
TEST_CASE("hi-pdlp-cuda-sparse-product-differentials",
          "[pdlp][hipdlp][phase5][cuda][primitives][sparse]") {
  struct Fixture {
    HighsLp lp;
    std::vector<double> x;
    std::vector<double> y;
  };
  const std::vector<Fixture> fixtures = {
      {makeEqualityLp(2, 3, {0, 0, 0, 0}, {}, {}),
       {2.0, -1.0, 4.0},
       {3.0, -2.0}},
      {makeEqualityLp(0, 0, {0}, {}, {}), {}, {}},
      {makeEqualityLp(3, 3, {0, 1, 2, 3}, {0, 1, 2}, {2.0, -3.0, 4.0}),
       {1.0, 2.0, -1.0},
       {1.0, 2.0, -1.0}},
      {makeEqualityLp(2, 3, {0, 2, 3, 4}, {0, 1, 1, 0}, {1.0, -1.0, 3.0, 2.0}),
       {2.0, -1.0, 4.0},
       {3.0, -2.0}},
      {makeEqualityLp(2, 2, {0, 1, 2}, {0, 1}, {1e-200, 1e200}),
       {1e200, 1e-200},
       {1e200, 1e-200}}};

  for (std::size_t fixture = 0; fixture < fixtures.size(); ++fixture) {
    CAPTURE(fixture);
    const ProductResult cpu =
        runProducts(hipdlp::backend::Kind::kCpu, fixtures[fixture].lp,
                    fixtures[fixture].x, fixtures[fixture].y);
    const ProductResult cuda =
        runProducts(hipdlp::backend::Kind::kCuda, fixtures[fixture].lp,
                    fixtures[fixture].x, fixtures[fixture].y);
    checkApproxVector(cuda.ax, cpu.ax, 1e-13);
    checkApproxVector(cuda.aty, cpu.aty, 1e-13);
  }
}

TEST_CASE("hi-pdlp-cuda-vector-fused-and-reduction-differentials",
          "[pdlp][hipdlp][phase5][cuda][primitives][vector][reduction]") {
  const PrimitiveResult cpu = runPrimitiveSuite(hipdlp::backend::Kind::kCpu);
  const PrimitiveResult cuda = runPrimitiveSuite(hipdlp::backend::Kind::kCuda);
  checkApproxVector(cuda.vector_update, cpu.vector_update);
  CHECK(cuda.dot == Approx(cpu.dot).epsilon(1e-14));
  CHECK(cuda.norm == Approx(cpu.norm).epsilon(1e-14));
  CHECK(cuda.difference_norm == Approx(cpu.difference_norm).epsilon(1e-14));
  CHECK(cuda.metrics.primal_residual ==
        Approx(cpu.metrics.primal_residual).epsilon(1e-14));
  CHECK(cuda.metrics.dual_residual ==
        Approx(cpu.metrics.dual_residual).epsilon(1e-14));
  CHECK(cuda.metrics.objective_gap == cpu.metrics.objective_gap);
  checkApproxVector(cuda.primal_update, cpu.primal_update);
  checkApproxVector(cuda.dual_update, cpu.dual_update);
  checkApproxVector(cuda.halpern_primal, cpu.halpern_primal);
  checkApproxVector(cuda.halpern_dual, cpu.halpern_dual);
  checkApproxVector(cuda.average_primal, cpu.average_primal);
  checkApproxVector(cuda.average_dual, cpu.average_dual);
}

TEST_CASE("hi-pdlp-cuda-unsupported-and-transactional-resources",
          "[pdlp][hipdlp][phase5][cuda][primitives][resources]") {
  std::unique_ptr<hipdlp::backend::Backend> backend = makeCudaBackend();
  CHECK(backend->requireCapability(hipdlp::backend::Capability::kBatchedVectors)
            .code == hipdlp::api::ErrorCode::kUnsupportedCapability);
  CHECK(
      backend
          ->requireCapability(hipdlp::backend::Capability::kSparseMatrixMatrix)
          .code == hipdlp::api::ErrorCode::kUnsupportedCapability);
  hipdlp::backend::VectorHandle vector;
  CHECK(backend->createBatchedVector(3, 2, vector).code ==
        hipdlp::api::ErrorCode::kUnsupportedCapability);
  CHECK_FALSE(vector.valid());
  hipdlp::backend::ProblemHandle problem;
  CHECK(backend->multiplyBatched(problem, vector, vector).code ==
        hipdlp::api::ErrorCode::kUnsupportedCapability);
  CHECK(backend->createVector(-1, vector).code ==
        hipdlp::api::ErrorCode::kInvalidArgument);
  CHECK_FALSE(vector.valid());
  hipdlp::model::Problem uninitialized;
  CHECK(backend->bindProblem(uninitialized, problem).code ==
        hipdlp::api::ErrorCode::kInvalidArgument);
  CHECK_FALSE(problem.valid());
#ifdef HIGHSINT64
  CHECK(backend->createVector(std::numeric_limits<HighsInt>::max(), vector)
            .code == hipdlp::api::ErrorCode::kAllocationFailure);
  CHECK_FALSE(vector.valid());
#endif
}

TEST_CASE("hi-pdlp-cuda-initialization-failures-are-transactional",
          "[pdlp][hipdlp][phase5][cuda][initialization][resources][failure]") {
  struct FailureSelection {
    hipdlp::backend::cuda::testing::FailurePoint failure_point;
    hipdlp::api::ErrorCode error_code;
  };
  const FailureSelection selections[] = {
      {hipdlp::backend::cuda::testing::FailurePoint::kBackendInitialization,
       hipdlp::api::ErrorCode::kBackendFailure},
      {hipdlp::backend::cuda::testing::FailurePoint::kDeviceAllocation,
       hipdlp::api::ErrorCode::kAllocationFailure}};

  for (const FailureSelection& selection : selections) {
    CAPTURE(static_cast<int>(selection.failure_point));
    hipdlp::backend::cuda::testing::resetResourceAudit();
    hipdlp::backend::cuda::testing::setFailurePoint(selection.failure_point);
    hipdlp::backend::FactoryResult explicit_cuda =
        hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda);
    hipdlp::backend::FactoryResult automatic =
        hipdlp::backend::createBackend(hipdlp::backend::Kind::kAuto);
    hipdlp::backend::cuda::testing::setFailurePoint(
        hipdlp::backend::cuda::testing::FailurePoint::kNone);

    CHECK_FALSE(explicit_cuda.ok());
    CHECK(explicit_cuda.error.code == selection.error_code);
    CHECK_FALSE(explicit_cuda.error.message.empty());
    CHECK(explicit_cuda.backend.get() == nullptr);
    REQUIRE(automatic.ok());
    CHECK(automatic.backend->descriptor().kind == hipdlp::backend::Kind::kCpu);
    const hipdlp::backend::cuda::testing::ResourceAudit audit =
        hipdlp::backend::cuda::testing::resourceAudit();
    CHECK(audit.stream_creations == audit.stream_destructions);
    CHECK(audit.cublas_creations == audit.cublas_destructions);
    CHECK(audit.cusparse_creations == audit.cusparse_destructions);
    CHECK(audit.allocation_calls == audit.free_calls);

    hipdlp::backend::FactoryResult recovered =
        hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda);
    INFO(recovered.error.message);
    REQUIRE(recovered.ok());
  }

  hipdlp::backend::FactoryResult factory =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda);
  REQUIRE(factory.ok());
  hipdlp::backend::cuda::testing::setFailurePoint(
      hipdlp::backend::cuda::testing::FailurePoint::kDeviceAllocation);
  hipdlp::backend::VectorHandle vector;
  const hipdlp::api::Error vector_error =
      factory.backend->createVector(3, vector);
  HighsLp lp = makeEqualityLp(1, 1, {0, 1}, {0}, {1.0});
  hipdlp::model::ProblemBuildResult built = hipdlp::model::buildProblem(lp);
  REQUIRE(built.ok());
  hipdlp::backend::ProblemHandle problem;
  const hipdlp::api::Error problem_error =
      factory.backend->bindProblem(built.problem, problem);
  hipdlp::backend::cuda::testing::setFailurePoint(
      hipdlp::backend::cuda::testing::FailurePoint::kNone);
  CHECK(vector_error.code == hipdlp::api::ErrorCode::kAllocationFailure);
  CHECK_FALSE(vector.valid());
  CHECK(problem_error.code == hipdlp::api::ErrorCode::kAllocationFailure);
  CHECK_FALSE(problem.valid());
}

TEST_CASE("hi-pdlp-cuda-resource-lifetimes-are-balanced",
          "[pdlp][hipdlp][phase5][cuda][resources][ownership]") {
  hipdlp::model::ProblemBuildResult built = hipdlp::model::buildProblem(
      makeEqualityLp(2, 3, {0, 2, 3, 4}, {0, 1, 1, 0}, {1.0, -1.0, 3.0, 2.0}));
  REQUIRE(built.ok());
  hipdlp::backend::cuda::testing::resetResourceAudit();
  for (int repetition = 0; repetition < 4; ++repetition) {
    hipdlp::backend::FactoryResult factory =
        hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda);
    INFO(factory.error.message);
    REQUIRE(factory.ok());
    hipdlp::backend::ProblemHandle problem;
    hipdlp::backend::WorkspaceHandle workspace;
    hipdlp::backend::VectorHandle primal;
    hipdlp::backend::VectorHandle dual;
    REQUIRE(factory.backend->bindProblem(built.problem, problem).ok());
    REQUIRE(factory.backend->createWorkspace(problem, workspace).ok());
    REQUIRE(factory.backend->createVector(3, primal).ok());
    REQUIRE(factory.backend->createVector(2, dual).ok());
  }
  const hipdlp::backend::cuda::testing::ResourceAudit audit =
      hipdlp::backend::cuda::testing::resourceAudit();
  CHECK(audit.stream_creations == 4);
  CHECK(audit.stream_creations == audit.stream_destructions);
  CHECK(audit.cublas_creations == audit.cublas_destructions);
  CHECK(audit.cusparse_creations == audit.cusparse_destructions);
  CHECK(audit.allocation_calls > 0);
  CHECK(audit.allocation_calls == audit.free_calls);
  CHECK(audit.descriptor_creations > 0);
  CHECK(audit.descriptor_creations == audit.descriptor_destructions);
}

TEST_CASE("hi-pdlp-cuda-fixed-trace-iteration-has-no-resource-traffic",
          "[pdlp][hipdlp][phase5][cuda][trace][resources][audit]") {
  HighsLp lp = makeEqualityLp(1, 1, {0, 1}, {0}, {1.0});
  lp.col_cost_ = {-1.0};
  lp.col_lower_ = {0.0};
  lp.col_upper_ = {2.0};
  lp.row_lower_ = {2.0};
  lp.row_upper_ = {2.0};
  hipdlp::model::ProblemBuildResult built = hipdlp::model::buildProblem(lp);
  REQUIRE(built.ok());
  hipdlp::backend::FactoryResult factory =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda);
  REQUIRE(factory.ok());
  hipdlp::backend::Backend& backend = *factory.backend;
  hipdlp::backend::ProblemHandle problem;
  hipdlp::backend::WorkspaceHandle workspace;
  REQUIRE(backend.bindProblem(built.problem, problem).ok());
  REQUIRE(backend.createWorkspace(problem, workspace).ok());

  hipdlp::backend::VectorHandle current_primal;
  hipdlp::backend::VectorHandle current_dual;
  hipdlp::backend::VectorHandle current_ax;
  hipdlp::backend::VectorHandle current_aty;
  hipdlp::backend::VectorHandle primal_anchor;
  hipdlp::backend::VectorHandle dual_anchor;
  hipdlp::backend::VectorHandle primal_candidate;
  hipdlp::backend::VectorHandle dual_candidate;
  hipdlp::backend::VectorHandle reflected_primal;
  hipdlp::backend::VectorHandle reflected_dual;
  hipdlp::backend::VectorHandle reflected_ax;
  hipdlp::backend::VectorHandle next_primal;
  hipdlp::backend::VectorHandle next_dual;
  hipdlp::backend::VectorHandle point_primal;
  hipdlp::backend::VectorHandle point_dual;
  hipdlp::backend::VectorHandle point_ax;
  hipdlp::backend::VectorHandle point_aty;
  for (hipdlp::backend::VectorHandle* vector :
       {&current_primal, &current_dual, &current_ax, &current_aty,
        &primal_anchor, &dual_anchor, &primal_candidate, &dual_candidate,
        &reflected_primal, &reflected_dual, &reflected_ax, &next_primal,
        &next_dual, &point_primal, &point_dual, &point_ax, &point_aty})
    REQUIRE(backend.createVector(1, *vector).ok());
  REQUIRE(backend.uploadVector({0.0}, current_primal).ok());
  REQUIRE(backend.uploadVector({0.0}, current_dual).ok());
  REQUIRE(backend.copy(current_primal, primal_anchor).ok());
  REQUIRE(backend.copy(current_dual, dual_anchor).ok());
  REQUIRE(backend.multiplyTranspose(problem, current_dual, current_aty).ok());

  hipdlp::backend::cuda::testing::resetResourceAudit();
  REQUIRE(backend
              .fusedProjectedPrimalUpdate(problem, workspace, current_primal,
                                          current_aty, 0.5, primal_candidate)
              .ok());
  REQUIRE(backend.copy(primal_candidate, reflected_primal).ok());
  REQUIRE(backend.scale(2.0, reflected_primal).ok());
  REQUIRE(backend.axpy(-1.0, current_primal, reflected_primal).ok());
  REQUIRE(backend.multiply(problem, reflected_primal, reflected_ax).ok());
  REQUIRE(backend
              .fusedProjectedDualUpdate(problem, workspace, current_dual,
                                        reflected_ax, 0.25, dual_candidate)
              .ok());
  REQUIRE(backend.copy(dual_candidate, reflected_dual).ok());
  REQUIRE(backend.scale(2.0, reflected_dual).ok());
  REQUIRE(backend.axpy(-1.0, current_dual, reflected_dual).ok());
  REQUIRE(backend
              .fusedVectorBlend(workspace, primal_anchor, dual_anchor,
                                reflected_primal, reflected_dual, 0.5,
                                next_primal, next_dual)
              .ok());
  REQUIRE(backend.multiply(problem, next_primal, current_ax).ok());
  REQUIRE(backend.multiplyTranspose(problem, next_dual, current_aty).ok());
  REQUIRE(backend.copy(next_primal, point_primal).ok());
  REQUIRE(backend.copy(next_dual, point_dual).ok());
  REQUIRE(backend.copy(current_ax, point_ax).ok());
  REQUIRE(backend.copy(current_aty, point_aty).ok());

  const hipdlp::backend::cuda::testing::ResourceAudit iteration_audit =
      hipdlp::backend::cuda::testing::resourceAudit();
  CHECK(iteration_audit.allocation_calls == 0);
  CHECK(iteration_audit.free_calls == 0);
  CHECK(iteration_audit.descriptor_creations == 0);
  CHECK(iteration_audit.descriptor_destructions == 0);
  CHECK(iteration_audit.host_to_device_transfers == 0);
  CHECK(iteration_audit.device_to_host_transfers == 0);
  CHECK(iteration_audit.device_to_device_transfers > 0);
  CHECK(iteration_audit.graph_capture_calls == 0);
  CHECK(iteration_audit.graph_instantiation_calls == 0);
  CHECK(iteration_audit.graph_launch_calls == 0);
  CHECK(iteration_audit.graph_update_calls == 0);
  CHECK(iteration_audit.graph_capture_primitive_enqueue_calls == 0);
  CHECK(iteration_audit.direct_primitive_enqueue_calls > 0);
  CHECK(iteration_audit.synchronization_calls == 0);
  REQUIRE(backend.synchronize().ok());
}

TEST_CASE("hi-pdlp-cuda-graph-replay-resource-audit",
          "[pdlp][hipdlp][phase5g][cuda][graph][resources][audit]") {
  CudaCommandBufferFixture fixture = makeCudaCommandBufferFixture();
  REQUIRE(recordCudaCommandBuffer(fixture).ok());
  hipdlp::backend::cuda::testing::resetResourceAudit();
  hipdlp::backend::CommandBufferStatistics statistics;
  REQUIRE(fixture.recorder->replay(fixture.command_buffer, statistics).ok());
  CHECK(statistics.replay_calls == 1);
  CHECK(statistics.direct_operation_enqueue_calls == 0);
  CHECK(statistics.allocation_calls == 0);
  CHECK(statistics.descriptor_creation_calls == 0);
  CHECK(statistics.host_to_device_transfer_calls == 0);
  CHECK(statistics.device_to_host_transfer_calls == 0);
  CHECK(statistics.synchronization_calls == 0);

  const hipdlp::backend::cuda::testing::ResourceAudit audit =
      hipdlp::backend::cuda::testing::resourceAudit();
  CHECK(audit.graph_capture_calls == 0);
  CHECK(audit.graph_instantiation_calls == 0);
  CHECK(audit.graph_launch_calls == 1);
  CHECK(audit.graph_update_calls == 0);
  CHECK(audit.allocation_calls == 0);
  CHECK(audit.free_calls == 0);
  CHECK(audit.descriptor_creations == 0);
  CHECK(audit.descriptor_destructions == 0);
  CHECK(audit.host_to_device_transfers == 0);
  CHECK(audit.device_to_host_transfers == 0);
  CHECK(audit.device_to_device_transfers == 0);
  CHECK(audit.direct_primitive_enqueue_calls == 0);
  CHECK(audit.synchronization_calls == 0);
  REQUIRE(fixture.backend->synchronize().ok());
  CHECK(download(*fixture.backend, fixture.destination) ==
        std::vector<double>({6.0}));
}

TEST_CASE("hi-pdlp-cuda-graph-resource-lifetimes-are-balanced",
          "[pdlp][hipdlp][phase5g][cuda][graph][resources][ownership]") {
  hipdlp::backend::cuda::testing::resetResourceAudit();
  for (int repetition = 0; repetition < 4; ++repetition) {
    CudaCommandBufferFixture fixture = makeCudaCommandBufferFixture();
    REQUIRE(recordCudaCommandBuffer(fixture).ok());
    hipdlp::backend::CommandBufferStatistics statistics;
    REQUIRE(fixture.recorder->replay(fixture.command_buffer, statistics).ok());
    REQUIRE(fixture.backend->synchronize().ok());
  }
  const hipdlp::backend::cuda::testing::ResourceAudit audit =
      hipdlp::backend::cuda::testing::resourceAudit();
  CHECK(audit.stream_creations == 4);
  CHECK(audit.stream_creations == audit.stream_destructions);
  CHECK(audit.cublas_creations == audit.cublas_destructions);
  CHECK(audit.cusparse_creations == audit.cusparse_destructions);
  CHECK(audit.allocation_calls > 0);
  CHECK(audit.allocation_calls == audit.free_calls);
  CHECK(audit.descriptor_creations > 0);
  CHECK(audit.descriptor_creations == audit.descriptor_destructions);
  CHECK(audit.graph_capture_calls == 4);
  CHECK(audit.graph_instantiation_calls == 4);
  CHECK(audit.graph_creations == 4);
  CHECK(audit.graph_creations == audit.graph_destructions);
  CHECK(audit.graph_exec_creations == 4);
  CHECK(audit.graph_exec_creations == audit.graph_exec_destructions);
  CHECK(audit.graph_launch_calls == 4);
  CHECK(audit.graph_update_calls == 0);
  CHECK(audit.graph_capture_primitive_enqueue_calls > 0);
}
#endif