#include <sys/stat.h>

#include <algorithm>
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
#include "pdlp/hipdlp/algorithm/halpern_pdhg_step.hpp"
#include "pdlp/hipdlp/backend/backend_factory.hpp"
#include "pdlp/hipdlp/backend/check.hpp"
#include "pdlp/hipdlp/backend/command_buffer.hpp"
#include "pdlp/hipdlp/backend/cpu/cpu_backend.hpp"
#include "pdlp/hipdlp/backend/cuda/certificate_screen.hpp"
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
#include "pdlp/hipdlp/backend/cuda/cuda_backend_testing.hpp"
#endif
#include "pdlp/hipdlp/model/problem.hpp"
#include "pdlp/hipdlp/model/scaling.hpp"
#include "pdlp/hipdlp/verification/certificate_verifier.hpp"

namespace {

TEST_CASE("hi-pdlp-step-only-computes-update-products",
          "[hipdlp][backend][step-fusion][cpu]") {
  class CountingBackend : public hipdlp::backend::cpu::CpuBackend {
   public:
    int forward_products = 0;
    int transpose_products = 0;
    int reflection_passes = 0;

    hipdlp::api::Error copy(const hipdlp::backend::VectorHandle& input,
                           hipdlp::backend::VectorHandle& output) override {
      ++reflection_passes;
      return CpuBackend::copy(input, output);
    }

    hipdlp::api::Error scale(double factor,
                            hipdlp::backend::VectorHandle& output) override {
      ++reflection_passes;
      return CpuBackend::scale(factor, output);
    }

    hipdlp::api::Error axpy(double factor,
                           const hipdlp::backend::VectorHandle& input,
                           hipdlp::backend::VectorHandle& output) override {
      ++reflection_passes;
      return CpuBackend::axpy(factor, input, output);
    }

    hipdlp::api::Error multiply(
        const hipdlp::backend::ProblemHandle& problem,
        const hipdlp::backend::VectorHandle& input,
        hipdlp::backend::VectorHandle& output) override {
      ++forward_products;
      return CpuBackend::multiply(problem, input, output);
    }

    hipdlp::api::Error multiplyTranspose(
        const hipdlp::backend::ProblemHandle& problem,
        const hipdlp::backend::VectorHandle& input,
        hipdlp::backend::VectorHandle& output) override {
      ++transpose_products;
      return CpuBackend::multiplyTranspose(problem, input, output);
    }
  } backend;
  HighsLp lp;
  lp.num_col_ = 1;
  lp.num_row_ = 1;
  lp.col_cost_ = {-2.0};
  lp.col_lower_ = {0.0};
  lp.col_upper_ = {2.0};
  lp.row_lower_ = {1.0};
  lp.row_upper_ = {kHighsInf};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 1;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1};
  lp.a_matrix_.index_ = {0};
  lp.a_matrix_.value_ = {1.0};
  auto built = hipdlp::model::buildProblem(lp);
  REQUIRE(built.ok());
  hipdlp::backend::ProblemHandle problem;
  hipdlp::backend::WorkspaceHandle workspace;
  REQUIRE(backend.bindProblem(built.problem, problem).ok());
  REQUIRE(backend.createWorkspace(problem, workspace).ok());
  hipdlp::backend::VectorHandle vectors[12];
  for (auto& vector : vectors) REQUIRE(backend.createVector(1, vector).ok());
  hipdlp::algorithm::HalpernPdhgState state;
  state.problem = &problem;
  state.workspace = &workspace;
  state.current_primal = &vectors[0];
  state.current_dual = &vectors[1];
  state.next_primal = &vectors[2];
  state.next_dual = &vectors[3];
  state.primal_anchor = &vectors[4];
  state.dual_anchor = &vectors[5];
  state.primal_candidate = &vectors[6];
  state.dual_candidate = &vectors[7];
  state.reflected_primal = &vectors[8];
  state.reflected_dual = &vectors[9];
  state.current_aty = &vectors[10];
  state.reflected_ax = &vectors[11];
  state.primal_step_size = 0.5;
  state.dual_step_size = 0.25;
  REQUIRE(hipdlp::algorithm::enqueueHalpernPdhgStep(backend, state).ok());
  CHECK(backend.forward_products == 1);
  CHECK(backend.transpose_products == 1);
  CHECK(backend.reflection_passes == 0);
  std::vector<double> values;
  REQUIRE(backend.downloadVector(*state.primal_candidate, values).ok());
  CHECK(values == std::vector<double>{1.0});
  REQUIRE(backend.downloadVector(*state.reflected_primal, values).ok());
  CHECK(values == std::vector<double>{2.0});
  REQUIRE(backend.downloadVector(*state.current_primal, values).ok());
  CHECK(values == std::vector<double>{1.0});
  REQUIRE(backend.downloadVector(*state.dual_candidate, values).ok());
  CHECK(values == std::vector<double>{0.0});
  REQUIRE(backend.downloadVector(*state.current_dual, values).ok());
  CHECK(values == std::vector<double>{0.0});
}

void checkBackendTiming(const hipdlp::backend::Kind kind) {
  hipdlp::backend::FactoryResult factory = hipdlp::backend::createBackend(kind);
  INFO(factory.error.message);
  REQUIRE(factory.ok());
  hipdlp::backend::Backend& backend = *factory.backend;
  CHECK(backend.timingUsesEvents() == (kind == hipdlp::backend::Kind::kCuda));
  double seconds = -7.0;
  CHECK(backend.endTiming(seconds).code ==
    hipdlp::api::ErrorCode::kInvalidArgument);
  CHECK(seconds == -7.0);
  CHECK(backend.beginTiming().code == hipdlp::api::ErrorCode::kInvalidArgument);
  REQUIRE(backend.prepareTiming().ok());
  REQUIRE(backend.prepareTiming().ok());
  hipdlp::backend::VectorHandle vector;
  REQUIRE(backend.createVector(4096, vector).ok());
  std::vector<double> values;
  for (int repetition = 0; repetition < 6; ++repetition) {
    REQUIRE(backend.beginTiming().ok());
    CHECK(backend.beginTiming().code == hipdlp::api::ErrorCode::kInvalidArgument);
    CHECK(backend.prepareTiming().code ==
      hipdlp::api::ErrorCode::kInvalidArgument);
    REQUIRE(backend.fill(vector, 3.0).ok());
    REQUIRE(backend.scale(2.0, vector).ok());
    REQUIRE(backend.endTiming(seconds).ok());
    CHECK(std::isfinite(seconds));
    CHECK(seconds >= 0.0);
    const double completed_seconds = seconds;
    CHECK(backend.endTiming(seconds).code ==
      hipdlp::api::ErrorCode::kInvalidArgument);
    CHECK(seconds == completed_seconds);
    REQUIRE(backend.downloadVector(vector, values).ok());
    REQUIRE(values.size() == 4096);
    CHECK(std::all_of(values.begin(), values.end(),
          [](double value) { return value == 6.0; }));
  }
}

TEST_CASE("hi-pdlp-cpu-opt-in-timing-lifecycle",
      "[hipdlp][backend][timing][cpu]") {
  checkBackendTiming(hipdlp::backend::Kind::kCpu);
}

TEST_CASE("hi-pdlp-cpu-timing-completion-failure-preserves-output",
      "[hipdlp][backend][timing][cpu]") {
  class CompletionFailureBackend : public hipdlp::backend::cpu::CpuBackend {
   public:
    bool fail_completion = false;
    hipdlp::api::Error synchronize() override {
  if (fail_completion)
    return hipdlp::api::Error(hipdlp::api::ErrorCode::kBackendFailure,
              "test completion failure");
  return CpuBackend::synchronize();
    }
  } backend;
  double seconds = 19.0;
  REQUIRE(backend.prepareTiming().ok());
  backend.fail_completion = true;
  CHECK(backend.beginTiming().code == hipdlp::api::ErrorCode::kBackendFailure);
  CHECK(backend.endTiming(seconds).code ==
    hipdlp::api::ErrorCode::kInvalidArgument);
  CHECK(seconds == 19.0);
  backend.fail_completion = false;
  REQUIRE(backend.beginTiming().ok());
  backend.fail_completion = true;
  CHECK(backend.endTiming(seconds).code ==
    hipdlp::api::ErrorCode::kBackendFailure);
  CHECK(seconds == 19.0);
  backend.fail_completion = false;
  CHECK(backend.endTiming(seconds).code ==
    hipdlp::api::ErrorCode::kInvalidArgument);
  REQUIRE(backend.beginTiming().ok());
  REQUIRE(backend.endTiming(seconds).ok());
  CHECK(std::isfinite(seconds));
  CHECK(seconds >= 0.0);
}

TEST_CASE("hi-pdlp-cpu-timers-are-instance-local",
      "[hipdlp][backend][timing][cpu]") {
  auto first = hipdlp::backend::createBackend(hipdlp::backend::Kind::kCpu);
  auto second = hipdlp::backend::createBackend(hipdlp::backend::Kind::kCpu);
  REQUIRE(first.ok());
  REQUIRE(second.ok());
  REQUIRE(first.backend->prepareTiming().ok());
  REQUIRE(second.backend->prepareTiming().ok());
  REQUIRE(first.backend->beginTiming().ok());
  REQUIRE(second.backend->beginTiming().ok());
  double seconds = -1.0;
  REQUIRE(first.backend->endTiming(seconds).ok());
  CHECK(seconds >= 0.0);
  REQUIRE(second.backend->endTiming(seconds).ok());
  CHECK(seconds >= 0.0);
}

#ifdef HIGHS_HIPDLP_CUDA_DETECTED
TEST_CASE("hi-pdlp-cuda-opt-in-timing-lifecycle",
      "[hipdlp][backend][timing][cuda]") {
  checkBackendTiming(hipdlp::backend::Kind::kCuda);
}

TEST_CASE("hi-pdlp-cuda-timing-and-recording-are-exclusive",
      "[hipdlp][backend][timing][cuda]") {
  auto factory = hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda);
  REQUIRE(factory.ok());
  auto* recorder = dynamic_cast<hipdlp::backend::CommandBufferRecorder*>(
  factory.backend.get());
  REQUIRE(recorder != nullptr);
  REQUIRE(factory.backend->prepareTiming().ok());
  REQUIRE(factory.backend->beginTiming().ok());
  CHECK(recorder->beginRecording().code ==
    hipdlp::api::ErrorCode::kInvalidArgument);
  double seconds = -13.0;
  REQUIRE(factory.backend->endTiming(seconds).ok());
  CHECK(seconds >= 0.0);
  REQUIRE(recorder->beginRecording().ok());
  CHECK(factory.backend->prepareTiming().code ==
    hipdlp::api::ErrorCode::kInvalidArgument);
  CHECK(factory.backend->beginTiming().code ==
    hipdlp::api::ErrorCode::kInvalidArgument);
  seconds = -13.0;
  CHECK(factory.backend->endTiming(seconds).code ==
    hipdlp::api::ErrorCode::kInvalidArgument);
  CHECK(seconds == -13.0);
  REQUIRE(recorder->cancelRecording().ok());
}
#endif

struct VectorDownloadFixture {
  std::unique_ptr<hipdlp::backend::Backend> backend;
  hipdlp::backend::VectorHandle vectors[4];
  std::vector<double> values[4];
  hipdlp::backend::VectorDownload downloads[4];
  std::vector<double> saved_values[4];
  const double* saved_data[4];
  std::size_t saved_capacity[4];

  explicit VectorDownloadFixture(const hipdlp::backend::Kind kind) {
    hipdlp::backend::FactoryResult factory =
        hipdlp::backend::createBackend(kind);
    INFO(factory.error.message);
    REQUIRE(factory.ok());
    backend = std::move(factory.backend);
    const HighsInt sizes[] = {3, 2, 2, 3};
    for (std::size_t index = 0; index < 4; ++index) {
      REQUIRE(backend->createVector(sizes[index], vectors[index]).ok());
      REQUIRE(backend->fill(vectors[index], double(index + 1)).ok());
      values[index].assign(1, -double(index + 1));
      downloads[index] = {&vectors[index], &values[index]};
    }
    rememberOutputs();
  }

  void rememberOutputs() {
    for (std::size_t index = 0; index < 4; ++index) {
      saved_values[index] = values[index];
      saved_data[index] = values[index].data();
      saved_capacity[index] = values[index].capacity();
    }
  }

  void checkStorage() const {
    for (std::size_t index = 0; index < 4; ++index) {
      CHECK(values[index].data() == saved_data[index]);
      CHECK(values[index].capacity() == saved_capacity[index]);
    }
  }

  void checkUnchanged() const {
    checkStorage();
    for (std::size_t index = 0; index < 4; ++index)
      CHECK(values[index] == saved_values[index]);
  }

  void checkDownloaded() const {
    CHECK(values[0] == std::vector<double>({1.0, 1.0, 1.0}));
    CHECK(values[1] == std::vector<double>({2.0, 2.0}));
    CHECK(values[2] == std::vector<double>({3.0, 3.0}));
    CHECK(values[3] == std::vector<double>({4.0, 4.0, 4.0}));
  }

  void prepare() {
    REQUIRE(backend->prepareVectorDownloads(downloads, 4).ok());
    for (std::size_t index = 0; index < 4; ++index)
      CHECK(values[index] == saved_values[index]);
    rememberOutputs();
  }
};

void checkRepeatedVectorDownloads(const hipdlp::backend::Kind kind) {
  VectorDownloadFixture fixture(kind);
  fixture.prepare();
  for (int repetition = 0; repetition < 8; ++repetition) {
    REQUIRE(fixture.backend->prepareVectorDownloads(fixture.downloads, 4).ok());
    REQUIRE(fixture.backend->downloadVectors(fixture.downloads, 4).ok());
    fixture.checkDownloaded();
    fixture.checkStorage();
  }
  REQUIRE(fixture.backend->downloadVector(fixture.vectors[0], fixture.values[0])
              .ok());
  REQUIRE(fixture.backend
              ->downloadSolution(fixture.vectors[0], fixture.vectors[1],
                                 fixture.values[0], fixture.values[1])
              .ok());
  fixture.checkDownloaded();
  fixture.checkStorage();
  REQUIRE(fixture.backend->fill(fixture.vectors[0], 9.0).ok());
  REQUIRE(fixture.backend->downloadVectors(fixture.downloads, 4).ok());
  CHECK(fixture.values[0] == std::vector<double>({9.0, 9.0, 9.0}));
  fixture.checkStorage();
}

void checkInvalidVectorDownloads(const hipdlp::backend::Kind kind) {
  VectorDownloadFixture fixture(kind);
  hipdlp::backend::FactoryResult foreign = hipdlp::backend::createBackend(kind);
  REQUIRE(foreign.ok());
  hipdlp::backend::VectorHandle foreign_vector;
  hipdlp::backend::VectorHandle empty_vector;
  hipdlp::backend::VectorHandle wrong_type(
      std::unique_ptr<hipdlp::backend::Resource>(new hipdlp::backend::Resource()));
  REQUIRE(foreign.backend->createVector(3, foreign_vector).ok());
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
  const auto before = hipdlp::backend::cuda::testing::resourceAudit();
#endif
  for (int invalid = 0; invalid < 7; ++invalid) {
    hipdlp::backend::VectorDownload downloads[4];
    std::copy_n(fixture.downloads, 4, downloads);
    if (invalid == 1) downloads[3].vector = nullptr;
    if (invalid == 2) downloads[3].vector = &empty_vector;
    if (invalid == 3) downloads[3].vector = &foreign_vector;
    if (invalid == 4) downloads[3].values = nullptr;
    if (invalid == 5) downloads[3].values = downloads[0].values;
    if (invalid == 6) downloads[3].vector = &wrong_type;
    for (int repetition = 0; repetition < 2; ++repetition) {
      const auto* request = invalid == 0 ? nullptr : downloads;
      CHECK(fixture.backend->prepareVectorDownloads(request, 4).code ==
            hipdlp::api::ErrorCode::kInvalidArgument);
      fixture.checkUnchanged();
      CHECK(fixture.backend->downloadVectors(request, 4).code ==
            hipdlp::api::ErrorCode::kInvalidArgument);
      fixture.checkUnchanged();
    }
  }
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
  if (kind == hipdlp::backend::Kind::kCuda) {
    const auto after = hipdlp::backend::cuda::testing::resourceAudit();
    CHECK(after.pinned_allocation_calls == before.pinned_allocation_calls);
    CHECK(after.allocation_calls == before.allocation_calls);
    CHECK(after.device_to_host_transfers == before.device_to_host_transfers);
    CHECK(after.device_to_host_bytes == before.device_to_host_bytes);
    CHECK(after.synchronization_calls == before.synchronization_calls);
  }
#endif
  REQUIRE(fixture.backend->downloadVectors(fixture.downloads, 4).ok());
  fixture.checkDownloaded();
}

void checkEmptyVectorDownloads(const hipdlp::backend::Kind kind) {
  VectorDownloadFixture fixture(kind);
  REQUIRE(fixture.backend->prepareVectorDownloads(nullptr, 0).ok());
  REQUIRE(fixture.backend->downloadVectors(nullptr, 0).ok());
  fixture.checkUnchanged();
  hipdlp::backend::VectorHandle empty;
  REQUIRE(fixture.backend->createVector(0, empty).ok());
  std::vector<double> values = {99.0};
  const double* const data = values.data();
  const std::size_t capacity = values.capacity();
  const hipdlp::backend::VectorDownload download = {&empty, &values};
  REQUIRE(fixture.backend->prepareVectorDownloads(&download, 1).ok());
  CHECK(values == std::vector<double>({99.0}));
  for (int repetition = 0; repetition < 3; ++repetition) {
    REQUIRE(fixture.backend->downloadVector(empty, values).ok());
    CHECK(values.empty());
    CHECK(values.data() == data);
    CHECK(values.capacity() == capacity);
  }
  const hipdlp::backend::VectorDownload repeated[] = {
      {&fixture.vectors[0], &fixture.values[0]},
      {&fixture.vectors[0], &fixture.values[1]},
      {&empty, &fixture.values[2]}};
  REQUIRE(fixture.backend->prepareVectorDownloads(repeated, 3).ok());
  fixture.rememberOutputs();
  REQUIRE(fixture.backend->downloadVectors(repeated, 3).ok());
  CHECK(fixture.values[0] == std::vector<double>({1.0, 1.0, 1.0}));
  CHECK(fixture.values[1] == fixture.values[0]);
  CHECK(fixture.values[2].empty());
  CHECK(fixture.values[3] == fixture.saved_values[3]);
  fixture.checkStorage();
}

struct DownloadReservationProbe : hipdlp::backend::Backend {
  static hipdlp::api::Error failSecondReservation(
      const hipdlp::backend::VectorDownload* downloads) {
    DownloadOutputs outputs(downloads, 2);
    hipdlp::api::Error error = outputs.reserve(0, 3);
    if (!error.ok()) return error;
    error = outputs.reserve(1, downloads[1].values->max_size() + 1);
    if (!error.ok()) return error;
    outputs.publishReservations();
    return hipdlp::api::Error();
  }
};

}  // namespace

TEST_CASE("hi-pdlp-cpu-batch-downloads-reuse-prepared-outputs",
          "[pdlp][hipdlp][phase8][backend][downloads][cpu]") {
  checkRepeatedVectorDownloads(hipdlp::backend::Kind::kCpu);
}

TEST_CASE("hi-pdlp-cpu-batch-downloads-validate-before-publication",
          "[pdlp][hipdlp][phase8][backend][downloads][cpu]") {
  checkInvalidVectorDownloads(hipdlp::backend::Kind::kCpu);
}

TEST_CASE("hi-pdlp-cpu-batch-downloads-empty-and-repeated-sources",
          "[pdlp][hipdlp][phase8][backend][downloads][cpu]") {
  checkEmptyVectorDownloads(hipdlp::backend::Kind::kCpu);
}

TEST_CASE("hi-pdlp-download-reservation-failure-is-transactional",
          "[pdlp][hipdlp][phase8][backend][downloads][cpu]") {
  VectorDownloadFixture fixture(hipdlp::backend::Kind::kCpu);
  for (int repetition = 0; repetition < 3; ++repetition) {
    CHECK(DownloadReservationProbe::failSecondReservation(fixture.downloads)
              .code == hipdlp::api::ErrorCode::kAllocationFailure);
    fixture.checkUnchanged();
  }
}

TEST_CASE("hi-pdlp-default-batch-download-is-transactional",
          "[pdlp][hipdlp][phase8][backend][downloads][cpu]") {
  using hipdlp::backend::Backend;
  VectorDownloadFixture fixture(hipdlp::backend::Kind::kCpu);
  hipdlp::backend::FactoryResult foreign =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kCpu);
  REQUIRE(foreign.ok());
  hipdlp::backend::VectorHandle foreign_vector;
  REQUIRE(foreign.backend->createVector(3, foreign_vector).ok());
  fixture.downloads[3].vector = &foreign_vector;
  CHECK(fixture.backend->Backend::downloadVectors(fixture.downloads, 4).code ==
        hipdlp::api::ErrorCode::kInvalidArgument);
  fixture.checkUnchanged();
  fixture.downloads[3].vector = &fixture.vectors[3];
  REQUIRE(fixture.backend->Backend::downloadVectors(fixture.downloads, 4).ok());
  fixture.checkDownloaded();
}

#ifdef HIGHS_HIPDLP_CUDA_DETECTED
TEST_CASE("hi-pdlp-cuda-batch-downloads-reuse-prepared-outputs",
          "[pdlp][hipdlp][phase8][backend][downloads][cuda]") {
  checkRepeatedVectorDownloads(hipdlp::backend::Kind::kCuda);
}

TEST_CASE("hi-pdlp-cuda-batch-downloads-validate-before-publication",
          "[pdlp][hipdlp][phase8][backend][downloads][cuda]") {
  checkInvalidVectorDownloads(hipdlp::backend::Kind::kCuda);
}

TEST_CASE("hi-pdlp-cuda-batch-downloads-empty-and-repeated-sources",
          "[pdlp][hipdlp][phase8][backend][downloads][cuda]") {
  checkEmptyVectorDownloads(hipdlp::backend::Kind::kCuda);
}

TEST_CASE("hi-pdlp-cuda-batch-downloads-one-wait-and-real-bytes",
      "[pdlp][hipdlp][phase8][backend][downloads][cuda]") {
  using hipdlp::backend::cuda::testing::resourceAudit;
  const auto initial = resourceAudit();
  {
    VectorDownloadFixture fixture(hipdlp::backend::Kind::kCuda);
    CHECK(resourceAudit().pinned_allocation_calls ==
      initial.pinned_allocation_calls);
    fixture.prepare();
    const auto prepared = resourceAudit();
    CHECK(prepared.pinned_allocation_calls - initial.pinned_allocation_calls == 5);
    CHECK(prepared.device_to_host_transfers == initial.device_to_host_transfers);
    for (int repetition = 0; repetition < 8; ++repetition) {
  REQUIRE(fixture.backend->prepareVectorDownloads(fixture.downloads, 4).ok());
  REQUIRE(fixture.backend->downloadVectors(fixture.downloads, 4).ok());
  fixture.checkDownloaded();
  fixture.checkStorage();
    }
    const auto repeated = resourceAudit();
    CHECK(repeated.device_to_host_transfers - prepared.device_to_host_transfers ==
      8 * 5);
    CHECK(repeated.device_to_host_bytes - prepared.device_to_host_bytes ==
      8 * (10 * sizeof(double) + sizeof(int)));
    CHECK(repeated.synchronization_calls - prepared.synchronization_calls == 8);
    CHECK(repeated.pinned_allocation_calls == prepared.pinned_allocation_calls);
    CHECK(repeated.pinned_free_calls == prepared.pinned_free_calls);
    CHECK(repeated.allocation_calls == prepared.allocation_calls);
    CHECK(repeated.free_calls == prepared.free_calls);
    CHECK(repeated.host_to_device_transfers == prepared.host_to_device_transfers);
    CHECK(repeated.descriptor_creations == prepared.descriptor_creations);
    CHECK(repeated.stream_creations == prepared.stream_creations);
    CHECK(repeated.cublas_creations == prepared.cublas_creations);
    CHECK(repeated.cusparse_creations == prepared.cusparse_creations);
    CHECK(repeated.direct_primitive_enqueue_calls ==
      prepared.direct_primitive_enqueue_calls);
    REQUIRE(fixture.backend->downloadVector(fixture.vectors[0], fixture.values[0])
        .ok());
    const auto single = resourceAudit();
    CHECK(single.device_to_host_transfers - repeated.device_to_host_transfers ==
      2);
    CHECK(single.device_to_host_bytes - repeated.device_to_host_bytes ==
      3 * sizeof(double) + sizeof(int));
    CHECK(single.synchronization_calls - repeated.synchronization_calls == 1);
    REQUIRE(fixture.backend
        ->downloadSolution(fixture.vectors[0], fixture.vectors[1],
               fixture.values[0], fixture.values[1])
        .ok());
    const auto solution = resourceAudit();
    CHECK(solution.device_to_host_transfers - single.device_to_host_transfers ==
      3);
    CHECK(solution.device_to_host_bytes - single.device_to_host_bytes ==
      5 * sizeof(double) + sizeof(int));
    CHECK(solution.synchronization_calls - single.synchronization_calls == 1);
    fixture.checkStorage();

    hipdlp::backend::VectorHandle empty;
    REQUIRE(fixture.backend->createVector(0, empty).ok());
    std::vector<double> empty_values = {99.0};
    const auto before_empty = resourceAudit();
    REQUIRE(fixture.backend->downloadVector(empty, empty_values).ok());
    CHECK(empty_values.empty());
    const auto after_empty = resourceAudit();
    CHECK(after_empty.device_to_host_transfers -
      before_empty.device_to_host_transfers == 1);
    CHECK(after_empty.device_to_host_bytes - before_empty.device_to_host_bytes ==
      sizeof(int));
    CHECK(after_empty.synchronization_calls - before_empty.synchronization_calls ==
      1);
    CHECK(after_empty.pinned_allocation_calls ==
      before_empty.pinned_allocation_calls);
    REQUIRE(fixture.backend->prepareVectorDownloads(nullptr, 0).ok());
    REQUIRE(fixture.backend->downloadVectors(nullptr, 0).ok());
    const auto empty_batch = resourceAudit();
    CHECK(empty_batch.device_to_host_transfers ==
      after_empty.device_to_host_transfers);
    CHECK(empty_batch.device_to_host_bytes == after_empty.device_to_host_bytes);
    CHECK(empty_batch.synchronization_calls == after_empty.synchronization_calls);
    CHECK(empty_batch.pinned_allocation_calls ==
      after_empty.pinned_allocation_calls);
  }
  const auto released = resourceAudit();
  CHECK(released.pinned_allocation_calls - initial.pinned_allocation_calls ==
    released.pinned_free_calls - initial.pinned_free_calls);
  CHECK(released.allocation_calls - initial.allocation_calls ==
    released.free_calls - initial.free_calls);
}

TEST_CASE("hi-pdlp-cuda-batch-downloads-preserve-outputs-on-sticky-error",
      "[pdlp][hipdlp][phase8][backend][downloads][cuda]") {
  using hipdlp::backend::cuda::testing::resourceAudit;
  VectorDownloadFixture fixture(hipdlp::backend::Kind::kCuda);
  fixture.prepare();
  REQUIRE(fixture.backend
      ->fill(fixture.vectors[3], std::numeric_limits<double>::max())
      .ok());
  REQUIRE(fixture.backend->scale(2.0, fixture.vectors[3]).ok());
  REQUIRE(fixture.backend->fill(fixture.vectors[3], 4.0).ok());
  const auto before = resourceAudit();
  for (int repetition = 0; repetition < 3; ++repetition) {
    CHECK(fixture.backend->downloadVectors(fixture.downloads, 4).code ==
      hipdlp::api::ErrorCode::kBackendFailure);
    fixture.checkUnchanged();
  }
  const auto after = resourceAudit();
  CHECK(after.device_to_host_transfers - before.device_to_host_transfers == 3 * 5);
  CHECK(after.device_to_host_bytes - before.device_to_host_bytes ==
    3 * (10 * sizeof(double) + sizeof(int)));
  CHECK(after.synchronization_calls - before.synchronization_calls == 3);
  CHECK(after.pinned_allocation_calls == before.pinned_allocation_calls);
  CHECK(after.pinned_free_calls == before.pinned_free_calls);

  std::vector<double> cold_values[4];
  hipdlp::backend::VectorDownload cold_downloads[4];
  const double* cold_data[4];
  std::size_t cold_capacity[4];
  for (std::size_t index = 0; index < 4; ++index) {
    cold_values[index].assign(1, -99.0);
    cold_data[index] = cold_values[index].data();
    cold_capacity[index] = cold_values[index].capacity();
    cold_downloads[index] = {&fixture.vectors[index], &cold_values[index]};
  }
  for (int repetition = 0; repetition < 3; ++repetition) {
    CHECK(fixture.backend->downloadVectors(cold_downloads, 4).code ==
      hipdlp::api::ErrorCode::kBackendFailure);
    for (std::size_t index = 0; index < 4; ++index) {
  CHECK(cold_values[index] == std::vector<double>({-99.0}));
  CHECK(cold_values[index].data() == cold_data[index]);
  CHECK(cold_values[index].capacity() == cold_capacity[index]);
    }
  }
}

TEST_CASE("hi-pdlp-cuda-download-preparation-failure-and-reuse",
      "[pdlp][hipdlp][phase8][backend][downloads][cuda]") {
  using namespace hipdlp::backend::cuda::testing;
  const auto initial = resourceAudit();
  {
    VectorDownloadFixture fixture(hipdlp::backend::Kind::kCuda);
    const auto before = resourceAudit();
    setFailurePoint(FailurePoint::kPinnedAllocation);
    const auto failed =
    fixture.backend->prepareVectorDownloads(fixture.downloads, 4);
    setFailurePoint(FailurePoint::kNone);
    CHECK(failed.code == hipdlp::api::ErrorCode::kAllocationFailure);
    fixture.checkUnchanged();
    CHECK(resourceAudit().pinned_allocation_calls == before.pinned_allocation_calls);
    CHECK(resourceAudit().device_to_host_transfers == before.device_to_host_transfers);
    REQUIRE(fixture.backend->prepareVectorDownloads(fixture.downloads, 1).ok());
    fixture.rememberOutputs();
    const auto partially_prepared = resourceAudit();
    setFailurePoint(FailurePoint::kPinnedAllocation);
    const auto failed_later =
    fixture.backend->prepareVectorDownloads(fixture.downloads, 4);
    setFailurePoint(FailurePoint::kNone);
    CHECK(failed_later.code == hipdlp::api::ErrorCode::kAllocationFailure);
    fixture.checkUnchanged();
    CHECK(resourceAudit().pinned_allocation_calls ==
      partially_prepared.pinned_allocation_calls);
    fixture.prepare();
    const auto prepared = resourceAudit();
    setFailurePoint(FailurePoint::kPinnedAllocation);
    const auto repeated =
    fixture.backend->prepareVectorDownloads(fixture.downloads, 4);
    const auto downloaded =
    fixture.backend->downloadVectors(fixture.downloads, 4);
    setFailurePoint(FailurePoint::kNone);
    REQUIRE(repeated.ok());
    REQUIRE(downloaded.ok());
    fixture.checkDownloaded();
    fixture.checkStorage();
    CHECK(resourceAudit().pinned_allocation_calls == prepared.pinned_allocation_calls);
  }
  const auto released = resourceAudit();
  CHECK(released.pinned_allocation_calls - initial.pinned_allocation_calls ==
    released.pinned_free_calls - initial.pinned_free_calls);
}

TEST_CASE("hi-pdlp-cuda-download-staging-follows-vector-lifetime",
      "[pdlp][hipdlp][phase8][backend][downloads][cuda]") {
  using hipdlp::backend::cuda::testing::resourceAudit;
  const auto initial = resourceAudit();
  hipdlp::backend::VectorHandle survivor;
  {
    VectorDownloadFixture fixture(hipdlp::backend::Kind::kCuda);
    fixture.prepare();
    survivor = std::move(fixture.vectors[0]);
    CHECK_FALSE(fixture.vectors[0].valid());
    CHECK(fixture.backend->downloadVectors(fixture.downloads, 4).code ==
      hipdlp::api::ErrorCode::kInvalidArgument);
    fixture.checkUnchanged();
    fixture.downloads[0].vector = &survivor;
    const auto moved = resourceAudit();
    REQUIRE(fixture.backend->downloadVectors(fixture.downloads, 4).ok());
    fixture.checkDownloaded();
    fixture.checkStorage();
    CHECK(resourceAudit().pinned_allocation_calls == moved.pinned_allocation_calls);
    fixture.backend.reset();
  }
  const auto retained = resourceAudit();
  CHECK((retained.pinned_allocation_calls - initial.pinned_allocation_calls) -
    (retained.pinned_free_calls - initial.pinned_free_calls) == 1);
  survivor.reset();
  const auto released = resourceAudit();
  CHECK(released.pinned_allocation_calls - initial.pinned_allocation_calls ==
    released.pinned_free_calls - initial.pinned_free_calls);
  CHECK(released.allocation_calls - initial.allocation_calls ==
    released.free_calls - initial.free_calls);
  CHECK(released.descriptor_creations - initial.descriptor_creations ==
    released.descriptor_destructions - initial.descriptor_destructions);
}

TEST_CASE("hi-pdlp-cuda-batch-downloads-reject-recording-before-allocation",
      "[pdlp][hipdlp][phase8][backend][downloads][cuda]") {
  using hipdlp::backend::cuda::testing::resourceAudit;
  VectorDownloadFixture fixture(hipdlp::backend::Kind::kCuda);
  auto* recorder =
  dynamic_cast<hipdlp::backend::CommandBufferRecorder*>(fixture.backend.get());
  REQUIRE(recorder != nullptr);
  REQUIRE(recorder->beginRecording().ok());
  const auto before = resourceAudit();
  const auto preparation =
  fixture.backend->prepareVectorDownloads(fixture.downloads, 4);
  const auto download = fixture.backend->downloadVectors(fixture.downloads, 4);
  const auto after = resourceAudit();
  REQUIRE(recorder->cancelRecording().ok());
  CHECK(preparation.code == hipdlp::api::ErrorCode::kInvalidArgument);
  CHECK(download.code == hipdlp::api::ErrorCode::kInvalidArgument);
  fixture.checkUnchanged();
  CHECK(after.pinned_allocation_calls == before.pinned_allocation_calls);
  CHECK(after.allocation_calls == before.allocation_calls);
  CHECK(after.device_to_host_transfers == before.device_to_host_transfers);
  CHECK(after.synchronization_calls == before.synchronization_calls);
}
#endif

TEST_CASE("hi-pdlp-cpu-solution-download-is-transactional",
          "[pdlp][hipdlp][phase8][backend][downloads][cpu]") {
  hipdlp::backend::FactoryResult factory =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kCpu);
  REQUIRE(factory.ok());
  hipdlp::backend::VectorHandle primal;
  hipdlp::backend::VectorHandle dual;
  REQUIRE(factory.backend->createVector(3, primal).ok());
  REQUIRE(factory.backend->uploadVector({1.0, 2.0, 3.0}, primal).ok());
  std::vector<double> primal_values = {-11.0, -12.0};
  std::vector<double> dual_values = {-21.0};
  primal_values.reserve(8);
  dual_values.reserve(8);
  const double* const primal_data = primal_values.data();
  const double* const dual_data = dual_values.data();
  const std::size_t primal_capacity = primal_values.capacity();
  const std::size_t dual_capacity = dual_values.capacity();
  for (int repetition = 0; repetition < 3; ++repetition) {
    const hipdlp::api::Error error = factory.backend->downloadSolution(
        primal, dual, primal_values, dual_values);
    CHECK(error.code == hipdlp::api::ErrorCode::kInvalidArgument);
    CHECK(primal_values == std::vector<double>({-11.0, -12.0}));
    CHECK(dual_values == std::vector<double>({-21.0}));
    CHECK(primal_values.data() == primal_data);
    CHECK(dual_values.data() == dual_data);
    CHECK(primal_values.capacity() == primal_capacity);
    CHECK(dual_values.capacity() == dual_capacity);
  }
}

TEST_CASE(
    "hi-pdlp-certificate-screen-defaults-allow-host-verification",
    "[pdlp][hipdlp][phase8][backend][check-packet][certificate-screen][cpu]") {
  const hipdlp::backend::CheckInput input;
  const hipdlp::backend::CheckResult result;
  CHECK(input.certificate_tolerance == 1e-4);
  CHECK(result.primal_ray_possible);
  CHECK(result.farkas_possible);
}

TEST_CASE(
    "hi-pdlp-certificate-intervals-enclose-both-host-normalizations",
    "[pdlp][hipdlp][phase8][backend][check-packet][certificate-screen][cpu]") {
  using namespace hipdlp::backend::cuda::certificate;
  for (const double scale : {1.0, 1e150, 1e-200}) {
    const std::vector<double> values = {3.0, -4.0, 0.0, 1e-100};
    std::vector<Interval> recovered;
    std::vector<double> host;
    NormBounds norm = {};
    double host_norm = 0.0;
    for (const double value : values) {
      recovered.push_back(recoverProduct(value, scale));
      host.push_back(
          static_cast<double>(static_cast<long double>(value) * scale));
      CHECK(recovered.back().lower <= host.back());
      CHECK(recovered.back().upper >= host.back());
      addNorm(norm, recovered.back());
      host_norm = std::hypot(host_norm, host.back());
    }
    const Interval normalization = hostNormBounds(norm, values.size());
    REQUIRE(positiveFinite(normalization));
    CHECK(normalization.lower <= host_norm);
    CHECK(normalization.upper >= host_norm);
    NormBounds normalized_norm = {};
    double verifier_norm = 0.0;
    for (std::size_t index = 0; index < values.size(); ++index) {
      const Interval normalized = divide(recovered[index], normalization);
      host[index] /= host_norm;
      CHECK(normalized.lower <= host[index]);
      CHECK(normalized.upper >= host[index]);
      addNorm(normalized_norm, normalized);
      verifier_norm = std::hypot(verifier_norm, host[index]);
    }
    const Interval verification =
        hostNormBounds(normalized_norm, values.size());
    REQUIRE(positiveFinite(verification));
    CHECK(verification.lower <= verifier_norm);
    CHECK(verification.upper >= verifier_norm);
  }
}

TEST_CASE(
    "hi-pdlp-certificate-intervals-preserve-tolerance-and-serial-sums",
    "[pdlp][hipdlp][phase8][backend][check-packet][certificate-screen][cpu]") {
  using namespace hipdlp::backend::cuda::certificate;
  const Interval unit = {1.0, 1.0};
  const double tolerance = 1e-4;
  CHECK_FALSE(rejectViolation(tolerance, unit, tolerance));
  CHECK(rejectViolation(2.0 * tolerance, unit, tolerance));
  CHECK_FALSE(rejectViolation(tolerance, {0.0, kHighsInf}, tolerance));
  CHECK(violationLower({-2.0, -1.0}, 0.0, kHighsInf, false) == 1.0);
  CHECK(violationLower({1.0, 2.0}, -1.0, 1.0, false) == 1.0);
  CHECK(violationLower({1.0, 2.0}, -1.0, 1.0, true) == 0.0);
  CHECK(violationLower({-1.0, 1.0}, -kHighsInf, kHighsInf, true) == 0.0);
  SumBounds sum = {};
  double serial = 0.0;
  for (const double value : {1e16, 1.0, -1e16, 0.25}) {
    addSum(sum, {value, value});
    serial += value;
  }
  const Interval enclosure = serialSumBounds(sum, 4.0);
  CHECK(enclosure.lower <= serial);
  CHECK(enclosure.upper >= serial);
  CHECK_FALSE(positiveFinite(hostNormBounds(NormBounds(), 0.0)));
  const double nan = std::numeric_limits<double>::quiet_NaN();
  CHECK_FALSE(finite(recoverProduct(nan, 1.0)));
}

TEST_CASE(
    "hi-pdlp-certificate-screen-decisions-are-necessary-conditions-only",
    "[pdlp][hipdlp][phase8][backend][check-packet][certificate-screen][cpu]") {
  using namespace hipdlp::backend::cuda::certificate;
  ScreenPartial partial = {};
  const Interval unit = {1.0, 1.0};
  partial.primal_violation = 0.5;
  partial.farkas_violation = 0.5;
  ScreenFlags flags =
      finishScreen(partial, unit, unit, 3.0, 2.0, 1.0, true, false, 1e-4);
  CHECK_FALSE(flags.primal_ray_possible);
  CHECK_FALSE(flags.farkas_possible);
  partial.farkas_violation = 0.0;
  addSum(partial.row_support, {-2.0, -2.0});
  flags = finishScreen(partial, unit, unit, 3.0, 2.0, 1.0, true, false, 1e-4);
  CHECK_FALSE(flags.farkas_possible);
  flags = finishScreen(partial, unit, unit, 3.0, 2.0, 1.0, false, false, 1e-4);
  CHECK(flags.farkas_possible);
  for (const double tolerance :
       {0.0, -1.0, kHighsInf, std::numeric_limits<double>::quiet_NaN()}) {
    flags = finishScreen(partial, unit, unit, 3.0, 2.0, 1.0, true, false,
                         tolerance);
    CHECK(flags.primal_ray_possible);
    CHECK(flags.farkas_possible);
  }
  flags = finishScreen(partial, unit, unit, 3.0, 2.0, 1.0, true, true, 1e-4);
  CHECK(flags.primal_ray_possible);
  CHECK(flags.farkas_possible);
  partial = ScreenPartial();
  addSum(partial.objective, {-2.0, -2.0});
  flags = finishScreen(partial, unit, unit, 3.0, 2.0, 1.0, false, false, 1e-4);
  CHECK(flags.primal_ray_possible);
  flags = finishScreen(partial, unit, unit, 3.0, 2.0, -1.0, false, false, 1e-4);
  CHECK_FALSE(flags.primal_ray_possible);
}

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

struct ResidentCheckFixture {
  std::unique_ptr<hipdlp::backend::Backend> backend;
  hipdlp::model::Problem problem;
  hipdlp::backend::ProblemHandle bound;
  hipdlp::backend::WorkspaceHandle workspace;
  hipdlp::backend::VectorHandle primal;
  hipdlp::backend::VectorHandle dual;
  hipdlp::backend::VectorHandle ax;
  hipdlp::backend::VectorHandle aty;
  hipdlp::backend::VectorHandle delta_primal;
  hipdlp::backend::VectorHandle delta_dual;
  hipdlp::backend::VectorHandle aty_delta;
  hipdlp::backend::VectorHandle primal_anchor;
  hipdlp::backend::VectorHandle dual_anchor;
  hipdlp::backend::CheckInput input;

  ResidentCheckFixture(const hipdlp::backend::Kind kind, const HighsLp& lp,
                       const bool scaled = false, const double diagonal = 0.0) {
    hipdlp::backend::FactoryResult factory =
        hipdlp::backend::createBackend(kind);
    INFO(factory.error.message);
    REQUIRE(factory.ok());
    backend = std::move(factory.backend);
    hipdlp::model::ProblemBuildResult built = hipdlp::model::buildProblem(lp);
    INFO(built.error.message);
    REQUIRE(built.ok());
    problem = std::move(built.problem);
    if (scaled) REQUIRE(hipdlp::model::applyCurtisReidScaling(problem, 4).ok());
    if (diagonal > 0.0)
      REQUIRE(problem
                  .applyDiagonalScaling(
                      std::vector<double>(problem.numRows(), 1.0 / diagonal),
                      std::vector<double>(problem.numCols(), diagonal), 0)
                  .ok());
    REQUIRE(backend->bindProblem(problem, bound).ok());
    REQUIRE(backend->createWorkspace(bound, workspace).ok());
    for (hipdlp::backend::VectorHandle* handle :
         {&primal, &aty, &delta_primal, &aty_delta, &primal_anchor})
      REQUIRE(backend->createVector(problem.numCols(), *handle).ok());
    for (hipdlp::backend::VectorHandle* handle :
         {&dual, &ax, &delta_dual, &dual_anchor})
      REQUIRE(backend->createVector(problem.numRows(), *handle).ok());
    input.primal = &primal;
    input.dual = &dual;
    input.ax = &ax;
    input.aty = &aty;
  }

  void prepare(const HighsLp& lp) {
    const hipdlp::api::Error error =
        backend->prepareCheckWorkspace(problem, lp, workspace);
    INFO(error.message);
    REQUIRE(error.ok());
  }
};

hipdlp::backend::cuda::certificate::ScreenFlags hostCertificateScreen(
    const hipdlp::model::Problem& problem, const HighsLp& lp,
    const std::vector<double>& primal, const std::vector<double>& dual,
    const std::vector<double>& aty, const double tolerance) {
  using namespace hipdlp::backend::cuda::certificate;
  std::vector<Interval> recovered_primal(lp.num_col_);
  std::vector<Interval> recovered_dual(lp.num_row_);
  std::vector<Interval> recovered_column(lp.num_col_);
  NormBounds primal_norms[2] = {};
  NormBounds dual_norms[2] = {};
  bool unsafe = false;
  bool straddle_zero = true;
  const auto& transform = problem.transformation();
  for (HighsInt col = 0; col < lp.num_col_; ++col) {
    recovered_primal[col] =
        recoverProduct(primal[col], problem.scaling().col_scale[col]);
    recovered_column[col] =
        recoverQuotient(-aty[col], problem.scaling().col_scale[col]);
    unsafe |= !safeRecovery(recovered_primal[col], primal[col]) ||
              !safeRecovery(recovered_column[col], -aty[col]);
    addNorm(primal_norms[col % 2], recovered_primal[col]);
    addNorm(dual_norms[col % 2], recovered_column[col]);
    straddle_zero &= lp.col_lower_[col] <= 0.0 && lp.col_upper_[col] >= 0.0;
  }
  for (HighsInt row = 0; row < lp.num_row_; ++row) {
    const HighsInt transformed = transform.original_to_transformed_row[row];
    recovered_dual[row] =
        multiply(recoverProduct(dual[transformed],
                                problem.scaling().row_scale[transformed]),
                 transform.original_row_sign[row]);
    unsafe |= !safeRecovery(recovered_dual[row], dual[transformed]);
    addNorm(dual_norms[row % 2], recovered_dual[row]);
  }
  mergeNorm(primal_norms[0], primal_norms[1]);
  mergeNorm(dual_norms[0], dual_norms[1]);
  const Interval primal_normalization =
      hostNormBounds(primal_norms[0], lp.num_col_);
  const Interval dual_normalization =
      hostNormBounds(dual_norms[0], double(lp.num_col_) + lp.num_row_ + 1.0);
  unsafe |= !positiveFinite(primal_normalization) ||
            !positiveFinite(dual_normalization);
  if (unsafe) return {true, true};
  NormBounds normalized_primal[2] = {};
  NormBounds normalized_dual[2] = {};
  for (HighsInt col = 0; col < lp.num_col_; ++col) {
    recovered_primal[col] = divide(recovered_primal[col], primal_normalization);
    recovered_column[col] = divide(recovered_column[col], dual_normalization);
    addNorm(normalized_primal[col % 2], recovered_primal[col]);
    addNorm(normalized_dual[col % 2], recovered_column[col]);
  }
  for (HighsInt row = 0; row < lp.num_row_; ++row) {
    recovered_dual[row] = divide(recovered_dual[row], dual_normalization);
    addNorm(normalized_dual[row % 2], recovered_dual[row]);
  }
  mergeNorm(normalized_primal[0], normalized_primal[1]);
  mergeNorm(normalized_dual[0], normalized_dual[1]);
  const Interval primal_verifier =
      hostNormBounds(normalized_primal[0], lp.num_col_);
  const Interval dual_verifier = hostNormBounds(
      normalized_dual[0], double(lp.num_col_) + lp.num_row_ + 1.0);
  std::vector<Interval> activity(lp.num_row_, Interval{0.0, 0.0});
  const auto& matrix = lp.a_matrix_;
  const bool colwise = matrix.format_ == MatrixFormat::kColwise;
  const HighsInt vectors = colwise ? lp.num_col_ : lp.num_row_;
  for (HighsInt vector = 0; vector < vectors; ++vector) {
    for (HighsInt entry = matrix.start_[vector];
         entry < matrix.start_[vector + 1]; ++entry) {
      const HighsInt row = colwise ? matrix.index_[entry] : vector;
      const HighsInt col = colwise ? vector : matrix.index_[entry];
      activity[row] = add(
          activity[row], multiply(recovered_primal[col], matrix.value_[entry]));
    }
  }
  ScreenPartial partials[2] = {};
  for (HighsInt col = 0; col < lp.num_col_; ++col) {
    ScreenPartial& partial = partials[col % 2];
    partial.primal_violation =
        maximum(partial.primal_violation,
                violationLower(recovered_primal[col], lp.col_lower_[col],
                               lp.col_upper_[col], false));
    partial.farkas_violation =
        maximum(partial.farkas_violation,
                violationLower(recovered_column[col], lp.col_lower_[col],
                               lp.col_upper_[col], true));
    addSum(partial.objective,
           multiply(recovered_primal[col], lp.col_cost_[col]));
  }
  for (HighsInt row = 0; row < lp.num_row_; ++row) {
    ScreenPartial& partial = partials[row % 2];
    partial.primal_violation =
        maximum(partial.primal_violation,
                violationLower(activity[row], lp.row_lower_[row],
                               lp.row_upper_[row], false));
    partial.farkas_violation =
        maximum(partial.farkas_violation,
                violationLower(recovered_dual[row], lp.row_lower_[row],
                               lp.row_upper_[row], true));
    if (straddle_zero)
      addSum(partial.row_support,
             rowSupport(recovered_dual[row], lp.row_lower_[row],
                        lp.row_upper_[row]));
  }
  mergeScreen(partials[0], partials[1]);
  return finishScreen(partials[0], primal_verifier, dual_verifier, lp.num_col_,
                      lp.num_row_,
                      lp.sense_ == ObjSense::kMinimize ? 1.0 : -1.0,
                      straddle_zero, unsafe, tolerance);
}

double certificateHostNorm(const std::vector<double>& values) {
  double norm = 0.0;
  for (const double value : values) norm = std::hypot(norm, value);
  return norm;
}

void checkCertificateCandidate(const hipdlp::backend::Kind kind,
                               const HighsLp& lp,
                               const std::vector<double>& original_primal,
                               const std::vector<double>& original_dual,
                               const std::vector<double>& original_column,
                               const bool scaled, const int expected_ray,
                               const int expected_farkas,
                               const bool exclude_both = false,
                               const double tolerance = 1e-4,
                               const double diagonal = 0.0) {
  ResidentCheckFixture fixture(kind, lp, scaled, diagonal);
  fixture.prepare(lp);
  const auto& transform = fixture.problem.transformation();
  const auto& scaling = fixture.problem.scaling();
  std::vector<double> primal(fixture.problem.numCols(), 0.0);
  std::vector<double> dual(fixture.problem.numRows(), 0.0);
  std::vector<double> aty(fixture.problem.numCols(), 0.0);
  for (HighsInt col = 0; col < lp.num_col_; ++col) {
    primal[col] = original_primal[col] / scaling.col_scale[col];
    aty[col] = -original_column[col] * scaling.col_scale[col];
  }
  for (HighsInt row = 0; row < lp.num_row_; ++row) {
    const HighsInt transformed = transform.original_to_transformed_row[row];
    dual[transformed] = original_dual[row] * transform.original_row_sign[row] /
                        scaling.row_scale[transformed];
  }
  std::vector<double> recovered_primal;
  std::vector<double> recovered_dual;
  std::vector<double> recovered_column;
  std::vector<double> negative_aty(aty.size());
  for (std::size_t col = 0; col < aty.size(); ++col)
    negative_aty[col] = -aty[col];
  REQUIRE(
      fixture.problem.recoverPrimalDirection(primal, recovered_primal).ok());
  REQUIRE(fixture.problem.recoverDualDirection(dual, recovered_dual).ok());
  REQUIRE(fixture.problem
              .recoverReducedCostDirection(negative_aty, recovered_column)
              .ok());
  const double primal_norm = certificateHostNorm(recovered_primal);
  const double dual_norm = std::hypot(certificateHostNorm(recovered_dual),
                                      certificateHostNorm(recovered_column));
  bool ray_pass = false;
  bool farkas_pass = false;
  if (std::isfinite(primal_norm) && primal_norm > 0.0) {
    for (double& value : recovered_primal) value /= primal_norm;
    ray_pass =
        hipdlp::verification::verifyPrimalRay(lp, recovered_primal, tolerance)
            .pass;
  }
  if (std::isfinite(dual_norm) && dual_norm > 0.0) {
    for (double& value : recovered_dual) value /= dual_norm;
    for (double& value : recovered_column) value /= dual_norm;
    farkas_pass = hipdlp::verification::verifyFarkasCertificate(
                      lp, recovered_dual, recovered_column, tolerance)
                      .pass;
  }
  if (expected_ray >= 0) REQUIRE(ray_pass == (expected_ray != 0));
  if (expected_farkas >= 0) REQUIRE(farkas_pass == (expected_farkas != 0));
  const auto screen =
      hostCertificateScreen(fixture.problem, lp, primal, dual, aty, tolerance);
  if (ray_pass) CHECK(screen.primal_ray_possible);
  if (farkas_pass) CHECK(screen.farkas_possible);
  if (exclude_both) {
    CHECK_FALSE(ray_pass);
    CHECK_FALSE(farkas_pass);
    CHECK_FALSE(screen.primal_ray_possible);
    CHECK_FALSE(screen.farkas_possible);
  }
  fixture.input.certificate_tolerance = tolerance;
  REQUIRE(fixture.backend->uploadVector(primal, fixture.primal).ok());
  REQUIRE(fixture.backend->uploadVector(dual, fixture.dual).ok());
  REQUIRE(fixture.backend->uploadVector(aty, fixture.aty).ok());
  REQUIRE(fixture.backend->multiply(fixture.bound, fixture.primal, fixture.ax)
              .ok());
  hipdlp::backend::CheckResult result;
  REQUIRE(fixture.backend
              ->evaluateCheck(fixture.problem, fixture.workspace, fixture.input,
                              result)
              .ok());
  CHECK((result.valid || result.requires_host_oracle));
  if (ray_pass) CHECK(result.primal_ray_possible);
  if (farkas_pass) CHECK(result.farkas_possible);
  if (kind == hipdlp::backend::Kind::kCpu) {
    CHECK(result.primal_ray_possible);
    CHECK(result.farkas_possible);
  } else {
    CHECK(result.device_to_host_transfers == 1);
    CHECK(result.device_to_host_bytes == 12 * sizeof(double) + 4 * sizeof(int));
    CHECK(result.synchronizations == 1);
    if (exclude_both) {
      CHECK_FALSE(result.primal_ray_possible);
      CHECK_FALSE(result.farkas_possible);
    }
  }
}

void setCertificateBounds(const int type, double& lower, double& upper) {
  lower = type == 1 || type == 4 ? -kHighsInf : 0.0;
  upper = type == 0 || type == 4 ? kHighsInf : (type == 3 ? 0.0 : 1.0);
}

void checkAcceptedCertificateCandidates(const hipdlp::backend::Kind kind) {
  for (const bool scaled : {false, true}) {
    for (const bool rowwise : {false, true}) {
      for (const ObjSense sense : {ObjSense::kMinimize, ObjSense::kMaximize}) {
        for (int col_type = 0; col_type < 5; ++col_type) {
          for (int row_type = 0; row_type < 5; ++row_type) {
            CAPTURE(scaled, rowwise, sense, col_type, row_type);
            HighsLp ray =
                makeEqualityLp(2, 3, {0, 0, 1, 2}, {0, 1}, {2.0, 0.25});
            ray.sense_ = sense;
            ray.offset_ = 17.0;
            ray.col_cost_[0] = sense == ObjSense::kMinimize ? -2.0 : 2.0;
            setCertificateBounds(col_type, ray.col_lower_[1],
                                 ray.col_upper_[1]);
            setCertificateBounds(row_type, ray.row_lower_[0],
                                 ray.row_upper_[0]);
            double direction = col_type == 1 ? -0.5 : 0.5;
            if (col_type == 2 || col_type == 3 || row_type == 2 ||
                row_type == 3 || (row_type == 0 && direction < 0.0) ||
                (row_type == 1 && direction > 0.0))
              direction = 0.0;
            if (rowwise) ray.a_matrix_.ensureRowwise();
            checkCertificateCandidate(kind, ray, {1.0, direction, 0.0},
                                      {1.0, 0.0}, {0.0, 0.0, 0.0}, scaled, 1,
                                      -1);

            const HighsInt active = row_type == 4 ? 1 : 0;
            const double sign = row_type == 1 ? -1.0 : 1.0;
            const double coefficient =
                col_type == 0 ? -2.0 : (col_type == 4 ? 0.0 : 2.0);
            const double rhs = col_type == 2 || col_type == 3 ? 3.0 : 1.0;
            HighsLp farkas = makeEqualityLp(2, 3, {0, 1, 1, 1}, {active},
                                            {sign * coefficient});
            farkas.sense_ = sense;
            setCertificateBounds(col_type, farkas.col_lower_[0],
                                 farkas.col_upper_[0]);
            if (col_type == 1) farkas.col_upper_[0] = 0.0;
            if (col_type == 3)
              farkas.col_lower_[0] = farkas.col_upper_[0] = 1.0;
            if (row_type == 4) {
              farkas.row_lower_[0] = -kHighsInf;
              farkas.row_upper_[0] = kHighsInf;
            }
            farkas.row_lower_[active] = row_type == 1 ? -kHighsInf : rhs;
            farkas.row_upper_[active] =
                row_type == 0 || row_type == 4 ? kHighsInf : sign * rhs;
            if (row_type == 2) farkas.row_upper_[active] = rhs + 1.0;
            std::vector<double> multiplier(2, 0.0);
            multiplier[active] = sign;
            if (rowwise) farkas.a_matrix_.ensureRowwise();
            checkCertificateCandidate(kind, farkas, {1.0, 0.5, 0.25},
                                      multiplier, {-coefficient, 0.0, 0.0},
                                      scaled, 0, 1);
          }
        }
      }
    }
  }
}

void checkFeasibleCertificateExclusions(const hipdlp::backend::Kind kind) {
  for (const bool scaled : {false, true}) {
    for (const ObjSense sense : {ObjSense::kMinimize, ObjSense::kMaximize}) {
      for (const double sign : {-1.0, 1.0}) {
        HighsLp lp = makeEqualityLp(2, 3, {0, 1, 2, 2}, {0, 1}, {1.0, 1.0});
        lp.sense_ = sense;
        lp.row_lower_ = lp.row_upper_ = {1.0, 1.0};
        lp.col_lower_ = {0.0, 0.0, 0.0};
        lp.col_cost_ = {2.0 * sign, sign, 0.0};
        checkCertificateCandidate(kind, lp, {1.0, 1.0, 1.0}, {2.0 * sign, sign},
                                  {-2.0 * sign, -sign, 0.0}, scaled, 0, 0,
                                  true);
      }
    }
  }
}

void checkCertificateRoundingCandidates(const hipdlp::backend::Kind kind) {
  HighsLp scaled_ray =
      makeEqualityLp(1, 2, {0, 1, 2}, {0, 0}, {2.0, -0.25});
  scaled_ray.col_cost_[0] = -1.0;
  checkCertificateCandidate(kind, scaled_ray, {1.0, 8.0}, {1.0}, {-2.0, 0.25},
                            true, 1, 0);
  HighsLp scaled_farkas =
      makeEqualityLp(2, 1, {0, 2}, {0, 1}, {1.0, 2.0});
  scaled_farkas.row_lower_ = scaled_farkas.row_upper_ = {1.0, 1.0};
  checkCertificateCandidate(kind, scaled_farkas, {1.0}, {2.0, -1.0}, {0.0},
                            true, 0, 1);

  HighsLp bound_scale = makeEqualityLp(1, 1, {0, 1}, {0}, {-1e-5});
  bound_scale.col_lower_ = bound_scale.col_upper_ = {1e5};
  for (const bool scaled : {false, true})
    checkCertificateCandidate(kind, bound_scale, {1e5}, {1.0}, {-1e-5}, scaled,
                              0, 1);

  HighsLp duplicate =
      makeEqualityLp(1, 1, {0, 3}, {0, 0, 0}, {1e16, 1.0, -1e16});
  duplicate.col_cost_[0] = -1.0;
  checkCertificateCandidate(kind, duplicate, {1.0}, {1.0}, {0.0}, false, 1, 0);

  HighsLp cancellation = makeEqualityLp(4, 1, {0, 0}, {}, {});
  cancellation.row_lower_ = cancellation.row_upper_ = {1e16, -1.0, -1e16, 0.5};
  checkCertificateCandidate(kind, cancellation, {1.0}, {1.0, 1.0, 1.0, 1.0},
                            {0.0}, false, 0, 1);

  for (const double tolerance : {1e-2, 1e-4, 1e-8}) {
    HighsLp near = makeEqualityLp(0, 2, {0, 0, 0}, {}, {});
    near.col_cost_ = {-1.0, 0.0};
    near.col_lower_[1] = -1.0;
    near.col_upper_[1] = 1.0;
    checkCertificateCandidate(kind, near, {1.0, std::nextafter(tolerance, 0.0)},
                              {}, {1.0, 1.0}, false, 1, 0, false, tolerance);
    near.col_cost_[0] = -std::nextafter(tolerance, kHighsInf);
    checkCertificateCandidate(kind, near, {1.0, 0.0}, {}, {1.0, 1.0}, false, 1,
                              0, false, tolerance);
    HighsLp near_sign = makeEqualityLp(1, 1, {0, 1}, {0}, {-0.25 * tolerance});
    near_sign.col_lower_[0] = 0.0;
    near_sign.row_lower_[0] = 1.0;
    near_sign.row_upper_[0] = kHighsInf;
    for (const bool scaled : {false, true})
      checkCertificateCandidate(kind, near_sign, {1.0}, {1.0}, {-0.5 * tolerance},
                                scaled, 0, 1, false, tolerance);
  }
  HighsLp extreme = makeEqualityLp(1, 2, {0, 0, 1}, {0}, {1.0});
  extreme.col_cost_[0] = -1.0;
  for (const double diagonal : {1e150, 1e-150, 1e300, 1e-300})
    checkCertificateCandidate(kind, extreme, {1.0, 0.0}, {1.0}, {0.0, -1.0},
                              false, 1, 0, false, 1e-4, diagonal);
}

void checkCertificateUncertainRecovery(const hipdlp::backend::Kind kind) {
  HighsLp lp = makeEqualityLp(0, 1, {0, 0}, {}, {});
  for (const double diagonal : {1e100, 1e-100}) {
    ResidentCheckFixture fixture(kind, lp, false, diagonal);
    fixture.prepare(lp);
    const double value = diagonal > 1.0 ? 1e300 : 1e-300;
    const std::vector<double> primal = {value};
    const std::vector<double> aty = {diagonal};
    const auto screen =
        hostCertificateScreen(fixture.problem, lp, primal, {}, aty, 1e-4);
    CHECK(screen.primal_ray_possible);
    CHECK(screen.farkas_possible);
    std::vector<double> recovered;
    CHECK_FALSE(fixture.problem.recoverPrimalDirection(primal, recovered).ok());
    REQUIRE(fixture.backend->uploadVector(primal, fixture.primal).ok());
    REQUIRE(fixture.backend->uploadVector(aty, fixture.aty).ok());
    hipdlp::backend::CheckResult result;
    REQUIRE(fixture.backend
                ->evaluateCheck(fixture.problem, fixture.workspace,
                                fixture.input, result)
                .ok());
    CHECK(result.valid);
    CHECK(result.primal_ray_possible);
    CHECK(result.farkas_possible);
  }
  ResidentCheckFixture fixture(kind, lp);
  fixture.prepare(lp);
  REQUIRE(fixture.backend->uploadVector({1.0}, fixture.primal).ok());
  REQUIRE(fixture.backend->uploadVector({1.0}, fixture.aty).ok());
  for (const double tolerance :
       {0.0, -1.0, kHighsInf, std::numeric_limits<double>::quiet_NaN()}) {
    fixture.input.certificate_tolerance = tolerance;
    hipdlp::backend::CheckResult result;
    REQUIRE(fixture.backend
                ->evaluateCheck(fixture.problem, fixture.workspace,
                                fixture.input, result)
                .ok());
    CHECK(result.valid);
    CHECK(result.primal_ray_possible);
    CHECK(result.farkas_possible);
  }
}

void checkRawMetrics(const hipdlp::algorithm::RawConvergenceMetrics& actual,
                     const hipdlp::algorithm::RawConvergenceMetrics& expected) {
  CHECK(actual.primal_residual == Approx(expected.primal_residual).epsilon(1e-12));
  CHECK(actual.dual_residual == Approx(expected.dual_residual).epsilon(1e-12));
  CHECK(actual.objective_gap == Approx(expected.objective_gap).epsilon(1e-12));
  CHECK(actual.primal_objective == Approx(expected.primal_objective).epsilon(1e-12));
  CHECK(actual.dual_objective == Approx(expected.dual_objective).epsilon(1e-12));
  CHECK(actual.unscaled_rhs_norm == expected.unscaled_rhs_norm);
  CHECK(actual.unscaled_cost_norm == expected.unscaled_cost_norm);
}

void checkEmptyPacket(const hipdlp::backend::CheckResult& result,
                      const bool requires_host_oracle = false) {
  CHECK_FALSE(result.valid);
  checkRawMetrics(result.raw, hipdlp::algorithm::RawConvergenceMetrics());
  CHECK(result.primal_delta_norm == 0.0);
  CHECK(result.dual_delta_norm == 0.0);
  CHECK(result.delta_dot == 0.0);
  CHECK(result.primal_anchor_distance == 0.0);
  CHECK(result.dual_anchor_distance == 0.0);
    CHECK(result.requires_host_oracle == requires_host_oracle);
  CHECK(result.primal_ray_possible);
  CHECK(result.farkas_possible);
  CHECK(result.scalar_count == 0);
    CHECK(result.device_to_host_transfers == (requires_host_oracle ? 1 : 0));
    CHECK(result.device_to_host_bytes ==
      (requires_host_oracle ? 12 * sizeof(double) + 4 * sizeof(int) : 0));
    CHECK(result.synchronizations == (requires_host_oracle ? 1 : 0));
}

void checkResidentRawOracle(const hipdlp::backend::Kind kind) {
  for (const bool scaled : {false, true}) {
    for (const ObjSense sense : {ObjSense::kMinimize, ObjSense::kMaximize}) {
      HighsLp lp = makeEqualityLp(4, 4, {0, 1, 2, 3, 4}, {0, 1, 2, 3},
                                  {1.0, 4.0, 0.25, 2.0});
      lp.sense_ = sense;
      lp.offset_ = -7.25;
      lp.col_cost_ = {2.0, -3.0, 5.0, -7.0};
      lp.col_lower_ = {-2.0, -kHighsInf, 0.0, -kHighsInf};
      lp.col_upper_ = {3.0, kHighsInf, kHighsInf, 4.0};
      lp.row_lower_ = {1.0, -2.0, -kHighsInf, -3.0};
      lp.row_upper_ = {1.0, kHighsInf, 5.0, 6.0};
      ResidentCheckFixture fixture(kind, lp, scaled);
      fixture.prepare(lp);
      const std::vector<double> primal(fixture.problem.numCols(), 0.75);
      const std::vector<double> dual(fixture.problem.numRows(), -0.5);
      const std::vector<double> ax(fixture.problem.numRows(), 0.25);
      const std::vector<double> aty(fixture.problem.numCols(), 1.5);
      REQUIRE(fixture.backend->uploadVector(primal, fixture.primal).ok());
      REQUIRE(fixture.backend->uploadVector(dual, fixture.dual).ok());
      REQUIRE(fixture.backend->uploadVector(ax, fixture.ax).ok());
      REQUIRE(fixture.backend->uploadVector(aty, fixture.aty).ok());
      hipdlp::algorithm::RawConvergenceMetrics expected;
      REQUIRE(hipdlp::algorithm::computeRawConvergenceMetrics(
                  fixture.problem, primal, dual, ax, aty, expected).ok());
      hipdlp::backend::CheckResult result;
      const hipdlp::api::Error error = fixture.backend->evaluateCheck(
          fixture.problem, fixture.workspace, fixture.input, result);
      INFO(error.message);
      REQUIRE(error.ok());
      REQUIRE(result.valid);
      CHECK_FALSE(result.requires_host_oracle);
      checkRawMetrics(result.raw, expected);
      CHECK(result.scalar_count == 7);
      CHECK(result.primal_delta_norm == 0.0);
      CHECK(result.dual_delta_norm == 0.0);
      CHECK(result.delta_dot == 0.0);
      CHECK(result.primal_anchor_distance == 0.0);
      CHECK(result.dual_anchor_distance == 0.0);
      if (kind == hipdlp::backend::Kind::kCpu) {
        CHECK(result.primal_ray_possible);
        CHECK(result.farkas_possible);
        CHECK(result.device_to_host_transfers == 0);
        CHECK(result.device_to_host_bytes == 0);
        CHECK(result.synchronizations == 0);
      }
    }
  }
}

void checkResidentOptionalPrimitives(const hipdlp::backend::Kind kind) {
  HighsLp lp = makeEqualityLp(2, 2, {0, 1, 2}, {0, 1}, {1.0, 1.0});
  ResidentCheckFixture fixture(kind, lp);
  fixture.prepare(lp);
  REQUIRE(fixture.backend->uploadVector({6.0, 8.0}, fixture.primal).ok());
  REQUIRE(fixture.backend->uploadVector({5.0, 12.0}, fixture.dual).ok());
  REQUIRE(fixture.backend->uploadVector({3.0, 4.0}, fixture.delta_primal).ok());
  REQUIRE(fixture.backend->uploadVector({0.0, 12.0}, fixture.delta_dual).ok());
  REQUIRE(fixture.backend->uploadVector({2.0, -1.0}, fixture.aty_delta).ok());
  REQUIRE(fixture.backend->uploadVector({3.0, 4.0}, fixture.primal_anchor).ok());
  REQUIRE(fixture.backend->uploadVector({0.0, 0.0}, fixture.dual_anchor).ok());
  for (const bool fixed_point : {false, true}) {
    for (const bool anchor_distances : {false, true}) {
      fixture.input.fixed_point = fixed_point;
      fixture.input.anchor_distances = anchor_distances;
      fixture.input.delta_primal = fixed_point ? &fixture.delta_primal : nullptr;
      fixture.input.delta_dual = fixed_point ? &fixture.delta_dual : nullptr;
      fixture.input.aty_delta = fixed_point ? &fixture.aty_delta : nullptr;
      fixture.input.primal_anchor = anchor_distances ? &fixture.primal_anchor : nullptr;
      fixture.input.dual_anchor = anchor_distances ? &fixture.dual_anchor : nullptr;
      hipdlp::backend::CheckResult result;
      const hipdlp::api::Error error = fixture.backend->evaluateCheck(
          fixture.problem, fixture.workspace, fixture.input, result);
      INFO(error.message);
      REQUIRE(error.ok());
      REQUIRE(result.valid);
      CHECK(result.primal_delta_norm == (fixed_point ? 5.0 : 0.0));
      CHECK(result.dual_delta_norm == (fixed_point ? 12.0 : 0.0));
      CHECK(result.delta_dot == (fixed_point ? 2.0 : 0.0));
      CHECK(result.primal_anchor_distance == (anchor_distances ? 5.0 : 0.0));
      CHECK(result.dual_anchor_distance == (anchor_distances ? 13.0 : 0.0));
      CHECK(result.scalar_count == 7 + (fixed_point ? 3 : 0) +
                                       (anchor_distances ? 2 : 0));
      if (kind == hipdlp::backend::Kind::kCuda) {
        CHECK(result.device_to_host_transfers == 1);
        CHECK(result.device_to_host_bytes == 12 * sizeof(double) + 4 * sizeof(int));
        CHECK(result.synchronizations == 1);
      }
    }
  }
  REQUIRE(fixture.backend->fusedProjectedPrimalUpdate(
              fixture.bound, fixture.workspace, fixture.primal, fixture.aty,
              0.0, fixture.delta_primal).ok());
}

void checkResidentFixedPointPreparation(const hipdlp::backend::Kind kind) {
  HighsLp lp = makeEqualityLp(2, 2, {0, 1, 2}, {0, 1}, {2.0, 3.0});
  ResidentCheckFixture fixture(kind, lp);
  fixture.prepare(lp);
  REQUIRE(fixture.backend->uploadVector({6.0, 8.0}, fixture.primal).ok());
  REQUIRE(fixture.backend->uploadVector({5.0, 12.0}, fixture.dual).ok());
  REQUIRE(fixture.backend->uploadVector({3.0, 4.0}, fixture.primal_anchor).ok());
  REQUIRE(fixture.backend->uploadVector({1.0, 2.0}, fixture.dual_anchor).ok());
  REQUIRE(fixture.backend->multiply(fixture.bound, fixture.primal, fixture.ax).ok());
  REQUIRE(fixture.backend->multiplyTranspose(fixture.bound, fixture.dual,
                                              fixture.aty).ok());
  fixture.input.fixed_point = true;
  fixture.input.anchor_distances = true;
  fixture.input.delta_primal = &fixture.delta_primal;
  fixture.input.delta_dual = &fixture.delta_dual;
  fixture.input.aty_delta = &fixture.aty_delta;
  fixture.input.primal_anchor = &fixture.primal_anchor;
  fixture.input.dual_anchor = &fixture.dual_anchor;
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
  const auto before = hipdlp::backend::cuda::testing::resourceAudit();
#endif
  for (int repetition = 0; repetition < 3; ++repetition) {
    REQUIRE(fixture.backend->prepareCheckFixedPoint(
                fixture.bound, fixture.workspace, fixture.primal, fixture.dual,
                fixture.primal_anchor, fixture.dual_anchor, fixture.delta_primal,
                fixture.delta_dual, fixture.aty_delta).ok());
    hipdlp::backend::CheckResult result;
    REQUIRE(fixture.backend->evaluateCheck(
                fixture.problem, fixture.workspace, fixture.input, result).ok());
    REQUIRE(result.valid);
    CHECK_FALSE(result.requires_host_oracle);
    CHECK(result.scalar_count == 12);
    CHECK(result.primal_delta_norm == 5.0);
    CHECK(result.dual_delta_norm == Approx(std::sqrt(116.0)));
    CHECK(result.delta_dot == 144.0);
    CHECK(result.primal_anchor_distance == 5.0);
    CHECK(result.dual_anchor_distance == Approx(std::sqrt(116.0)));
  }
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
  if (kind == hipdlp::backend::Kind::kCuda) {
    const auto after = hipdlp::backend::cuda::testing::resourceAudit();
    CHECK(after.allocation_calls == before.allocation_calls);
    CHECK(after.free_calls == before.free_calls);
    CHECK(after.pinned_allocation_calls == before.pinned_allocation_calls);
    CHECK(after.pinned_free_calls == before.pinned_free_calls);
    CHECK(after.host_to_device_transfers == before.host_to_device_transfers);
    CHECK(after.device_to_host_transfers - before.device_to_host_transfers == 3);
    CHECK(after.device_to_host_bytes - before.device_to_host_bytes == 3 * 112);
    CHECK(after.synchronization_calls - before.synchronization_calls == 3);
  }
#endif
  checkApproxVector(download(*fixture.backend, fixture.primal), {6.0, 8.0}, 0.0);
  checkApproxVector(download(*fixture.backend, fixture.dual), {5.0, 12.0}, 0.0);
  checkApproxVector(download(*fixture.backend, fixture.primal_anchor), {3.0, 4.0}, 0.0);
  checkApproxVector(download(*fixture.backend, fixture.dual_anchor), {1.0, 2.0}, 0.0);
  checkApproxVector(download(*fixture.backend, fixture.delta_primal), {3.0, 4.0}, 0.0);
  checkApproxVector(download(*fixture.backend, fixture.delta_dual), {4.0, 10.0}, 0.0);
  checkApproxVector(download(*fixture.backend, fixture.aty_delta), {8.0, 30.0}, 0.0);
}

void checkResidentEmptyShapes(const hipdlp::backend::Kind kind) {
  for (const HighsInt rows : {0, 3}) {
    for (const HighsInt cols : {0, 2}) {
      HighsLp lp = makeEqualityLp(rows, cols,
                                  std::vector<HighsInt>(cols + 1, 0), {}, {});
      lp.offset_ = 2.75;
      ResidentCheckFixture fixture(kind, lp);
      fixture.prepare(lp);
      fixture.input.fixed_point = true;
      fixture.input.anchor_distances = true;
      fixture.input.delta_primal = &fixture.delta_primal;
      fixture.input.delta_dual = &fixture.delta_dual;
      fixture.input.aty_delta = &fixture.aty_delta;
      fixture.input.primal_anchor = &fixture.primal_anchor;
      fixture.input.dual_anchor = &fixture.dual_anchor;
      hipdlp::backend::CheckResult result;
      REQUIRE(fixture.backend->evaluateCheck(
                  fixture.problem, fixture.workspace, fixture.input, result).ok());
      REQUIRE(result.valid);
      hipdlp::algorithm::RawConvergenceMetrics expected;
      REQUIRE(hipdlp::algorithm::computeRawConvergenceMetrics(
                  fixture.problem, std::vector<double>(cols, 0.0),
                  std::vector<double>(rows, 0.0), std::vector<double>(rows, 0.0),
                  std::vector<double>(cols, 0.0), expected).ok());
      checkRawMetrics(result.raw, expected);
      CHECK(result.raw.primal_objective == 2.75);
      CHECK(result.raw.dual_objective == 2.75);
      CHECK(result.scalar_count == 12);
      CHECK(result.primal_delta_norm == 0.0);
      CHECK(result.dual_delta_norm == 0.0);
      CHECK(result.delta_dot == 0.0);
      CHECK(result.primal_anchor_distance == 0.0);
      CHECK(result.dual_anchor_distance == 0.0);
      CHECK(result.primal_ray_possible);
      CHECK(result.farkas_possible);
    }
  }
}

void checkResidentOracleOverflow(const hipdlp::backend::Kind kind,
                                 const bool objective_product) {
  const double maximum = std::numeric_limits<double>::max();
  const double epsilon = std::numeric_limits<double>::epsilon();
  const HighsInt size = objective_product ? 1 : 17;
  std::vector<HighsInt> starts(size + 1);
  std::vector<HighsInt> indices(size);
  for (HighsInt index = 0; index < size; ++index) {
    starts[index] = index;
    indices[index] = index;
  }
  starts[size] = size;
  HighsLp lp = makeEqualityLp(size, size, starts, indices,
                              std::vector<double>(size, 1.0));
  lp.col_lower_.assign(size, -kHighsInf);
  lp.col_upper_.assign(size, kHighsInf);
  lp.col_cost_.assign(size, objective_product ? maximum : 0.0);
  lp.row_lower_.assign(size, objective_product ? 1.0 : 0.0);
  lp.row_upper_ = lp.row_lower_;
  lp.offset_ = objective_product ? -maximum : 0.0;
  ResidentCheckFixture fixture(kind, lp);
  fixture.prepare(lp);
  std::vector<double> primal(size, std::ldexp(maximum, -27));
  primal[0] = objective_product ? 1.0 + epsilon : maximum;
  const std::vector<double> dual(size, objective_product ? maximum : 0.0);
  const std::vector<double> ax = primal;
  const std::vector<double> aty = dual;
  REQUIRE(fixture.backend->uploadVector(primal, fixture.primal).ok());
  REQUIRE(fixture.backend->uploadVector(dual, fixture.dual).ok());
  REQUIRE(fixture.backend->uploadVector(ax, fixture.ax).ok());
  REQUIRE(fixture.backend->uploadVector(aty, fixture.aty).ok());
  hipdlp::algorithm::RawConvergenceMetrics oracle;
  const auto oracle_error = hipdlp::algorithm::computeRawConvergenceMetrics(
      fixture.problem, primal, dual, ax, aty, oracle);
  INFO(oracle_error.message);
  REQUIRE(oracle_error.ok());
  CHECK(oracle.primal_residual == (objective_product ? epsilon : maximum));
  CHECK(oracle.dual_residual == 0.0);
  CHECK(oracle.primal_objective ==
        (objective_product ? std::ldexp(maximum, -52) : 0.0));
  CHECK(oracle.dual_objective == 0.0);
  CHECK(oracle.objective_gap == oracle.primal_objective);
  hipdlp::backend::CheckResult result;
  result.valid = true;
  result.raw.primal_residual = 42.0;
  const auto error = fixture.backend->evaluateCheck(
      fixture.problem, fixture.workspace, fixture.input, result);
  INFO(error.message);
  CHECK(error.ok());
  if (kind == hipdlp::backend::Kind::kCuda) {
    checkEmptyPacket(result, true);
  } else {
    REQUIRE(result.valid);
    checkRawMetrics(result.raw, oracle);
  }
  std::vector<double> downloaded_primal;
  std::vector<double> downloaded_dual;
  std::vector<double> downloaded_ax;
  std::vector<double> downloaded_aty;
  REQUIRE(fixture.backend->downloadVector(fixture.primal, downloaded_primal).ok());
  REQUIRE(fixture.backend->downloadVector(fixture.dual, downloaded_dual).ok());
  REQUIRE(fixture.backend->downloadVector(fixture.ax, downloaded_ax).ok());
  REQUIRE(fixture.backend->downloadVector(fixture.aty, downloaded_aty).ok());
  CHECK(downloaded_primal == primal);
  CHECK(downloaded_dual == dual);
  CHECK(downloaded_ax == ax);
  CHECK(downloaded_aty == aty);
  hipdlp::algorithm::RawConvergenceMetrics downloaded_oracle;
  REQUIRE(hipdlp::algorithm::computeRawConvergenceMetrics(
              fixture.problem, downloaded_primal, downloaded_dual,
              downloaded_ax, downloaded_aty, downloaded_oracle).ok());
  checkRawMetrics(downloaded_oracle, oracle);
}

void checkResidentStoppingRounding(const hipdlp::backend::Kind kind) {
  for (const int metric : {0, 1, 2, 3, 4}) {
    for (const HighsInt size : {17, 257}) {
      INFO("metric " << metric << ", size " << size);
      const double magnitude = metric == 4
          ? std::numeric_limits<double>::denorm_min() : std::ldexp(1.0, -14);
      std::vector<HighsInt> starts(size + 1);
      std::vector<HighsInt> indices(size);
      for (HighsInt index = 0; index < size; ++index) {
        starts[index] = index;
        indices[index] = index;
      }
      starts[size] = size;
      HighsLp lp = makeEqualityLp(size, size, starts, indices,
                                  std::vector<double>(size, 1.0));
      lp.col_lower_.assign(size, -kHighsInf);
      lp.col_upper_.assign(size, kHighsInf);
      lp.col_cost_.assign(size, metric == 3 ? 0.1 : 0.0);
      lp.row_lower_.assign(size, metric == 3 ? 1.0 : 0.0);
      if (metric == 1) lp.row_lower_[0] = 1.0;
      if (metric == 2) lp.col_cost_[0] = 1.0;
      lp.row_upper_ = lp.row_lower_;
      lp.offset_ = metric == 3 ? 0.3 : 0.0;
      ResidentCheckFixture fixture(kind, lp);
      fixture.prepare(lp);
      std::vector<double> primal(size, std::ldexp(magnitude, -27));
      primal[0] = magnitude;
      std::vector<double> dual(size, 0.0);
      if (metric == 1) primal[0] += 1.0;
      if (metric == 2) {
        for (HighsInt index = 0; index < size; ++index)
          dual[index] = lp.col_cost_[index] - primal[index];
        primal.assign(size, 0.0);
      }
      if (metric == 3) primal.assign(size, 1.0);
      REQUIRE(fixture.backend->uploadVector(primal, fixture.primal).ok());
      REQUIRE(fixture.backend->uploadVector(dual, fixture.dual).ok());
      REQUIRE(fixture.backend->uploadVector(primal, fixture.ax).ok());
      REQUIRE(fixture.backend->uploadVector(dual, fixture.aty).ok());
      hipdlp::algorithm::StoppingInput expected;
      expected.tolerances.primal = 1.0;
      expected.tolerances.dual = 1.0;
      expected.tolerances.gap = 1.0;
      REQUIRE(hipdlp::algorithm::computeRawConvergenceMetrics(
                  fixture.problem, primal, dual, primal, dual, expected.raw).ok());
      if (metric != 3)
        CHECK((metric == 2 ? expected.raw.dual_residual
                           : expected.raw.primal_residual) == magnitude);
      const auto normalized = hipdlp::algorithm::evaluateStopping(expected);
      REQUIRE(normalized.ok());
      const double threshold = metric == 3 ? normalized.metrics.relative_gap
          : metric == 2 ? normalized.metrics.relative_dual
                        : normalized.metrics.relative_primal;
      for (const double tolerance :
           {threshold, std::nextafter(threshold, kHighsInf),
            std::nextafter(threshold, 0.0)}) {
        if (tolerance == 0.0) continue;
        if (metric == 3) {
          expected.tolerances.gap = tolerance;
          fixture.input.gap_tolerance = tolerance;
        } else if (metric == 2) {
          expected.tolerances.dual = tolerance;
          fixture.input.dual_tolerance = tolerance;
        } else {
          expected.tolerances.primal = tolerance;
          fixture.input.primal_tolerance = tolerance;
        }
        const auto oracle_stopping = hipdlp::algorithm::evaluateStopping(expected);
        REQUIRE(oracle_stopping.ok());
        REQUIRE(oracle_stopping.decision ==
          (tolerance < threshold
               ? hipdlp::algorithm::StoppingDecision::kContinue
               : hipdlp::algorithm::StoppingDecision::kOptimal));
        hipdlp::backend::CheckResult result;
        REQUIRE(fixture.backend->evaluateCheck(
                    fixture.problem, fixture.workspace, fixture.input, result).ok());
        hipdlp::algorithm::StoppingInput actual = expected;
        if (kind == hipdlp::backend::Kind::kCuda) {
          checkEmptyPacket(result, true);
          const auto downloaded_primal = download(*fixture.backend, fixture.primal);
          const auto downloaded_dual = download(*fixture.backend, fixture.dual);
          const auto downloaded_ax = download(*fixture.backend, fixture.ax);
          const auto downloaded_aty = download(*fixture.backend, fixture.aty);
          REQUIRE(hipdlp::algorithm::computeRawConvergenceMetrics(
                      fixture.problem, downloaded_primal, downloaded_dual,
                      downloaded_ax, downloaded_aty, actual.raw).ok());
        } else {
          REQUIRE(result.valid);
          actual.raw = result.raw;
        }
        const auto stopping = hipdlp::algorithm::evaluateStopping(actual);
        REQUIRE(stopping.ok());
        CHECK(stopping.decision == oracle_stopping.decision);
      }
      if (metric == 3) {
        fixture.input.dual_tolerance = 1e-30;
        hipdlp::backend::CheckResult result;
        REQUIRE(fixture.backend->evaluateCheck(
                    fixture.problem, fixture.workspace, fixture.input, result).ok());
        if (kind == hipdlp::backend::Kind::kCuda) checkEmptyPacket(result, true);
      }
    }
  }
}

void checkResidentDualProductOverflow(const hipdlp::backend::Kind kind) {
  const double maximum = std::numeric_limits<double>::max();
  const double multiplier = 1.0 + std::numeric_limits<double>::epsilon();
  for (const int support : {0, 1, 2}) {
    INFO("support " << support);
    HighsLp lp = makeEqualityLp(1, 1, {0, 1}, {0}, {1.0});
    lp.offset_ = -maximum;
    lp.col_cost_[0] = support == 0 ? 1.0
                                 : support == 1 ? multiplier : -multiplier;
    lp.col_lower_[0] = support == 1 ? maximum : -kHighsInf;
    lp.col_upper_[0] = support == 2 ? -maximum : kHighsInf;
    lp.row_lower_[0] = lp.row_upper_[0] = support == 2 ? -maximum : maximum;
    ResidentCheckFixture fixture(kind, lp);
    fixture.prepare(lp);
    const std::vector<double> primal = {lp.row_lower_[0]};
    const std::vector<double> dual = {support == 0 ? multiplier : 0.0};
    REQUIRE(fixture.backend->uploadVector(primal, fixture.primal).ok());
    REQUIRE(fixture.backend->uploadVector(dual, fixture.dual).ok());
    REQUIRE(fixture.backend->uploadVector(primal, fixture.ax).ok());
    REQUIRE(fixture.backend->uploadVector(dual, fixture.aty).ok());
    hipdlp::algorithm::RawConvergenceMetrics oracle;
    REQUIRE(hipdlp::algorithm::computeRawConvergenceMetrics(
                fixture.problem, primal, dual, primal, dual, oracle).ok());
    CHECK(oracle.primal_objective ==
          (support == 0 ? 0.0 : std::ldexp(maximum, -52)));
    CHECK(oracle.dual_objective == std::ldexp(maximum, -52));
    CHECK(oracle.primal_residual == 0.0);
    hipdlp::backend::CheckResult result;
    REQUIRE(fixture.backend->evaluateCheck(
                fixture.problem, fixture.workspace, fixture.input, result).ok());
    if (kind == hipdlp::backend::Kind::kCuda) {
      checkEmptyPacket(result, true);
    } else {
      REQUIRE(result.valid);
      checkRawMetrics(result.raw, oracle);
    }
    CHECK(download(*fixture.backend, fixture.primal) == primal);
    CHECK(download(*fixture.backend, fixture.dual) == dual);
    CHECK(download(*fixture.backend, fixture.ax) == primal);
    CHECK(download(*fixture.backend, fixture.aty) == dual);
  }
}

void checkResidentReducedCostOverflow(const hipdlp::backend::Kind kind) {
  const double maximum = std::numeric_limits<double>::max();
  for (const bool bounded : {false, true}) {
    HighsLp lp = makeEqualityLp(1, 1, {0, 1}, {0}, {1.0});
    lp.col_cost_[0] = maximum;
    lp.col_lower_[0] = bounded ? 0.0 : -kHighsInf;
    lp.col_upper_[0] = kHighsInf;
    ResidentCheckFixture fixture(kind, lp);
    fixture.prepare(lp);
    REQUIRE(fixture.backend->uploadVector({-maximum}, fixture.dual).ok());
    REQUIRE(fixture.backend->uploadVector({-maximum}, fixture.aty).ok());
    hipdlp::algorithm::RawConvergenceMetrics oracle;
    CHECK_FALSE(hipdlp::algorithm::computeRawConvergenceMetrics(
                    fixture.problem, {0.0}, {-maximum}, {0.0}, {-maximum},
                    oracle).ok());
    if (bounded) CHECK(oracle.dual_residual == 0.0);
    hipdlp::backend::CheckResult result;
    const auto error = fixture.backend->evaluateCheck(
        fixture.problem, fixture.workspace, fixture.input, result);
    const bool deferred = kind == hipdlp::backend::Kind::kCuda;
    CHECK(error.ok() == deferred);
    checkEmptyPacket(result, deferred);
    CHECK(download(*fixture.backend, fixture.primal) == std::vector<double>{0.0});
    CHECK(download(*fixture.backend, fixture.dual) == std::vector<double>{-maximum});
    CHECK(download(*fixture.backend, fixture.ax) == std::vector<double>{0.0});
    CHECK(download(*fixture.backend, fixture.aty) == std::vector<double>{-maximum});
    REQUIRE(fixture.backend->fill(fixture.aty, maximum).ok());
    REQUIRE(fixture.backend->evaluateCheck(
                fixture.problem, fixture.workspace, fixture.input, result).ok());
    CHECK(result.valid);
  }
}

void checkResidentToleranceValidation(const hipdlp::backend::Kind kind) {
  HighsLp lp = makeEqualityLp(1, 1, {0, 1}, {0}, {1.0});
  ResidentCheckFixture fixture(kind, lp);
  fixture.prepare(lp);
  using hipdlp::backend::CheckInput;
  double CheckInput::* const fields[] = {
      &CheckInput::primal_tolerance, &CheckInput::dual_tolerance,
      &CheckInput::gap_tolerance};
  for (const auto field : fields) {
    for (const double tolerance :
         {0.0, 1e-4, std::numeric_limits<double>::max()}) {
      fixture.input.*field = tolerance;
      hipdlp::backend::CheckResult result;
      REQUIRE(fixture.backend->evaluateCheck(
                  fixture.problem, fixture.workspace, fixture.input, result).ok());
      REQUIRE(result.valid);
      CHECK_FALSE(result.requires_host_oracle);
      CHECK(result.raw.primal_residual == 0.0);
      CHECK(result.raw.dual_residual == 0.0);
      CHECK(result.raw.objective_gap == 0.0);
    }
    for (const double tolerance :
         {-1.0, kHighsInf, std::numeric_limits<double>::quiet_NaN()}) {
      fixture.input.*field = tolerance;
      hipdlp::backend::CheckResult result;
      result.valid = true;
      result.raw.primal_residual = 42.0;
      const auto error = fixture.backend->evaluateCheck(
          fixture.problem, fixture.workspace, fixture.input, result);
      CHECK(error.code == hipdlp::api::ErrorCode::kInvalidArgument);
      checkEmptyPacket(result);
    }
    fixture.input.*field = 0.0;
  }
  hipdlp::backend::CheckResult result;
  REQUIRE(fixture.backend->evaluateCheck(
              fixture.problem, fixture.workspace, fixture.input, result).ok());
  CHECK(result.valid);
  CHECK(download(*fixture.backend, fixture.primal) == std::vector<double>{0.0});
}

void checkResidentStableNorms(const hipdlp::backend::Kind kind) {
  HighsLp lp = makeEqualityLp(2, 2, {0, 1, 2}, {0, 1}, {1.0, 1.0});
  ResidentCheckFixture fixture(kind, lp);
  fixture.prepare(lp);
  const std::vector<double> ax = {1e300, -1e300};
  const std::vector<double> aty = {-1e300, 1e300};
  REQUIRE(fixture.backend->uploadVector(ax, fixture.ax).ok());
  REQUIRE(fixture.backend->uploadVector(aty, fixture.aty).ok());
  REQUIRE(fixture.backend->uploadVector({3e150, 4e150}, fixture.delta_primal).ok());
  REQUIRE(fixture.backend->uploadVector({3e-200, 4e-200}, fixture.delta_dual).ok());
  REQUIRE(fixture.backend->uploadVector({1e-150, 1e-150}, fixture.aty_delta).ok());
  fixture.input.fixed_point = true;
  fixture.input.delta_primal = &fixture.delta_primal;
  fixture.input.delta_dual = &fixture.delta_dual;
  fixture.input.aty_delta = &fixture.aty_delta;
  hipdlp::backend::CheckResult result;
  REQUIRE(fixture.backend->evaluateCheck(
              fixture.problem, fixture.workspace, fixture.input, result).ok());
  REQUIRE(result.valid);
  hipdlp::algorithm::RawConvergenceMetrics expected;
  REQUIRE(hipdlp::algorithm::computeRawConvergenceMetrics(
              fixture.problem, {0.0, 0.0}, {0.0, 0.0}, ax, aty, expected).ok());
  checkRawMetrics(result.raw, expected);
  CHECK(result.raw.primal_residual / 1e300 == Approx(std::sqrt(2.0)));
  CHECK(result.raw.dual_residual / 1e300 == Approx(std::sqrt(2.0)));
  CHECK(result.primal_delta_norm / 1e150 == Approx(5.0));
  CHECK(result.dual_delta_norm / 1e-200 == Approx(5.0));
  CHECK(result.delta_dot == Approx(7.0));
}

void checkResidentNumericalFailures(const hipdlp::backend::Kind kind) {
  for (const bool make_nan : {false, true}) {
    HighsLp lp = makeEqualityLp(1, 1, {0, 1}, {0}, {1e200});
    ResidentCheckFixture fixture(kind, lp);
    fixture.prepare(lp);
    for (const double invalid : {std::numeric_limits<double>::infinity(),
                                 std::numeric_limits<double>::quiet_NaN()})
      CHECK_FALSE(fixture.backend->uploadVector({invalid}, fixture.primal).ok());
    hipdlp::backend::CheckResult result;
    REQUIRE(fixture.backend->evaluateCheck(
                fixture.problem, fixture.workspace, fixture.input, result).ok());
    REQUIRE(result.valid);
    REQUIRE(fixture.backend->uploadVector({1e200}, fixture.primal).ok());
    const hipdlp::api::Error multiply_error = fixture.backend->multiply(
        fixture.bound, fixture.primal, fixture.ax);
    CHECK(multiply_error.ok() == (kind == hipdlp::backend::Kind::kCuda));
    if (make_nan) {
      const hipdlp::api::Error blend_error = fixture.backend->fusedVectorBlend(
          fixture.workspace, fixture.primal_anchor, fixture.dual_anchor,
          fixture.ax, fixture.ax, 1.0, fixture.aty, fixture.dual);
      CHECK(blend_error.ok() == (kind == hipdlp::backend::Kind::kCuda));
    }
    CHECK_FALSE(fixture.backend->evaluateCheck(
                    fixture.problem, fixture.workspace, fixture.input, result).ok());
    checkEmptyPacket(result);
  }

  for (const int failure : {0, 1, 2}) {
    HighsLp lp = makeEqualityLp(2, 2, {0, 1, 2}, {0, 1}, {1.0, 1.0});
    ResidentCheckFixture fixture(kind, lp);
    fixture.prepare(lp);
    const double largest = std::numeric_limits<double>::max();
    if (failure == 0) {
      REQUIRE(fixture.backend->uploadVector({largest, largest}, fixture.ax).ok());
    } else if (failure == 1) {
      REQUIRE(fixture.backend->uploadVector({largest, 0.0}, fixture.delta_primal).ok());
      REQUIRE(fixture.backend->uploadVector({2.0, 0.0}, fixture.aty_delta).ok());
      fixture.input.fixed_point = true;
      fixture.input.delta_primal = &fixture.delta_primal;
      fixture.input.delta_dual = &fixture.delta_dual;
      fixture.input.aty_delta = &fixture.aty_delta;
    } else {
      REQUIRE(fixture.backend->uploadVector({largest, 0.0}, fixture.primal).ok());
      REQUIRE(fixture.backend->uploadVector({-largest, 0.0}, fixture.primal_anchor).ok());
      fixture.input.anchor_distances = true;
      fixture.input.primal_anchor = &fixture.primal_anchor;
      fixture.input.dual_anchor = &fixture.dual_anchor;
    }
    hipdlp::backend::CheckResult result;
    const bool deferred = kind == hipdlp::backend::Kind::kCuda;
    const hipdlp::api::Error error = fixture.backend->evaluateCheck(
        fixture.problem, fixture.workspace, fixture.input, result);
    CHECK(error.ok() == deferred);
    checkEmptyPacket(result, deferred);
    if (deferred) {
      std::vector<double> mandatory;
      REQUIRE(fixture.backend->downloadVector(fixture.primal, mandatory).ok());
      if (failure == 0) {
        const auto ax = download(*fixture.backend, fixture.ax);
        hipdlp::algorithm::RawConvergenceMetrics oracle;
        CHECK_FALSE(hipdlp::algorithm::computeRawConvergenceMetrics(
                        fixture.problem, mandatory, {0.0, 0.0}, ax,
                        {0.0, 0.0}, oracle).ok());
        REQUIRE(fixture.backend->fill(fixture.ax, 0.0).ok());
      }
      fixture.input.fixed_point = false;
      fixture.input.anchor_distances = false;
      REQUIRE(fixture.backend->evaluateCheck(
                  fixture.problem, fixture.workspace, fixture.input, result).ok());
      CHECK(result.valid);
      if (failure == 0) continue;
      double optional_value = 0.0;
      const hipdlp::api::Error late_error = failure == 1
          ? fixture.backend->dot(fixture.delta_primal, fixture.aty_delta,
                                  optional_value)
          : fixture.backend->differenceNorm2(fixture.primal,
                                              fixture.primal_anchor,
                                              optional_value);
      CHECK_FALSE(late_error.ok());
    }
  }
}

void checkResidentInvalidHandles(const hipdlp::backend::Kind kind) {
  HighsLp lp = makeEqualityLp(1, 2, {0, 1, 1}, {0}, {1.0});
  ResidentCheckFixture fixture(kind, lp);
  hipdlp::backend::CheckResult result;
  CHECK_FALSE(fixture.backend->evaluateCheck(
                  fixture.problem, fixture.workspace, fixture.input, result).ok());
  checkEmptyPacket(result);
  fixture.prepare(lp);
  fixture.input.fixed_point = true;
  fixture.input.anchor_distances = true;
  fixture.input.delta_primal = &fixture.delta_primal;
  fixture.input.delta_dual = &fixture.delta_dual;
  fixture.input.aty_delta = &fixture.aty_delta;
  fixture.input.primal_anchor = &fixture.primal_anchor;
  fixture.input.dual_anchor = &fixture.dual_anchor;
  hipdlp::backend::CheckResult valid;
  REQUIRE(fixture.backend->evaluateCheck(
              fixture.problem, fixture.workspace, fixture.input, valid).ok());
  REQUIRE(valid.valid);
  using hipdlp::backend::CheckInput;
  using hipdlp::backend::VectorHandle;
  const VectorHandle* CheckInput::* fields[] = {
      &CheckInput::primal, &CheckInput::dual, &CheckInput::ax, &CheckInput::aty,
      &CheckInput::delta_primal, &CheckInput::delta_dual, &CheckInput::aty_delta,
      &CheckInput::primal_anchor, &CheckInput::dual_anchor};
  const bool row_field[] = {false, true, true, false, false, true, false, false, true};
  hipdlp::backend::FactoryResult other = hipdlp::backend::createBackend(kind);
  REQUIRE(other.ok());
  VectorHandle foreign_primal;
  VectorHandle foreign_dual;
  REQUIRE(other.backend->createVector(2, foreign_primal).ok());
  REQUIRE(other.backend->createVector(1, foreign_dual).ok());
  for (std::size_t index = 0; index < 9; ++index) {
    for (const int failure : {0, 1, 2}) {
      CheckInput broken = fixture.input;
      if (failure == 0) broken.*fields[index] = nullptr;
      if (failure == 1)
        broken.*fields[index] = row_field[index] ? &fixture.primal : &fixture.dual;
      if (failure == 2)
        broken.*fields[index] = row_field[index] ? &foreign_dual : &foreign_primal;
      result = valid;
      CHECK_FALSE(fixture.backend->evaluateCheck(
                      fixture.problem, fixture.workspace, broken, result).ok());
      checkEmptyPacket(result);
    }
  }
  hipdlp::backend::WorkspaceHandle empty_workspace;
  result = valid;
  CHECK_FALSE(fixture.backend->prepareCheckWorkspace(
                  fixture.problem, lp, empty_workspace).ok());
  CHECK_FALSE(fixture.backend->evaluateCheck(
                  fixture.problem, empty_workspace, fixture.input, result).ok());
  checkEmptyPacket(result);

  hipdlp::backend::ProblemHandle other_bound;
  hipdlp::backend::WorkspaceHandle other_workspace;
  REQUIRE(other.backend->bindProblem(fixture.problem, other_bound).ok());
  REQUIRE(other.backend->createWorkspace(other_bound, other_workspace).ok());
  REQUIRE(other.backend->prepareCheckWorkspace(fixture.problem, lp, other_workspace).ok());
  result = valid;
  CHECK_FALSE(fixture.backend->prepareCheckWorkspace(
                  fixture.problem, lp, other_workspace).ok());
  CHECK_FALSE(fixture.backend->evaluateCheck(
                  fixture.problem, other_workspace, fixture.input, result).ok());
  checkEmptyPacket(result);

  const HighsLp wrong_lp = makeEqualityLp(2, 3, {0, 0, 0, 0}, {}, {});
  hipdlp::model::ProblemBuildResult wrong = hipdlp::model::buildProblem(wrong_lp);
  REQUIRE(wrong.ok());
  result = valid;
  CHECK_FALSE(fixture.backend->evaluateCheck(
                  wrong.problem, fixture.workspace, fixture.input, result).ok());
  checkEmptyPacket(result);
  CHECK_FALSE(fixture.backend->prepareCheckWorkspace(
                  wrong.problem, wrong_lp, fixture.workspace).ok());
  CHECK_FALSE(fixture.backend->evaluateCheck(
                  fixture.problem, fixture.workspace, fixture.input, result).ok());
  checkEmptyPacket(result);
  fixture.prepare(lp);
  fixture.bound.reset();
  REQUIRE(fixture.backend->evaluateCheck(
              fixture.problem, fixture.workspace, fixture.input, result).ok());
  CHECK(result.valid);
}

void checkResidentObjectiveExtremes(const hipdlp::backend::Kind kind) {
  const std::vector<std::vector<double>> values = {
      {1e16, 1.0, -1e16, 0.0},
      {1e308, 1e308, -1e308, -1e308},
      {1.0, 1.0, 1.0, 1.0}};
  for (const bool primal_objective : {false, true}) {
    for (std::size_t index = 0; index < values.size(); ++index) {
      const HighsInt rows = primal_objective ? 0 : 4;
      const HighsInt cols = primal_objective ? 4 : 0;
      HighsLp lp = makeEqualityLp(rows, cols,
                                  std::vector<HighsInt>(cols + 1, 0), {}, {});
      lp.offset_ = index == 2 ? 1e16 : 0.0;
      lp.col_cost_.assign(cols, 1.0);
      lp.row_lower_.assign(rows, 1.0);
      lp.row_upper_.assign(rows, 1.0);
      ResidentCheckFixture fixture(kind, lp);
      fixture.prepare(lp);
      REQUIRE(fixture.backend->uploadVector(
                  values[index], primal_objective ? fixture.primal : fixture.dual).ok());
      const std::vector<double> primal = primal_objective ? values[index] : std::vector<double>();
      const std::vector<double> dual = primal_objective ? std::vector<double>() : values[index];
      hipdlp::algorithm::RawConvergenceMetrics oracle;
      const hipdlp::api::Error oracle_error = hipdlp::algorithm::computeRawConvergenceMetrics(
          fixture.problem, primal, dual, std::vector<double>(rows, 0.0),
          std::vector<double>(cols, 0.0), oracle);
      CHECK(oracle_error.ok() == (index != 1));
      if (oracle_error.ok()) {
        CHECK(oracle.primal_objective == lp.offset_);
        CHECK(oracle.dual_objective == lp.offset_);
        CHECK(oracle.objective_gap == 0.0);
      }
      hipdlp::backend::CheckResult result;
      const hipdlp::api::Error error = fixture.backend->evaluateCheck(
          fixture.problem, fixture.workspace, fixture.input, result);
      if (kind == hipdlp::backend::Kind::kCpu) {
        CHECK(error.ok() == oracle_error.ok());
        if (oracle_error.ok()) {
          REQUIRE(result.valid);
          checkRawMetrics(result.raw, oracle);
        } else {
          checkEmptyPacket(result);
        }
      } else {
        REQUIRE(error.ok());
        CHECK_FALSE(result.valid);
        CHECK(result.requires_host_oracle);
        CHECK(result.scalar_count == 0);
        checkRawMetrics(result.raw, hipdlp::algorithm::RawConvergenceMetrics());
        CHECK(result.primal_delta_norm == 0.0);
        CHECK(result.dual_delta_norm == 0.0);
        CHECK(result.delta_dot == 0.0);
        CHECK(result.primal_anchor_distance == 0.0);
        CHECK(result.dual_anchor_distance == 0.0);
        CHECK(result.device_to_host_transfers == 1);
        CHECK(result.device_to_host_bytes == 12 * sizeof(double) + 4 * sizeof(int));
        CHECK(result.synchronizations == 1);
      }
    }
  }
}

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

TEST_CASE(
    "hi-pdlp-cpu-certificate-screen-preserves-verifier-accepted-candidates",
    "[pdlp][hipdlp][phase8][backend][check-packet][certificate-screen][cpu]") {
  checkAcceptedCertificateCandidates(hipdlp::backend::Kind::kCpu);
}

TEST_CASE(
    "hi-pdlp-cpu-certificate-screen-excludes-positive-and-negative-cost-"
    "feasible-points",
    "[pdlp][hipdlp][phase8][backend][check-packet][certificate-screen][cpu]") {
  checkFeasibleCertificateExclusions(hipdlp::backend::Kind::kCpu);
}

TEST_CASE(
    "hi-pdlp-cpu-certificate-screen-preserves-rounding-sensitive-candidates",
    "[pdlp][hipdlp][phase8][backend][check-packet][certificate-screen][cpu]") {
  checkCertificateRoundingCandidates(hipdlp::backend::Kind::kCpu);
}

TEST_CASE(
    "hi-pdlp-cpu-certificate-screen-uncertainty-keeps-downloads-possible",
    "[pdlp][hipdlp][phase8][backend][check-packet][certificate-screen][cpu]") {
  checkCertificateUncertainRecovery(hipdlp::backend::Kind::kCpu);
}

TEST_CASE("hi-pdlp-check-packet-defaults-are-invalid",
          "[pdlp][hipdlp][phase8][backend][check-packet][cpu]") {
  hipdlp::backend::CheckInput input;
  CHECK(input.primal == nullptr);
  CHECK(input.dual == nullptr);
  CHECK(input.ax == nullptr);
  CHECK(input.aty == nullptr);
  CHECK(input.delta_primal == nullptr);
  CHECK(input.delta_dual == nullptr);
  CHECK(input.aty_delta == nullptr);
  CHECK(input.primal_anchor == nullptr);
  CHECK(input.dual_anchor == nullptr);
  CHECK_FALSE(input.fixed_point);
  CHECK_FALSE(input.anchor_distances);
  CHECK(input.primal_tolerance == 0.0);
  CHECK(input.dual_tolerance == 0.0);
  CHECK(input.gap_tolerance == 0.0);
  checkEmptyPacket(hipdlp::backend::CheckResult());
}

TEST_CASE("hi-pdlp-cpu-check-packet-matches-raw-oracle",
          "[pdlp][hipdlp][phase8][backend][check-packet][cpu]") {
  checkResidentRawOracle(hipdlp::backend::Kind::kCpu);
}

TEST_CASE("hi-pdlp-cpu-check-packet-optional-primitives",
          "[pdlp][hipdlp][phase8][backend][check-packet][cpu]") {
  checkResidentOptionalPrimitives(hipdlp::backend::Kind::kCpu);
}

TEST_CASE("hi-pdlp-cpu-check-fixed-point-preparation",
          "[pdlp][hipdlp][phase8][backend][check-packet][optional-error][cpu]") {
  checkResidentFixedPointPreparation(hipdlp::backend::Kind::kCpu);
}

TEST_CASE("hi-pdlp-cpu-check-packet-empty-rectangular-shapes",
          "[pdlp][hipdlp][phase8][backend][check-packet][cpu]") {
  checkResidentEmptyShapes(hipdlp::backend::Kind::kCpu);
}

TEST_CASE("hi-pdlp-cpu-check-packet-overflow-safe-norms",
          "[pdlp][hipdlp][phase8][backend][check-packet][cpu]") {
  checkResidentStableNorms(hipdlp::backend::Kind::kCpu);
}

TEST_CASE("hi-pdlp-cpu-check-packet-oracle-norm-overflow-boundary",
          "[pdlp][hipdlp][phase8][backend][check-packet][oracle-rounding][cpu]") {
  checkResidentOracleOverflow(hipdlp::backend::Kind::kCpu, false);
}

TEST_CASE("hi-pdlp-cpu-check-packet-oracle-fused-product-overflow",
          "[pdlp][hipdlp][phase8][backend][check-packet][oracle-rounding][cpu]") {
  checkResidentOracleOverflow(hipdlp::backend::Kind::kCpu, true);
}

TEST_CASE("hi-pdlp-cpu-check-packet-oracle-dual-product-overflow",
          "[pdlp][hipdlp][phase8][backend][check-packet][oracle-rounding][cpu]") {
  checkResidentDualProductOverflow(hipdlp::backend::Kind::kCpu);
}

TEST_CASE("hi-pdlp-cpu-check-packet-reduced-cost-overflow-defers-to-oracle",
          "[pdlp][hipdlp][phase8][backend][check-packet][oracle-rounding][cpu]") {
  checkResidentReducedCostOverflow(hipdlp::backend::Kind::kCpu);
}

TEST_CASE("hi-pdlp-cpu-check-packet-inclusive-norm-threshold",
          "[pdlp][hipdlp][phase8][backend][check-packet][stopping-rounding][cpu]") {
  checkResidentStoppingRounding(hipdlp::backend::Kind::kCpu);
}

TEST_CASE("hi-pdlp-cpu-check-packet-validates-stopping-tolerances",
          "[pdlp][hipdlp][phase8][backend][check-packet][check-tolerances][cpu]") {
  checkResidentToleranceValidation(hipdlp::backend::Kind::kCpu);
}

TEST_CASE("hi-pdlp-cpu-check-packet-numerical-failures",
          "[pdlp][hipdlp][phase8][backend][check-packet][cpu]") {
  checkResidentNumericalFailures(hipdlp::backend::Kind::kCpu);
}

TEST_CASE("hi-pdlp-cpu-check-packet-preserves-serial-objective-semantics",
          "[pdlp][hipdlp][phase8][backend][check-packet][cpu]") {
  checkResidentObjectiveExtremes(hipdlp::backend::Kind::kCpu);
}

TEST_CASE("hi-pdlp-cpu-check-packet-rejects-invalid-handles-transactionally",
          "[pdlp][hipdlp][phase8][backend][check-packet][cpu]") {
  checkResidentInvalidHandles(hipdlp::backend::Kind::kCpu);
}

#ifdef HIGHS_HIPDLP_CUDA_DETECTED
TEST_CASE(
    "hi-pdlp-cuda-certificate-screen-preserves-verifier-accepted-candidates",
    "[pdlp][hipdlp][phase8][backend][check-packet][certificate-screen][cuda]") {
  checkAcceptedCertificateCandidates(hipdlp::backend::Kind::kCuda);
}

TEST_CASE(
    "hi-pdlp-cuda-certificate-screen-excludes-positive-and-negative-cost-"
    "feasible-points",
    "[pdlp][hipdlp][phase8][backend][check-packet][certificate-screen][cuda]") {
  checkFeasibleCertificateExclusions(hipdlp::backend::Kind::kCuda);
}

TEST_CASE(
    "hi-pdlp-cuda-certificate-screen-preserves-rounding-sensitive-candidates",
    "[pdlp][hipdlp][phase8][backend][check-packet][certificate-screen][cuda]") {
  checkCertificateRoundingCandidates(hipdlp::backend::Kind::kCuda);
}

TEST_CASE(
    "hi-pdlp-cuda-certificate-screen-uncertainty-keeps-downloads-possible",
    "[pdlp][hipdlp][phase8][backend][check-packet][certificate-screen][cuda]") {
  checkCertificateUncertainRecovery(hipdlp::backend::Kind::kCuda);
}

TEST_CASE(
    "hi-pdlp-cuda-certificate-screen-owns-original-lp-data",
    "[pdlp][hipdlp][phase8][backend][check-packet][certificate-screen][cuda]") {
  HighsLp lp = makeEqualityLp(1, 1, {0, 1}, {0}, {1.0});
  lp.col_cost_[0] = -1.0;
  lp.col_lower_[0] = 0.0;
  lp.row_lower_[0] = lp.row_upper_[0] = 1.0;
  ResidentCheckFixture fixture(hipdlp::backend::Kind::kCuda, lp);
  fixture.prepare(lp);
  REQUIRE(fixture.backend->uploadVector({1.0}, fixture.primal).ok());
  REQUIRE(fixture.backend->uploadVector({-1.0}, fixture.dual).ok());
  REQUIRE(fixture.backend->uploadVector({1.0}, fixture.ax).ok());
  REQUIRE(fixture.backend->uploadVector({-1.0}, fixture.aty).ok());
  lp = HighsLp();
  fixture.bound.reset();
  hipdlp::backend::CheckResult result;
  REQUIRE(fixture.backend
              ->evaluateCheck(fixture.problem, fixture.workspace, fixture.input,
                              result)
              .ok());
  CHECK_FALSE(result.primal_ray_possible);
  CHECK_FALSE(result.farkas_possible);
  CHECK(result.device_to_host_transfers == 1);
  CHECK(result.device_to_host_bytes == 12 * sizeof(double) + 4 * sizeof(int));
  CHECK(result.synchronizations == 1);
}

TEST_CASE("hi-pdlp-cuda-check-packet-rejects-invalid-handles-transactionally",
          "[pdlp][hipdlp][phase8][backend][check-packet][cuda]") {
  checkResidentInvalidHandles(hipdlp::backend::Kind::kCuda);
}

TEST_CASE("hi-pdlp-cuda-check-packet-matches-raw-oracle",
          "[pdlp][hipdlp][phase8][backend][check-packet][cuda]") {
  checkResidentRawOracle(hipdlp::backend::Kind::kCuda);
}

TEST_CASE("hi-pdlp-cuda-check-packet-optional-primitives",
          "[pdlp][hipdlp][phase8][backend][check-packet][cuda]") {
  checkResidentOptionalPrimitives(hipdlp::backend::Kind::kCuda);
}

TEST_CASE("hi-pdlp-cuda-check-fixed-point-preparation",
          "[pdlp][hipdlp][phase8][backend][check-packet][optional-error][cuda]") {
  checkResidentFixedPointPreparation(hipdlp::backend::Kind::kCuda);
}

TEST_CASE("hi-pdlp-cuda-check-preparation-validates-before-writing",
          "[pdlp][hipdlp][phase8][backend][check-packet][optional-error][cuda]") {
  HighsLp lp = makeEqualityLp(1, 1, {0, 1}, {0}, {1.0});
  ResidentCheckFixture fixture(hipdlp::backend::Kind::kCuda, lp);
  ResidentCheckFixture foreign(hipdlp::backend::Kind::kCuda, lp);
  hipdlp::backend::VectorHandle invalid;
  hipdlp::backend::VectorHandle wrong_size;
  REQUIRE(fixture.backend->createVector(2, wrong_size).ok());
  CHECK_FALSE(fixture.backend->prepareCheckFixedPoint(
                  fixture.bound, fixture.workspace, fixture.primal, fixture.dual,
                  fixture.primal_anchor, fixture.dual_anchor, fixture.delta_primal,
                  fixture.delta_dual, fixture.aty_delta).ok());
  fixture.prepare(lp);
  REQUIRE(fixture.backend->fill(fixture.delta_primal, 7.0).ok());
  REQUIRE(fixture.backend->fill(fixture.delta_dual, 8.0).ok());
  REQUIRE(fixture.backend->fill(fixture.aty_delta, 9.0).ok());
  hipdlp::backend::VectorHandle* handles[] = {
      &fixture.primal, &fixture.dual, &fixture.primal_anchor,
      &fixture.dual_anchor, &fixture.delta_primal, &fixture.delta_dual,
      &fixture.aty_delta};
  const auto prepare = [&]() {
    return fixture.backend->prepareCheckFixedPoint(
        fixture.bound, fixture.workspace, *handles[0], *handles[1], *handles[2],
        *handles[3], *handles[4], *handles[5], *handles[6]);
  };
  const auto before = hipdlp::backend::cuda::testing::resourceAudit();
  for (std::size_t index = 0; index < 7; ++index) {
    hipdlp::backend::VectorHandle* original = handles[index];
    for (hipdlp::backend::VectorHandle* replacement :
         {&invalid, &wrong_size, &foreign.primal}) {
      handles[index] = replacement;
      CHECK(prepare().code == hipdlp::api::ErrorCode::kInvalidArgument);
    }
    handles[index] = original;
  }
  handles[4] = handles[0];
  CHECK(prepare().code == hipdlp::api::ErrorCode::kInvalidArgument);
  CHECK(fixture.backend->prepareCheckFixedPoint(
              foreign.bound, fixture.workspace, fixture.primal, fixture.dual,
              fixture.primal_anchor, fixture.dual_anchor, fixture.delta_primal,
              fixture.delta_dual, fixture.aty_delta).code ==
        hipdlp::api::ErrorCode::kInvalidArgument);
  CHECK(fixture.backend->prepareCheckFixedPoint(
              fixture.bound, foreign.workspace, fixture.primal, fixture.dual,
              fixture.primal_anchor, fixture.dual_anchor, fixture.delta_primal,
              fixture.delta_dual, fixture.aty_delta).code ==
        hipdlp::api::ErrorCode::kInvalidArgument);
  const auto after = hipdlp::backend::cuda::testing::resourceAudit();
  CHECK(after.direct_primitive_enqueue_calls == before.direct_primitive_enqueue_calls);
  CHECK(after.device_to_device_transfers == before.device_to_device_transfers);
  CHECK(after.device_to_host_transfers == before.device_to_host_transfers);
  CHECK(after.synchronization_calls == before.synchronization_calls);
  CHECK(download(*fixture.backend, fixture.delta_primal) == std::vector<double>{7.0});
  CHECK(download(*fixture.backend, fixture.delta_dual) == std::vector<double>{8.0});
  CHECK(download(*fixture.backend, fixture.aty_delta) == std::vector<double>{9.0});
  hipdlp::backend::CheckResult result;
  REQUIRE(fixture.backend->evaluateCheck(
              fixture.problem, fixture.workspace, fixture.input, result).ok());
  CHECK(result.valid);
}

TEST_CASE("hi-pdlp-cuda-check-packet-empty-rectangular-shapes",
          "[pdlp][hipdlp][phase8][backend][check-packet][cuda]") {
  checkResidentEmptyShapes(hipdlp::backend::Kind::kCuda);
}

TEST_CASE("hi-pdlp-cuda-check-packet-overflow-safe-norms",
          "[pdlp][hipdlp][phase8][backend][check-packet][cuda]") {
  checkResidentStableNorms(hipdlp::backend::Kind::kCuda);
}

TEST_CASE("hi-pdlp-cuda-check-packet-oracle-norm-overflow-boundary",
          "[pdlp][hipdlp][phase8][backend][check-packet][oracle-rounding][cuda]") {
  checkResidentOracleOverflow(hipdlp::backend::Kind::kCuda, false);
}

TEST_CASE("hi-pdlp-cuda-check-packet-oracle-fused-product-overflow",
          "[pdlp][hipdlp][phase8][backend][check-packet][oracle-rounding][cuda]") {
  checkResidentOracleOverflow(hipdlp::backend::Kind::kCuda, true);
}

TEST_CASE("hi-pdlp-cuda-check-packet-oracle-dual-product-overflow",
          "[pdlp][hipdlp][phase8][backend][check-packet][oracle-rounding][cuda]") {
  checkResidentDualProductOverflow(hipdlp::backend::Kind::kCuda);
}

TEST_CASE("hi-pdlp-cuda-check-packet-reduced-cost-overflow-defers-to-oracle",
          "[pdlp][hipdlp][phase8][backend][check-packet][oracle-rounding][cuda]") {
  checkResidentReducedCostOverflow(hipdlp::backend::Kind::kCuda);
}

TEST_CASE("hi-pdlp-cuda-check-packet-inclusive-norm-threshold",
          "[pdlp][hipdlp][phase8][backend][check-packet][stopping-rounding][cuda]") {
  checkResidentStoppingRounding(hipdlp::backend::Kind::kCuda);
}

TEST_CASE("hi-pdlp-cuda-check-packet-validates-stopping-tolerances",
          "[pdlp][hipdlp][phase8][backend][check-packet][check-tolerances][cuda]") {
  checkResidentToleranceValidation(hipdlp::backend::Kind::kCuda);
}

TEST_CASE("hi-pdlp-cuda-check-packet-numerical-failures",
          "[pdlp][hipdlp][phase8][backend][check-packet][cuda]") {
  checkResidentNumericalFailures(hipdlp::backend::Kind::kCuda);
}

TEST_CASE("hi-pdlp-cuda-check-optional-scratch-error-is-deferred",
          "[pdlp][hipdlp][phase8][backend][check-packet][optional-error][cuda]") {
  const double maximum = std::numeric_limits<double>::max();
  HighsLp lp = makeEqualityLp(1, 1, {0, 1}, {0}, {1.0});
  ResidentCheckFixture fixture(hipdlp::backend::Kind::kCuda, lp);
  fixture.prepare(lp);
  REQUIRE(fixture.backend->fill(fixture.primal, maximum).ok());
  REQUIRE(fixture.backend->multiply(fixture.bound, fixture.primal,
                                     fixture.ax).ok());
  REQUIRE(fixture.backend->fill(fixture.primal_anchor, -maximum).ok());
  hipdlp::backend::CheckResult raw;
  REQUIRE(fixture.backend->evaluateCheck(fixture.problem, fixture.workspace,
                                         fixture.input, raw).ok());
  REQUIRE(raw.valid);
  CHECK(raw.raw.primal_residual == maximum);

  bool prior_error = false;
  SECTION("optional scratch overflow leaves mandatory inputs downloadable") {}
  SECTION("pending mandatory error survives optional preparation") {
    prior_error = true;
    REQUIRE(fixture.backend->scale(2.0, fixture.ax).ok());
    REQUIRE(fixture.backend->fill(fixture.ax, maximum).ok());
  }

  fixture.input.fixed_point = true;
  fixture.input.delta_primal = &fixture.delta_primal;
  fixture.input.delta_dual = &fixture.delta_dual;
  fixture.input.aty_delta = &fixture.aty_delta;
  const auto before = hipdlp::backend::cuda::testing::resourceAudit();
  REQUIRE(fixture.backend->prepareCheckFixedPoint(
              fixture.bound, fixture.workspace, fixture.primal, fixture.dual,
              fixture.primal_anchor, fixture.dual_anchor, fixture.delta_primal,
              fixture.delta_dual, fixture.aty_delta).ok());
  hipdlp::backend::CheckResult result;
  const hipdlp::api::Error error = fixture.backend->evaluateCheck(
      fixture.problem, fixture.workspace, fixture.input, result);
  const auto after = hipdlp::backend::cuda::testing::resourceAudit();
  CHECK(after.allocation_calls == before.allocation_calls);
  CHECK(after.pinned_allocation_calls == before.pinned_allocation_calls);
  CHECK(after.host_to_device_transfers == before.host_to_device_transfers);
  CHECK(after.device_to_host_transfers - before.device_to_host_transfers == 1);
  CHECK(after.device_to_host_bytes - before.device_to_host_bytes ==
        12 * sizeof(double) + 4 * sizeof(int));
  CHECK(after.synchronization_calls - before.synchronization_calls == 1);
  std::vector<double> downloaded = {42.0};
  if (prior_error) {
    CHECK_FALSE(error.ok());
    checkEmptyPacket(result);
    const hipdlp::api::Error preparation_error =
        fixture.backend->prepareCheckFixedPoint(
            fixture.bound, fixture.workspace, fixture.primal, fixture.dual,
            fixture.primal_anchor, fixture.dual_anchor, fixture.delta_primal,
            fixture.delta_dual, fixture.aty_delta);
    CHECK(preparation_error.code == error.code);
    CHECK(preparation_error.message == error.message);
    const hipdlp::api::Error repeated = fixture.backend->evaluateCheck(
        fixture.problem, fixture.workspace, fixture.input, result);
    CHECK(repeated.code == error.code);
    CHECK(repeated.message == error.message);
    checkEmptyPacket(result);
    CHECK_FALSE(fixture.backend->downloadVector(fixture.primal, downloaded).ok());
    CHECK(downloaded == std::vector<double>{42.0});
  } else {
    INFO(error.message);
    CHECK(error.ok());
    CHECK_FALSE(result.valid);
    CHECK(result.requires_host_oracle);
    CHECK(result.primal_ray_possible);
    CHECK(result.farkas_possible);
    CHECK(result.scalar_count == 0);
    CHECK(result.raw.primal_residual == 0.0);
    CHECK(result.raw.dual_residual == 0.0);
    CHECK(result.raw.primal_objective == 0.0);
    CHECK(result.raw.dual_objective == 0.0);
    CHECK(result.raw.objective_gap == 0.0);
    CHECK(result.raw.unscaled_rhs_norm == 0.0);
    CHECK(result.raw.unscaled_cost_norm == 0.0);
    CHECK(result.primal_delta_norm == 0.0);
    CHECK(result.dual_delta_norm == 0.0);
    CHECK(result.delta_dot == 0.0);
    CHECK(result.primal_anchor_distance == 0.0);
    CHECK(result.dual_anchor_distance == 0.0);
    CHECK(result.device_to_host_transfers == 1);
    CHECK(result.device_to_host_bytes == 12 * sizeof(double) + 4 * sizeof(int));
    CHECK(result.synchronizations == 1);
    CHECK(fixture.backend->downloadVector(fixture.primal, downloaded).ok());
    CHECK(downloaded == std::vector<double>{maximum});
    REQUIRE(fixture.backend->copy(fixture.primal, fixture.primal_anchor).ok());
    REQUIRE(fixture.backend->prepareCheckFixedPoint(
                fixture.bound, fixture.workspace, fixture.primal, fixture.dual,
                fixture.primal_anchor, fixture.dual_anchor, fixture.delta_primal,
                fixture.delta_dual, fixture.aty_delta).ok());
    REQUIRE(fixture.backend->evaluateCheck(
                fixture.problem, fixture.workspace, fixture.input, result).ok());
    CHECK(result.valid);
    CHECK_FALSE(result.requires_host_oracle);
    CHECK(result.raw.primal_residual == maximum);
    CHECK(result.primal_delta_norm == 0.0);
  }
}

TEST_CASE("hi-pdlp-cuda-check-optional-reduction-overflow-is-deferred",
          "[pdlp][hipdlp][phase8][backend][check-packet][optional-error][cuda]") {
  const double maximum = std::numeric_limits<double>::max();
  const char* scenarios[] = {"primal delta norm", "dual delta norm", "dot sum",
                             "primal anchor norm", "dual anchor norm",
                             "dual anchor subtraction"};
  for (int scenario = 0; scenario < 6; ++scenario) {
    for (const bool mandatory_overflow : {false, true}) {
      INFO(scenarios[scenario]);
      INFO(mandatory_overflow);
      HighsLp lp = makeEqualityLp(2, 2, {0, 1, 2}, {0, 1}, {1.0, 1.0});
      ResidentCheckFixture fixture(hipdlp::backend::Kind::kCuda, lp);
      fixture.prepare(lp);
      fixture.input.fixed_point = true;
      fixture.input.anchor_distances = true;
      fixture.input.delta_primal = &fixture.delta_primal;
      fixture.input.delta_dual = &fixture.delta_dual;
      fixture.input.aty_delta = &fixture.aty_delta;
      fixture.input.primal_anchor = &fixture.primal_anchor;
      fixture.input.dual_anchor = &fixture.dual_anchor;
      if (scenario == 0 || scenario == 1) {
        REQUIRE(fixture.backend->fill(scenario == 0 ? fixture.delta_primal
                                                    : fixture.delta_dual,
                                        maximum).ok());
      } else if (scenario == 2) {
        REQUIRE(fixture.backend->fill(fixture.delta_primal, 1.0).ok());
        REQUIRE(fixture.backend->fill(fixture.aty_delta, maximum).ok());
      } else if (scenario == 3 || scenario == 4) {
        REQUIRE(fixture.backend->fill(scenario == 3 ? fixture.primal_anchor
                                                    : fixture.dual_anchor,
                                        maximum).ok());
      } else {
        REQUIRE(fixture.backend->fill(fixture.dual, maximum).ok());
        REQUIRE(fixture.backend->fill(fixture.dual_anchor, -maximum).ok());
      }
      if (mandatory_overflow)
        REQUIRE(fixture.backend->fill(fixture.ax, maximum).ok());
      hipdlp::backend::CheckResult result;
      const hipdlp::api::Error error = fixture.backend->evaluateCheck(
          fixture.problem, fixture.workspace, fixture.input, result);
      REQUIRE(error.ok());
      checkEmptyPacket(result, true);
      REQUIRE(fixture.backend->synchronize().ok());
      const auto primal = download(*fixture.backend, fixture.primal);
      const auto dual = download(*fixture.backend, fixture.dual);
      const auto ax = download(*fixture.backend, fixture.ax);
      const auto aty = download(*fixture.backend, fixture.aty);
      hipdlp::algorithm::RawConvergenceMetrics oracle;
      CHECK(hipdlp::algorithm::computeRawConvergenceMetrics(
                fixture.problem, primal, dual, ax, aty, oracle).ok() ==
            !mandatory_overflow);
      if (mandatory_overflow)
        REQUIRE(fixture.backend->fill(fixture.ax, 0.0).ok());
      fixture.input.fixed_point = false;
      fixture.input.anchor_distances = false;
      REQUIRE(fixture.backend->evaluateCheck(
                  fixture.problem, fixture.workspace, fixture.input, result).ok());
      REQUIRE(result.valid);
      double optional_value = 0.0;
      hipdlp::api::Error late_error;
      if (scenario == 0 || scenario == 1)
        late_error = fixture.backend->norm2(
            scenario == 0 ? fixture.delta_primal : fixture.delta_dual,
            optional_value);
      else if (scenario == 2)
        late_error = fixture.backend->dot(fixture.delta_primal,
                                          fixture.aty_delta, optional_value);
      else if (scenario == 3)
        late_error = fixture.backend->differenceNorm2(
            fixture.primal, fixture.primal_anchor, optional_value);
      else
        late_error = fixture.backend->differenceNorm2(
            fixture.dual, fixture.dual_anchor, optional_value);
      CHECK_FALSE(late_error.ok());
    }
  }
}

TEST_CASE("hi-pdlp-cuda-check-optional-spmv-error-is-local-to-preparation",
          "[pdlp][hipdlp][phase8][backend][check-packet][optional-error][cuda]") {
  HighsLp lp = makeEqualityLp(1, 1, {0, 1}, {0}, {2.0});
  ResidentCheckFixture fixture(hipdlp::backend::Kind::kCuda, lp);
  fixture.prepare(lp);
  REQUIRE(fixture.backend->fill(fixture.primal_anchor, -1.0).ok());
  REQUIRE(fixture.backend->fill(fixture.dual_anchor,
                                 -std::numeric_limits<double>::max()).ok());
  REQUIRE(fixture.backend->prepareCheckFixedPoint(
              fixture.bound, fixture.workspace, fixture.primal, fixture.dual,
              fixture.primal_anchor, fixture.dual_anchor, fixture.delta_primal,
              fixture.delta_dual, fixture.aty_delta).ok());
  REQUIRE(fixture.backend->fill(fixture.delta_dual, 0.0).ok());
  REQUIRE(fixture.backend->fill(fixture.aty_delta, 0.0).ok());
  fixture.input.fixed_point = true;
  fixture.input.delta_primal = &fixture.delta_primal;
  fixture.input.delta_dual = &fixture.delta_dual;
  fixture.input.aty_delta = &fixture.aty_delta;
  hipdlp::backend::CheckResult result;
  REQUIRE(fixture.backend->evaluateCheck(
              fixture.problem, fixture.workspace, fixture.input, result).ok());
  checkEmptyPacket(result, true);
  CHECK(download(*fixture.backend, fixture.primal) == std::vector<double>{0.0});
  fixture.input.fixed_point = false;
  REQUIRE(fixture.backend->evaluateCheck(
              fixture.problem, fixture.workspace, fixture.input, result).ok());
  CHECK(result.valid);
  REQUIRE(fixture.backend->fill(fixture.dual_anchor, 0.0).ok());
  REQUIRE(fixture.backend->prepareCheckFixedPoint(
              fixture.bound, fixture.workspace, fixture.primal, fixture.dual,
              fixture.primal_anchor, fixture.dual_anchor, fixture.delta_primal,
              fixture.delta_dual, fixture.aty_delta).ok());
  fixture.input.fixed_point = true;
  REQUIRE(fixture.backend->evaluateCheck(
              fixture.problem, fixture.workspace, fixture.input, result).ok());
  CHECK(result.valid);
  CHECK_FALSE(result.requires_host_oracle);
  CHECK(result.primal_delta_norm == 1.0);
  CHECK(result.dual_delta_norm == 0.0);
  CHECK(result.delta_dot == 0.0);
}

TEST_CASE("hi-pdlp-cuda-check-packet-one-readback-no-new-resources",
          "[pdlp][hipdlp][phase8][backend][check-packet][cuda]") {
  using hipdlp::backend::cuda::testing::ResourceAudit;
  using hipdlp::backend::cuda::testing::resourceAudit;
  hipdlp::backend::cuda::testing::resetResourceAudit();
  for (const HighsInt size : {1, 4099, 262147}) {
    HighsLp lp = makeEqualityLp(size, size,
                                std::vector<HighsInt>(size + 1, 0), {}, {});
    ResidentCheckFixture fixture(hipdlp::backend::Kind::kCuda, lp);
    fixture.prepare(lp);
    REQUIRE(fixture.backend->fill(fixture.primal, 1.0).ok());
    REQUIRE(fixture.backend->fill(fixture.dual, 2.0).ok());
    REQUIRE(fixture.backend->fill(fixture.ax, 3.0).ok());
    REQUIRE(fixture.backend->fill(fixture.aty, 4.0).ok());
    hipdlp::backend::CheckResult first;
    const ResourceAudit before = resourceAudit();
    for (int repetition = 0; repetition < 5; ++repetition) {
      hipdlp::backend::CheckResult result;
      REQUIRE(fixture.backend->evaluateCheck(
                  fixture.problem, fixture.workspace, fixture.input, result).ok());
      REQUIRE(result.valid);
      CHECK(result.raw.primal_residual == Approx(3.0 * std::sqrt(size)));
      CHECK(result.raw.dual_residual == Approx(4.0 * std::sqrt(size)));
      CHECK_FALSE(result.primal_ray_possible);
      CHECK_FALSE(result.farkas_possible);
      CHECK(result.scalar_count == 7);
      CHECK(result.device_to_host_transfers == 1);
      CHECK(result.synchronizations == 1);
      CHECK(result.device_to_host_bytes >= 12 * sizeof(double) + sizeof(int));
      CHECK(result.device_to_host_bytes <= 128);
      if (repetition == 0) first = result;
      CHECK(result.raw.primal_residual == first.raw.primal_residual);
      CHECK(result.raw.dual_residual == first.raw.dual_residual);
      CHECK(result.primal_ray_possible == first.primal_ray_possible);
      CHECK(result.farkas_possible == first.farkas_possible);
      CHECK(result.device_to_host_bytes == first.device_to_host_bytes);
    }
    const ResourceAudit after = resourceAudit();
    CHECK(after.device_to_host_transfers - before.device_to_host_transfers == 5);
    CHECK(after.device_to_host_bytes - before.device_to_host_bytes ==
          5 * first.device_to_host_bytes);
    CHECK(after.synchronization_calls - before.synchronization_calls == 5);
    CHECK(after.allocation_calls == before.allocation_calls);
    CHECK(after.free_calls == before.free_calls);
    CHECK(after.pinned_allocation_calls == before.pinned_allocation_calls);
    CHECK(after.pinned_free_calls == before.pinned_free_calls);
    CHECK(after.host_to_device_transfers == before.host_to_device_transfers);
    CHECK(after.stream_creations == before.stream_creations);
    CHECK(after.cublas_creations == before.cublas_creations);
    CHECK(after.cusparse_creations == before.cusparse_creations);
    CHECK(after.descriptor_creations == before.descriptor_creations);
    CHECK(after.graph_capture_calls == before.graph_capture_calls);
  }
  const ResourceAudit released = resourceAudit();
  CHECK(released.pinned_allocation_calls > 0);
  CHECK(released.pinned_allocation_calls == released.pinned_free_calls);
  CHECK(released.allocation_calls == released.free_calls);
}

TEST_CASE("hi-pdlp-cuda-check-packet-preparation-failure-releases-resources",
          "[pdlp][hipdlp][phase8][backend][check-packet][cuda]") {
  using hipdlp::backend::cuda::testing::FailurePoint;
  using hipdlp::backend::cuda::testing::resourceAudit;
  using hipdlp::backend::cuda::testing::setFailurePoint;
  for (const FailurePoint failure : {FailurePoint::kDeviceAllocation,
                                     FailurePoint::kPinnedAllocation}) {
    HighsLp lp = makeEqualityLp(1, 1, {0, 1}, {0}, {1.0});
    ResidentCheckFixture fixture(hipdlp::backend::Kind::kCuda, lp);
    const auto before = resourceAudit();
    setFailurePoint(failure);
    const hipdlp::api::Error error = fixture.backend->prepareCheckWorkspace(
        fixture.problem, lp, fixture.workspace);
    setFailurePoint(FailurePoint::kNone);
    CHECK(error.code == hipdlp::api::ErrorCode::kAllocationFailure);
    const auto after = resourceAudit();
    CHECK(after.allocation_calls - before.allocation_calls ==
          after.free_calls - before.free_calls);
    CHECK(after.pinned_allocation_calls - before.pinned_allocation_calls ==
          after.pinned_free_calls - before.pinned_free_calls);
    hipdlp::backend::CheckResult result;
    CHECK_FALSE(fixture.backend->evaluateCheck(
                    fixture.problem, fixture.workspace, fixture.input, result).ok());
    checkEmptyPacket(result);
    fixture.prepare(lp);
    REQUIRE(fixture.backend->evaluateCheck(
                fixture.problem, fixture.workspace, fixture.input, result).ok());
    CHECK(result.valid);
  }
}

TEST_CASE("hi-pdlp-cuda-check-packet-preserves-pending-sticky-error",
          "[pdlp][hipdlp][phase8][backend][check-packet][cuda]") {
  HighsLp lp = makeEqualityLp(1, 1, {0, 1}, {0}, {1.0});
  ResidentCheckFixture fixture(hipdlp::backend::Kind::kCuda, lp);
  fixture.prepare(lp);
  REQUIRE(fixture.backend->fill(fixture.primal,
                                 std::numeric_limits<double>::max()).ok());
  REQUIRE(fixture.backend->scale(2.0, fixture.primal).ok());
  REQUIRE(fixture.backend->fill(fixture.primal, 0.0).ok());
  const auto before = hipdlp::backend::cuda::testing::resourceAudit();
  hipdlp::backend::CheckResult result;
  const hipdlp::api::Error first = fixture.backend->evaluateCheck(
      fixture.problem, fixture.workspace, fixture.input, result);
  CHECK_FALSE(first.ok());
  checkEmptyPacket(result);
  const auto after = hipdlp::backend::cuda::testing::resourceAudit();
  CHECK(after.device_to_host_transfers - before.device_to_host_transfers == 1);
  CHECK(after.synchronization_calls - before.synchronization_calls == 1);
  const hipdlp::api::Error repeated = fixture.backend->evaluateCheck(
      fixture.problem, fixture.workspace, fixture.input, result);
  CHECK(repeated.code == first.code);
  CHECK(repeated.message == first.message);
  checkEmptyPacket(result);
}

TEST_CASE("hi-pdlp-cuda-check-packet-requests-serial-objective-oracle",
          "[pdlp][hipdlp][phase8][backend][check-packet][cuda]") {
  checkResidentObjectiveExtremes(hipdlp::backend::Kind::kCuda);
}
#endif

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

void checkProjectionReflection(const hipdlp::backend::Kind kind) {
  HighsLp lp = makeEqualityLp(2, 2, {0, 1, 2}, {0, 1}, {1.0, 1.0});
  lp.col_cost_ = {1.0, -2.0};
  lp.col_lower_ = {0.0, -kHighsInf};
  lp.col_upper_ = {2.0, kHighsInf};
  lp.row_lower_ = {1.0, 2.0};
  lp.row_upper_ = {1.0, kHighsInf};
  auto built = hipdlp::model::buildProblem(lp);
  REQUIRE(built.ok());
  auto factory = hipdlp::backend::createBackend(kind);
  REQUIRE(factory.ok());
  auto& backend = *factory.backend;
  hipdlp::backend::ProblemHandle problem;
  hipdlp::backend::WorkspaceHandle workspace;
  REQUIRE(backend.bindProblem(built.problem, problem).ok());
  REQUIRE(backend.createWorkspace(problem, workspace).ok());
  hipdlp::backend::VectorHandle current, product, candidate, reflected;
  for (auto* vector : {&current, &product, &candidate, &reflected})
    REQUIRE(backend.createVector(2, *vector).ok());
  REQUIRE(backend.uploadVector({1.0, 3.0}, current).ok());
  REQUIRE(backend.uploadVector({4.0, -1.0}, product).ok());

  SECTION("projected candidates and reflections") {
    REQUIRE(backend.fusedProjectedPrimalUpdateAndReflect(
        problem, workspace, current, product, 0.5, candidate, reflected).ok());
    CHECK(download(backend, candidate) == std::vector<double>({2.0, 3.5}));
    CHECK(download(backend, reflected) == std::vector<double>({3.0, 4.0}));
    CHECK(download(backend, current) == std::vector<double>({1.0, 3.0}));
    REQUIRE(backend.uploadVector({0.5, -1.0}, current).ok());
    REQUIRE(backend.uploadVector({0.0, 3.0}, product).ok());
    REQUIRE(backend.fusedProjectedDualUpdateAndReflect(
        problem, workspace, current, product, 2.0, candidate, reflected).ok());
    CHECK(download(backend, candidate) == std::vector<double>({2.5, 0.0}));
    CHECK(download(backend, reflected) == std::vector<double>({4.5, 1.0}));
  }

  SECTION("invalid arguments preserve inputs and outputs") {
    REQUIRE(backend.fill(candidate, 13.0).ok());
    REQUIRE(backend.fill(reflected, 17.0).ok());
    using hipdlp::api::ErrorCode;
    for (const double step : {-1.0, kHighsInf,
                             std::numeric_limits<double>::quiet_NaN()}) {
      CHECK(backend.fusedProjectedPrimalUpdateAndReflect(
          problem, workspace, current, product, step, candidate, reflected).code ==
            ErrorCode::kInvalidArgument);
      CHECK(backend.fusedProjectedDualUpdateAndReflect(
          problem, workspace, current, product, step, candidate, reflected).code ==
            ErrorCode::kInvalidArgument);
    }
    for (auto* aliased : {&current, &product, &candidate}) {
      CHECK(backend.fusedProjectedPrimalUpdateAndReflect(
          problem, workspace, current, product, 0.5, candidate, *aliased).code ==
            ErrorCode::kInvalidArgument);
    }
    for (auto* aliased : {&current, &product, &reflected}) {
      CHECK(backend.fusedProjectedDualUpdateAndReflect(
          problem, workspace, current, product, 0.5, *aliased, reflected).code ==
            ErrorCode::kInvalidArgument);
    }
    hipdlp::backend::VectorHandle short_output;
    REQUIRE(backend.createVector(1, short_output).ok());
    CHECK(backend.fusedProjectedPrimalUpdateAndReflect(
        problem, workspace, current, product, 0.5, short_output, reflected).code ==
          ErrorCode::kInvalidArgument);
    CHECK(backend.fusedProjectedDualUpdateAndReflect(
        problem, workspace, current, product, 0.5, short_output, reflected).code ==
          ErrorCode::kInvalidArgument);
    CHECK(download(backend, candidate) == std::vector<double>({13.0, 13.0}));
    CHECK(download(backend, reflected) == std::vector<double>({17.0, 17.0}));
    CHECK(download(backend, current) == std::vector<double>({1.0, 3.0}));
    CHECK(download(backend, product) == std::vector<double>({4.0, -1.0}));
  }

  SECTION("primal intermediate doubling overflow is rejected") {
    REQUIRE(backend.uploadVector({1.0, std::numeric_limits<double>::max()}, current).ok());
    REQUIRE(backend.uploadVector({1.0, -2.0}, product).ok());
    const auto error = backend.fusedProjectedPrimalUpdateAndReflect(
        problem, workspace, current, product, 0.5, candidate, reflected);
    CHECK((!error.ok() || !backend.synchronize().ok()));
  }

  SECTION("dual intermediate doubling overflow is rejected") {
    REQUIRE(backend.uploadVector({std::numeric_limits<double>::max(), 0.0}, current).ok());
    REQUIRE(backend.uploadVector({1.0, 2.0}, product).ok());
    const auto error = backend.fusedProjectedDualUpdateAndReflect(
        problem, workspace, current, product, 0.5, candidate, reflected);
    CHECK((!error.ok() || !backend.synchronize().ok()));
  }
}

TEST_CASE("hi-pdlp-cpu-projection-reflection-contract",
          "[hipdlp][backend][step-fusion][cpu]") {
  checkProjectionReflection(hipdlp::backend::Kind::kCpu);
}

TEST_CASE("hi-pdlp-cpu-projection-reflection-matches-composed-primitives",
          "[hipdlp][backend][step-fusion][cpu]") {
  for (const HighsInt size : {0, 1, 257}) {
    CAPTURE(size);
    std::vector<HighsInt> starts(size + 1), indices(size);
    std::vector<double> values(size, 1.0), current(size), product(size);
    for (HighsInt index = 0; index < size; ++index) {
      starts[index] = index;
      indices[index] = index;
      current[index] = (index % 11 - 5) * 0.137;
      product[index] = (index % 7 - 3) * 0.971;
    }
    starts[size] = size;
    auto lp = makeEqualityLp(size, size, starts, indices, values);
    auto built = hipdlp::model::buildProblem(lp);
    REQUIRE(built.ok());
    hipdlp::backend::cpu::CpuBackend backend;
    hipdlp::backend::ProblemHandle problem;
    hipdlp::backend::WorkspaceHandle workspace;
    REQUIRE(backend.bindProblem(built.problem, problem).ok());
    REQUIRE(backend.createWorkspace(problem, workspace).ok());
    hipdlp::backend::VectorHandle input, matrix_product, projected, reflected,
        reference_projected, reference_reflected;
    for (auto* vector : {&input, &matrix_product, &projected, &reflected,
                         &reference_projected, &reference_reflected})
      REQUIRE(backend.createVector(size, *vector).ok());
    REQUIRE(backend.uploadVector(current, input).ok());
    REQUIRE(backend.uploadVector(product, matrix_product).ok());
    for (const double step : {0.0, 0.5, 0.731}) {
      CAPTURE(step);
      REQUIRE(backend.fusedProjectedPrimalUpdateAndReflect(problem, workspace,
          input, matrix_product, step, projected, reflected).ok());
      REQUIRE(backend.Backend::fusedProjectedPrimalUpdateAndReflect(problem,
          workspace, input, matrix_product, step, reference_projected,
          reference_reflected).ok());
      CHECK(download(backend, projected) == download(backend, reference_projected));
      CHECK(download(backend, reflected) == download(backend, reference_reflected));
      REQUIRE(backend.fusedProjectedDualUpdateAndReflect(problem, workspace,
          input, matrix_product, step, projected, reflected).ok());
      REQUIRE(backend.Backend::fusedProjectedDualUpdateAndReflect(problem,
          workspace, input, matrix_product, step, reference_projected,
          reference_reflected).ok());
      CHECK(download(backend, projected) == download(backend, reference_projected));
      CHECK(download(backend, reflected) == download(backend, reference_reflected));
    }
  }
}

#ifdef HIGHS_HIPDLP_CUDA_DETECTED
TEST_CASE("hi-pdlp-cuda-projection-reflection-contract",
          "[hipdlp][backend][step-fusion][cuda]") {
  checkProjectionReflection(hipdlp::backend::Kind::kCuda);
}
#endif

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