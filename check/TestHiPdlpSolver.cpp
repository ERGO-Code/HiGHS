#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <limits>
#include <new>
#include <string>
#include <vector>

#include "HCheckConfig.h"
#include "HConfig.h"
#include "Highs.h"
#include "catch.hpp"
#include "lp_data/HConst.h"
#include "lp_data/HighsLp.h"
#include "pdlp/hipdlp/algorithm/solver.hpp"
#include "pdlp/hipdlp/api/options.hpp"
#include "pdlp/hipdlp/api/solver.hpp"
#include "pdlp/hipdlp/backend/backend_factory.hpp"
#include "pdlp/hipdlp/verification/kkt_verifier.hpp"

#if defined(HIGHS_HIPDLP_CUDA_DETECTED) && defined(HIGHS_HIPDLP_TESTING)
#include "pdlp/hipdlp/backend/cuda/cuda_backend_testing.hpp"
#endif

namespace solve_allocation_testing {
#if defined(__linux__) && defined(__GLIBC__)
const bool available = true;
#else
const bool available = false;
#endif
thread_local bool enabled = false;
thread_local std::size_t calls = 0;
thread_local std::size_t bytes = 0;
thread_local bool fail_captured_allocation = false;
thread_local std::size_t injected_failures = 0;
}

#if defined(__linux__) && defined(__GLIBC__)
void* operator new(const std::size_t size) {
#if defined(HIGHS_HIPDLP_CUDA_DETECTED) && defined(HIGHS_HIPDLP_TESTING)
  if (solve_allocation_testing::fail_captured_allocation) {
    const auto audit = hipdlp::backend::cuda::testing::resourceAudit();
    if (audit.graph_capture_primitive_enqueue_calls > 0 &&
        audit.graph_instantiation_calls == 0) {
      solve_allocation_testing::fail_captured_allocation = false;
      ++solve_allocation_testing::injected_failures;
      throw std::bad_alloc();
    }
  }
#endif
  void* allocation = std::malloc(size == 0 ? 1 : size);
  if (allocation == nullptr) throw std::bad_alloc();
  if (solve_allocation_testing::enabled) {
    ++solve_allocation_testing::calls;
    solve_allocation_testing::bytes += size;
  }
  return allocation;
}

void* operator new[](const std::size_t size) { return ::operator new(size); }
void operator delete(void* allocation) noexcept { std::free(allocation); }
void operator delete[](void* allocation) noexcept { std::free(allocation); }
#if defined(__cpp_sized_deallocation)
void operator delete(void* allocation, std::size_t) noexcept {
  std::free(allocation);
}
void operator delete[](void* allocation, std::size_t) noexcept {
  std::free(allocation);
}
#endif
#endif

namespace {

const double kSolverTolerance = 1e-4;
const double kCudaTolerance = 1e-9;

bool cudaBackendUsable() {
  return hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda).ok();
}

HighsLp makeTinyFeasibleLp() {
  HighsLp lp;
  lp.num_col_ = 1;
  lp.num_row_ = 1;
  lp.col_cost_ = {1.0};
  lp.col_lower_ = {0.0};
  lp.col_upper_ = {kHighsInf};
  lp.row_lower_ = {1.0};
  lp.row_upper_ = {1.0};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 1;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1};
  lp.a_matrix_.index_ = {0};
  lp.a_matrix_.value_ = {1.0};
  lp.sense_ = ObjSense::kMinimize;
  return lp;
}

HighsLp makeScaledFeasibleLp() {
  HighsLp lp = makeTinyFeasibleLp();
  lp.col_cost_[0] = 3.0;
  lp.row_lower_[0] = 8.0;
  lp.row_upper_[0] = 8.0;
  lp.a_matrix_.value_[0] = 4.0;
  return lp;
}

HighsLp makeMaximizationLp() {
  HighsLp lp = makeTinyFeasibleLp();
  lp.col_cost_[0] = 3.0;
  lp.row_lower_[0] = -kHighsInf;
  lp.row_upper_[0] = 4.0;
  lp.a_matrix_.value_[0] = 2.0;
  lp.sense_ = ObjSense::kMaximize;
  lp.offset_ = 5.0;
  return lp;
}

HighsLp makeAlreadyOptimalLp() {
  HighsLp lp = makeTinyFeasibleLp();
  lp.col_cost_ = {0.0};
  lp.row_lower_ = {0.0};
  lp.row_upper_ = {kHighsInf};
  return lp;
}

// Every column is boxed, so the relative dual residual is identically zero and
// the PID primal-weight ratio is always undefined.
HighsLp makeBoxedColumnLp() {
  HighsLp lp = makeTinyFeasibleLp();
  lp.col_upper_ = {5.0};
  return lp;
}

HighsLp makeResidentDiagonalLp(const HighsInt dimension,
                               const bool ill_conditioned) {
  HighsLp lp;
  lp.num_col_ = dimension;
  lp.num_row_ = dimension;
  lp.col_lower_.assign(dimension, 0.0);
  lp.col_upper_.assign(dimension, 5.0);
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = dimension;
  lp.a_matrix_.num_row_ = dimension;
  lp.a_matrix_.start_.clear();
  for (HighsInt coordinate = 0; coordinate < dimension; ++coordinate) {
    const double coefficient =
        ill_conditioned && coordinate % 2 != 0 ? 0.01 : 1.0;
    lp.col_cost_.push_back(-coefficient);
    lp.row_lower_.push_back(coefficient);
    lp.row_upper_.push_back(coefficient);
    lp.a_matrix_.start_.push_back(coordinate);
    lp.a_matrix_.index_.push_back(coordinate);
    lp.a_matrix_.value_.push_back(coefficient);
  }
  lp.a_matrix_.start_.push_back(dimension);
  lp.sense_ = ObjSense::kMinimize;
  return lp;
}

class RecordingLogSink : public hipdlp::algorithm::SolveLogSink {
 public:
  void info(const std::string& message) override {
    messages.push_back(message);
  }
  void warning(const std::string& message) override {
    messages.push_back(message);
    warnings.push_back(message);
  }

  std::size_t countContaining(const std::string& needle) const {
    std::size_t count = 0;
    for (std::size_t index = 0; index < messages.size(); ++index) {
      if (messages[index].find(needle) != std::string::npos) ++count;
    }
    return count;
  }

  std::vector<std::string> messages;
  std::vector<std::string> warnings;
};

hipdlp::api::SolverOptions makeCpuFixedOptions() {
  hipdlp::api::SolverOptions options;
  options.implementation = hipdlp::api::SolverImplementation::kRewrite;
  options.backend = hipdlp::backend::Kind::kCpu;
  options.time_limit = kHighsInf;
  options.kkt_tolerance = kSolverTolerance;
  options.features_off = kPdlpAllFeaturesOn;
  options.iteration_limit = 10000;
  options.scaling.mode = hipdlp::model::ScalingMode::kOff;
  options.policy.step_size = hipdlp::api::StepSizePolicy::kFixed;
  options.policy.restart = hipdlp::api::RestartPolicy::kHalpern;
  options.optimality_tolerance = kSolverTolerance;
  return options;
}

#if defined(HIGHS_HIPDLP_CUDA_DETECTED) && defined(HIGHS_HIPDLP_TESTING)
TEST_CASE("hi-pdlp-production-command-buffer-matches-direct",
      "[hipdlp][solver][cuda][production-command-buffer]") {
  const HighsLp lp = makeResidentDiagonalLp(17, true);
  for (const bool pid : {false, true}) {
    for (const HighsInt limit : {0, 1, 39, 40, 41, 79, 80, 81, 121, 400}) {
  CAPTURE(pid, limit);
  auto options = makeCpuFixedOptions();
  options.backend = hipdlp::backend::Kind::kCuda;
  options.iteration_limit = limit;
  options.optimality_tolerance = 1e-14;
  options.kkt_tolerance = 1e-14;
  if (pid) options.policy.step_size = hipdlp::api::StepSizePolicy::kPid;
  hipdlp::algorithm::SolveDiagnostics direct_diagnostics;
  hipdlp::algorithm::SolveContext direct_context;
  direct_context.use_command_buffers = false;
  direct_context.diagnostics = &direct_diagnostics;
  const auto direct = hipdlp::algorithm::solve(
      options, hipdlp::model::OriginalModel(lp), direct_context);
  INFO(direct.error.message);
  REQUIRE(direct.ok());
  REQUIRE(direct.iteration_count == limit);
  CHECK(direct_diagnostics.command_buffer_replays == 0);

  hipdlp::algorithm::SolveDiagnostics recorded_diagnostics;
  hipdlp::algorithm::SolveContext recorded_context;
  recorded_context.diagnostics = &recorded_diagnostics;
  hipdlp::backend::cuda::testing::resetResourceAudit();
  const auto recorded = hipdlp::algorithm::solve(
      options, hipdlp::model::OriginalModel(lp), recorded_context);
  const auto audit = hipdlp::backend::cuda::testing::resourceAudit();
  INFO(recorded.error.message);
  REQUIRE(recorded.ok());
  CHECK(recorded.status == direct.status);
  CHECK(recorded.iteration_count == direct.iteration_count);
  CHECK(recorded.solution.col_value == direct.solution.col_value);
  CHECK(recorded.solution.row_value == direct.solution.row_value);
  CHECK(recorded.solution.col_dual == direct.solution.col_dual);
  CHECK(recorded.solution.row_dual == direct.solution.row_dual);
  CHECK(recorded_diagnostics.restarts == direct_diagnostics.restarts);
  CHECK(recorded_diagnostics.primal_weight_updates ==
    direct_diagnostics.primal_weight_updates);
  CHECK(recorded_diagnostics.restart_reference_checks ==
    direct_diagnostics.restart_reference_checks);
  CHECK(recorded_diagnostics.command_buffer_replays ==
    static_cast<std::size_t>(limit / 40));
  CHECK(recorded_diagnostics.recorded_iterations ==
    static_cast<std::size_t>(39 * (limit / 40)));
  CHECK(audit.graph_capture_calls == (limit >= 40 ? 1 : 0));
  CHECK(audit.graph_instantiation_calls == (limit >= 40 ? 1 : 0));
  CHECK(audit.graph_launch_calls == static_cast<std::size_t>(limit / 40));
  CHECK(recorded_diagnostics.resource_creation_calls_during_iterations == 0);
    }
  }
}

TEST_CASE("hi-pdlp-production-command-buffer-zero-time-and-failures",
          "[hipdlp][solver][cuda][production-command-buffer]") {
  const HighsLp lp = makeResidentDiagonalLp(17, true);
  auto options = makeCpuFixedOptions();
  options.backend = hipdlp::backend::Kind::kCuda;
  options.iteration_limit = 80;
  options.optimality_tolerance = 1e-14;
  options.kkt_tolerance = 1e-14;
  SECTION("zero time never records or replays") {
    options.time_limit = 0.0;
    hipdlp::backend::cuda::testing::resetResourceAudit();
    const auto result = hipdlp::algorithm::solve(
        options, hipdlp::model::OriginalModel(lp));
    REQUIRE(result.ok());
    CHECK(result.status == hipdlp::api::SolverStatus::kTimeLimit);
    CHECK(result.iteration_count == 0);
    const auto audit = hipdlp::backend::cuda::testing::resourceAudit();
    CHECK(audit.graph_capture_calls == 0);
    CHECK(audit.graph_launch_calls == 0);
  }
  SECTION("capture finalization launch and asynchronous errors propagate") {
    using hipdlp::backend::cuda::testing::FailurePoint;
    for (const auto failure : {FailurePoint::kGraphCapture,
                               FailurePoint::kGraphInstantiation,
                               FailurePoint::kGraphLaunch,
                               FailurePoint::kGraphAsynchronousExecution}) {
      CAPTURE(static_cast<int>(failure));
      hipdlp::backend::cuda::testing::resetResourceAudit();
      hipdlp::backend::cuda::testing::setFailurePoint(failure);
      const auto result = hipdlp::algorithm::solve(
          options, hipdlp::model::OriginalModel(lp));
      hipdlp::backend::cuda::testing::setFailurePoint(FailurePoint::kNone);
      CHECK_FALSE(result.ok());
      CHECK_FALSE(result.has_solution);
      CHECK(result.error.code == hipdlp::api::ErrorCode::kBackendFailure);
      const auto audit = hipdlp::backend::cuda::testing::resourceAudit();
      CHECK(audit.allocation_calls == audit.free_calls);
      CHECK(audit.pinned_allocation_calls == audit.pinned_free_calls);
      CHECK(audit.graph_creations == audit.graph_destructions);
      CHECK(audit.graph_exec_creations == audit.graph_exec_destructions);
      const auto recovered = hipdlp::algorithm::solve(
          options, hipdlp::model::OriginalModel(lp));
      REQUIRE(recovered.ok());
      CHECK(recovered.iteration_count == 80);
    }
  }
  SECTION("host allocation failure during capture cancels before cleanup") {
    if (!solve_allocation_testing::available) return;
    hipdlp::backend::cuda::testing::resetResourceAudit();
    solve_allocation_testing::injected_failures = 0;
    solve_allocation_testing::fail_captured_allocation = true;
    const auto failed = hipdlp::algorithm::solve(
        options, hipdlp::model::OriginalModel(lp));
    solve_allocation_testing::fail_captured_allocation = false;
    CHECK(solve_allocation_testing::injected_failures == 1);
    CHECK_FALSE(failed.ok());
    CHECK(failed.error.code == hipdlp::api::ErrorCode::kAllocationFailure);
    const auto audit = hipdlp::backend::cuda::testing::resourceAudit();
    CHECK(audit.allocation_calls == audit.free_calls);
    CHECK(audit.pinned_allocation_calls == audit.pinned_free_calls);
    const auto recovered = hipdlp::algorithm::solve(
        options, hipdlp::model::OriginalModel(lp));
    INFO(recovered.error.message);
    REQUIRE(recovered.ok());
    CHECK(recovered.iteration_count == 80);
  }
}

TEST_CASE("hi-pdlp-production-command-buffer-stops-on-expired-positive-budget",
          "[hipdlp][solver][cuda][production-command-buffer]") {
  class ExpiringBudget : public hipdlp::algorithm::SolveObserver {
   public:
    hipdlp::api::SolverOptions* options;
    HighsInt stop_iteration;
    bool expired = false;

    void onEvent(const hipdlp::algorithm::SolveObservation& event) noexcept override {
      if (event.phase == hipdlp::algorithm::SolvePhase::kStep && !event.begin &&
          event.iteration == stop_iteration) {
        options->time_limit = std::numeric_limits<double>::min();
        expired = true;
      }
    }
  } observer;
  const HighsLp lp = makeResidentDiagonalLp(17, true);
  for (const HighsInt stop_iteration : {1, 40}) {
    CAPTURE(stop_iteration);
    auto options = makeCpuFixedOptions();
    options.backend = hipdlp::backend::Kind::kCuda;
    options.iteration_limit = 400;
    options.time_limit = 60.0;
    options.optimality_tolerance = 1e-14;
    options.kkt_tolerance = 1e-14;
    observer.options = &options;
    observer.stop_iteration = stop_iteration;
    observer.expired = false;
    hipdlp::algorithm::SolveContext context;
    context.observer = &observer;
    hipdlp::algorithm::SolveDiagnostics diagnostics;
    context.diagnostics = &diagnostics;
    const auto result = hipdlp::algorithm::solve(
        options, hipdlp::model::OriginalModel(lp), context);
    REQUIRE(observer.expired);
    REQUIRE(result.ok());
    CHECK(result.status == hipdlp::api::SolverStatus::kTimeLimit);
    CHECK(result.iteration_count == stop_iteration);
    CHECK(diagnostics.command_buffer_replays == (stop_iteration == 40 ? 1 : 0));
    options.time_limit = kHighsInf;
    options.iteration_limit = stop_iteration;
    hipdlp::algorithm::SolveContext direct_context;
    direct_context.use_command_buffers = false;
    const auto direct = hipdlp::algorithm::solve(
        options, hipdlp::model::OriginalModel(lp), direct_context);
    REQUIRE(direct.ok());
    CHECK(result.solution.col_value == direct.solution.col_value);
    CHECK(result.solution.row_value == direct.solution.row_value);
    CHECK(result.solution.col_dual == direct.solution.col_dual);
    CHECK(result.solution.row_dual == direct.solution.row_dual);
  }
}
#endif

class HostPreparationObserver : public hipdlp::algorithm::SolveObserver {
 public:
  bool iterating = false;
  std::size_t preparations = 0;
  std::size_t allocations = 0;
  std::size_t allocated_bytes = 0;

  void onEvent(const hipdlp::algorithm::SolveObservation& event) noexcept
      override {
    using hipdlp::algorithm::SolvePhase;
    if (event.phase == SolvePhase::kIterations) iterating = event.begin;
    if (!iterating || event.phase != SolvePhase::kHostPreparation) return;
    if (event.begin) {
      solve_allocation_testing::calls = 0;
      solve_allocation_testing::bytes = 0;
      solve_allocation_testing::enabled = true;
    } else {
      solve_allocation_testing::enabled = false;
      ++preparations;
      allocations += solve_allocation_testing::calls;
      allocated_bytes += solve_allocation_testing::bytes;
    }
  }
};

#if defined(HIGHS_HIPDLP_CUDA_DETECTED) && defined(HIGHS_HIPDLP_TESTING)
using SolveCudaAudit = hipdlp::backend::cuda::testing::ResourceAudit;

std::size_t SolveCudaAudit::* const kSolveCudaFields[] = {
    &SolveCudaAudit::allocation_calls,
    &SolveCudaAudit::allocation_bytes,
    &SolveCudaAudit::free_calls,
    &SolveCudaAudit::free_bytes,
    &SolveCudaAudit::pinned_allocation_calls,
    &SolveCudaAudit::pinned_allocation_bytes,
    &SolveCudaAudit::pinned_free_calls,
    &SolveCudaAudit::pinned_free_bytes,
    &SolveCudaAudit::host_to_device_transfers,
    &SolveCudaAudit::host_to_device_bytes,
    &SolveCudaAudit::device_to_host_transfers,
    &SolveCudaAudit::device_to_host_bytes,
    &SolveCudaAudit::synchronization_calls};

SolveCudaAudit solveCudaDelta(const SolveCudaAudit& before,
                             const SolveCudaAudit& after) {
  SolveCudaAudit delta;
  for (const auto field : kSolveCudaFields)
    delta.*field = after.*field - before.*field;
  return delta;
}

bool noCudaResourcesOrUploads(const SolveCudaAudit& delta) {
  return delta.allocation_calls == 0 && delta.free_calls == 0 &&
         delta.pinned_allocation_calls == 0 && delta.pinned_free_calls == 0 &&
         delta.host_to_device_transfers == 0 && delta.host_to_device_bytes == 0;
}
#endif

struct ObservedPhase {
  std::size_t begins = 0;
  std::size_t ends = 0;
  std::size_t allocations = 0;
  std::size_t bytes = 0;
  std::size_t start_allocations = 0;
  std::size_t start_bytes = 0;
  HighsInt start_iteration = 0;
#if defined(HIGHS_HIPDLP_CUDA_DETECTED) && defined(HIGHS_HIPDLP_TESTING)
  SolveCudaAudit cuda_start;
  SolveCudaAudit cuda_total;
#endif
};

class ProductionSolveObserver : public hipdlp::algorithm::SolveObserver {
 public:
  ObservedPhase phases[static_cast<std::size_t>(
      hipdlp::algorithm::SolvePhase::kCount)]{};
  hipdlp::algorithm::SolvePhase stack[32]{};
  std::size_t depth = 0;
  bool balanced = true;
  bool iterating = false;
  bool boundary_active = false;
  bool boundary_special = false;
  std::size_t ordinary_boundaries = 0;
  std::size_t special_boundaries = 0;
  std::size_t terminal_boundaries = 0;
  std::size_t ordinary_allocations = 0;
  std::size_t complete_blocks = 0;
  std::size_t candidate_downloads = 0;
  std::size_t first_candidate_download_allocations = 0;
  bool audit_cuda = false;
  bool ordinary_cuda_valid = true;
#if defined(HIGHS_HIPDLP_CUDA_DETECTED) && defined(HIGHS_HIPDLP_TESTING)
  SolveCudaAudit first_candidate_cuda;
#endif

  const ObservedPhase& phase(const hipdlp::algorithm::SolvePhase value) const {
    return phases[static_cast<std::size_t>(value)];
  }

  void onEvent(const hipdlp::algorithm::SolveObservation& event) noexcept
      override {
    using hipdlp::algorithm::SolvePhase;
    ObservedPhase& current = phases[static_cast<std::size_t>(event.phase)];
    if (event.begin) {
      if (event.phase == SolvePhase::kSolve)
        solve_allocation_testing::enabled = true;
      if (depth < 32)
        stack[depth++] = event.phase;
      else
        balanced = false;
      ++current.begins;
      current.start_allocations = solve_allocation_testing::calls;
      current.start_bytes = solve_allocation_testing::bytes;
      current.start_iteration = event.iteration;
    #if defined(HIGHS_HIPDLP_CUDA_DETECTED) && defined(HIGHS_HIPDLP_TESTING)
      if (audit_cuda)
        current.cuda_start = hipdlp::backend::cuda::testing::resourceAudit();
    #endif
      if (event.phase == SolvePhase::kIterations) iterating = true;
      if (event.phase == SolvePhase::kBoundary) {
        boundary_active = true;
        boundary_special = false;
      }
      if (boundary_active &&
          (event.phase == SolvePhase::kDownload ||
           event.phase == SolvePhase::kHostVerification ||
           event.phase == SolvePhase::kFixedPointFallback))
        boundary_special = true;
      return;
    }

    const std::size_t allocations =
        solve_allocation_testing::calls - current.start_allocations;
  #if defined(HIGHS_HIPDLP_CUDA_DETECTED) && defined(HIGHS_HIPDLP_TESTING)
    SolveCudaAudit cuda_delta;
    if (audit_cuda) {
      cuda_delta = solveCudaDelta(
        current.cuda_start, hipdlp::backend::cuda::testing::resourceAudit());
      for (const auto field : kSolveCudaFields)
      current.cuda_total.*field += cuda_delta.*field;
    }
  #endif
    ++current.ends;
    current.allocations += allocations;
    current.bytes += solve_allocation_testing::bytes - current.start_bytes;
    if (depth == 0 || stack[depth - 1] != event.phase)
      balanced = false;
    else
      --depth;
    if (event.phase == SolvePhase::kDownload && iterating) {
      if (candidate_downloads++ == 0) {
        first_candidate_download_allocations = allocations;
#if defined(HIGHS_HIPDLP_CUDA_DETECTED) && defined(HIGHS_HIPDLP_TESTING)
        first_candidate_cuda = cuda_delta;
#endif
      }
    }
    if (event.phase == SolvePhase::kBoundary) {
      if (iterating) {
        if (!event.continued)
          ++terminal_boundaries;
        else if (boundary_special)
          ++special_boundaries;
        else {
          ++ordinary_boundaries;
          ordinary_allocations += allocations;
#if defined(HIGHS_HIPDLP_CUDA_DETECTED) && defined(HIGHS_HIPDLP_TESTING)
          if (audit_cuda)
            ordinary_cuda_valid = ordinary_cuda_valid &&
                noCudaResourcesOrUploads(cuda_delta) &&
                cuda_delta.device_to_host_transfers == 1 &&
                cuda_delta.device_to_host_bytes == 112 &&
                cuda_delta.synchronization_calls == 1;
#endif
        }
      }
      boundary_active = false;
    }
    if (event.phase == SolvePhase::kIterationBlock &&
        event.iteration - current.start_iteration == 40)
      ++complete_blocks;
    if (event.phase == SolvePhase::kIterations) iterating = false;
    if (event.phase == SolvePhase::kSolve)
      solve_allocation_testing::enabled = false;
  }
};

}  // namespace

TEST_CASE("hi-pdlp-early-kkt-option-is-explicitly-opt-in",
    "[pdlp][hipdlp][solver][options][early-kkt-option][algorithm-study]") {
  Highs highs;
  highs.setOptionValue("output_flag", false);
  bool enabled = true;
  REQUIRE(highs.getOptionValue("hipdlp_early_kkt_check", enabled) ==
    HighsStatus::kOk);
  CHECK_FALSE(enabled);
    const hipdlp::api::OptionTranslationResult defaults =
        hipdlp::api::translateOptions(highs.getOptions(),
              hipdlp::api::defaultMigrationOptions());
    REQUIRE(defaults.ok());
    CHECK_FALSE(defaults.options.early_kkt_check);
  REQUIRE(highs.setOptionValue("hipdlp_early_kkt_check", true) ==
    HighsStatus::kOk);
  REQUIRE(highs.getOptionValue("hipdlp_early_kkt_check", enabled) ==
    HighsStatus::kOk);
  CHECK(enabled);
    const hipdlp::api::OptionTranslationResult opted_in =
        hipdlp::api::translateOptions(highs.getOptions(),
              hipdlp::api::defaultMigrationOptions());
    REQUIRE(opted_in.ok());
    CHECK(opted_in.options.early_kkt_check);
}

TEST_CASE("hi-pdlp-periodic-kkt-probe-rejects-offset-hidden-defect",
          "[pdlp][hipdlp][solver][periodic-kkt][algorithm-study]") {
  HighsLp lp = makeTinyFeasibleLp();
  lp.num_col_ = 3;
  lp.num_row_ = 2;
  lp.col_cost_ = {0.0, 0.0, 0.0};
  lp.col_lower_ = {10000.0, 0.0, 0.0};
  lp.col_upper_ = {10000.0, kHighsInf, kHighsInf};
  lp.row_lower_ = {0.0, 1.0};
  lp.row_upper_ = {0.0, 1.0};
  lp.offset_ = 1e6;
  lp.a_matrix_.num_col_ = 3;
  lp.a_matrix_.num_row_ = 2;
  lp.a_matrix_.start_ = {0, 0, 1, 2};
  lp.a_matrix_.index_ = {0, 1};
  lp.a_matrix_.value_ = {1.0, 0.03};
  class ProbeObserver : public hipdlp::algorithm::SolveObserver {
   public:
    std::size_t preparations = 0;
    void onEvent(const hipdlp::algorithm::SolveObservation& event) noexcept
        override {
      if (event.begin && event.iteration == 400 &&
          event.phase == hipdlp::algorithm::SolvePhase::kHostPreparation)
        ++preparations;
    }
  };
  std::vector<hipdlp::backend::Kind> backends = {hipdlp::backend::Kind::kCpu};
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
  if (cudaBackendUsable()) backends.push_back(hipdlp::backend::Kind::kCuda);
#endif
  for (const hipdlp::backend::Kind backend : backends) {
    for (const double tolerance : {1e-4, 1e-10}) {
      for (const bool enabled : {false, true}) {
      CAPTURE(static_cast<int>(backend), tolerance, enabled);
      hipdlp::api::SolverOptions options = makeCpuFixedOptions();
      options.backend = backend;
      options.early_kkt_check = enabled;
      options.iteration_limit = 401;
      options.kkt_tolerance = tolerance;
      options.time_limit = 2.0;
      ProbeObserver observer;
      hipdlp::algorithm::SolveContext context;
      context.observer = &observer;
      const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
          options, hipdlp::model::OriginalModel(lp), context);
      INFO(result.error.message);
      REQUIRE(result.ok());
      CHECK(observer.preparations == (enabled ? 1 : 0));
      CHECK(result.status == hipdlp::api::SolverStatus::kIterationLimit);
      CHECK(result.iteration_count == 401);
      REQUIRE(result.solution.col_value.size() == 3);
      const double internal_primal =
          std::abs(0.03 * result.solution.col_value[2] - 1.0) / 2.0;
      CHECK(internal_primal > options.optimality_tolerance);
      CHECK(internal_primal < 0.01);
      const hipdlp::verification::KktVerificationResult checked =
          hipdlp::verification::verifyKkt(lp, result.solution, tolerance);
      REQUIRE(checked.input_valid);
      REQUIRE(checked.metrics_finite);
      CHECK_FALSE(checked.pass);
      CHECK_FALSE(checked.objective_gap_ok);
      HighsLp unshifted = lp;
      unshifted.offset_ = 0.0;
      const auto unshifted_check = hipdlp::verification::verifyKkt(
          unshifted, result.solution, tolerance);
      CHECK(unshifted_check.pass == checked.pass);
      CHECK(unshifted_check.metrics.relative_objective_gap ==
        checked.metrics.relative_objective_gap);
      }
    }
  }
}

TEST_CASE("hi-pdlp-original-rejection-preserves-initial-limits",
          "[pdlp][hipdlp][solver][original-stopping]") {
  HighsLp lp = makeTinyFeasibleLp();
  lp.col_cost_ = {-1e9};
  lp.col_upper_ = {0.005};
  lp.row_lower_ = {0.0};
  lp.row_upper_ = {kHighsInf};
  lp.offset_ = 1e12;
  std::vector<hipdlp::backend::Kind> backends = {hipdlp::backend::Kind::kCpu};
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
  if (cudaBackendUsable()) backends.push_back(hipdlp::backend::Kind::kCuda);
#endif
  for (const auto backend : backends) {
    for (const HighsInt limit : {0, 1, 40}) {
      CAPTURE(static_cast<int>(backend), limit);
      auto options = makeCpuFixedOptions();
      options.backend = backend;
      options.iteration_limit = limit;
      const auto result = hipdlp::algorithm::solve(
          options, hipdlp::model::OriginalModel(lp));
      INFO(result.error.message);
      REQUIRE(result.ok());
      REQUIRE(result.has_solution);
      CHECK(result.iteration_count <= limit);
      const auto checked = hipdlp::verification::verifyKkt(lp, result.solution, 1e-4);
      REQUIRE(checked.input_valid);
      REQUIRE(checked.metrics_finite);
      if (limit == 0) {
        CHECK(result.status == hipdlp::api::SolverStatus::kIterationLimit);
        CHECK(result.iteration_count == 0);
        CHECK_FALSE(checked.pass);
      } else {
        CHECK(result.status == hipdlp::api::SolverStatus::kOptimal);
        CHECK(checked.pass);
        CHECK(checked.complementarity_ok);
        CHECK(result.solution.col_value == std::vector<double>{0.005});
      }
    }
    auto options = makeCpuFixedOptions();
    options.backend = backend;
    options.time_limit = 0.0;
    const auto result = hipdlp::algorithm::solve(
        options, hipdlp::model::OriginalModel(lp));
    REQUIRE(result.ok());
    CHECK(result.status == hipdlp::api::SolverStatus::kTimeLimit);
    CHECK(result.iteration_count == 0);
    options.iteration_limit = 0;
    const auto both_limits = hipdlp::algorithm::solve(
      options, hipdlp::model::OriginalModel(lp));
    REQUIRE(both_limits.ok());
    CHECK(both_limits.status == hipdlp::api::SolverStatus::kIterationLimit);
    CHECK(both_limits.iteration_count == 0);
    const auto already_optimal = hipdlp::algorithm::solve(
      options, hipdlp::model::OriginalModel(makeAlreadyOptimalLp()));
    REQUIRE(already_optimal.ok());
    CHECK(already_optimal.status == hipdlp::api::SolverStatus::kOptimal);
    CHECK(already_optimal.iteration_count == 0);
  }
}

TEST_CASE("hi-pdlp-solve-equivalent-lps-recover-original-kkt",
          "[pdlp][hipdlp][solver][algorithm-study]") {
  HighsLp original = makeTinyFeasibleLp();
  original.num_col_ = 3;
  original.num_row_ = 2;
  original.col_cost_ = {1.0, 2.0, 3.0};
  original.col_lower_ = {0.0, 0.0, 0.0};
  original.col_upper_ = {kHighsInf, kHighsInf, 5.0};
  original.row_lower_ = {3.0, 2.0};
  original.row_upper_ = {3.0, 4.0};
  original.offset_ = 7.0;
  original.a_matrix_.num_col_ = 3;
  original.a_matrix_.num_row_ = 2;
  original.a_matrix_.start_ = {0, 1, 3, 4};
  original.a_matrix_.index_ = {0, 0, 1, 1};
  original.a_matrix_.value_ = {1.0, 1.0, 1.0, 1.0};

  struct EquivalentCase {
    bool row_scaling;
    bool col_scaling;
    bool permutation;
    double objective_scale;
  };
  const EquivalentCase cases[] = {
      {false, false, false, 1.0}, {true, false, false, 1.0},
      {false, true, false, 1.0}, {false, false, true, 1.0},
      {false, false, false, 8.0}, {true, true, true, 0.125}};
  std::vector<hipdlp::backend::Kind> backends = {hipdlp::backend::Kind::kCpu};
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
  if (cudaBackendUsable()) backends.push_back(hipdlp::backend::Kind::kCuda);
#endif
  for (const hipdlp::backend::Kind backend : backends) {
    for (const EquivalentCase& test_case : cases) {
      CAPTURE(static_cast<int>(backend), test_case.row_scaling,
              test_case.col_scaling, test_case.permutation,
              test_case.objective_scale);
      const std::vector<double> row_scale = test_case.row_scaling
          ? std::vector<double>{0.01, 100.0} : std::vector<double>{1.0, 1.0};
      const std::vector<double> col_scale = test_case.col_scaling
          ? std::vector<double>{10.0, 0.1, 2.0}
          : std::vector<double>{1.0, 1.0, 1.0};
      const std::vector<HighsInt> row_order = test_case.permutation
          ? std::vector<HighsInt>{1, 0} : std::vector<HighsInt>{0, 1};
      const std::vector<HighsInt> col_order = test_case.permutation
          ? std::vector<HighsInt>{2, 0, 1} : std::vector<HighsInt>{0, 1, 2};
      HighsLp equivalent = original;
      equivalent.offset_ *= test_case.objective_scale;
      equivalent.a_matrix_.start_.clear();
      equivalent.a_matrix_.index_.clear();
      equivalent.a_matrix_.value_.clear();
      for (HighsInt row = 0; row < original.num_row_; ++row) {
        const HighsInt source_row = row_order[row];
        equivalent.row_lower_[row] =
            row_scale[source_row] * original.row_lower_[source_row];
        equivalent.row_upper_[row] =
            row_scale[source_row] * original.row_upper_[source_row];
      }
      for (HighsInt col = 0; col < original.num_col_; ++col) {
        const HighsInt source_col = col_order[col];
        equivalent.col_cost_[col] = test_case.objective_scale *
            col_scale[source_col] * original.col_cost_[source_col];
        equivalent.col_lower_[col] =
            original.col_lower_[source_col] / col_scale[source_col];
        equivalent.col_upper_[col] =
            original.col_upper_[source_col] / col_scale[source_col];
        equivalent.a_matrix_.start_.push_back(
            static_cast<HighsInt>(equivalent.a_matrix_.value_.size()));
        for (HighsInt row = 0; row < original.num_row_; ++row) {
          for (HighsInt entry = original.a_matrix_.start_[source_col];
               entry < original.a_matrix_.start_[source_col + 1]; ++entry) {
            if (original.a_matrix_.index_[entry] != row_order[row]) continue;
            equivalent.a_matrix_.index_.push_back(row);
            equivalent.a_matrix_.value_.push_back(
                row_scale[row_order[row]] * original.a_matrix_.value_[entry] *
                col_scale[source_col]);
          }
        }
      }
      equivalent.a_matrix_.start_.push_back(
          static_cast<HighsInt>(equivalent.a_matrix_.value_.size()));
      hipdlp::api::SolverOptions options = makeCpuFixedOptions();
      options.backend = backend;
      options.scaling.mode = hipdlp::model::ScalingMode::kLegacy;
      options.early_kkt_check = true;
      options.policy.step_size = hipdlp::api::StepSizePolicy::kPid;
      options.kkt_tolerance = options.optimality_tolerance = 1e-8;
      options.time_limit = 2.0;
      const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
          options, hipdlp::model::OriginalModel(equivalent));
      INFO(result.error.message);
      REQUIRE(result.ok());
      REQUIRE(result.status == hipdlp::api::SolverStatus::kOptimal);
      REQUIRE(hipdlp::verification::verifyKkt(equivalent, result.solution,
                                             kSolverTolerance).pass);
      HighsSolution recovered = result.solution;
      for (HighsInt col = 0; col < original.num_col_; ++col) {
        const HighsInt source_col = col_order[col];
        recovered.col_value[source_col] =
            col_scale[source_col] * result.solution.col_value[col];
        recovered.col_dual[source_col] = result.solution.col_dual[col] /
            (col_scale[source_col] * test_case.objective_scale);
      }
      for (HighsInt row = 0; row < original.num_row_; ++row) {
        const HighsInt source_row = row_order[row];
        recovered.row_value[source_row] =
            result.solution.row_value[row] / row_scale[source_row];
        recovered.row_dual[source_row] = row_scale[source_row] *
            result.solution.row_dual[row] / test_case.objective_scale;
      }
      const hipdlp::verification::KktVerificationResult checked =
          hipdlp::verification::verifyKkt(original, recovered, kSolverTolerance);
      REQUIRE(checked.input_valid);
      REQUIRE(checked.metrics_finite);
      CHECK(checked.pass);
      CHECK(recovered.col_value[0] == Approx(1.0).margin(kSolverTolerance));
      CHECK(recovered.col_value[1] == Approx(2.0).margin(kSolverTolerance));
      CHECK(recovered.col_value[2] == Approx(0.0).margin(kSolverTolerance));
      CHECK(checked.metrics.primal_objective ==
            Approx(12.0).margin(kSolverTolerance));
    }
  }
}

TEST_CASE("hi-pdlp-solve-preallocates-first-postsolve-work-buffer",
          "[pdlp][hipdlp][phase8][allocation-audit][cpu]") {
  hipdlp::api::SolverOptions options = makeCpuFixedOptions();
  options.iteration_limit = 1;
  for (const HighsInt dimension : {1, 257}) {
    INFO(dimension);
    const HighsLp lp = makeResidentDiagonalLp(dimension, false);
    HostPreparationObserver observer;
    hipdlp::algorithm::SolveContext context;
    context.observer = &observer;
    const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
        options, hipdlp::model::OriginalModel(lp), context);
    INFO(result.error.message);
    REQUIRE(result.ok());
    CHECK(result.iteration_count == 1);
    CHECK(observer.preparations == 1);
    if (solve_allocation_testing::available) {
      CHECK(observer.allocations == 0);
      CHECK(observer.allocated_bytes == 0);
    }
    CHECK_FALSE(observer.iterating);
    CHECK_FALSE(solve_allocation_testing::enabled);
  }
}

TEST_CASE("hi-pdlp-solve-observes-production-phases-without-step-allocation",
          "[pdlp][hipdlp][phase8][allocation-audit][cpu]") {
  using hipdlp::algorithm::SolvePhase;
  hipdlp::api::SolverOptions options = makeCpuFixedOptions();
  options.iteration_limit = 400;
  options.optimality_tolerance = 1e-12;
  options.kkt_tolerance = 1e-12;
  SECTION("fixed") {}
  SECTION("PID") {
    options.policy.step_size = hipdlp::api::StepSizePolicy::kPid;
  }
  for (const HighsInt dimension : {1, 257, 4097}) {
    INFO(dimension);
    const HighsLp lp = makeResidentDiagonalLp(dimension, dimension > 1);
    ProductionSolveObserver observer;
    hipdlp::algorithm::SolveDiagnostics diagnostics;
    hipdlp::algorithm::SolveContext context;
    context.observer = &observer;
    context.diagnostics = &diagnostics;
    const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
        options, hipdlp::model::OriginalModel(lp), context);
    solve_allocation_testing::enabled = false;
    INFO(result.error.message);
    REQUIRE(result.ok());
    CHECK(observer.balanced);
    CHECK(observer.depth == 0);
    CHECK(observer.phase(SolvePhase::kSolve).ends == 1);
    CHECK(observer.phase(SolvePhase::kPreprocessing).ends == 1);
    CHECK(observer.phase(SolvePhase::kBackendSetup).ends >= 1);
    CHECK(observer.phase(SolvePhase::kOperatorNorm).ends == 1);
    CHECK(observer.phase(SolvePhase::kIterations).ends == 1);
    CHECK(observer.phase(SolvePhase::kStep).ends ==
          static_cast<std::size_t>(result.iteration_count));
    if (solve_allocation_testing::available) {
      CHECK(observer.phase(SolvePhase::kStep).allocations == 0);
      CHECK(observer.phase(SolvePhase::kDownload).allocations > 0);
      CHECK(observer.first_candidate_download_allocations == 0);
    }
    CHECK(observer.candidate_downloads == diagnostics.candidate_downloads);
    CHECK(observer.phase(SolvePhase::kBoundary).ends ==
          diagnostics.check_boundaries);
    CHECK(observer.phase(SolvePhase::kReferenceFixedPoint).ends ==
          diagnostics.restart_reference_checks);
    CHECK(observer.complete_blocks ==
          static_cast<std::size_t>(result.iteration_count / 40));
    CHECK(observer.complete_blocks >= 3);
    CHECK(diagnostics.restarts > 0);
    CHECK(observer.terminal_boundaries == 1);
    CHECK(observer.ordinary_boundaries == 0);
    CHECK(observer.special_boundaries + observer.terminal_boundaries + 1 ==
          diagnostics.check_boundaries);
    if (solve_allocation_testing::available)
      CHECK(observer.ordinary_allocations == 0);
    CHECK_FALSE(observer.iterating);
  }
}

TEST_CASE("hi-pdlp-solve-observer-balances-early-returns",
          "[pdlp][hipdlp][phase8][allocation-audit][cpu]") {
  using hipdlp::algorithm::SolvePhase;
  hipdlp::api::SolverOptions options = makeCpuFixedOptions();
  HighsLp lp = makeTinyFeasibleLp();
  bool expected_success = true;
  SECTION("invalid options") {
    options.iteration_limit = -1;
    expected_success = false;
  }
  SECTION("invalid model") {
    lp.a_matrix_.value_[0] = std::numeric_limits<double>::quiet_NaN();
    expected_success = false;
  }
  SECTION("zero iteration budget") { options.iteration_limit = 0; }
  SECTION("zero time budget") { options.time_limit = 0.0; }
  SECTION("optimal initial point") { lp = makeAlreadyOptimalLp(); }
  ProductionSolveObserver observer;
  hipdlp::algorithm::SolveContext context;
  context.observer = &observer;
  const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
      options, hipdlp::model::OriginalModel(lp), context);
  solve_allocation_testing::enabled = false;
  CHECK(result.ok() == expected_success);
  CHECK(observer.balanced);
  CHECK(observer.depth == 0);
  CHECK(observer.phase(SolvePhase::kSolve).begins == 1);
  CHECK(observer.phase(SolvePhase::kSolve).ends == 1);
  CHECK(observer.phase(SolvePhase::kIterations).begins == 0);
  CHECK(observer.phase(SolvePhase::kStep).begins == 0);
}

TEST_CASE("hi-pdlp-solve-timings-are-opt-in-and-preserve-results",
          "[pdlp][hipdlp][phase8][allocation-audit][cpu]") {
  const HighsLp lp = makeResidentDiagonalLp(257, true);
  hipdlp::api::SolverOptions options = makeCpuFixedOptions();
  options.iteration_limit = 400;
  options.optimality_tolerance = 1e-12;
  options.kkt_tolerance = 1e-12;
  options.policy.step_size = hipdlp::api::StepSizePolicy::kPid;
  const hipdlp::api::SolveResult plain = hipdlp::algorithm::solve(
      options, hipdlp::model::OriginalModel(lp));
  REQUIRE(plain.ok());

  ProductionSolveObserver observer;
  hipdlp::algorithm::SolveContext context;
  CHECK(context.diagnostics == nullptr);
  CHECK(context.observer == nullptr);
  context.observer = &observer;
  const hipdlp::api::SolveResult observed = hipdlp::algorithm::solve(
      options, hipdlp::model::OriginalModel(lp), context);
  solve_allocation_testing::enabled = false;
  REQUIRE(observed.ok());

  hipdlp::algorithm::SolveDiagnostics diagnostics;
  context.diagnostics = &diagnostics;
  context.observer = nullptr;
  const hipdlp::api::SolveResult timed = hipdlp::algorithm::solve(
      options, hipdlp::model::OriginalModel(lp), context);
  REQUIRE(timed.ok());
  for (const hipdlp::api::SolveResult* result : {&observed, &timed}) {
    CHECK(result->status == plain.status);
    CHECK(result->iteration_count == plain.iteration_count);
    CHECK(result->solution.col_value == plain.solution.col_value);
    CHECK(result->solution.row_value == plain.solution.row_value);
    CHECK(result->solution.col_dual == plain.solution.col_dual);
    CHECK(result->solution.row_dual == plain.solution.row_dual);
  }
  CHECK(diagnostics.total_wall_seconds > 0.0);
  CHECK(diagnostics.preprocessing_seconds > 0.0);
  CHECK(diagnostics.backend_setup_seconds > 0.0);
  CHECK(diagnostics.operator_norm_seconds > 0.0);
  CHECK(diagnostics.iterations_seconds > 0.0);
  CHECK(diagnostics.iteration_block_seconds > 0.0);
  CHECK(diagnostics.step_seconds > 0.0);
  CHECK(diagnostics.boundary_seconds > 0.0);
  CHECK(diagnostics.download_seconds > 0.0);
  CHECK(diagnostics.host_preparation_seconds > 0.0);
  CHECK(diagnostics.host_verification_seconds > 0.0);
  CHECK(diagnostics.reference_fixed_point_seconds > 0.0);
  CHECK(diagnostics.fixed_point_fallback_seconds == 0.0);
  CHECK(diagnostics.iteration_block_seconds <= diagnostics.iterations_seconds);
  CHECK(diagnostics.step_seconds <= diagnostics.iterations_seconds);
  CHECK(diagnostics.total_wall_seconds >=
        diagnostics.preprocessing_seconds + diagnostics.backend_setup_seconds +
            diagnostics.operator_norm_seconds + diagnostics.iterations_seconds);
}

#if defined(HIGHS_HIPDLP_CUDA_DETECTED) && defined(HIGHS_HIPDLP_TESTING)
TEST_CASE("hi-pdlp-solve-prewarms-first-cuda-candidate-download",
          "[pdlp][hipdlp][phase8][allocation-audit][cuda]") {
  REQUIRE(cudaBackendUsable());
  hipdlp::api::SolverOptions options = makeCpuFixedOptions();
  options.backend = hipdlp::backend::Kind::kCuda;
  options.iteration_limit = 1;
  for (const HighsInt dimension : {1, 257, 4097}) {
    INFO(dimension);
    const HighsLp lp = makeResidentDiagonalLp(dimension, false);
    ProductionSolveObserver observer;
    observer.audit_cuda = true;
    hipdlp::algorithm::SolveContext context;
    context.observer = &observer;
    const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
        options, hipdlp::model::OriginalModel(lp), context);
    solve_allocation_testing::enabled = false;
    INFO(result.error.message);
    REQUIRE(result.ok());
    CHECK(result.iteration_count == 1);
    CHECK(observer.balanced);
    CHECK(observer.candidate_downloads == 1);
    if (solve_allocation_testing::available)
      CHECK(observer.first_candidate_download_allocations == 0);
    CHECK(noCudaResourcesOrUploads(observer.first_candidate_cuda));
    CHECK(observer.first_candidate_cuda.device_to_host_transfers == 5);
    CHECK(observer.first_candidate_cuda.device_to_host_bytes ==
          32 * static_cast<std::size_t>(dimension) + 4);
    CHECK(observer.first_candidate_cuda.synchronization_calls == 1);
    CHECK(observer.ordinary_boundaries == 0);
    CHECK(observer.terminal_boundaries == 1);
  }
}

TEST_CASE("hi-pdlp-solve-audits-each-ordinary-cuda-boundary",
          "[pdlp][hipdlp][phase8][allocation-audit][cuda]") {
  using hipdlp::algorithm::SolvePhase;
  REQUIRE(cudaBackendUsable());
  hipdlp::api::SolverOptions options = makeCpuFixedOptions();
  options.backend = hipdlp::backend::Kind::kCuda;
  options.iteration_limit = 400;
  options.optimality_tolerance = 1e-12;
  options.kkt_tolerance = 1e-12;
  SECTION("fixed") {}
  SECTION("PID") {
    options.policy.step_size = hipdlp::api::StepSizePolicy::kPid;
  }
  for (const HighsInt dimension : {1, 257, 4097}) {
    INFO(dimension);
    const HighsLp lp = makeResidentDiagonalLp(dimension, false);
    ProductionSolveObserver observer;
    observer.audit_cuda = true;
    hipdlp::algorithm::SolveDiagnostics diagnostics;
    hipdlp::algorithm::SolveContext context;
    context.diagnostics = &diagnostics;
    context.observer = &observer;
    const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
        options, hipdlp::model::OriginalModel(lp), context);
    solve_allocation_testing::enabled = false;
    INFO(result.error.message);
    REQUIRE(result.ok());
    CHECK(observer.balanced);
    CHECK(observer.depth == 0);
    REQUIRE(observer.ordinary_boundaries > 0);
    CHECK(observer.ordinary_cuda_valid);
    if (solve_allocation_testing::available)
      CHECK(observer.ordinary_allocations == 0);
    CHECK(observer.terminal_boundaries == 1);
    CHECK(observer.special_boundaries + observer.ordinary_boundaries +
              observer.terminal_boundaries + 1 == diagnostics.check_boundaries);
    const ObservedPhase& steps = observer.phase(SolvePhase::kStep);
          if (solve_allocation_testing::available) CHECK(steps.allocations == 0);
        CHECK(steps.cuda_total.allocation_calls == 0);
        CHECK(steps.cuda_total.free_calls == 0);
        CHECK(steps.cuda_total.pinned_allocation_calls == 0);
        CHECK(steps.cuda_total.pinned_free_calls == 0);
        CHECK(steps.cuda_total.host_to_device_transfers == diagnostics.command_buffer_replays);
        CHECK(steps.cuda_total.host_to_device_bytes ==
          diagnostics.command_buffer_replays * 117 * sizeof(double));
    CHECK(steps.cuda_total.device_to_host_transfers == 0);
    CHECK(steps.cuda_total.synchronization_calls == 0);
    CHECK(diagnostics.restarts > 0);
    CHECK(diagnostics.primal_weight_updates ==
          (options.policy.step_size == hipdlp::api::StepSizePolicy::kPid
               ? diagnostics.restarts : 0));
    CHECK(observer.phase(SolvePhase::kReferenceFixedPoint).ends ==
          diagnostics.restart_reference_checks);
    if (solve_allocation_testing::available)
      CHECK(observer.first_candidate_download_allocations == 0);
    CHECK(noCudaResourcesOrUploads(observer.first_candidate_cuda));

    const hipdlp::api::SolveResult plain = hipdlp::algorithm::solve(
        options, hipdlp::model::OriginalModel(lp));
    REQUIRE(plain.ok());
    CHECK(result.status == plain.status);
    CHECK(result.iteration_count == plain.iteration_count);
    CHECK(result.solution.col_value == plain.solution.col_value);
    CHECK(result.solution.row_value == plain.solution.row_value);
    CHECK(result.solution.col_dual == plain.solution.col_dual);
    CHECK(result.solution.row_dual == plain.solution.row_dual);
  }
}
#endif

TEST_CASE("hi-pdlp-resident-check-production-fixed-and-pid",
          "[pdlp][hipdlp][phase8][solver][cpu]") {
  hipdlp::api::SolverOptions options = makeCpuFixedOptions();
  options.optimality_tolerance = 1e-12;
  options.kkt_tolerance = 1e-12;
  options.iteration_limit = 200000;
  // Production scaling: like legacy, unscaled PID can drive the primal weight
  // to extremes once the primal iterate is pinned at its bounds.
  options.scaling.mode = hipdlp::model::ScalingMode::kAuto;
  SECTION("fixed") {}
  SECTION("PID") {
    options.policy.step_size = hipdlp::api::StepSizePolicy::kPid;
  }

  for (const HighsInt dimension : {1, 2, 33}) {
    INFO(dimension);
    const HighsLp lp = makeResidentDiagonalLp(dimension, dimension > 1);
    hipdlp::algorithm::SolveDiagnostics diagnostics;
    hipdlp::algorithm::SolveContext context;
    context.diagnostics = &diagnostics;
    const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
        options, hipdlp::model::OriginalModel(lp), context);
    INFO(result.error.message);
    REQUIRE(result.ok());
    REQUIRE(result.status == hipdlp::api::SolverStatus::kOptimal);
    REQUIRE(result.has_solution);
    CHECK(result.iteration_count % 40 == 0);
    CHECK(diagnostics.ordinary_check_boundaries >= 3);
    CHECK(diagnostics.resident_check_calls ==
          diagnostics.ordinary_check_boundaries);
    CHECK(diagnostics.check_boundaries ==
          diagnostics.resident_check_calls + 1);
    CHECK(diagnostics.ordinary_check_d2h_bytes == 0);
    CHECK(diagnostics.ordinary_check_transfers == 0);
    CHECK(diagnostics.ordinary_check_synchronizations == 0);
    CHECK(diagnostics.candidate_downloads == diagnostics.resident_check_calls);
    CHECK(diagnostics.host_oracle_checks == 1);
    CHECK(diagnostics.certificate_checks == diagnostics.resident_check_calls);
    CHECK(diagnostics.restarts > 0);
    CHECK(diagnostics.restart_reference_checks == diagnostics.restarts);
    CHECK(diagnostics.fixed_point_evaluations + 1 ==
          diagnostics.ordinary_check_boundaries +
              diagnostics.restart_reference_checks);
    CHECK(diagnostics.resource_creation_calls_before_iterations == 18);
    CHECK(diagnostics.resource_creation_calls_during_iterations == 0);
    CHECK(diagnostics.primal_weight_updates ==
          (options.policy.step_size == hipdlp::api::StepSizePolicy::kPid
               ? diagnostics.restarts
               : 0));
    const hipdlp::verification::KktVerificationResult verification =
        hipdlp::verification::verifyKkt(lp, result.solution,
                                      options.kkt_tolerance);
    INFO(verification.error);
    REQUIRE(verification.pass);
    REQUIRE(result.solution.col_value.size() ==
            static_cast<std::size_t>(dimension));
    for (const double value : result.solution.col_value)
      CHECK(value == Approx(1.0).margin(1e-8));
  }
}

TEST_CASE("hi-pdlp-resident-check-kkt-rejection-continues-to-next-candidate",
          "[pdlp][hipdlp][phase8][solver][cpu][kkt-gate]") {
  const HighsLp lp = makeTinyFeasibleLp();
  hipdlp::api::SolverOptions options = makeCpuFixedOptions();
  options.iteration_limit = 41;
  options.optimality_tolerance = 0.4;
  options.kkt_tolerance = 1.0;
  hipdlp::algorithm::SolveDiagnostics accepted_diagnostics;
  hipdlp::algorithm::SolveContext context;
  context.diagnostics = &accepted_diagnostics;
  const hipdlp::api::SolveResult accepted = hipdlp::algorithm::solve(
      options, hipdlp::model::OriginalModel(lp), context);
  REQUIRE(accepted.ok());
  REQUIRE(accepted.status == hipdlp::api::SolverStatus::kOptimal);
  REQUIRE(accepted.iteration_count == 40);
  CHECK(accepted_diagnostics.resident_check_calls == 1);
  CHECK(accepted_diagnostics.fixed_point_evaluations == 0);

  options.kkt_tolerance = 1e-12;
  hipdlp::algorithm::SolveDiagnostics rejected_diagnostics;
  context.diagnostics = &rejected_diagnostics;
  const hipdlp::api::SolveResult rejected = hipdlp::algorithm::solve(
      options, hipdlp::model::OriginalModel(lp), context);
  INFO(rejected.error.message);
  REQUIRE(rejected.ok());
  CHECK(rejected.status == hipdlp::api::SolverStatus::kIterationLimit);
  CHECK(rejected.iteration_count == 41);
  CHECK(rejected.solution.col_value != accepted.solution.col_value);
  CHECK(rejected.solution.row_dual != accepted.solution.row_dual);
  CHECK(rejected_diagnostics.resident_check_calls == 2);
  CHECK(rejected_diagnostics.ordinary_check_boundaries == 1);
  CHECK(rejected_diagnostics.candidate_downloads == 2);
  CHECK(rejected_diagnostics.certificate_checks == 3);
  CHECK(rejected_diagnostics.fixed_point_evaluations == 1);
  CHECK(rejected_diagnostics.restart_reference_checks == 0);
  const hipdlp::verification::KktVerificationResult verification =
      hipdlp::verification::verifyKkt(lp, rejected.solution,
                                    options.kkt_tolerance);
  REQUIRE(verification.input_valid);
  REQUIRE(verification.metrics_finite);
  CHECK_FALSE(verification.pass);
}

TEST_CASE("hi-pdlp-rejected-kkt-verifications-back-off",
          "[pdlp][hipdlp][solver][cpu][kkt-gate]") {
  // The loose transformed tolerance passes at every boundary while the
  // original one is out of reach, so every due verification is rejected. The
  // CPU verifies every boundary; a device backend spaces the checks 1, 2, 4,
  // 8 and then 16 boundaries apart. The limit boundary is always checked.
  const HighsLp lp = makeResidentDiagonalLp(17, true);
  std::vector<hipdlp::backend::Kind> backends = {hipdlp::backend::Kind::kCpu};
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
  if (cudaBackendUsable()) backends.push_back(hipdlp::backend::Kind::kCuda);
#endif
  for (const hipdlp::backend::Kind backend : backends) {
    CAPTURE(static_cast<int>(backend));
    hipdlp::api::SolverOptions options = makeCpuFixedOptions();
    options.backend = backend;
    options.optimality_tolerance = 0.4;
    options.kkt_tolerance = 1e-12;
    options.iteration_limit = 1600;
    hipdlp::algorithm::SolveDiagnostics diagnostics;
    hipdlp::algorithm::SolveContext context;
    context.diagnostics = &diagnostics;
    const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
        options, hipdlp::model::OriginalModel(lp), context);
    INFO(result.error.message);
    REQUIRE(result.ok());
    CHECK(result.status == hipdlp::api::SolverStatus::kIterationLimit);
    CHECK(result.iteration_count == 1600);
    // Every boundary on the CPU; on a device, iterations 40, 80, 160, 320,
    // 640, 1280 and the limit at 1600.
    const std::size_t expected =
        backend == hipdlp::backend::Kind::kCpu ? 40 : 7;
    CHECK(diagnostics.original_kkt_verifications == expected);
    CHECK(diagnostics.original_kkt_rejections == expected);
  }
}

TEST_CASE("hi-pdlp-resident-check-limit-does-not-finalize-unused-fpe",
          "[pdlp][hipdlp][phase8][solver][cpu][iteration-limit]") {
  const HighsLp lp = makeResidentDiagonalLp(2, true);
  hipdlp::api::SolverOptions options = makeCpuFixedOptions();
  options.iteration_limit = 40;
  options.optimality_tolerance = 1e-12;
  hipdlp::algorithm::SolveDiagnostics diagnostics;
  hipdlp::algorithm::SolveContext context;
  context.diagnostics = &diagnostics;
  const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
      options, hipdlp::model::OriginalModel(lp), context);
  INFO(result.error.message);
  REQUIRE(result.ok());
  CHECK(result.status == hipdlp::api::SolverStatus::kIterationLimit);
  CHECK(result.iteration_count == 40);
  CHECK(diagnostics.resident_check_calls == 1);
  CHECK(diagnostics.ordinary_check_boundaries == 1);
  CHECK(diagnostics.fixed_point_evaluations == 0);
  CHECK(diagnostics.restart_reference_checks == 0);
  CHECK(diagnostics.candidate_downloads == 1);
  CHECK(diagnostics.certificate_checks == 2);
}

namespace {

void checkResidentTerminalOverflow(const hipdlp::backend::Kind kind) {
  HighsLp lp = makeBoxedColumnLp();
  lp.row_lower_[0] = 1e160;
  lp.row_upper_[0] = 1e160;
  lp.col_upper_[0] = 5e160;
  hipdlp::api::SolverOptions options = makeCpuFixedOptions();
  options.backend = kind;
  options.iteration_limit = 40;
  options.optimality_tolerance = 1e-12;
  hipdlp::api::SolverStatus expected_status =
      hipdlp::api::SolverStatus::kIterationLimit;
  bool expect_error = false;
  SECTION("iteration limit precedes unused norm-square overflow") {}
  SECTION("independent optimality precedes unused norm-square overflow") {
    options.iteration_limit = 41;
    options.optimality_tolerance = 0.4;
    options.kkt_tolerance = 1.0;
    expected_status = hipdlp::api::SolverStatus::kOptimal;
  }
  SECTION("a continuing boundary still fails on norm-square overflow") {
    options.iteration_limit = 41;
    expect_error = true;
  }
  hipdlp::algorithm::SolveDiagnostics diagnostics;
  hipdlp::algorithm::SolveContext context;
  context.diagnostics = &diagnostics;
  const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
      options, hipdlp::model::OriginalModel(lp), context);
  INFO(result.error.message);
  CHECK(diagnostics.resident_check_calls == 1);
  CHECK(diagnostics.ordinary_check_boundaries == 1);
  CHECK(diagnostics.restart_reference_checks == 0);
  CHECK(diagnostics.resource_creation_calls_during_iterations == 0);
  if (kind == hipdlp::backend::Kind::kCuda) {
    CHECK(diagnostics.ordinary_check_transfers == 1);
    CHECK(diagnostics.ordinary_check_d2h_bytes == 112);
    CHECK(diagnostics.ordinary_check_synchronizations == 1);
  }
  if (expect_error) {
    CHECK_FALSE(result.ok());
    CHECK(result.status == hipdlp::api::SolverStatus::kError);
    CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidArgument);
    CHECK(result.error.message.find("norm square") != std::string::npos);
    CHECK_FALSE(result.has_solution);
    CHECK_FALSE(result.solution.value_valid);
    CHECK_FALSE(result.solution.dual_valid);
    CHECK(result.solution.col_value.empty());
    CHECK(result.solution.col_dual.empty());
    CHECK(result.solution.row_value.empty());
    CHECK(result.solution.row_dual.empty());
    CHECK(diagnostics.fixed_point_evaluations == 1);
  } else {
    REQUIRE(result.ok());
    CHECK(result.status == expected_status);
    CHECK(result.iteration_count == 40);
    REQUIRE(result.has_solution);
    CHECK(diagnostics.fixed_point_evaluations == 0);
    const hipdlp::verification::KktVerificationResult verification =
        hipdlp::verification::verifyKkt(lp, result.solution,
                                      options.kkt_tolerance);
    REQUIRE(verification.input_valid);
    REQUIRE(verification.metrics_finite);
    if (expected_status == hipdlp::api::SolverStatus::kOptimal)
      CHECK(verification.pass);
  }
}

}

TEST_CASE("hi-pdlp-resident-check-terminal-precedence-over-fpe-overflow",
          "[pdlp][hipdlp][phase8][solver][cpu][overflow][optional-error]") {
  checkResidentTerminalOverflow(hipdlp::backend::Kind::kCpu);
}

#ifdef HIGHS_HIPDLP_CUDA_DETECTED
TEST_CASE("hi-pdlp-resident-check-cuda-terminal-precedence-over-fpe-overflow",
          "[pdlp][hipdlp][phase8][solver][cuda][overflow][optional-error]") {
  checkResidentTerminalOverflow(hipdlp::backend::Kind::kCuda);
}

TEST_CASE("hi-pdlp-resident-check-production-cuda-host-oracle-fallback",
          "[pdlp][hipdlp][phase8][solver][cuda][host-oracle]") {
  REQUIRE(cudaBackendUsable());
  HighsLp lp = makeResidentDiagonalLp(1, false);
  // sum|a_i| above DBL_MAX / 4 trips the overflow guard on every check.
  lp.offset_ = 5e307;
  hipdlp::api::SolverOptions options = makeCpuFixedOptions();
  options.optimality_tolerance = 1e-12;
  options.kkt_tolerance = 1e-12;
  options.policy.step_size = hipdlp::api::StepSizePolicy::kPid;
  hipdlp::algorithm::SolveDiagnostics cpu_diagnostics;
  hipdlp::algorithm::SolveContext context;
  context.diagnostics = &cpu_diagnostics;
  const hipdlp::api::SolveResult cpu = hipdlp::algorithm::solve(
      options, hipdlp::model::OriginalModel(lp), context);
  INFO(cpu.error.message);
  REQUIRE(cpu.ok());
  REQUIRE(cpu.status == hipdlp::api::SolverStatus::kOptimal);
  hipdlp::algorithm::SolveDiagnostics gpu_diagnostics;
  context.diagnostics = &gpu_diagnostics;
  options.backend = hipdlp::backend::Kind::kCuda;
  const hipdlp::api::SolveResult gpu = hipdlp::algorithm::solve(
      options, hipdlp::model::OriginalModel(lp), context);
  INFO(gpu.error.message);
  REQUIRE(gpu.ok());
  CHECK(gpu.status == cpu.status);
  CHECK(gpu.iteration_count == cpu.iteration_count);
  CHECK(gpu_diagnostics.ordinary_check_boundaries >= 3);
  CHECK(gpu_diagnostics.host_oracle_checks ==
        gpu_diagnostics.resident_check_calls + 1);
  CHECK(gpu_diagnostics.candidate_downloads ==
        gpu_diagnostics.resident_check_calls);
  CHECK(gpu_diagnostics.fixed_point_evaluations + 1 ==
        gpu_diagnostics.ordinary_check_boundaries +
            gpu_diagnostics.restart_reference_checks);
  CHECK(gpu_diagnostics.restarts == cpu_diagnostics.restarts);
  CHECK(gpu_diagnostics.primal_weight_updates ==
        cpu_diagnostics.primal_weight_updates);
  CHECK(gpu_diagnostics.ordinary_check_transfers ==
        gpu_diagnostics.ordinary_check_boundaries);
  CHECK(gpu_diagnostics.ordinary_check_d2h_bytes ==
        112 * gpu_diagnostics.ordinary_check_boundaries);
  CHECK(gpu_diagnostics.ordinary_check_synchronizations ==
        gpu_diagnostics.ordinary_check_boundaries);
  const hipdlp::verification::KktVerificationResult verification =
      hipdlp::verification::verifyKkt(lp, gpu.solution, options.kkt_tolerance);
  INFO(verification.error);
  REQUIRE(verification.pass);
}

TEST_CASE("hi-pdlp-resident-check-production-cuda-absorbed-offset-stays-on-device",
          "[pdlp][hipdlp][phase8][solver][cuda][absorbed-offset]") {
  REQUIRE(cudaBackendUsable());
  HighsLp lp = makeResidentDiagonalLp(1, false);
  lp.offset_ = 1e150;
  hipdlp::api::SolverOptions options = makeCpuFixedOptions();
  options.optimality_tolerance = 1e-12;
  options.kkt_tolerance = 1e-12;
  options.policy.step_size = hipdlp::api::StepSizePolicy::kPid;
  hipdlp::algorithm::SolveDiagnostics cpu_diagnostics;
  hipdlp::algorithm::SolveContext context;
  context.diagnostics = &cpu_diagnostics;
  const hipdlp::api::SolveResult cpu = hipdlp::algorithm::solve(
      options, hipdlp::model::OriginalModel(lp), context);
  INFO(cpu.error.message);
  REQUIRE(cpu.ok());
  REQUIRE(cpu.status == hipdlp::api::SolverStatus::kOptimal);
  hipdlp::algorithm::SolveDiagnostics gpu_diagnostics;
  context.diagnostics = &gpu_diagnostics;
  options.backend = hipdlp::backend::Kind::kCuda;
  const hipdlp::api::SolveResult gpu = hipdlp::algorithm::solve(
      options, hipdlp::model::OriginalModel(lp), context);
  INFO(gpu.error.message);
  REQUIRE(gpu.ok());
  CHECK(gpu.status == cpu.status);
  CHECK(gpu.iteration_count == cpu.iteration_count);
  CHECK(gpu_diagnostics.ordinary_check_boundaries >= 3);
  // Every small term is absorbed by the 1e150 offset, but the interval in
  // guardCheckStopping accepts the device packet, so only the terminal
  // boundary ever reaches the host oracle.
  CHECK(gpu_diagnostics.host_oracle_checks == 1);
  CHECK(gpu_diagnostics.candidate_downloads == 1);
  CHECK(gpu_diagnostics.restarts == cpu_diagnostics.restarts);
  CHECK(gpu_diagnostics.primal_weight_updates ==
        cpu_diagnostics.primal_weight_updates);
  CHECK(gpu_diagnostics.ordinary_check_transfers ==
        gpu_diagnostics.ordinary_check_boundaries);
  CHECK(gpu_diagnostics.ordinary_check_d2h_bytes ==
        112 * gpu_diagnostics.ordinary_check_boundaries);
  CHECK(gpu_diagnostics.ordinary_check_synchronizations ==
        gpu_diagnostics.ordinary_check_boundaries);
  const hipdlp::verification::KktVerificationResult verification =
      hipdlp::verification::verifyKkt(lp, gpu.solution, options.kkt_tolerance);
  INFO(verification.error);
  REQUIRE(verification.pass);
}

TEST_CASE("hi-pdlp-resident-check-production-cuda-scalar-transport",
          "[pdlp][hipdlp][phase8][solver][cuda]") {
  REQUIRE(cudaBackendUsable());
  hipdlp::api::SolverOptions options = makeCpuFixedOptions();
  options.optimality_tolerance = 1e-12;
  options.kkt_tolerance = 1e-12;
  SECTION("fixed") {}
  SECTION("PID") {
    options.policy.step_size = hipdlp::api::StepSizePolicy::kPid;
  }
  for (const HighsInt dimension : {1, 257}) {
    INFO(dimension);
    const HighsLp lp = makeResidentDiagonalLp(dimension, false);
    hipdlp::algorithm::SolveDiagnostics cpu_diagnostics;
    hipdlp::algorithm::SolveContext context;
    context.diagnostics = &cpu_diagnostics;
    options.backend = hipdlp::backend::Kind::kCpu;
    const hipdlp::api::SolveResult cpu = hipdlp::algorithm::solve(
        options, hipdlp::model::OriginalModel(lp), context);
    INFO(cpu.error.message);
    REQUIRE(cpu.ok());
    REQUIRE(cpu.status == hipdlp::api::SolverStatus::kOptimal);

    hipdlp::algorithm::SolveDiagnostics gpu_diagnostics;
    context.diagnostics = &gpu_diagnostics;
    options.backend = hipdlp::backend::Kind::kCuda;
    const hipdlp::api::SolveResult gpu = hipdlp::algorithm::solve(
        options, hipdlp::model::OriginalModel(lp), context);
    INFO(gpu.error.message);
    REQUIRE(gpu.ok());
    REQUIRE(gpu.status == cpu.status);
    CHECK(gpu.iteration_count == cpu.iteration_count);
    CHECK(gpu_diagnostics.restarts == cpu_diagnostics.restarts);
    CHECK(gpu_diagnostics.primal_weight_updates ==
          cpu_diagnostics.primal_weight_updates);
    CHECK(gpu_diagnostics.ordinary_check_boundaries >= 3);
    CHECK(gpu_diagnostics.ordinary_check_transfers ==
          gpu_diagnostics.ordinary_check_boundaries);
    CHECK(gpu_diagnostics.ordinary_check_d2h_bytes ==
          112 * gpu_diagnostics.ordinary_check_boundaries);
    CHECK(gpu_diagnostics.ordinary_check_synchronizations ==
          gpu_diagnostics.ordinary_check_boundaries);
    CHECK(gpu_diagnostics.host_oracle_checks == 1);
    CHECK(gpu_diagnostics.candidate_downloads == 1);
    CHECK(gpu_diagnostics.certificate_checks == 1);
    CHECK(gpu_diagnostics.resource_creation_calls_during_iterations == 0);
    const hipdlp::verification::KktVerificationResult verification =
        hipdlp::verification::verifyKkt(lp, gpu.solution,
                                      options.kkt_tolerance);
    INFO(verification.error);
    REQUIRE(verification.pass);
    REQUIRE(gpu.solution.col_value.size() == cpu.solution.col_value.size());
    for (std::size_t coordinate = 0;
         coordinate < cpu.solution.col_value.size(); ++coordinate) {
      CHECK(gpu.solution.col_value[coordinate] ==
            Approx(cpu.solution.col_value[coordinate]).margin(kCudaTolerance));
      CHECK(gpu.solution.row_dual[coordinate] ==
            Approx(cpu.solution.row_dual[coordinate]).margin(kCudaTolerance));
    }
  }
}
#endif

TEST_CASE("hi-pdlp-rewrite-solver-solves-tiny-feasible-lp",
          "[pdlp][hipdlp][phase7][solver][optimal]") {
  const HighsLp lp = makeTinyFeasibleLp();
  const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
      makeCpuFixedOptions(), hipdlp::model::OriginalModel(lp));

  INFO(result.error.message);
  REQUIRE(result.ok());
  CHECK(result.status == hipdlp::api::SolverStatus::kOptimal);
  REQUIRE(result.has_solution);
  REQUIRE(result.solution.value_valid);
  REQUIRE(result.solution.dual_valid);
  REQUIRE(result.solution.col_value.size() == 1);
  REQUIRE(result.solution.row_value.size() == 1);
  REQUIRE(result.solution.col_dual.size() == 1);
  REQUIRE(result.solution.row_dual.size() == 1);
  CHECK(result.solution.col_value[0] == Approx(1.0).margin(kSolverTolerance));
  CHECK(result.solution.row_value[0] == Approx(1.0).margin(kSolverTolerance));

  const hipdlp::verification::KktVerificationResult verification =
      hipdlp::verification::verifyKkt(lp, result.solution, kSolverTolerance);
  INFO(verification.error);
  CHECK(verification.pass);
}

TEST_CASE("hi-pdlp-rewrite-solver-returns-final-candidate-at-iteration-limit",
          "[pdlp][hipdlp][phase7][solver][iteration-limit]") {
  const HighsLp lp = makeTinyFeasibleLp();
  hipdlp::api::SolverOptions options = makeCpuFixedOptions();
  options.iteration_limit = 3;
  const hipdlp::api::SolveResult result =
      hipdlp::algorithm::solve(options, hipdlp::model::OriginalModel(lp));

  INFO(result.error.message);
  REQUIRE(result.ok());
  CHECK(result.status == hipdlp::api::SolverStatus::kIterationLimit);
  CHECK(result.iteration_count == 3);
  REQUIRE(result.has_solution);
  REQUIRE(result.solution.value_valid);
  REQUIRE(result.solution.dual_valid);
  REQUIRE(result.solution.col_value.size() == 1);
  REQUIRE(result.solution.row_dual.size() == 1);
  CHECK(result.solution.col_value[0] == Approx(0.994008).margin(1e-12));
  CHECK(result.solution.row_value[0] == Approx(0.994008).margin(1e-12));
  CHECK(result.solution.row_dual[0] == Approx(1.009960032).margin(1e-12));
  CHECK(result.solution.col_dual[0] == Approx(-0.009960032).margin(1e-12));
}

TEST_CASE("hi-pdlp-rewrite-solver-honors-zero-time-limit-without-iteration",
          "[pdlp][hipdlp][phase7][solver][time-limit]") {
  const HighsLp lp = makeTinyFeasibleLp();
  hipdlp::api::SolverOptions options = makeCpuFixedOptions();
  options.time_limit = 0.0;
  const hipdlp::api::SolveResult result =
      hipdlp::algorithm::solve(options, hipdlp::model::OriginalModel(lp));

  INFO(result.error.message);
  REQUIRE(result.ok());
  CHECK(result.status == hipdlp::api::SolverStatus::kTimeLimit);
  CHECK(result.iteration_count == 0);
  REQUIRE(result.has_solution);
  REQUIRE(result.solution.value_valid);
  REQUIRE(result.solution.dual_valid);
  REQUIRE(result.solution.col_value.size() == 1);
  REQUIRE(result.solution.row_dual.size() == 1);
  CHECK(result.solution.col_value[0] == 0.0);
  CHECK(result.solution.row_value[0] == 0.0);
  CHECK(result.solution.row_dual[0] == 0.0);
  CHECK(result.solution.col_dual[0] == 1.0);

  const hipdlp::verification::KktVerificationResult verification =
      hipdlp::verification::verifyKkt(lp, result.solution, kSolverTolerance);
  CHECK_FALSE(verification.pass);
}

TEST_CASE("hi-pdlp-rewrite-solver-postsolves-scaled-solution",
          "[pdlp][hipdlp][phase7][solver][scaling][postsolve]") {
  const HighsLp lp = makeScaledFeasibleLp();
  hipdlp::api::SolverOptions options = makeCpuFixedOptions();
  options.scaling.mode = hipdlp::model::ScalingMode::kCurtisReid;
  const hipdlp::api::SolveResult result =
      hipdlp::algorithm::solve(options, hipdlp::model::OriginalModel(lp));

  INFO(result.error.message);
  REQUIRE(result.ok());
  CHECK(result.status == hipdlp::api::SolverStatus::kOptimal);
  REQUIRE(result.has_solution);
  REQUIRE(result.solution.value_valid);
  REQUIRE(result.solution.dual_valid);
  CHECK(result.solution.col_value[0] == Approx(2.0).margin(kSolverTolerance));
  CHECK(result.solution.row_value[0] ==
        Approx(8.0).margin(4.0 * kSolverTolerance));
  CHECK(result.solution.row_dual[0] == Approx(0.75).margin(kSolverTolerance));

  const hipdlp::verification::KktVerificationResult verification =
      hipdlp::verification::verifyKkt(lp, result.solution, kSolverTolerance);
  INFO(verification.error);
  CHECK(verification.pass);
}

TEST_CASE("hi-pdlp-rewrite-scaling-stages-solve-paired-inequalities",
          "[pdlp][hipdlp][phase7][solver][scaling][sctest]") {
  Highs highs;
  REQUIRE(highs.readModel(std::string(HIGHS_DIR) +
                          "/check/instances/sctest.mps") ==
          HighsStatus::kOk);
  const HighsLp lp = highs.getLp();
  hipdlp::api::SolverOptions options = makeCpuFixedOptions();
  options.iteration_limit = 100000;
  options.policy.step_size = hipdlp::api::StepSizePolicy::kPid;
  options.scaling.mode = hipdlp::model::ScalingMode::kLegacy;

  SECTION("Ruiz") {
    options.scaling.legacy_scaling_mode = kPdlpScalingRuiz;
  }
  SECTION("Pock-Chambolle") {
    options.scaling.legacy_scaling_mode = kPdlpScalingPC;
  }

  const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
      options, hipdlp::model::OriginalModel(lp));
  INFO(result.error.message);
  REQUIRE(result.ok());
  REQUIRE(result.status == hipdlp::api::SolverStatus::kOptimal);
  REQUIRE(result.iteration_count < options.iteration_limit);
  const hipdlp::verification::KktVerificationResult verification =
      hipdlp::verification::verifyKkt(lp, result.solution,
                                      options.kkt_tolerance);
  INFO(verification.error);
  CHECK(verification.pass);
}

TEST_CASE("hi-pdlp-rewrite-solver-restores-maximization-sense-and-offset",
          "[pdlp][hipdlp][phase7][solver][postsolve][maximization]") {
  const HighsLp lp = makeMaximizationLp();
  const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
      makeCpuFixedOptions(), hipdlp::model::OriginalModel(lp));

  INFO(result.error.message);
  REQUIRE(result.ok());
  CHECK(result.status == hipdlp::api::SolverStatus::kOptimal);
  REQUIRE(result.has_solution);
  REQUIRE(result.solution.value_valid);
  REQUIRE(result.solution.dual_valid);
  CHECK(result.solution.col_value[0] == Approx(2.0).margin(kSolverTolerance));
  CHECK(result.solution.row_value[0] ==
        Approx(4.0).margin(2.0 * kSolverTolerance));
  CHECK(result.solution.row_dual[0] == Approx(1.5).margin(kSolverTolerance));
  const double objective =
      lp.offset_ + lp.col_cost_[0] * result.solution.col_value[0];
  CHECK(objective == Approx(11.0).margin(3.0 * kSolverTolerance));

  const hipdlp::verification::KktVerificationResult verification =
      hipdlp::verification::verifyKkt(lp, result.solution, kSolverTolerance);
  INFO(verification.error);
  CHECK(verification.pass);
}

TEST_CASE("hi-pdlp-rewrite-solver-requires-independent-kkt-acceptance",
          "[pdlp][hipdlp][phase7][solver][kkt-gate]") {
  const HighsLp lp = makeTinyFeasibleLp();
  hipdlp::api::SolverOptions loose = makeCpuFixedOptions();
  loose.iteration_limit = 41;
  loose.optimality_tolerance = 0.4;
  loose.kkt_tolerance = 1.0;
  const hipdlp::api::SolveResult accepted =
      hipdlp::algorithm::solve(loose, hipdlp::model::OriginalModel(lp));
  INFO(accepted.error.message);
  REQUIRE(accepted.ok());
  REQUIRE(accepted.status == hipdlp::api::SolverStatus::kOptimal);
  REQUIRE(accepted.iteration_count == 40);

  hipdlp::api::SolverOptions strict = loose;
  strict.kkt_tolerance = 1e-12;
  const hipdlp::api::SolveResult rejected =
      hipdlp::algorithm::solve(strict, hipdlp::model::OriginalModel(lp));
  INFO(rejected.error.message);
  REQUIRE(rejected.ok());
  CHECK(rejected.status == hipdlp::api::SolverStatus::kIterationLimit);
  CHECK(rejected.iteration_count == 41);
  REQUIRE(rejected.has_solution);
  REQUIRE(rejected.solution.value_valid);
  REQUIRE(rejected.solution.dual_valid);
  CHECK(rejected.solution.col_value != accepted.solution.col_value);
  CHECK(rejected.solution.row_dual != accepted.solution.row_dual);

  const hipdlp::verification::KktVerificationResult verification =
      hipdlp::verification::verifyKkt(lp, rejected.solution,
                                      strict.kkt_tolerance);
  REQUIRE(verification.input_valid);
  REQUIRE(verification.metrics_finite);
  CHECK_FALSE(verification.pass);
}

TEST_CASE("hi-pdlp-rewrite-solver-rejects-unsupported-execution-scope",
          "[pdlp][hipdlp][phase7][solver][options][errors]") {
  HighsLp malformed = makeTinyFeasibleLp();
  malformed.a_matrix_.start_.clear();

  SECTION("an unavailable backend is rejected before model processing") {
    hipdlp::api::SolverOptions options = makeCpuFixedOptions();
    options.backend = hipdlp::backend::Kind::kOpenmp;
    const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
        options, hipdlp::model::OriginalModel(malformed));
    CHECK_FALSE(result.ok());
    CHECK(result.status == hipdlp::api::SolverStatus::kError);
    CHECK(result.error.code == hipdlp::api::ErrorCode::kBackendUnavailable);
    CHECK(result.iteration_count == 0);
    CHECK_FALSE(result.has_solution);
  }

  SECTION("an invalid backend value is rejected before model processing") {
    hipdlp::api::SolverOptions options = makeCpuFixedOptions();
    options.backend = static_cast<hipdlp::backend::Kind>(-1);
    const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
        options, hipdlp::model::OriginalModel(malformed));
    CHECK_FALSE(result.ok());
    CHECK(result.status == hipdlp::api::SolverStatus::kError);
    CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidOption);
    CHECK(result.error.message ==
          "HiPDLP rewrite solve received an invalid backend");
    CHECK(result.iteration_count == 0);
    CHECK_FALSE(result.has_solution);
  }

  SECTION("an invalid step-size policy is rejected before model processing") {
    hipdlp::api::SolverOptions options = makeCpuFixedOptions();
    options.policy.step_size = static_cast<hipdlp::api::StepSizePolicy>(7);
    const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
        options, hipdlp::model::OriginalModel(malformed));
    CHECK_FALSE(result.ok());
    CHECK(result.status == hipdlp::api::SolverStatus::kError);
    CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidOption);
    CHECK(result.error.message ==
          "HiPDLP rewrite solve received an invalid step-size policy");
    CHECK(result.iteration_count == 0);
    CHECK_FALSE(result.has_solution);
  }
}

TEST_CASE("hi-pdlp-rewrite-solver-returns-structured-model-errors",
          "[pdlp][hipdlp][phase7][solver][model][errors]") {
  SECTION("malformed model preserves the model error") {
    HighsLp malformed = makeTinyFeasibleLp();
    malformed.col_cost_.clear();
    const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
        makeCpuFixedOptions(), hipdlp::model::OriginalModel(malformed));
    CHECK_FALSE(result.ok());
    CHECK(result.status == hipdlp::api::SolverStatus::kError);
    CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidModel);
    CHECK_FALSE(result.error.message.empty());
    CHECK(result.iteration_count == 0);
    CHECK_FALSE(result.has_solution);
    CHECK_FALSE(result.solution.value_valid);
    CHECK_FALSE(result.solution.dual_valid);
  }

  SECTION("zero operator reports the deferred trivial-LP path") {
    HighsLp zero_operator = makeTinyFeasibleLp();
    zero_operator.a_matrix_.start_ = {0, 0};
    zero_operator.a_matrix_.index_.clear();
    zero_operator.a_matrix_.value_.clear();
    const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
        makeCpuFixedOptions(), hipdlp::model::OriginalModel(zero_operator));
    CHECK_FALSE(result.ok());
    CHECK(result.status == hipdlp::api::SolverStatus::kError);
    CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidModel);
    CHECK(result.error.message ==
          "HiPDLP rewrite does not yet support zero-operator LPs");
    CHECK(result.iteration_count == 0);
    CHECK_FALSE(result.has_solution);
    CHECK_FALSE(result.solution.value_valid);
    CHECK_FALSE(result.solution.dual_valid);
  }
}

TEST_CASE("hi-pdlp-public-options-register-defaults-and-round-trip",
          "[pdlp][hipdlp][phase7][options][public]") {
  Highs highs;
  REQUIRE(highs.setOptionValue("output_flag", false) == HighsStatus::kOk);

  std::string backend;
  REQUIRE(highs.getOptionValue(kHiPdlpBackendString, backend) ==
          HighsStatus::kOk);
  CHECK(backend == kHiPdlpAutoString);
  std::string scaling;
  REQUIRE(highs.getOptionValue(kHiPdlpScalingString, scaling) ==
          HighsStatus::kOk);
  CHECK(scaling == kHiPdlpAutoString);
  HighsInt sweeps = 0;
  REQUIRE(highs.getOptionValue(kHiPdlpCurtisReidIterationsString, sweeps) ==
          HighsStatus::kOk);
  CHECK(sweeps == kHiPdlpDefaultCurtisReidIterations);

  SECTION("every documented backend string is accepted") {
    const std::string values[] = {
        kHiPdlpAutoString, kHiPdlpBackendCpuString, kHiPdlpBackendCudaString,
        kHiPdlpBackendOpenmpString, kHiPdlpBackendRocmString};
    const hipdlp::backend::Kind kinds[] = {
        hipdlp::backend::Kind::kAuto, hipdlp::backend::Kind::kCpu,
        hipdlp::backend::Kind::kCuda, hipdlp::backend::Kind::kOpenmp,
        hipdlp::backend::Kind::kRocm};
    for (std::size_t index = 0; index < 5; ++index) {
      REQUIRE(highs.setOptionValue(kHiPdlpBackendString, values[index]) ==
              HighsStatus::kOk);
      std::string read_back;
      REQUIRE(highs.getOptionValue(kHiPdlpBackendString, read_back) ==
              HighsStatus::kOk);
      CHECK(read_back == values[index]);
      hipdlp::backend::Kind kind = hipdlp::backend::Kind::kAuto;
      REQUIRE(hipdlp::api::backendKindFromString(read_back, kind));
      CHECK(kind == kinds[index]);
    }
  }

  SECTION("every documented scaling string is accepted") {
    const std::string values[] = {kHiPdlpAutoString, kHighsOffString,
                                  kHiPdlpScalingCurtisReidString,
                                  kHiPdlpScalingLegacyString};
    const hipdlp::model::ScalingMode modes[] = {
        hipdlp::model::ScalingMode::kAuto, hipdlp::model::ScalingMode::kOff,
        hipdlp::model::ScalingMode::kCurtisReid,
        hipdlp::model::ScalingMode::kLegacy};
    for (std::size_t index = 0; index < 4; ++index) {
      REQUIRE(highs.setOptionValue(kHiPdlpScalingString, values[index]) ==
              HighsStatus::kOk);
      std::string read_back;
      REQUIRE(highs.getOptionValue(kHiPdlpScalingString, read_back) ==
              HighsStatus::kOk);
      CHECK(read_back == values[index]);
      hipdlp::model::ScalingMode mode = hipdlp::model::ScalingMode::kAuto;
      REQUIRE(hipdlp::api::scalingModeFromString(read_back, mode));
      CHECK(mode == modes[index]);
    }
  }

  SECTION("illegal option values are rejected and leave the option alone") {
    CHECK(highs.setOptionValue(kHiPdlpBackendString, "gpu") ==
          HighsStatus::kError);
    CHECK(highs.setOptionValue(kHiPdlpScalingString, "ruiz") ==
          HighsStatus::kError);
    CHECK(highs.setOptionValue(kHiPdlpCurtisReidIterationsString,
                               HighsInt{0}) == HighsStatus::kError);
    CHECK(highs.setOptionValue(kHiPdlpCurtisReidIterationsString,
                               HighsInt{-1}) == HighsStatus::kError);
    CHECK(highs.getOptions().hipdlp_backend == kHiPdlpAutoString);
    CHECK(highs.getOptions().hipdlp_scaling == kHiPdlpAutoString);
    CHECK(highs.getOptions().hipdlp_curtis_reid_iterations ==
          kHiPdlpDefaultCurtisReidIterations);
    hipdlp::backend::Kind kind = hipdlp::backend::Kind::kCpu;
    CHECK_FALSE(hipdlp::api::backendKindFromString("gpu", kind));
    hipdlp::model::ScalingMode mode = hipdlp::model::ScalingMode::kOff;
    CHECK_FALSE(hipdlp::api::scalingModeFromString("ruiz", mode));
  }

  SECTION("a positive Curtis-Reid count reaches the rewrite scaling options") {
    REQUIRE(highs.setOptionValue(kHiPdlpCurtisReidIterationsString,
                                 HighsInt{7}) == HighsStatus::kOk);
    REQUIRE(highs.setOptionValue(kHiPdlpScalingString,
                                 kHiPdlpScalingCurtisReidString) ==
            HighsStatus::kOk);
    const hipdlp::api::MigrationOptions rewrite(
        hipdlp::api::SolverImplementation::kRewrite,
        hipdlp::backend::Kind::kAuto);
    const hipdlp::api::OptionTranslationResult translated =
        hipdlp::api::translateOptions(highs.getOptions(), rewrite);
    INFO(translated.error.message);
    REQUIRE(translated.ok());
    CHECK(translated.options.scaling.mode ==
          hipdlp::model::ScalingMode::kCurtisReid);
    CHECK(translated.options.scaling.curtis_reid_sweeps == 7);
  }
}

TEST_CASE("hi-pdlp-effective-policy-translation",
          "[pdlp][hipdlp][phase7][options][migration][policy]") {
  Highs highs;
  REQUIRE(highs.setOptionValue("output_flag", false) == HighsStatus::kOk);
  const hipdlp::api::MigrationOptions rewrite(
      hipdlp::api::SolverImplementation::kRewrite,
      hipdlp::backend::Kind::kAuto);

  SECTION("the public default migration selects rewrite") {
    const hipdlp::api::MigrationOptions defaults =
        hipdlp::api::defaultMigrationOptions();
    CHECK(defaults.implementation ==
          hipdlp::api::SolverImplementation::kRewrite);
    CHECK(defaults.backend == hipdlp::backend::Kind::kAuto);
    const hipdlp::api::OptionTranslationResult translated =
        hipdlp::api::translateOptions(highs.getOptions(), defaults);
    INFO(translated.error.message);
    REQUIRE(translated.ok());
    CHECK(translated.options.implementation ==
          hipdlp::api::SolverImplementation::kRewrite);
  }

  SECTION("default HiGHS options translate for the rewrite selector") {
    REQUIRE(highs.getOptions().pdlp_step_size_strategy ==
            kPdlpStepSizeStrategyAdaptive);
    REQUIRE(highs.getOptions().pdlp_restart_strategy ==
            kPdlpRestartStrategyAdaptive);
    const hipdlp::api::OptionTranslationResult translated =
        hipdlp::api::translateOptions(highs.getOptions(), rewrite);
    INFO(translated.error.message);
    REQUIRE(translated.ok());
    CHECK(translated.options.implementation ==
          hipdlp::api::SolverImplementation::kRewrite);
    CHECK(translated.options.backend == hipdlp::backend::Kind::kAuto);
    CHECK(translated.options.scaling.mode == hipdlp::model::ScalingMode::kAuto);
    CHECK(translated.options.scaling.curtis_reid_sweeps ==
          kHiPdlpDefaultCurtisReidIterations);
    CHECK(translated.options.policy.step_size ==
          hipdlp::api::StepSizePolicy::kFixed);
    CHECK_FALSE(translated.options.policy.step_size ==
                hipdlp::api::StepSizePolicy::kPid);
    CHECK(translated.options.policy.step_size_substituted);
    CHECK(translated.options.policy.requested_step_size_strategy ==
          kPdlpStepSizeStrategyAdaptive);
    CHECK(translated.options.policy.restart ==
          hipdlp::api::RestartPolicy::kHalpern);
    CHECK(translated.options.policy.restart_substituted);
    CHECK(translated.options.policy.requested_restart_strategy ==
          kPdlpRestartStrategyAdaptive);
  }

  SECTION("explicit fixed and PID step sizes are preserved") {
    REQUIRE(highs.setOptionValue("pdlp_step_size_strategy",
                                 kPdlpStepSizeStrategyFixed) ==
            HighsStatus::kOk);
    const hipdlp::api::OptionTranslationResult fixed =
        hipdlp::api::translateOptions(highs.getOptions(), rewrite);
    REQUIRE(fixed.ok());
    CHECK(fixed.options.policy.step_size ==
          hipdlp::api::StepSizePolicy::kFixed);
    CHECK_FALSE(fixed.options.policy.step_size_substituted);

    REQUIRE(highs.setOptionValue("pdlp_step_size_strategy",
                                 kPdlpStepSizeStrategyPid) == HighsStatus::kOk);
    const hipdlp::api::OptionTranslationResult pid =
        hipdlp::api::translateOptions(highs.getOptions(), rewrite);
    REQUIRE(pid.ok());
    CHECK(pid.options.policy.step_size == hipdlp::api::StepSizePolicy::kPid);
    CHECK_FALSE(pid.options.policy.step_size_substituted);
  }

  SECTION("Malitsky-Pock is an explicit unsupported-option error") {
    REQUIRE(highs.setOptionValue("pdlp_step_size_strategy",
                                 kPdlpStepSizeStrategyMalitskyPock) ==
            HighsStatus::kOk);
    const hipdlp::api::OptionTranslationResult translated =
        hipdlp::api::translateOptions(highs.getOptions(), rewrite);
    CHECK_FALSE(translated.ok());
    CHECK(translated.error.code == hipdlp::api::ErrorCode::kInvalidOption);
    CHECK(translated.error.message.find("Malitsky-Pock") != std::string::npos);

    // The retained internal legacy selector keeps its own semantics.
    const hipdlp::api::MigrationOptions legacy_selector(
        hipdlp::api::SolverImplementation::kLegacy,
        hipdlp::backend::Kind::kAuto);
    const hipdlp::api::OptionTranslationResult legacy =
        hipdlp::api::translateOptions(highs.getOptions(), legacy_selector);
    INFO(legacy.error.message);
    CHECK(legacy.ok());
  }

  SECTION("restart requests all resolve to Halpern") {
    const HighsInt requested[] = {
        kPdlpRestartStrategyOff, kPdlpRestartStrategyFixed,
        kPdlpRestartStrategyAdaptive, kPdlpRestartStrategyHalpern};
    for (std::size_t index = 0; index < 4; ++index) {
      REQUIRE(highs.setOptionValue("pdlp_restart_strategy", requested[index]) ==
              HighsStatus::kOk);
      const hipdlp::api::OptionTranslationResult translated =
          hipdlp::api::translateOptions(highs.getOptions(), rewrite);
      INFO(translated.error.message);
      REQUIRE(translated.ok());
      CHECK(translated.options.policy.restart ==
            hipdlp::api::RestartPolicy::kHalpern);
      CHECK(translated.options.policy.restart_substituted ==
            (requested[index] != kPdlpRestartStrategyHalpern));
    }
    REQUIRE(highs.setOptionValue("pdlp_restart_strategy",
                                 kPdlpRestartStrategyHalpern) ==
            HighsStatus::kOk);
    HighsOptions restart_off_options = highs.getOptions();
    restart_off_options.pdlp_features_off = kPdlpRestartOff;
    const hipdlp::api::OptionTranslationResult restart_off =
        hipdlp::api::translateOptions(restart_off_options, rewrite);
    REQUIRE(restart_off.ok());
    CHECK(restart_off.options.policy.restart ==
          hipdlp::api::RestartPolicy::kHalpern);
    CHECK(restart_off.options.policy.restart_substituted);
  }

  SECTION("a migration backend override outranks the public option") {
    REQUIRE(highs.setOptionValue(kHiPdlpBackendString,
                                 kHiPdlpBackendCudaString) == HighsStatus::kOk);
    const hipdlp::api::OptionTranslationResult from_public =
        hipdlp::api::translateOptions(highs.getOptions(), rewrite);
    REQUIRE(from_public.ok());
    CHECK(from_public.options.backend == hipdlp::backend::Kind::kCuda);

    const hipdlp::api::MigrationOptions override_cpu(
        hipdlp::api::SolverImplementation::kRewrite,
        hipdlp::backend::Kind::kCpu);
    const hipdlp::api::OptionTranslationResult from_override =
        hipdlp::api::translateOptions(highs.getOptions(), override_cpu);
    REQUIRE(from_override.ok());
    CHECK(from_override.options.backend == hipdlp::backend::Kind::kCpu);
  }

  SECTION("an invalid stored option string is rejected") {
    HighsOptions options = highs.getOptions();
    options.hipdlp_backend = "gpu";
    CHECK_FALSE(hipdlp::api::translateOptions(options, rewrite).ok());
    options = highs.getOptions();
    options.hipdlp_scaling = "ruiz";
    CHECK_FALSE(hipdlp::api::translateOptions(options, rewrite).ok());
    options = highs.getOptions();
    options.hipdlp_curtis_reid_iterations = 0;
    CHECK_FALSE(hipdlp::api::translateOptions(options, rewrite).ok());
  }
}

TEST_CASE("hi-pdlp-rewrite-solves-with-default-style-highs-options",
          "[pdlp][hipdlp][phase7][solver][options][defaults]") {
  Highs highs;
  REQUIRE(highs.setOptionValue("output_flag", false) == HighsStatus::kOk);
  // Nothing else is changed: adaptive step size and auto scaling are defaults.
  const hipdlp::api::MigrationOptions rewrite(
      hipdlp::api::SolverImplementation::kRewrite, hipdlp::backend::Kind::kCpu);
  const hipdlp::api::OptionTranslationResult translated =
      hipdlp::api::translateOptions(highs.getOptions(), rewrite);
  INFO(translated.error.message);
  REQUIRE(translated.ok());
  REQUIRE(translated.options.policy.step_size ==
          hipdlp::api::StepSizePolicy::kFixed);

  const HighsLp lp = makeTinyFeasibleLp();
  RecordingLogSink sink;
  hipdlp::algorithm::SolveContext context;
  context.log = &sink;
  const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
      translated.options, hipdlp::model::OriginalModel(lp), context);
  INFO(result.error.message);
  REQUIRE(result.ok());
  CHECK(result.status == hipdlp::api::SolverStatus::kOptimal);
  REQUIRE(result.has_solution);
  CHECK(result.solution.col_value[0] == Approx(1.0).margin(1e-5));

  const hipdlp::verification::KktVerificationResult verification =
      hipdlp::verification::verifyKkt(lp, result.solution,
                                      translated.options.kkt_tolerance);
  INFO(verification.error);
  CHECK(verification.pass);

  // The adaptive substitution is reported once and is never PID.
  CHECK(sink.countContaining("adaptive") == 1);
  CHECK(sink.countContaining("effective fixed") == 1);
  CHECK(sink.countContaining("PID") == 0);
  CHECK(sink.countContaining("HiPDLP rewrite backend") == 1);
  CHECK(sink.countContaining("HiPDLP rewrite restart") == 1);
  CHECK(sink.countContaining("HiPDLP rewrite scaling") == 1);
  CHECK(sink.countContaining("HiPDLP rewrite finished") == 1);
}

TEST_CASE("hi-pdlp-rewrite-logs-only-every-fortieth-convergence-boundary",
          "[pdlp][hipdlp][phase7][solver][diagnostics]") {
  Highs highs;
  REQUIRE(highs.readModel(std::string(HIGHS_DIR) +
                          "/check/instances/sctest.mps") ==
          HighsStatus::kOk);
  const HighsLp lp = highs.getLp();
  hipdlp::api::SolverOptions options = makeCpuFixedOptions();
  options.iteration_limit = 1601;
  options.optimality_tolerance = std::numeric_limits<double>::denorm_min();
  options.kkt_tolerance = std::numeric_limits<double>::denorm_min();
  options.policy.step_size = hipdlp::api::StepSizePolicy::kPid;
  options.scaling.mode = hipdlp::model::ScalingMode::kAuto;
  RecordingLogSink sink;
  hipdlp::algorithm::SolveContext context;
  context.log = &sink;

  const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
      options, hipdlp::model::OriginalModel(lp), context);

  INFO(result.error.message);
  REQUIRE(result.ok());
  REQUIRE(result.status == hipdlp::api::SolverStatus::kIterationLimit);
  REQUIRE(result.iteration_count == 1601);
  REQUIRE(sink.countContaining("HiPDLP rewrite boundary") == 1);
  std::string boundary;
  for (std::size_t index = 0; index < sink.messages.size(); ++index) {
    if (sink.messages[index].find("HiPDLP rewrite boundary") !=
        std::string::npos) {
      boundary = sink.messages[index];
    }
  }
  CHECK(boundary.find("iteration 1600") != std::string::npos);
  CHECK(boundary.find("relative primal ") != std::string::npos);
  CHECK(boundary.find("dual ") != std::string::npos);
  CHECK(boundary.find("gap ") != std::string::npos);
  CHECK(boundary.find("tau ") != std::string::npos);
  CHECK(boundary.find("sigma ") != std::string::npos);
  CHECK(boundary.find("omega ") != std::string::npos);
  CHECK(boundary.find("eta ") != std::string::npos);
  CHECK(boundary.find("restart condition ") != std::string::npos);
  CHECK(boundary.find("PID outcome ") != std::string::npos);
}

TEST_CASE("hi-pdlp-rewrite-accepts-an-already-optimal-starting-point",
          "[pdlp][hipdlp][phase7][solver][iteration-zero]") {
  const HighsLp lp = makeAlreadyOptimalLp();
  hipdlp::api::SolverOptions options = makeCpuFixedOptions();
  hipdlp::algorithm::SolveDiagnostics diagnostics;
  hipdlp::algorithm::SolveContext context;
  context.diagnostics = &diagnostics;
  const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
      options, hipdlp::model::OriginalModel(lp), context);

  INFO(result.error.message);
  REQUIRE(result.ok());
  CHECK(result.status == hipdlp::api::SolverStatus::kOptimal);
  CHECK(result.iteration_count == 0);
  REQUIRE(result.has_solution);
  REQUIRE(result.solution.value_valid);
  REQUIRE(result.solution.dual_valid);
  CHECK(result.solution.col_value[0] == 0.0);

  // The initial point still had to survive postsolve and the KKT verifier.
  const hipdlp::verification::KktVerificationResult verification =
      hipdlp::verification::verifyKkt(lp, result.solution,
                                      options.kkt_tolerance);
  INFO(verification.error);
  CHECK(verification.pass);
  CHECK(diagnostics.restarts == 0);
  CHECK(diagnostics.resource_creation_calls_during_iterations == 0);
}

TEST_CASE("hi-pdlp-rewrite-allocates-every-backend-vector-before-iterating",
          "[pdlp][hipdlp][phase7][solver][resources]") {
  const HighsLp lp = makeTinyFeasibleLp();
  hipdlp::api::SolverOptions options = makeCpuFixedOptions();
  options.iteration_limit = 200;
  options.optimality_tolerance = 1e-9;
  hipdlp::algorithm::SolveDiagnostics diagnostics;
  hipdlp::algorithm::SolveContext context;
  context.diagnostics = &diagnostics;
  const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
      options, hipdlp::model::OriginalModel(lp), context);

  INFO(result.error.message);
  REQUIRE(result.ok());
  REQUIRE(result.iteration_count > 40);
  // One problem binding, one workspace and sixteen iteration vectors.
  CHECK(diagnostics.resource_creation_calls_before_iterations == 18);
  CHECK(diagnostics.resource_creation_calls_during_iterations == 0);
  CHECK(diagnostics.restarts > 0);
  CHECK(diagnostics.primal_weight_updates == 0);
}

TEST_CASE("hi-pdlp-rewrite-consumes-the-pid-primal-weight-update",
          "[pdlp][hipdlp][phase7][solver][pid][restart]") {
  hipdlp::api::SolverOptions fixed = makeCpuFixedOptions();
  fixed.iteration_limit = 200;
  fixed.optimality_tolerance = 1e-9;
  hipdlp::api::SolverOptions pid = fixed;
  pid.policy.step_size = hipdlp::api::StepSizePolicy::kPid;

  SECTION("a defined residual ratio changes the later trajectory") {
    const HighsLp lp = makeTinyFeasibleLp();
    hipdlp::algorithm::SolveDiagnostics fixed_diagnostics;
    hipdlp::algorithm::SolveContext fixed_context;
    fixed_context.diagnostics = &fixed_diagnostics;
    const hipdlp::api::SolveResult fixed_result = hipdlp::algorithm::solve(
        fixed, hipdlp::model::OriginalModel(lp), fixed_context);
    hipdlp::algorithm::SolveDiagnostics pid_diagnostics;
    hipdlp::algorithm::SolveContext pid_context;
    pid_context.diagnostics = &pid_diagnostics;
    const hipdlp::api::SolveResult pid_result = hipdlp::algorithm::solve(
        pid, hipdlp::model::OriginalModel(lp), pid_context);

    INFO(fixed_result.error.message);
    REQUIRE(fixed_result.ok());
    INFO(pid_result.error.message);
    REQUIRE(pid_result.ok());
    CHECK(fixed_diagnostics.restarts == pid_diagnostics.restarts);
    CHECK(fixed_diagnostics.primal_weight_updates == 0);
    CHECK(pid_diagnostics.primal_weight_updates == pid_diagnostics.restarts);
    CHECK(pid_diagnostics.primal_weight_updates > 0);
    CHECK(pid_result.iteration_count == fixed_result.iteration_count);
    CHECK(pid_result.solution.col_value != fixed_result.solution.col_value);
  }

  SECTION("a boxed column's zero dual residual does not disable the PID") {
    const HighsLp lp = makeBoxedColumnLp();
    const hipdlp::api::SolveResult fixed_result =
        hipdlp::algorithm::solve(fixed, hipdlp::model::OriginalModel(lp));
    hipdlp::algorithm::SolveDiagnostics pid_diagnostics;
    hipdlp::algorithm::SolveContext pid_context;
    pid_context.diagnostics = &pid_diagnostics;
    const hipdlp::api::SolveResult pid_result = hipdlp::algorithm::solve(
        pid, hipdlp::model::OriginalModel(lp), pid_context);

    INFO(fixed_result.error.message);
    REQUIRE(fixed_result.ok());
    INFO(pid_result.error.message);
    REQUIRE(pid_result.ok());
    CHECK(pid_diagnostics.primal_weight_updates > 0);
    // The bound projection makes a boxed column's dual residual zero; the
    // controller must still move the weight, and with it the trajectory.
    CHECK(pid_result.solution.col_value != fixed_result.solution.col_value);
  }
}

TEST_CASE("hi-pdlp-rewrite-postsolves-every-scaling-mode",
          "[pdlp][hipdlp][phase7][solver][scaling][postsolve]") {
  const HighsLp lp = makeScaledFeasibleLp();
  const hipdlp::model::ScalingMode modes[] = {
      hipdlp::model::ScalingMode::kAuto, hipdlp::model::ScalingMode::kOff,
      hipdlp::model::ScalingMode::kCurtisReid,
      hipdlp::model::ScalingMode::kLegacy};
  for (std::size_t index = 0; index < 4; ++index) {
    hipdlp::api::SolverOptions options = makeCpuFixedOptions();
    options.scaling.mode = modes[index];
    const hipdlp::api::SolveResult result =
        hipdlp::algorithm::solve(options, hipdlp::model::OriginalModel(lp));
    INFO(result.error.message);
    REQUIRE(result.ok());
    CHECK(result.status == hipdlp::api::SolverStatus::kOptimal);
    REQUIRE(result.has_solution);
    CHECK(result.solution.col_value[0] ==
          Approx(2.0).margin(4.0 * kSolverTolerance));
    const hipdlp::verification::KktVerificationResult verification =
        hipdlp::verification::verifyKkt(lp, result.solution,
                                        options.kkt_tolerance);
    INFO(verification.error);
    CHECK(verification.pass);
  }
}

TEST_CASE("hi-pdlp-rewrite-uses-only-the-remaining-global-time-budget",
          "[pdlp][hipdlp][phase7][solver][time-limit][api]") {
  Highs highs;
  REQUIRE(highs.setOptionValue("output_flag", false) == HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("pdlp_step_size_strategy",
                               kPdlpStepSizeStrategyFixed) == HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("pdlp_optimality_tolerance", kSolverTolerance) ==
          HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("kkt_tolerance", kSolverTolerance) ==
          HighsStatus::kOk);
  const HighsLp lp = makeTinyFeasibleLp();
  const hipdlp::api::MigrationOptions rewrite(
      hipdlp::api::SolverImplementation::kRewrite, hipdlp::backend::Kind::kCpu);

  SECTION("an exhausted global budget leaves nothing for the rewrite") {
    REQUIRE(highs.setOptionValue("time_limit", 0.005) == HighsStatus::kOk);
    HighsTimer timer;
    timer.setPrintfFlag(false, false);
    timer.start();
    while (timer.read() < 0.01) {
    }
    HighsBasis basis;
    basis.valid = true;
    HighsSolution solution;
    HighsModelStatus model_status = HighsModelStatus::kNotset;
    HighsInfo info;
    HighsCallback callback(&highs);
    REQUIRE(hipdlp::api::solve(highs.getOptions(), timer, lp, basis, solution,
                               model_status, info, callback,
                               rewrite) == HighsStatus::kOk);
    timer.stop();
    CHECK(model_status == HighsModelStatus::kTimeLimit);
    CHECK(info.pdlp_iteration_count == 0);
    REQUIRE(solution.value_valid);
    REQUIRE(solution.dual_valid);
  }

  SECTION("the infinite sentinel is preserved") {
    REQUIRE(highs.getOptions().time_limit == kHighsInf);
    HighsTimer timer;
    timer.setPrintfFlag(false, false);
    timer.start();
    while (timer.read() < 0.001) {
    }
    HighsBasis basis;
    basis.valid = true;
    HighsSolution solution;
    HighsModelStatus model_status = HighsModelStatus::kNotset;
    HighsInfo info;
    HighsCallback callback(&highs);
    REQUIRE(hipdlp::api::solve(highs.getOptions(), timer, lp, basis, solution,
                               model_status, info, callback,
                               rewrite) == HighsStatus::kOk);
    timer.stop();
    CHECK(model_status == HighsModelStatus::kOptimal);
    CHECK(info.pdlp_iteration_count > 0);
  }
}

TEST_CASE("hi-pdlp-rewrite-emits-no-user-log-when-output-is-off",
          "[pdlp][hipdlp][phase7][solver][logging]") {
  const HighsLp lp = makeTinyFeasibleLp();
  const hipdlp::api::MigrationOptions rewrite(
      hipdlp::api::SolverImplementation::kRewrite, hipdlp::backend::Kind::kCpu);
  const bool output_flags[] = {false, true};
  const char* log_files[] = {"hipdlp-rewrite-log-off.log",
                             "hipdlp-rewrite-log-on.log"};
  for (std::size_t index = 0; index < 2; ++index) {
    std::remove(log_files[index]);
    Highs highs;
    REQUIRE(highs.setOptionValue("log_to_console", false) == HighsStatus::kOk);
    REQUIRE(highs.setOptionValue(kLogFileString, log_files[index]) ==
            HighsStatus::kOk);
    REQUIRE(highs.setOptionValue("output_flag", output_flags[index]) ==
            HighsStatus::kOk);
    REQUIRE(highs.setOptionValue("pdlp_optimality_tolerance",
                                 kSolverTolerance) == HighsStatus::kOk);
    REQUIRE(highs.setOptionValue("kkt_tolerance", kSolverTolerance) ==
            HighsStatus::kOk);
    HighsTimer timer;
    timer.setPrintfFlag(false, false);
    HighsBasis basis;
    basis.valid = true;
    HighsSolution solution;
    HighsModelStatus model_status = HighsModelStatus::kNotset;
    HighsInfo info;
    HighsCallback callback(&highs);
    REQUIRE(hipdlp::api::solve(highs.getOptions(), timer, lp, basis, solution,
                               model_status, info, callback,
                               rewrite) == HighsStatus::kOk);
    CHECK(model_status == HighsModelStatus::kOptimal);
    REQUIRE(highs.setOptionValue(kLogFileString, "") == HighsStatus::kOk);

    std::ifstream stream(log_files[index]);
    const std::string contents((std::istreambuf_iterator<char>(stream)),
                               std::istreambuf_iterator<char>());
    stream.close();
    if (output_flags[index]) {
      CHECK(contents.find("HiPDLP rewrite") != std::string::npos);
    } else {
      CHECK(contents.empty());
    }
    std::remove(log_files[index]);
  }
}

TEST_CASE("hi-pdlp-public-path-never-enters-pdhg-for-a-zero-operator-lp",
          "[pdlp][hipdlp][phase7][solver][trivial]") {
  Highs highs;
  REQUIRE(highs.setOptionValue("output_flag", false) == HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("presolve", "off") == HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("solver", kHiPdlpString) == HighsStatus::kOk);

  SECTION("a zero constraint matrix is solved directly") {
    HighsLp lp = makeTinyFeasibleLp();
    lp.row_lower_ = {-kHighsInf};
    lp.row_upper_ = {kHighsInf};
    lp.a_matrix_.start_ = {0, 0};
    lp.a_matrix_.index_.clear();
    lp.a_matrix_.value_.clear();
    REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
    REQUIRE(highs.run() == HighsStatus::kOk);
    CHECK(highs.getModelStatus() == HighsModelStatus::kOptimal);
    CHECK(highs.getInfo().pdlp_iteration_count == 0);
    CHECK(highs.getSolution().col_value[0] == Approx(0.0).margin(1e-12));
  }

  SECTION("a row-free model is solved directly") {
    HighsLp lp;
    lp.num_col_ = 1;
    lp.num_row_ = 0;
    lp.col_cost_ = {1.0};
    lp.col_lower_ = {0.0};
    lp.col_upper_ = {kHighsInf};
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    lp.a_matrix_.num_col_ = 1;
    lp.a_matrix_.num_row_ = 0;
    lp.a_matrix_.start_ = {0, 0};
    lp.sense_ = ObjSense::kMinimize;
    REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
    REQUIRE(highs.run() == HighsStatus::kOk);
    CHECK(highs.getModelStatus() == HighsModelStatus::kOptimal);
    CHECK(highs.getInfo().pdlp_iteration_count == 0);
  }
}

TEST_CASE("hi-pdlp-public-path-routes-explicit-cuda-pid-to-rewrite",
          "[pdlp][hipdlp][phase7][solver][public][cuda][pid]") {
  const char* log_file = "hipdlp-public-cuda-pid.log";
  std::remove(log_file);
  Highs highs;
  REQUIRE(highs.setOptionValue("log_to_console", false) == HighsStatus::kOk);
  REQUIRE(highs.setOptionValue(kLogFileString, log_file) == HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("output_flag", true) == HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("solver", kHiPdlpString) == HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("presolve", "off") == HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("hipdlp_backend", "cuda") ==
          HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("hipdlp_scaling", "auto") ==
          HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("pdlp_step_size_strategy",
                               kPdlpStepSizeStrategyPid) == HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("pdlp_optimality_tolerance", kSolverTolerance) ==
          HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("kkt_tolerance", kSolverTolerance) ==
          HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("time_limit", 60.0) == HighsStatus::kOk);
  REQUIRE(highs.readModel(std::string(HIGHS_DIR) +
                          "/check/instances/sctest.mps") ==
          HighsStatus::kOk);

  const bool cuda_usable = cudaBackendUsable();
  const HighsStatus run_status = highs.run();
  REQUIRE(highs.setOptionValue(kLogFileString, "") == HighsStatus::kOk);
  std::ifstream stream(log_file);
  const std::string contents((std::istreambuf_iterator<char>(stream)),
                             std::istreambuf_iterator<char>());
  stream.close();
  std::remove(log_file);

  if (!cuda_usable) {
    CHECK(run_status == HighsStatus::kError);
    CHECK(highs.getModelStatus() == HighsModelStatus::kSolveError);
    CHECK(contents.find("requested HiPDLP backend is unavailable") !=
          std::string::npos);
    CHECK(contents.find("effective cpu") == std::string::npos);
    return;
  }

  REQUIRE(run_status == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
  const HighsLp& lp = highs.getLp();
  const HighsSolution& solution = highs.getSolution();
  REQUIRE(solution.value_valid);
  REQUIRE(solution.dual_valid);
  CHECK(solution.col_value.size() == static_cast<std::size_t>(lp.num_col_));
  CHECK(solution.row_value.size() == static_cast<std::size_t>(lp.num_row_));
  CHECK(solution.col_dual.size() == static_cast<std::size_t>(lp.num_col_));
  CHECK(solution.row_dual.size() == static_cast<std::size_t>(lp.num_row_));
  const hipdlp::verification::KktVerificationResult verification =
      hipdlp::verification::verifyKkt(lp, solution, kSolverTolerance);
  INFO(verification.error);
  CHECK(verification.pass);
  CHECK(contents.find(
            "HiPDLP rewrite backend: requested cuda, effective cuda") !=
        std::string::npos);
  CHECK(contents.find("effective PID") != std::string::npos);
}

TEST_CASE("hi-pdlp-rewrite-runs-one-control-loop-on-every-backend",
          "[pdlp][hipdlp][phase7][solver][backend][cuda]") {
  const bool cuda_usable = cudaBackendUsable();
  const HighsLp lp = makeScaledFeasibleLp();

  hipdlp::api::SolverOptions cpu_options = makeCpuFixedOptions();
  cpu_options.scaling.mode = hipdlp::model::ScalingMode::kCurtisReid;
  const hipdlp::api::SolveResult cpu_result =
      hipdlp::algorithm::solve(cpu_options, hipdlp::model::OriginalModel(lp));
  INFO(cpu_result.error.message);
  REQUIRE(cpu_result.ok());
  REQUIRE(cpu_result.status == hipdlp::api::SolverStatus::kOptimal);

  SECTION("an unimplemented backend is never a silent fallback") {
    const hipdlp::backend::Kind unavailable[] = {hipdlp::backend::Kind::kOpenmp,
                                                 hipdlp::backend::Kind::kRocm};
    for (std::size_t index = 0; index < 2; ++index) {
      hipdlp::api::SolverOptions options = cpu_options;
      options.backend = unavailable[index];
      const hipdlp::api::SolveResult result =
          hipdlp::algorithm::solve(options, hipdlp::model::OriginalModel(lp));
      CHECK_FALSE(result.ok());
      CHECK(result.status == hipdlp::api::SolverStatus::kError);
      CHECK(result.error.code == hipdlp::api::ErrorCode::kBackendUnavailable);
      CHECK_FALSE(result.has_solution);
      CHECK(result.iteration_count == 0);
    }
  }

  SECTION("auto selects CUDA when it is usable and CPU otherwise") {
    hipdlp::api::SolverOptions options = cpu_options;
    options.backend = hipdlp::backend::Kind::kAuto;
    RecordingLogSink sink;
    hipdlp::algorithm::SolveContext context;
    context.log = &sink;
    const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
        options, hipdlp::model::OriginalModel(lp), context);
    INFO(result.error.message);
    REQUIRE(result.ok());
    CHECK(result.status == cpu_result.status);
    CHECK(result.iteration_count == cpu_result.iteration_count);
    REQUIRE(result.has_solution);
    CHECK(result.solution.col_value[0] ==
          Approx(cpu_result.solution.col_value[0]).margin(kCudaTolerance));
    CHECK(result.solution.row_dual[0] ==
          Approx(cpu_result.solution.row_dual[0]).margin(kCudaTolerance));
    CHECK(sink.countContaining("requested auto, effective ") == 1);
    CHECK(sink.countContaining(cuda_usable ? "effective cuda"
                                           : "effective cpu") == 1);
    const hipdlp::verification::KktVerificationResult verification =
        hipdlp::verification::verifyKkt(lp, result.solution,
                                        options.kkt_tolerance);
    INFO(verification.error);
    CHECK(verification.pass);
  }

  SECTION("explicit CUDA either runs the same loop or reports unavailable") {
    hipdlp::api::SolverOptions options = cpu_options;
    options.backend = hipdlp::backend::Kind::kCuda;
    const hipdlp::api::SolveResult result =
        hipdlp::algorithm::solve(options, hipdlp::model::OriginalModel(lp));
    if (!cuda_usable) {
      CHECK_FALSE(result.ok());
      CHECK(result.status == hipdlp::api::SolverStatus::kError);
      CHECK(result.error.code == hipdlp::api::ErrorCode::kBackendUnavailable);
      CHECK_FALSE(result.has_solution);
      CHECK(result.iteration_count == 0);
      return;
    }
    INFO(result.error.message);
    REQUIRE(result.ok());
    CHECK(result.status == cpu_result.status);
    CHECK(result.iteration_count == cpu_result.iteration_count);
    REQUIRE(result.has_solution);
    CHECK(result.solution.col_value[0] ==
          Approx(cpu_result.solution.col_value[0]).margin(kCudaTolerance));
    CHECK(result.solution.row_value[0] ==
          Approx(cpu_result.solution.row_value[0]).margin(kCudaTolerance));
    CHECK(result.solution.row_dual[0] ==
          Approx(cpu_result.solution.row_dual[0]).margin(kCudaTolerance));
    CHECK(result.solution.col_dual[0] ==
          Approx(cpu_result.solution.col_dual[0]).margin(kCudaTolerance));
    const hipdlp::verification::KktVerificationResult verification =
        hipdlp::verification::verifyKkt(lp, result.solution,
                                        options.kkt_tolerance);
    INFO(verification.error);
    CHECK(verification.pass);
  }

  SECTION("explicit CUDA takes the same limit and restart decisions") {
    if (!cuda_usable) return;
    const HighsLp tiny = makeTinyFeasibleLp();

    hipdlp::api::SolverOptions cpu_limit = makeCpuFixedOptions();
    cpu_limit.iteration_limit = 137;
    cpu_limit.optimality_tolerance = 1e-16;
    hipdlp::algorithm::SolveDiagnostics cpu_diagnostics;
    hipdlp::algorithm::SolveContext cpu_context;
    cpu_context.diagnostics = &cpu_diagnostics;
    const hipdlp::api::SolveResult cpu_limited = hipdlp::algorithm::solve(
        cpu_limit, hipdlp::model::OriginalModel(tiny), cpu_context);

    hipdlp::api::SolverOptions cuda_limit = cpu_limit;
    cuda_limit.backend = hipdlp::backend::Kind::kCuda;
    hipdlp::algorithm::SolveDiagnostics cuda_diagnostics;
    hipdlp::algorithm::SolveContext cuda_context;
    cuda_context.diagnostics = &cuda_diagnostics;
    const hipdlp::api::SolveResult cuda_limited = hipdlp::algorithm::solve(
        cuda_limit, hipdlp::model::OriginalModel(tiny), cuda_context);

    INFO(cpu_limited.error.message);
    REQUIRE(cpu_limited.ok());
    INFO(cuda_limited.error.message);
    REQUIRE(cuda_limited.ok());
    // A non-multiple-of-40 limit is honoured exactly on both backends.
    CHECK(cpu_limited.status == hipdlp::api::SolverStatus::kIterationLimit);
    CHECK(cuda_limited.status == cpu_limited.status);
    CHECK(cpu_limited.iteration_count == 137);
    CHECK(cuda_limited.iteration_count == 137);
    CHECK(cuda_diagnostics.restarts == cpu_diagnostics.restarts);
    CHECK(cpu_diagnostics.restarts > 0);
    CHECK(cuda_diagnostics.resource_creation_calls_during_iterations == 0);
    CHECK(cuda_limited.solution.col_value[0] ==
          Approx(cpu_limited.solution.col_value[0]).margin(kCudaTolerance));

    hipdlp::api::SolverOptions cpu_pid = cpu_limit;
    cpu_pid.policy.step_size = hipdlp::api::StepSizePolicy::kPid;
    hipdlp::api::SolverOptions cuda_pid = cuda_limit;
    cuda_pid.policy.step_size = hipdlp::api::StepSizePolicy::kPid;
    hipdlp::algorithm::SolveDiagnostics cpu_pid_diagnostics;
    hipdlp::algorithm::SolveContext cpu_pid_context;
    cpu_pid_context.diagnostics = &cpu_pid_diagnostics;
    const hipdlp::api::SolveResult cpu_pid_result = hipdlp::algorithm::solve(
        cpu_pid, hipdlp::model::OriginalModel(tiny), cpu_pid_context);
    hipdlp::algorithm::SolveDiagnostics cuda_pid_diagnostics;
    hipdlp::algorithm::SolveContext cuda_pid_context;
    cuda_pid_context.diagnostics = &cuda_pid_diagnostics;
    const hipdlp::api::SolveResult cuda_pid_result = hipdlp::algorithm::solve(
        cuda_pid, hipdlp::model::OriginalModel(tiny), cuda_pid_context);
    INFO(cpu_pid_result.error.message);
    REQUIRE(cpu_pid_result.ok());
    INFO(cuda_pid_result.error.message);
    REQUIRE(cuda_pid_result.ok());
    CHECK(cpu_pid_diagnostics.primal_weight_updates > 0);
    CHECK(cuda_pid_diagnostics.primal_weight_updates ==
          cpu_pid_diagnostics.primal_weight_updates);
    CHECK(cuda_pid_result.status == cpu_pid_result.status);
    CHECK(cuda_pid_result.iteration_count == cpu_pid_result.iteration_count);
    CHECK(cuda_pid_result.solution.col_value[0] ==
          Approx(cpu_pid_result.solution.col_value[0]).margin(kCudaTolerance));
  }

  SECTION("explicit CUDA honours a zero time budget without iterating") {
    if (!cuda_usable) return;
    hipdlp::api::SolverOptions options = cpu_options;
    options.backend = hipdlp::backend::Kind::kCuda;
    options.time_limit = 0.0;
    const hipdlp::api::SolveResult result =
        hipdlp::algorithm::solve(options, hipdlp::model::OriginalModel(lp));
    INFO(result.error.message);
    REQUIRE(result.ok());
    CHECK(result.status == hipdlp::api::SolverStatus::kTimeLimit);
    CHECK(result.iteration_count == 0);
    REQUIRE(result.has_solution);
    REQUIRE(result.solution.value_valid);
    REQUIRE(result.solution.dual_valid);
  }
}
