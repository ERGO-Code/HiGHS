#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>

#include "HCheckConfig.h"
#include "catch.hpp"
#include "pdlp/hipdlp/algorithm/averaging_policy.hpp"

namespace {

hipdlp::algorithm::AveragingInput makeValidInput() {
  hipdlp::algorithm::AveragingInput input;
  input.primal_step_size = 4.0;
  input.dual_step_size = 9.0;
  input.completed_halpern_iterate = true;
  return input;
}

bool bitwiseEqual(const double left, const double right) {
  return std::memcmp(&left, &right, sizeof(double)) == 0;
}

void checkInvalidResult(const hipdlp::algorithm::AveragingResult& result) {
  CHECK_FALSE(result.ok());
  CHECK_FALSE(result.error.ok());
  CHECK_FALSE(result.next.history_valid);
  CHECK(result.next.sample_count == 0);
  CHECK(result.next.total_weight == 0.0);
  CHECK_FALSE(result.next.hasAverage());
  CHECK(result.sample_weight == 0.0);
  CHECK(result.normalization == 0.0);
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

}  // namespace

TEST_CASE("hi-pdlp-averaging-policy-empty-and-first-sample",
          "[pdlp][hipdlp][phase7b][averaging][policy][weights]") {
  const hipdlp::algorithm::AveragingState empty;
  CHECK(empty.history_valid);
  CHECK(empty.sample_count == 0);
  CHECK(empty.total_weight == 0.0);
  CHECK_FALSE(empty.hasAverage());

  const hipdlp::algorithm::AveragingResult result =
      hipdlp::algorithm::updateAveraging(makeValidInput());
  INFO(result.error.message);
  REQUIRE(result.ok());
  CHECK(result.sample_weight == 6.0);
  CHECK(result.next.sample_count == 1);
  CHECK(result.next.total_weight == 6.0);
  CHECK(result.normalization == 1.0 / 6.0);
  CHECK(result.next.history_valid);
  CHECK(result.next.hasAverage());
}

TEST_CASE("hi-pdlp-averaging-policy-variable-step-sequence",
          "[pdlp][hipdlp][phase7b][averaging][policy][weights]"
          "[normalization]") {
  hipdlp::algorithm::AveragingInput input = makeValidInput();
  const hipdlp::algorithm::AveragingResult first =
      hipdlp::algorithm::updateAveraging(input);
  REQUIRE(first.ok());

  input.previous = first.next;
  input.primal_step_size = 1.0;
  input.dual_step_size = 16.0;
  const hipdlp::algorithm::AveragingResult second =
      hipdlp::algorithm::updateAveraging(input);
  REQUIRE(second.ok());
  CHECK(second.sample_weight == 4.0);
  CHECK(second.next.sample_count == 2);
  CHECK(second.next.total_weight == 10.0);
  CHECK(second.normalization == 0.1);

  input.previous = second.next;
  input.primal_step_size = 0.25;
  input.dual_step_size = 4.0;
  const hipdlp::algorithm::AveragingResult third =
      hipdlp::algorithm::updateAveraging(input);
  REQUIRE(third.ok());
  CHECK(third.sample_weight == 1.0);
  CHECK(third.next.sample_count == 3);
  CHECK(third.next.total_weight == 11.0);
  CHECK(third.normalization == 1.0 / 11.0);
}

TEST_CASE("hi-pdlp-averaging-policy-requires-completed-halpern-iterate",
          "[pdlp][hipdlp][phase7b][averaging][policy][eligibility]"
          "[errors]") {
  hipdlp::algorithm::AveragingInput input = makeValidInput();
  input.completed_halpern_iterate = false;
  checkInvalidResult(hipdlp::algorithm::updateAveraging(input));
  CHECK(input.previous.history_valid);
  CHECK(input.previous.sample_count == 0);
  CHECK(input.previous.total_weight == 0.0);
}

TEST_CASE("hi-pdlp-averaging-policy-reset-is-exact-and-reusable",
          "[pdlp][hipdlp][phase7b][averaging][policy][reset]"
          "[determinism]") {
  const hipdlp::algorithm::AveragingResult populated =
      hipdlp::algorithm::updateAveraging(makeValidInput());
  REQUIRE(populated.ok());

  const hipdlp::algorithm::AveragingState first_reset =
      hipdlp::algorithm::resetAveraging();
  const hipdlp::algorithm::AveragingState second_reset =
      hipdlp::algorithm::resetAveraging();
  CHECK(first_reset.history_valid);
  CHECK(first_reset.sample_count == 0);
  CHECK_FALSE(std::signbit(first_reset.total_weight));
  CHECK(first_reset.total_weight == 0.0);
  CHECK_FALSE(first_reset.hasAverage());
  CHECK(second_reset.history_valid == first_reset.history_valid);
  CHECK(second_reset.sample_count == first_reset.sample_count);
  CHECK(bitwiseEqual(second_reset.total_weight, first_reset.total_weight));

  hipdlp::algorithm::AveragingInput reused = makeValidInput();
  reused.previous = first_reset;
  const hipdlp::algorithm::AveragingResult after_reset =
      hipdlp::algorithm::updateAveraging(reused);
  REQUIRE(after_reset.ok());
  CHECK(after_reset.next.sample_count == 1);
  CHECK(after_reset.sample_weight == populated.sample_weight);
  CHECK(after_reset.next.total_weight == populated.next.total_weight);
  CHECK(after_reset.normalization == populated.normalization);
}

TEST_CASE("hi-pdlp-averaging-policy-preserves-tiny-valid-weights",
          "[pdlp][hipdlp][phase7b][averaging][policy][normalization]"
          "[underflow]") {
  hipdlp::algorithm::AveragingInput input = makeValidInput();
  input.primal_step_size = std::numeric_limits<double>::denorm_min();
  input.dual_step_size = 1.0;
  const hipdlp::algorithm::AveragingResult result =
      hipdlp::algorithm::updateAveraging(input);
  INFO(result.error.message);
  REQUIRE(result.ok());
  CHECK(result.sample_weight ==
        std::sqrt(std::numeric_limits<double>::denorm_min()));
  CHECK(result.next.sample_count == 1);
  CHECK(result.next.total_weight == result.sample_weight);
  CHECK(result.normalization == 1.0 / result.sample_weight);
  CHECK(result.normalization != 1.0);
  CHECK(result.next.hasAverage());
}

TEST_CASE("hi-pdlp-averaging-policy-rejects-invalid-step-sizes",
          "[pdlp][hipdlp][phase7b][averaging][policy][errors]") {
  typedef double hipdlp::algorithm::AveragingInput::* StepField;
  const StepField fields[] = {
      &hipdlp::algorithm::AveragingInput::primal_step_size,
      &hipdlp::algorithm::AveragingInput::dual_step_size};
  const double invalid_values[] = {0.0, -1.0,
                                   std::numeric_limits<double>::quiet_NaN(),
                                   std::numeric_limits<double>::infinity(),
                                   -std::numeric_limits<double>::infinity()};
  for (const StepField field : fields) {
    for (const double value : invalid_values) {
      hipdlp::algorithm::AveragingInput input = makeValidInput();
      input.*field = value;
      checkInvalidResult(hipdlp::algorithm::updateAveraging(input));
    }
  }
}

TEST_CASE("hi-pdlp-averaging-policy-rejects-invalid-prior-state",
          "[pdlp][hipdlp][phase7b][averaging][policy][errors]"
          "[transactional]") {
  hipdlp::algorithm::AveragingInput input = makeValidInput();
  input.previous.history_valid = false;
  checkInvalidResult(hipdlp::algorithm::updateAveraging(input));

  input = makeValidInput();
  input.previous.sample_count = -1;
  checkInvalidResult(hipdlp::algorithm::updateAveraging(input));

  input = makeValidInput();
  input.previous.total_weight = 1.0;
  checkInvalidResult(hipdlp::algorithm::updateAveraging(input));

  input = makeValidInput();
  input.previous.sample_count = 1;
  checkInvalidResult(hipdlp::algorithm::updateAveraging(input));

  const double invalid_totals[] = {-1.0,
                                   std::numeric_limits<double>::quiet_NaN(),
                                   std::numeric_limits<double>::infinity()};
  for (const double total : invalid_totals) {
    input = makeValidInput();
    input.previous.sample_count = 1;
    input.previous.total_weight = total;
    checkInvalidResult(hipdlp::algorithm::updateAveraging(input));
  }

  input = makeValidInput();
  input.previous.total_weight = -0.0;
  checkInvalidResult(hipdlp::algorithm::updateAveraging(input));
}

TEST_CASE("hi-pdlp-averaging-policy-rejects-arithmetic-failures",
          "[pdlp][hipdlp][phase7b][averaging][policy][errors]"
          "[transactional]") {
  SECTION("step-size product underflow") {
    hipdlp::algorithm::AveragingInput input = makeValidInput();
    input.primal_step_size = std::numeric_limits<double>::denorm_min();
    input.dual_step_size = 0.5;
    checkInvalidResult(hipdlp::algorithm::updateAveraging(input));
  }

  SECTION("step-size product overflow") {
    hipdlp::algorithm::AveragingInput input = makeValidInput();
    input.primal_step_size = std::numeric_limits<double>::max();
    input.dual_step_size = 2.0;
    checkInvalidResult(hipdlp::algorithm::updateAveraging(input));
  }

  SECTION("accumulated-weight overflow") {
    hipdlp::algorithm::AveragingInput input = makeValidInput();
    input.previous.sample_count = 1;
    input.previous.total_weight = std::numeric_limits<double>::max();
    input.primal_step_size = 1.0;
    input.dual_step_size = 1.0;
    checkInvalidResult(hipdlp::algorithm::updateAveraging(input));
    CHECK(input.previous.sample_count == 1);
    CHECK(input.previous.total_weight == std::numeric_limits<double>::max());
  }

  SECTION("sample-count overflow") {
    hipdlp::algorithm::AveragingInput input = makeValidInput();
    input.previous.sample_count = std::numeric_limits<HighsInt>::max();
    input.previous.total_weight = 1.0;
    checkInvalidResult(hipdlp::algorithm::updateAveraging(input));
    CHECK(input.previous.sample_count == std::numeric_limits<HighsInt>::max());
    CHECK(input.previous.total_weight == 1.0);
  }
}

TEST_CASE("hi-pdlp-averaging-policy-is-bitwise-deterministic",
          "[pdlp][hipdlp][phase7b][averaging][policy][determinism]") {
  hipdlp::algorithm::AveragingInput input = makeValidInput();
  input.previous.sample_count = 3;
  input.previous.total_weight = 7.25;
  input.primal_step_size = 0.125;
  input.dual_step_size = 32.0;

  const hipdlp::algorithm::AveragingResult first =
      hipdlp::algorithm::updateAveraging(input);
  const hipdlp::algorithm::AveragingResult second =
      hipdlp::algorithm::updateAveraging(input);
  REQUIRE(first.ok());
  REQUIRE(second.ok());
  CHECK(first.error.code == second.error.code);
  CHECK(first.error.message == second.error.message);
  CHECK(first.next.history_valid == second.next.history_valid);
  CHECK(first.next.sample_count == second.next.sample_count);
  CHECK(bitwiseEqual(first.sample_weight, second.sample_weight));
  CHECK(bitwiseEqual(first.next.total_weight, second.next.total_weight));
  CHECK(bitwiseEqual(first.normalization, second.normalization));
}

TEST_CASE("hi-pdlp-averaging-policy-has-no-resource-or-decision-surface",
          "[pdlp][hipdlp][phase7b][averaging][policy][structure]") {
  const char* policy_files[] = {
      "highs/pdlp/hipdlp/algorithm/averaging_policy.hpp",
      "highs/pdlp/hipdlp/algorithm/averaging_policy.cc"};
  const char* forbidden_terms[] = {"backend/",
                                   "backend::",
                                   "ProblemHandle",
                                   "VectorHandle",
                                   "WorkspaceHandle",
                                   "CommandBufferHandle",
                                   "DeviceInfo",
                                   "std::vector",
                                   "cuda",
                                   "CUDA",
                                   "rocm",
                                   "ROCm",
                                   "#include <hip/",
                                   "hip_runtime",
                                   "hipblas",
                                   "hipsparse",
                                   "synchronize",
                                   "malloc(",
                                   "calloc(",
                                   "realloc(",
                                   "new ",
                                   "fusedAverageUpdate",
                                   "CommandBufferRecorder",
                                   "backend_factory",
                                   "evaluateStopping",
                                   "convergence",
                                   "verifyKkt",
                                   "kkt_verifier",
                                   "restart",
                                   "Restart"};
  for (const char* policy_file : policy_files) {
    CAPTURE(policy_file);
    const std::string contents = readSourceFile(policy_file);
    for (const char* forbidden_term : forbidden_terms) {
      CAPTURE(forbidden_term);
      CHECK(contents.find(forbidden_term) == std::string::npos);
    }
  }

  const hipdlp::algorithm::AveragingResult result =
      hipdlp::algorithm::updateAveraging(makeValidInput());
  REQUIRE(result.ok());
  CHECK(result.sample_weight == 6.0);
}