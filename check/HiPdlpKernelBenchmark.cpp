#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "HConfig.h"
#include "Highs.h"
#include "pdlp/hipdlp/algorithm/convergence.hpp"
#include "pdlp/hipdlp/algorithm/fixed_point_error.hpp"
#include "pdlp/hipdlp/algorithm/halpern_pdhg_step.hpp"
#include "pdlp/hipdlp/algorithm/operator_norm.hpp"
#include "pdlp/hipdlp/algorithm/step_size_policy.hpp"
#include "pdlp/hipdlp/backend/backend_factory.hpp"
#include "pdlp/hipdlp/model/scaling.hpp"
#include "pdlp/hipdlp/verification/certificate_verifier.hpp"
#include "pdlp/hipdlp/verification/kkt_verifier.hpp"

namespace {

namespace algorithm = hipdlp::algorithm;
namespace backend = hipdlp::backend;
namespace model = hipdlp::model;
namespace verification = hipdlp::verification;
using Clock = std::chrono::steady_clock;
using hipdlp::api::Error;

const HighsInt kSpmvIterations = 400;
const double kComparisonTolerance = 1e-9;
const double kVerificationTolerance = 1e-4;

#ifdef HIPDLP_KERNEL_LEGACY_CHECK
const bool kLegacyCheck = true;
#else
const bool kLegacyCheck = false;
#endif

void require(const Error& error, const char* operation) {
  if (!error.ok())
    throw std::runtime_error(std::string(operation) + ": " + error.message);
}

void expect(const bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

double secondsSince(const Clock::time_point start) {
  const double seconds =
      std::chrono::duration<double>(Clock::now() - start).count();
  expect(std::isfinite(seconds) && seconds >= 0.0,
         "invalid wall-clock duration");
  return seconds;
}

std::string quoted(const std::string& text) {
  std::string result = "\"";
  for (const unsigned char character : text) {
    if (character == '"' || character == '\\') {
      result += '\\';
      result += static_cast<char>(character);
    } else if (character < 0x20) {
      char encoded[7];
      std::snprintf(encoded, sizeof(encoded), "\\u%04x", unsigned(character));
      result += encoded;
    } else {
      result += static_cast<char>(character);
    }
  }
  return result + '"';
}

class JsonLine {
 public:
  JsonLine() {
    stream_.imbue(std::locale::classic());
    stream_ << std::setprecision(17) << std::boolalpha << '{';
  }

  template <typename Value>
  void add(const char* key, const Value value) {
    prefix(key);
    stream_ << value;
  }

  void add(const char* key, const double value) {
    expect(std::isfinite(value), std::string("non-finite JSON field: ") + key);
    prefix(key);
    stream_ << value;
  }

  void add(const char* key, const std::string& value) {
    prefix(key);
    stream_ << quoted(value);
  }

  void add(const char* key, const char* value) { add(key, std::string(value)); }

  void null(const char* key) {
    prefix(key);
    stream_ << "null";
  }

  std::string finish() { return stream_.str() + '}'; }

 private:
  void prefix(const char* key) {
    if (!first_) stream_ << ',';
    first_ = false;
    stream_ << quoted(key) << ':';
  }

  bool first_ = true;
  std::ostringstream stream_;
};

class NewFile {
 public:
  explicit NewFile(const std::string& filename)
      : file_(std::fopen(filename.c_str(), "wx")) {
    expect(file_ != nullptr, "cannot create new output file: " + filename);
  }
  ~NewFile() {
    if (file_ != nullptr) std::fclose(file_);
  }
  NewFile(const NewFile&) = delete;
  NewFile& operator=(const NewFile&) = delete;

  void write(const std::string& text) {
    expect(std::fwrite(text.data(), 1, text.size(), file_) == text.size(),
           "failed to write benchmark output");
  }
  FILE* get() const { return file_; }
  void close() {
    const bool failed = std::ferror(file_) != 0;
    const int result = std::fclose(file_);
    file_ = nullptr;
    expect(!failed && result == 0, "failed to finish benchmark output");
  }

 private:
  FILE* file_;
};

HighsInt positiveInteger(const std::string& text) {
  expect(!text.empty() && std::all_of(text.begin(), text.end(),
                                      [](char character) {
                                        return character >= '0' &&
                                               character <= '9';
                                      }),
         "invalid positive integer: " + text);
  const unsigned long long value = std::stoull(text);
  expect(value > 0 && value <= std::numeric_limits<HighsInt>::max(),
         "positive integer is out of range: " + text);
  return static_cast<HighsInt>(value);
}

bool boolean(const std::string& text) {
  expect(text == "true" || text == "false", "expected true or false: " + text);
  return text == "true";
}

struct Options {
  std::string instance;
  std::string output_directory;
  std::string reference_state;
  backend::Kind backend_kind;
  HighsInt warmups;
  HighsInt repetitions;
  bool pid_distances;
  bool certificates;
};

Options parseOptions(const int argc, char** argv) {
  std::map<std::string, std::string> arguments = {
      {"--instance", ""},          {"--backend", "cpu"},
      {"--output-dir", ""},        {"--reference-state", ""},
      {"--warmups", "1"},          {"--repetitions", "5"},
      {"--pid-distances", "true"}, {"--certificates", "true"}};
  std::set<std::string> seen;
  for (int index = 1; index < argc; index += 2) {
    const std::string flag = argv[index];
    expect(
        index + 1 < argc && arguments.count(flag) && seen.insert(flag).second,
        "unknown, duplicate or valueless option: " + flag);
    arguments[flag] = argv[index + 1];
  }
  Options options;
  options.instance = arguments.at("--instance");
  options.output_directory = arguments.at("--output-dir");
  options.reference_state = arguments.at("--reference-state");
  expect(
      !options.instance.empty() && !options.output_directory.empty(),
      "--instance and --output-dir are required; output directory must exist");
  const std::string kind = arguments.at("--backend");
  expect(kind == "cpu" || kind == "cuda", "--backend must be cpu or cuda");
#ifndef HIGHS_HIPDLP_CUDA_DETECTED
  expect(kind == "cpu", "this benchmark was compiled without CUDA");
#endif
  options.backend_kind =
      kind == "cpu" ? backend::Kind::kCpu : backend::Kind::kCuda;
  options.warmups = positiveInteger(arguments.at("--warmups"));
  options.repetitions = positiveInteger(arguments.at("--repetitions"));
  expect(options.repetitions >= 5 &&
             options.warmups <=
                 std::numeric_limits<HighsInt>::max() - options.repetitions,
         "at least one warmup and five measurements are required without "
         "overflow");
  options.pid_distances = boolean(arguments.at("--pid-distances"));
  options.certificates = boolean(arguments.at("--certificates"));
  return options;
}

struct Measurement {
  double wall = 0.0;
  double backend_seconds = 0.0;
  bool backend_timed = false;
  bool events = false;
};

template <typename Work>
Measurement measureHost(Work work) {
  const Clock::time_point start = Clock::now();
  work();
  Measurement result;
  result.wall = secondsSince(start);
  return result;
}

template <typename Work>
Measurement measureBackend(backend::Backend& selected, Work work) {
  require(selected.synchronize(), "timing start boundary");
  const Clock::time_point start = Clock::now();
#ifndef HIPDLP_KERNEL_LEGACY_CHECK
  require(selected.beginTiming(), "begin timing");
#endif
  work();
  Measurement result;
#ifndef HIPDLP_KERNEL_LEGACY_CHECK
  require(selected.endTiming(result.backend_seconds), "end timing");
  result.backend_timed = true;
  result.events = selected.timingUsesEvents();
#endif
  require(selected.synchronize(), "timing completion and finite check");
  result.wall = secondsSince(start);
  return result;
}

struct Vectors {
  backend::VectorHandle current_primal, next_primal, primal_anchor;
  backend::VectorHandle primal_candidate, reflected_primal, current_aty;
  backend::VectorHandle candidate_aty, aty_delta;
  backend::VectorHandle current_dual, next_dual, dual_anchor, dual_candidate;
  backend::VectorHandle reflected_dual, current_ax, reflected_ax, candidate_ax;

  void allocate(backend::Backend& selected, const model::Problem& problem) {
    backend::VectorHandle* columns[] = {
        &current_primal,   &next_primal, &primal_anchor, &primal_candidate,
        &reflected_primal, &current_aty, &candidate_aty, &aty_delta};
    backend::VectorHandle* rows[] = {
        &current_dual,   &next_dual,  &dual_anchor,  &dual_candidate,
        &reflected_dual, &current_ax, &reflected_ax, &candidate_ax};
    for (std::size_t index = 0; index < 8; ++index) {
      require(selected.createVector(problem.numCols(), *columns[index]),
              "column vector");
      require(selected.createVector(problem.numRows(), *rows[index]),
              "row vector");
    }
  }
};

struct CheckObservation {
  algorithm::RawConvergenceMetrics raw;
  double fixed_point_error = 0.0;
  double primal_distance = 0.0;
  double dual_distance = 0.0;
  double raw_oracle_error = 0.0;
  bool fixed_point = false;
  bool distances = false;
  bool packet_valid = false;
  bool host_fallback = false;
  bool primal_ray_possible = true;
  bool farkas_possible = true;
  bool kkt_pass = false;
  bool primal_ray_pass = false;
  bool farkas_pass = false;
  std::size_t packet_transfers = 0;
  std::size_t packet_bytes = 0;
};

double relativeError(const double left, const double right) {
  expect(std::isfinite(left) && std::isfinite(right),
         "non-finite comparison input");
  const double scale = std::max(1.0, std::max(std::abs(left), std::abs(right)));
  return std::abs(left / scale - right / scale);
}

double compareRaw(const algorithm::RawConvergenceMetrics& actual,
                  const algorithm::RawConvergenceMetrics& expected) {
  const double observed[] = {actual.primal_residual,   actual.dual_residual,
                             actual.objective_gap,     actual.primal_objective,
                             actual.dual_objective,    actual.unscaled_rhs_norm,
                             actual.unscaled_cost_norm};
  const double reference[] = {
      expected.primal_residual,   expected.dual_residual,
      expected.objective_gap,     expected.primal_objective,
      expected.dual_objective,    expected.unscaled_rhs_norm,
      expected.unscaled_cost_norm};
  double maximum = 0.0;
  for (std::size_t index = 0; index < 7; ++index)
    maximum =
        std::max(maximum, relativeError(observed[index], reference[index]));
  expect(maximum <= kComparisonTolerance,
         "raw check disagrees with serial host oracle");
  return maximum;
}

double norm(const std::vector<double>& values) {
  double result = 0.0;
  for (const double value : values) result = std::hypot(result, value);
  return result;
}

bool normalize(std::vector<double>& values, const double length) {
  if (!std::isfinite(length) || length <= 0.0) return false;
  for (double& value : values) {
    const double normalized = value / length;
    if (!std::isfinite(normalized) || (value != 0.0 && normalized == 0.0))
      return false;
    value = normalized;
  }
  return true;
}

class Workload {
 public:
  Workload(const model::Problem& problem, const HighsLp& lp)
      : problem(problem), lp(lp) {}

  void prepare(const backend::Kind kind) {
    backend::FactoryResult factory = backend::createBackend(kind);
    require(factory.error, "backend creation");
    expect(factory.backend != nullptr, "backend factory returned no backend");
    selected = std::move(factory.backend);
    require(selected->bindProblem(problem, problem_handle),
            "one-time problem upload");
    require(selected->createWorkspace(problem_handle, workspace),
            "workspace creation");
    vectors.allocate(*selected, problem);
    host[0].resize(problem.numCols());
    host[1].resize(problem.numRows());
    host[2].resize(problem.numRows());
    host[3].resize(problem.numCols());
    reduced.resize(problem.numCols());
    primal_direction.resize(lp.num_col_);
    dual_direction.resize(lp.num_row_);
    column_direction.resize(lp.num_col_);
#ifndef HIPDLP_KERNEL_LEGACY_CHECK
    require(selected->prepareCheckWorkspace(problem, lp, workspace),
            "check workspace preparation");
    require(selected->prepareTiming(), "optional timer preparation");
#endif
    reset();
#ifndef HIPDLP_KERNEL_LEGACY_CHECK
    const backend::VectorDownload candidate[] = {
        {state.primal_candidate, &host[0]},
        {state.dual_candidate, &host[1]},
        {&vectors.candidate_ax, &host[2]},
        {&vectors.candidate_aty, &host[3]}};
    const backend::VectorDownload current[] = {{state.current_primal, &host[0]},
                                               {state.current_dual, &host[1]},
                                               {&vectors.current_ax, &host[2]},
                                               {state.current_aty, &host[3]}};
    require(selected->prepareVectorDownloads(candidate, 4),
            "candidate download preparation");
    require(selected->prepareVectorDownloads(current, 4),
            "final-state download preparation");
#endif
    require(selected->synchronize(), "initialization and upload completion");
  }

  void initializeSteps() {
    estimate = algorithm::estimateOperatorNormSquared(
        *selected, problem_handle, problem.numRows(), problem.numCols());
    require(estimate.error, "operator-norm estimation");
    expect(!estimate.empty_matrix &&
               std::isfinite(estimate.operator_norm_squared) &&
               estimate.operator_norm_squared > 0.0,
           "nonpositive operator norm");
    const std::vector<model::ScalingStage>& stages = problem.scaling().stages;
    if (!stages.empty() &&
        stages.back() == model::ScalingStage::kPockChambolle) {
      upper_bound_squared = 1.0;
    } else {
      require(algorithm::computeOperatorNormSquaredUpperBound(
                  problem.matrix(), upper_bound_squared),
              "operator-norm upper bound");
    }
    algorithm::StepSizeInitializationInput input;
    input.operator_norm_squared =
        std::max(estimate.operator_norm_squared, upper_bound_squared);
    input.unscaled_rhs_norm = problem.unscaledNorms().row_bounds_l2;
    input.unscaled_cost_norm = problem.unscaledNorms().objective_l2;
    const algorithm::StepSizeInitializationResult initialized =
        algorithm::initializeStepSizes(input);
    require(initialized.error, "safe fixed-step initialization");
    steps = initialized.step_sizes;
    expect(steps.base_step * steps.base_step * upper_bound_squared < 1.0,
           "initialized steps do not satisfy the operator-norm safety bound");
    require(selected->synchronize(), "operator-norm completion");
  }

  void reset() {
    require(selected->synchronize(), "pre-reset finite check");
    state.problem = &problem_handle;
    state.workspace = &workspace;
    state.current_primal = &vectors.current_primal;
    state.current_dual = &vectors.current_dual;
    state.next_primal = &vectors.next_primal;
    state.next_dual = &vectors.next_dual;
    state.primal_anchor = &vectors.primal_anchor;
    state.dual_anchor = &vectors.dual_anchor;
    state.primal_candidate = &vectors.primal_candidate;
    state.dual_candidate = &vectors.dual_candidate;
    state.reflected_primal = &vectors.reflected_primal;
    state.reflected_dual = &vectors.reflected_dual;
#ifdef HIPDLP_KERNEL_LEGACY_CHECK
    state.current_ax = &vectors.current_ax;
#endif
    state.current_aty = &vectors.current_aty;
    state.reflected_ax = &vectors.reflected_ax;
    state.primal_step_size = steps.primal_step;
    state.dual_step_size = steps.dual_step;
    state.iteration = 0;
    require(selected->initialize(*state.current_primal, *state.current_dual),
            "zero start");
    require(selected->fill(*state.current_aty, 0.0), "zero transpose product");
    require(selected->fusedProjectedPrimalUpdate(
                problem_handle, workspace, *state.current_primal,
                *state.current_aty, 0.0, *state.next_primal),
            "start projection");
    std::swap(state.current_primal, state.next_primal);
    require(selected->copy(*state.current_primal, *state.primal_anchor),
            "primal anchor reset");
    require(selected->copy(*state.current_dual, *state.dual_anchor),
            "dual anchor reset");
    require(selected->multiply(problem_handle, *state.current_primal,
                               vectors.current_ax),
            "initial Ax");
    require(selected->multiplyTranspose(problem_handle, *state.current_dual,
                                        *state.current_aty),
            "initial Aty");
    require(selected->synchronize(), "reset completion");
  }

  void enqueueSteps(const HighsInt count) {
    for (HighsInt iteration = 0; iteration < count; ++iteration)
      require(algorithm::enqueueHalpernPdhgStep(*selected, state),
              "shared Halpern-PDHG step");
    expect(state.iteration == count,
           "incorrect fixed-workload iteration count");
  }

  void download(const bool candidate) {
#ifndef HIPDLP_KERNEL_LEGACY_CHECK
    if (!candidate)
      require(selected->multiply(problem_handle, *state.current_primal,
                                 vectors.current_ax),
              "snapshot Ax");
#endif
    const backend::VectorHandle* handles[] = {
        candidate ? state.primal_candidate : state.current_primal,
        candidate ? state.dual_candidate : state.current_dual,
        candidate ? &vectors.candidate_ax : &vectors.current_ax,
        candidate ? &vectors.candidate_aty : state.current_aty};
#ifdef HIPDLP_KERNEL_LEGACY_CHECK
    for (std::size_t index = 0; index < 4; ++index)
      require(selected->downloadVector(*handles[index], host[index]),
              "legacy four-vector download");
#else
    const backend::VectorDownload downloads[] = {{handles[0], &host[0]},
                                                 {handles[1], &host[1]},
                                                 {handles[2], &host[2]},
                                                 {handles[3], &host[3]}};
    require(selected->downloadVectors(downloads, 4),
            "prepared four-vector download");
#endif
  }

  algorithm::RawConvergenceMetrics rawOracle() const {
    algorithm::RawConvergenceMetrics raw;
    require(algorithm::computeRawConvergenceMetrics(problem, host[0], host[1],
                                                    host[2], host[3], raw),
            "serial raw CPU oracle");
    return raw;
  }

  algorithm::FixedPointErrorInput fixedPointInput() {
    algorithm::FixedPointErrorInput input;
    input.next_primal = state.primal_candidate;
    input.next_dual = state.dual_candidate;
    input.reflected_primal = state.reflected_primal;
    input.reflected_dual = state.reflected_dual;
    input.delta_primal = state.next_primal;
    input.delta_dual = state.next_dual;
    input.aty_delta = &vectors.aty_delta;
    input.omega = steps.primal_weight;
    input.eta = steps.base_step;
    return input;
  }

  void legacyCheck(CheckObservation& observation) {
    download(true);
    observation.raw = rawOracle();
    if (observation.fixed_point) {
      const algorithm::FixedPointErrorResult result =
          algorithm::evaluateFixedPointError(*selected, problem_handle,
                                             fixedPointInput());
      require(result.error, "legacy FPE evaluation");
      expect(result.ok(), "legacy FPE was not computed");
      observation.fixed_point_error = result.fixed_point_error;
    }
    if (observation.distances) {
      require(selected->differenceNorm2(*state.primal_candidate,
                                        *state.primal_anchor,
                                        observation.primal_distance),
              "primal anchor distance");
      require(
          selected->differenceNorm2(*state.dual_candidate, *state.dual_anchor,
                                    observation.dual_distance),
          "dual anchor distance");
    }
  }

  CheckObservation check(const bool fixed_point, const bool distances) {
    require(selected->multiply(problem_handle, *state.primal_candidate,
                               vectors.candidate_ax),
            "candidate Ax");
    require(selected->multiplyTranspose(problem_handle, *state.dual_candidate,
                                        vectors.candidate_aty),
            "candidate Aty");
    CheckObservation observation;
    observation.fixed_point = fixed_point;
    observation.distances = distances;
#ifdef HIPDLP_KERNEL_LEGACY_CHECK
    legacyCheck(observation);
#else
    backend::CheckInput input;
    input.primal = state.primal_candidate;
    input.dual = state.dual_candidate;
    input.ax = &vectors.candidate_ax;
    input.aty = &vectors.candidate_aty;
    input.fixed_point = fixed_point;
    input.anchor_distances = distances;
    input.primal_anchor = state.primal_anchor;
    input.dual_anchor = state.dual_anchor;
    input.certificate_tolerance = kVerificationTolerance;
    if (fixed_point) {
      const algorithm::FixedPointErrorInput fixed = fixedPointInput();
      require(
          selected->prepareCheckFixedPoint(
              problem_handle, workspace, *fixed.next_primal, *fixed.next_dual,
              *fixed.reflected_primal, *fixed.reflected_dual,
              *fixed.delta_primal, *fixed.delta_dual, *fixed.aty_delta),
          "resident FPE preparation");
      input.delta_primal = fixed.delta_primal;
      input.delta_dual = fixed.delta_dual;
      input.aty_delta = fixed.aty_delta;
    }
    backend::CheckResult packet;
    require(selected->evaluateCheck(problem, workspace, input, packet),
            "resident check packet");
    expect(packet.valid || packet.requires_host_oracle,
           "check packet supplied neither metrics nor fallback");
    observation.packet_valid = packet.valid;
    observation.host_fallback = packet.requires_host_oracle;
    observation.primal_ray_possible = packet.primal_ray_possible;
    observation.farkas_possible = packet.farkas_possible;
    observation.packet_transfers = packet.device_to_host_transfers;
    observation.packet_bytes = packet.device_to_host_bytes;
    if (packet.requires_host_oracle) {
      legacyCheck(observation);
    } else {
      observation.raw = packet.raw;
      observation.primal_distance = packet.primal_anchor_distance;
      observation.dual_distance = packet.dual_anchor_distance;
      if (fixed_point) {
        const algorithm::FixedPointErrorResult result =
            algorithm::computeFixedPointError(
                packet.primal_delta_norm, packet.dual_delta_norm,
                packet.delta_dot, steps.primal_weight, steps.base_step);
        require(result.error, "packet FPE scalar evaluation");
        expect(result.ok(), "packet FPE was not computed");
        observation.fixed_point_error = result.fixed_point_error;
      }
    }
#endif
    return observation;
  }

  void verifyKkt(CheckObservation& observation) {
    for (std::size_t index = 0; index < reduced.size(); ++index)
      reduced[index] = problem.costs()[index] - host[3][index];
    model::PostsolveResult recovered =
        model::postsolve(problem, host[0], host[1], reduced);
    require(recovered.error, "host postsolve");
    const verification::KktVerificationResult verified =
        verification::verifyKkt(lp, recovered.solution, kVerificationTolerance);
    expect(verified.input_valid && verified.metrics_finite,
           "invalid host KKT verification");
    observation.kkt_pass = verified.pass;
  }

  void verifyCertificates(CheckObservation& observation) {
    if (problem.recoverPrimalDirection(host[0], primal_direction).ok() &&
        normalize(primal_direction, norm(primal_direction)))
      observation.primal_ray_pass =
          verification::verifyPrimalRay(lp, primal_direction,
                                        kVerificationTolerance)
              .pass;
    for (std::size_t index = 0; index < reduced.size(); ++index)
      reduced[index] = -host[3][index];
    if (problem.recoverDualDirection(host[1], dual_direction).ok() &&
        problem.recoverReducedCostDirection(reduced, column_direction).ok()) {
      const double length =
          std::hypot(norm(dual_direction), norm(column_direction));
      if (normalize(dual_direction, length) &&
          normalize(column_direction, length))
        observation.farkas_pass =
            verification::verifyFarkasCertificate(
                lp, dual_direction, column_direction, kVerificationTolerance)
                .pass;
    }
    expect(!observation.primal_ray_pass || observation.primal_ray_possible,
           "certificate screen rejected a host-verified primal ray");
    expect(!observation.farkas_pass || observation.farkas_possible,
           "certificate screen rejected a host-verified Farkas certificate");
  }

  const model::Problem& problem;
  const HighsLp& lp;
  std::unique_ptr<backend::Backend> selected;
  backend::ProblemHandle problem_handle;
  backend::WorkspaceHandle workspace;
  Vectors vectors;
  algorithm::HalpernPdhgState state;
  algorithm::InitialStepSizes steps;
  algorithm::OperatorNormResult estimate;
  double upper_bound_squared = 0.0;
  std::vector<double> host[4];
  std::vector<double> reduced, primal_direction, dual_direction,
      column_direction;
};

class Report {
 public:
  explicit Report(const Options& options)
      : options_(options), file_(options.output_directory + "/samples.jsonl") {}

  void common(JsonLine& row) const {
    row.add("schema", "hipdlp-kernel-benchmark/1");
    row.add("instance", options_.instance);
    row.add("backend",
            options_.backend_kind == backend::Kind::kCpu ? "cpu" : "cuda");
    row.add("legacy_check", kLegacyCheck);
  }

  void emit(JsonLine& row) {
    const std::string line = row.finish() + '\n';
    file_.write(line);
    std::cout << line;
    expect(std::cout.good(), "failed to write benchmark stdout");
  }

  void sample(const char* field, const Measurement& measured,
              const HighsInt sequence, const HighsInt steps,
              const HighsInt work_items, const char* check_kind = "none",
              const CheckObservation* check = nullptr) {
    JsonLine row;
    common(row);
    const bool warmup = sequence >= 0 && sequence < options_.warmups;
    row.add("event", "sample");
    row.add("field", field);
    row.add("kind", sequence < 0 ? "setup" : warmup ? "warmup" : "measurement");
    row.add("index", sequence < 0 ? 0
                     : warmup     ? sequence
                                  : sequence - options_.warmups);
    row.add("sequence", sequence);
    row.add("steps", steps);
    row.add("samples", 1);
    row.add("work_items", work_items);
    row.add("check_kind", check_kind);
    row.add("wall_seconds", measured.wall);
    if (measured.backend_timed)
      row.add("backend_seconds", measured.backend_seconds);
    else
      row.null("backend_seconds");
    if (measured.events)
      row.add("event_seconds", measured.backend_seconds);
    else
      row.null("event_seconds");
    if (check != nullptr) {
      row.add("packet_valid", check->packet_valid);
      row.add("host_fallback", check->host_fallback);
      row.add("packet_d2h_transfers", check->packet_transfers);
      row.add("packet_d2h_bytes", check->packet_bytes);
      row.add("raw_oracle_matches", true);
      row.add("raw_oracle_max_relative_error", check->raw_oracle_error);
      row.add("primal_residual", check->raw.primal_residual);
      row.add("dual_residual", check->raw.dual_residual);
      row.add("objective_gap", check->raw.objective_gap);
      row.add("primal_objective", check->raw.primal_objective);
      row.add("dual_objective", check->raw.dual_objective);
      if (check->fixed_point)
        row.add("fixed_point_error", check->fixed_point_error);
      else
        row.null("fixed_point_error");
      if (check->distances) {
        row.add("primal_anchor_distance", check->primal_distance);
        row.add("dual_anchor_distance", check->dual_distance);
      }
      row.add("kkt_pass", check->kkt_pass);
      row.add("certificate_checked", options_.certificates);
      row.add("primal_ray_possible", check->primal_ray_possible);
      row.add("farkas_possible", check->farkas_possible);
      row.add("primal_ray_pass", check->primal_ray_pass);
      row.add("farkas_pass", check->farkas_pass);
    }
    emit(row);
  }

  void close() { file_.close(); }

 private:
  const Options& options_;
  NewFile file_;
};

const char* const kStateLabels[] = {"primal", "dual", "ax", "aty"};

void writeState(const std::string& filename, const Workload& workload) {
  NewFile output(filename);
  FILE* file = output.get();
  std::fprintf(file, "HIPDLP_KERNEL_STATE 1 %s\n400 %lld %lld\n%.17g %.17g\n",
               workload.selected->descriptor().name,
               static_cast<long long>(workload.problem.numRows()),
               static_cast<long long>(workload.problem.numCols()),
               workload.steps.primal_step, workload.steps.dual_step);
  for (std::size_t index = 0; index < 4; ++index) {
    const std::vector<double>& values = workload.host[index];
    std::fprintf(file, "%s %llu\n", kStateLabels[index],
                 static_cast<unsigned long long>(values.size()));
    for (const double value : values) {
      expect(std::isfinite(value), "non-finite final-state value");
      std::fprintf(file, "%.17g\n", value);
    }
  }
  output.close();
}

double compareState(const std::string& filename, const Workload& workload,
                    const bool require_cpu, const double tolerance) {
  std::ifstream input(filename);
  input.imbue(std::locale::classic());
  std::string magic, kind;
  int version = 0;
  HighsInt iterations = 0, rows = 0, columns = 0;
  double primal_step = 0.0, dual_step = 0.0;
  input >> magic >> version >> kind >> iterations >> rows >> columns >>
      primal_step >> dual_step;
  expect(input.good() && magic == "HIPDLP_KERNEL_STATE" && version == 1 &&
             iterations == 400 && rows == workload.problem.numRows() &&
             columns == workload.problem.numCols() &&
             (!require_cpu || kind == "cpu"),
         "invalid or non-CPU reference state: " + filename);
  double maximum =
      std::max(relativeError(primal_step, workload.steps.primal_step),
               relativeError(dual_step, workload.steps.dual_step));
  for (std::size_t index = 0; index < 4; ++index) {
    std::string label;
    std::size_t size = 0;
    input >> label >> size;
    expect(input.good() && label == kStateLabels[index] &&
               size == workload.host[index].size(),
           "reference-state vector dimensions do not match");
    for (const double actual : workload.host[index]) {
      double expected = 0.0;
      expect(static_cast<bool>(input >> expected),
             "truncated reference-state vector");
      maximum = std::max(maximum, relativeError(actual, expected));
    }
  }
  std::string trailing;
  expect(!(input >> trailing) && input.eof(),
         "unexpected reference-state trailing data");
  expect(maximum <= tolerance,
         "fixed 400-step state disagrees with reference: " + filename);
  return maximum;
}

void benchmarkCheck(Workload& workload, Report& report, const Options& options,
                    const HighsInt sequence, const HighsInt steps,
                    const bool fixed_point, const bool distances) {
  const char* kind = distances     ? "raw_fpe_pid_distances"
                     : fixed_point ? "raw_fpe"
                                   : "raw";
  CheckObservation observation;
  const Measurement boundary = measureBackend(*workload.selected, [&]() {
    observation = workload.check(fixed_point, distances);
  });
  const Measurement download =
      measureBackend(*workload.selected, [&]() { workload.download(true); });
  algorithm::RawConvergenceMetrics oracle;
  const Measurement raw = measureHost([&]() { oracle = workload.rawOracle(); });
  observation.raw_oracle_error = compareRaw(observation.raw, oracle);
  const Measurement kkt =
      measureHost([&]() { workload.verifyKkt(observation); });
  Measurement certificate;
  if (options.certificates)
    certificate =
        measureHost([&]() { workload.verifyCertificates(observation); });
  report.sample("check_boundary", boundary, sequence, steps, 1, kind,
                &observation);
  report.sample("download_candidate_four", download, sequence, steps, 4, kind,
                &observation);
  report.sample("raw_cpu_oracle", raw, sequence, steps, 1, kind, &observation);
  report.sample("host_postsolve_kkt_verify", kkt, sequence, steps, 1, kind,
                &observation);
  if (options.certificates)
    report.sample("host_certificate_verify", certificate, sequence, steps, 1,
                  kind, &observation);
}

int run(const Options& options) {
  Report report(options);
  Highs highs;
  expect(highs.setOptionValue("output_flag", false) == HighsStatus::kOk,
         "output option failed");
  expect(highs.setOptionValue("log_to_console", false) == HighsStatus::kOk,
         "logging option failed");
  expect(highs.setOptionValue("threads", HighsInt(1)) == HighsStatus::kOk,
         "serial option failed");
  const Measurement read = measureHost([&]() {
    expect(highs.readModel(options.instance) == HighsStatus::kOk,
           "failed to read model");
  });
  model::ProblemBuildResult built;
  const Measurement scaling = measureHost([&]() {
    built = model::buildProblem(highs.getLp());
    require(built.error, "model canonicalization");
    model::ScalingOptions scaling_options;
    scaling_options.mode = model::ScalingMode::kAuto;
    scaling_options.curtis_reid_sweeps = 20;
    require(model::applyScalingPipeline(built.problem, scaling_options),
            "auto scaling");
    expect(built.problem.matrix().numNonzeros() > 0,
           "zero-operator workload is unsupported");
  });
  Workload workload(built.problem, highs.getLp());
  const Measurement initialization =
      measureHost([&]() { workload.prepare(options.backend_kind); });
  const Measurement operator_norm =
      measureHost([&]() { workload.initializeSteps(); });
  JsonLine metadata;
  report.common(metadata);
  metadata.add("event", "metadata");
  metadata.add("step_implementation", "enqueueHalpernPdhgStep");
  metadata.add("scaling", "auto");
  metadata.add("curtis_reid_sweeps", 20);
  metadata.add("rows", built.problem.numRows());
  metadata.add("columns", built.problem.numCols());
  metadata.add("nonzeros", built.problem.matrix().numNonzeros());
  metadata.add("steady_state_vectors", 16);
  metadata.add("problem_uploads", 1);
  metadata.add("state_history", false);
  metadata.add("restarts", false);
  metadata.add("pid_updates", false);
  metadata.add("pid_distance_checks", options.pid_distances);
  metadata.add("host_certificate_checks", options.certificates);
  metadata.add("resident_certificate_screen_included",
               !kLegacyCheck && options.backend_kind == backend::Kind::kCuda);
  metadata.add("warmups", options.warmups);
  metadata.add("repetitions", options.repetitions);
  metadata.add("setup_samples", 1);
  metadata.add("spmv_iterations", kSpmvIterations);
  metadata.add("spmv_includes_finite_check", true);
  metadata.null("vendor_spmv_only_event_seconds");
  metadata.add("candidate_products_in_check_boundary", true);
  metadata.add("host_fallback_in_check_boundary", true);
  metadata.add("host_validation_outside_check_boundary", true);
  metadata.add("operator_norm_squared",
               workload.estimate.operator_norm_squared);
  metadata.add("operator_norm_iterations",
               workload.estimate.performed_iterations);
  metadata.add("operator_norm_upper_bound_squared",
               workload.upper_bound_squared);
  metadata.add("primal_step", workload.steps.primal_step);
  metadata.add("dual_step", workload.steps.dual_step);
  metadata.add("base_step", workload.steps.base_step);
  metadata.add("primal_weight", workload.steps.primal_weight);
  metadata.add("comparison_tolerance", kComparisonTolerance);
  metadata.add("verification_tolerance", kVerificationTolerance);
  metadata.add("reference_state", options.reference_state);
  metadata.add("timing", kLegacyCheck ? "synchronized_wall_only"
                         : options.backend_kind == backend::Kind::kCuda
                             ? "cuda_events_and_synchronized_wall"
                             : "cpu_steady_clock");
  report.emit(metadata);
  report.sample("model_read", read, -1, 0, 1);
  report.sample("model_build_scale", scaling, -1, 0, 1);
  report.sample("initialization_upload", initialization, -1, 0, 1);
  report.sample("operator_norm_and_steps", operator_norm, -1, 0, 1);
  const std::string state_file = options.output_directory + "/state-400.txt";
  double reference_error = 0.0;
  for (HighsInt sequence = 0; sequence < options.warmups + options.repetitions;
       ++sequence) {
    workload.reset();
    require(workload.selected->fill(workload.vectors.reflected_primal, 1.0),
            "SpMV primal input");
    require(workload.selected->fill(workload.vectors.reflected_dual, 1.0),
            "SpMV dual input");
    const Measurement ax = measureBackend(*workload.selected, [&]() {
      for (HighsInt iteration = 0; iteration < kSpmvIterations; ++iteration)
        require(workload.selected->multiply(workload.problem_handle,
                                            workload.vectors.reflected_primal,
                                            workload.vectors.candidate_ax),
                "SpMV Ax");
    });
    report.sample("spmv_ax_finite_check", ax, sequence, 0, kSpmvIterations);
    const Measurement aty = measureBackend(*workload.selected, [&]() {
      for (HighsInt iteration = 0; iteration < kSpmvIterations; ++iteration)
        require(workload.selected->multiplyTranspose(
                    workload.problem_handle, workload.vectors.reflected_dual,
                    workload.vectors.candidate_aty),
                "SpMV Aty");
    });
    report.sample("spmv_aty_finite_check", aty, sequence, 0, kSpmvIterations);
    for (const HighsInt count : {HighsInt(40), HighsInt(400)}) {
      workload.reset();
      const Measurement steps = measureBackend(
          *workload.selected, [&]() { workload.enqueueSteps(count); });
      report.sample("halpern_steps", steps, sequence, count, count);
      benchmarkCheck(workload, report, options, sequence, count, false, false);
      benchmarkCheck(workload, report, options, sequence, count, true, false);
      if (options.pid_distances)
        benchmarkCheck(workload, report, options, sequence, count, true, true);
      if (count == 400) {
        const Measurement final_download = measureBackend(
            *workload.selected, [&]() { workload.download(false); });
        report.sample("download_final_current_four", final_download, sequence,
                      count, 4);
        if (sequence == 0)
          writeState(state_file, workload);
        else
          compareState(state_file, workload, false, 0.0);
        if (!options.reference_state.empty())
          reference_error = std::max(
              reference_error, compareState(options.reference_state, workload,
                                            true, kComparisonTolerance));
      }
    }
  }
  JsonLine complete;
  report.common(complete);
  complete.add("event", "complete");
  complete.add("repeat_state_identical", true);
  complete.add("state_file", state_file);
  complete.add("cpu_reference_checked", !options.reference_state.empty());
  if (!options.reference_state.empty())
    complete.add("cpu_reference_max_relative_error", reference_error);
  else
    complete.null("cpu_reference_max_relative_error");
  complete.add("solver_optimality_claimed", false);
  report.emit(complete);
  report.close();
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  static_assert(
      sizeof(double) == 8 && std::numeric_limits<double>::digits == 53,
      "HiPDLP kernel benchmark requires FP64");
  try {
    return run(parseOptions(argc, argv));
  } catch (const std::exception& error) {
    std::cerr << "HiPDLP kernel benchmark: " << error.what() << '\n';
    return 1;
  }
}