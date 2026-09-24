#include <unistd.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

#include "HConfig.h"
#include "Highs.h"
#include "lp_data/HighsModelUtils.h"
#include "pdlp/HiPdlpWrapper.h"

namespace {

using Clock = std::chrono::steady_clock;

#ifdef HIPDLP_PERFORMANCE_LEGACY
const char* const kImplementation = "legacy";
const char* const kScaling = "legacy-native-default";
#ifdef CUPDLP_GPU
const bool kCompiledCuda = true;
const bool kCompiledCpu = false;
#else
const bool kCompiledCuda = false;
const bool kCompiledCpu = true;
#endif
#else
const char* const kImplementation = "rewrite";
const char* const kScaling = "auto";
const bool kCompiledCpu = true;
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
const bool kCompiledCuda = true;
#else
const bool kCompiledCuda = false;
#endif
#endif

double secondsSince(const Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}

HighsInt integer(const std::string& text) {
  std::size_t consumed = 0;
  const long long value = std::stoll(text, &consumed);
  if (consumed != text.size() || value < 0 ||
      value > std::numeric_limits<HighsInt>::max())
    throw std::invalid_argument("invalid nonnegative integer: " + text);
  return static_cast<HighsInt>(value);
}

std::string quoted(const std::string& text) {
  std::string result = "\"";
  for (const unsigned char character : text) {
    if (character == '"' || character == '\\') {
      result += '\\';
      result += static_cast<char>(character);
    } else if (character < 0x20) {
      char encoded[7];
      std::snprintf(encoded, sizeof(encoded), "\\u%04x",
                    static_cast<unsigned>(character));
      result += encoded;
    } else {
      result += static_cast<char>(character);
    }
  }
  return result + '"';
}

FILE* exclusiveFile(const std::string& filename) {
  FILE* file = std::fopen(filename.c_str(), "wx");
  if (!file) throw std::runtime_error("cannot create new file: " + filename);
  return file;
}

void closeFile(FILE* file) {
  const bool failed = std::ferror(file) != 0;
  if (std::fclose(file) != 0 || failed)
    throw std::runtime_error("failed to finish output file");
}

const char* statusName(const HighsModelStatus status) {
  switch (status) {
    case HighsModelStatus::kOptimal:
      return "OPTIMAL";
    case HighsModelStatus::kTimeLimit:
      return "TIME_LIMIT";
    case HighsModelStatus::kIterationLimit:
      return "ITERATION_LIMIT";
    case HighsModelStatus::kInfeasible:
      return "INFEASIBLE";
    case HighsModelStatus::kUnbounded:
      return "UNBOUNDED";
    case HighsModelStatus::kUnboundedOrInfeasible:
      return "UNBOUNDED_OR_INFEASIBLE";
    default:
      return "ERROR";
  }
}

}  // namespace

int main(int argc, char** argv) {
  static_assert(
      sizeof(double) == 8 && std::numeric_limits<double>::digits == 53,
      "HiPDLP performance driver requires FP64");
  const Clock::time_point process_start = Clock::now();
  try {
    std::map<std::string, std::string> arguments = {
        {"--instance", ""},
        {"--backend", "cpu"},
        {"--step-size", "pid"},
        {"--time-limit", "300"},
        {"--iteration-limit", "2147483647"},
        {"--output-dir", ""},
        {"--warmups", "1"},
        {"--repetitions", "5"},
        {"--seed", "20260908"}};
    std::set<std::string> seen;
    for (int index = 1; index < argc; index += 2) {
      const std::string flag = argv[index];
      if (index + 1 >= argc || !arguments.count(flag) ||
          !seen.insert(flag).second)
        throw std::invalid_argument("unknown, duplicate or valueless option: " +
                                    flag);
      arguments[flag] = argv[index + 1];
    }
    const std::string instance = arguments["--instance"];
    const std::string directory = arguments["--output-dir"];
    const std::string backend = arguments["--backend"];
    const std::string requested_step = arguments["--step-size"];
    if (instance.empty() || directory.empty())
      throw std::invalid_argument(
          "--instance and --output-dir are required; create a new output "
          "directory first");
    if ((backend != "cpu" && backend != "cuda") ||
        (backend == "cpu" && !kCompiledCpu) ||
        (backend == "cuda" && !kCompiledCuda))
      throw std::invalid_argument(
          "backend is not supported by this compiled driver: " + backend);
    const std::map<std::string, HighsInt> strategies = {
        {"fixed", kPdlpStepSizeStrategyFixed},
        {"pid", kPdlpStepSizeStrategyPid},
        {"adaptive", kPdlpStepSizeStrategyAdaptive}};
    if (!strategies.count(requested_step))
      throw std::invalid_argument("--step-size must be fixed, pid or adaptive");
    const HighsInt warmups = integer(arguments["--warmups"]);
    const HighsInt repetitions = integer(arguments["--repetitions"]);
    const HighsInt iteration_limit = integer(arguments["--iteration-limit"]);
    const HighsInt seed = integer(arguments["--seed"]);
    if (repetitions < 1 ||
        warmups > std::numeric_limits<HighsInt>::max() - repetitions)
      throw std::invalid_argument("invalid repetition count");
    std::size_t consumed = 0;
    const double time_limit = std::stod(arguments["--time-limit"], &consumed);
    if (consumed != arguments["--time-limit"].size() ||
        !std::isfinite(time_limit) || time_limit < 0)
      throw std::invalid_argument(
          "--time-limit must be finite and nonnegative");
    std::string effective_step = requested_step;
    if (requested_step == "adaptive") {
#ifdef HIPDLP_PERFORMANCE_LEGACY
      effective_step = "pid";
#else
      effective_step = "fixed";
#endif
    }
    FILE* guard = exclusiveFile(directory + "/.driver.lock");
    std::fprintf(guard, "%ld\n", static_cast<long>(getpid()));
    closeFile(guard);
    Highs highs;
    const auto requireOption = [](const HighsStatus status) {
      if (status != HighsStatus::kOk)
        throw std::runtime_error("failed to set solver option");
    };
    requireOption(highs.setOptionValue("output_flag", false));
    requireOption(highs.setOptionValue("log_to_console", false));
    requireOption(highs.setOptionValue("solver", "hipdlp"));
    requireOption(highs.setOptionValue("presolve", "off"));
    requireOption(highs.setOptionValue("parallel", "off"));
    requireOption(highs.setOptionValue("threads", HighsInt{1}));
    requireOption(highs.setOptionValue("random_seed", seed));
    requireOption(highs.setOptionValue("time_limit", time_limit));
    requireOption(
        highs.setOptionValue("pdlp_iteration_limit", iteration_limit));
    requireOption(highs.setOptionValue("kkt_tolerance", 1e-4));
    requireOption(highs.setOptionValue("primal_feasibility_tolerance", 1e-4));
    requireOption(highs.setOptionValue("dual_feasibility_tolerance", 1e-4));
    requireOption(highs.setOptionValue("pdlp_optimality_tolerance", 1e-4));
    requireOption(highs.setOptionValue("pdlp_step_size_strategy",
                                       strategies.at(requested_step)));
    requireOption(highs.setOptionValue("pdlp_restart_strategy",
                                       HighsInt{kPdlpRestartStrategyHalpern}));
#ifndef HIPDLP_PERFORMANCE_LEGACY
    requireOption(highs.setOptionValue("hipdlp_backend", backend));
    requireOption(highs.setOptionValue("hipdlp_scaling", "auto"));
#endif
    const Clock::time_point read_start = Clock::now();
    if (highs.readModel(instance) != HighsStatus::kOk)
      throw std::runtime_error("failed to read model: " + instance);
    const double read_time = secondsSince(read_start);
    const HighsLp& lp = highs.getLp();
    const HighsOptions& options = highs.getOptions();
    for (HighsInt sequence = 0; sequence < warmups + repetitions; ++sequence) {
      const Clock::time_point repeat_start = Clock::now();
      const bool warmup = sequence < warmups;
      const HighsInt index = warmup ? sequence : sequence - warmups;
      std::ostringstream basename;
      basename << (warmup ? "warmup-" : "measurement-") << std::setfill('0')
               << std::setw(3) << index;
      const std::string solution_name = basename.str() + ".sol";
      FILE* solution_file = exclusiveFile(directory + "/" + solution_name);
      HighsBasis basis;
      HighsSolution solution;
      HighsInfo info;
      HighsModelStatus model_status = HighsModelStatus::kNotset;
      HighsCallback callback(nullptr);
      HighsTimer timer;
      timer.setPrintfFlag(false, false);
      std::printf("{\"event\":\"begin\",\"sequence\":%lld}\n",
                  static_cast<long long>(sequence));
      std::fflush(stdout);
      timer.start();
      const Clock::time_point solve_start = Clock::now();
      const HighsStatus status = solveLpHiPdlp(
          options, timer, lp, basis, solution, model_status, info, callback);
      const double solver_time = secondsSince(solve_start);
      timer.stop();
      const HighsInt stored_status = model_status == HighsModelStatus::kOptimal
                                         ? kSolutionStatusFeasible
                                         : kSolutionStatusInfeasible;
      info.primal_solution_status =
          solution.value_valid ? stored_status : kSolutionStatusNone;
      info.dual_solution_status =
          solution.dual_valid ? stored_status : kSolutionStatusNone;
      writeSolutionFile(solution_file, options, highs.getModel(), basis,
                        solution, info, model_status, kSolutionStyleRaw);
      closeFile(solution_file);
      std::ostringstream json;
      json.imbue(std::locale::classic());
      json << std::setprecision(17) << std::boolalpha
           << "{\"schema\":\"hipdlp-performance/1\",\"event\":\"result\""
           << ",\"instance\":" << quoted(instance)
           << ",\"implementation\":" << quoted(kImplementation)
           << ",\"backend\":" << quoted(backend)
           << ",\"compiled_cpu\":" << kCompiledCpu
           << ",\"compiled_cuda\":" << kCompiledCuda
           << ",\"kind\":" << quoted(warmup ? "warmup" : "measurement")
           << ",\"index\":" << index << ",\"sequence\":" << sequence
           << ",\"cold\":" << (sequence == 0) << ",\"pid\":" << getpid()
           << ",\"seed\":" << seed << ",\"warmups\":" << warmups
           << ",\"repetitions\":" << repetitions
           << ",\"status\":" << quoted(statusName(model_status))
           << ",\"return_status\":"
           << quoted(status == HighsStatus::kOk        ? "ok"
                     : status == HighsStatus::kWarning ? "warning"
                                                       : "error")
           << ",\"iterations\":" << info.pdlp_iteration_count
           << ",\"read_time_sec\":" << read_time
           << ",\"solver_time_sec\":" << solver_time
           << ",\"repeat_wall_sec\":" << secondsSince(repeat_start)
           << ",\"end_to_end_sec\":" << secondsSince(process_start)
           << ",\"time_limit_sec\":" << time_limit
           << ",\"iteration_limit\":" << iteration_limit
           << ",\"requested_step_size\":" << quoted(requested_step)
           << ",\"requested_step_size_enum\":" << strategies.at(requested_step)
           << ",\"effective_step_size\":" << quoted(effective_step)
           << ",\"scaling\":" << quoted(kScaling)
           << ",\"pdlp_scaling_mode\":" << options.pdlp_scaling_mode
           << ",\"pdlp_ruiz_iterations\":" << options.pdlp_ruiz_iterations
           << ",\"kkt_tolerance\":" << options.kkt_tolerance
           << ",\"pdlp_optimality_tolerance\":"
           << options.pdlp_optimality_tolerance
           << ",\"primal_feasibility_tolerance\":"
           << options.primal_feasibility_tolerance
           << ",\"dual_feasibility_tolerance\":"
           << options.dual_feasibility_tolerance
           << ",\"presolve\":\"off\",\"threads\":1,\"precision\":\"fp64\""
           << ",\"solution_file\":" << quoted(solution_name) << "}\n";
      FILE* result_file =
          exclusiveFile(directory + "/" + basename.str() + ".json");
      std::fputs(json.str().c_str(), result_file);
      closeFile(result_file);
      std::fputs(json.str().c_str(), stdout);
      std::fflush(stdout);
    }
    return 0;
  } catch (const std::exception& error) {
    std::fprintf(stderr, "HiPdlpPerformance: %s\n", error.what());
    return 2;
  }
}