#include "ExperimentalPslpHiPdlp.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <sstream>

#include "lp_data/HighsLpUtils.h"
#include "PSLP/PSLP_sol.h"
#include "PSLP/PSLP_status.h"

namespace {

const bool kExperimentalOutputFlag = false;
const double kExperimentalKktTolerance = 1e-4;
const HighsInt kExperimentalPdlpRestartStrategy = 3;
const HighsInt kExperimentalPdlpScalingMode = 5;
const HighsInt kExperimentalPdlpStepSizeStrategy = 3;
const HighsInt kExperimentalPdlpIterationLimit = 859400;
const HighsInt kExperimentalLogDevLevel = 1;

template <typename T>
const T* dataOrNull(const std::vector<T>& values) {
  return values.empty() ? nullptr : values.data();
}

template <typename T>
T* dataOrNull(std::vector<T>& values) {
  return values.empty() ? nullptr : values.data();
}

bool checkedHighsIntToInt(const HighsInt value, int& converted,
                          std::string& error_message) {
  if (value < 0 || value > std::numeric_limits<int>::max()) {
    std::ostringstream stream;
    stream << "Value " << value << " does not fit in int for PSLP CSR data";
    error_message = stream.str();
    return false;
  }
  converted = static_cast<int>(value);
  return true;
}

bool checkedSizeTToHighsInt(const size_t value, HighsInt& converted,
                            std::string& error_message) {
  if (value > static_cast<size_t>(std::numeric_limits<HighsInt>::max())) {
    std::ostringstream stream;
    stream << "Value " << value << " does not fit in HighsInt";
    error_message = stream.str();
    return false;
  }
  converted = static_cast<HighsInt>(value);
  return true;
}

std::string basenameOf(const std::string& path) {
  const size_t slash = path.find_last_of("/\\");
  return slash == std::string::npos ? path : path.substr(slash + 1);
}

bool isModeValid(const std::string& mode) {
  return mode == "none" || mode == "plain" || mode == "highs" ||
         mode == "pslp";
}

std::string canonicalMode(const std::string& mode) {
  return mode == "plain" ? "none" : mode;
}

bool modelStatusHasReliableObjective(const HighsModelStatus model_status) {
  return model_status == HighsModelStatus::kOptimal ||
         model_status == HighsModelStatus::kModelEmpty;
}

HighsStatus setCommonHiPdlpOptions(Highs& highs, const bool use_presolve,
                                   std::string& error_message) {
  HighsStatus status = highs.setOptionValue("output_flag", kExperimentalOutputFlag);
  if (status != HighsStatus::kOk) {
    error_message = "Failed to set output_flag";
    return HighsStatus::kError;
  }
  status = highs.setOptionValue("solver", kHiPdlpString);
  if (status != HighsStatus::kOk) {
    error_message = "Failed to set solver=hipdlp";
    return HighsStatus::kError;
  }
  status = highs.setOptionValue("presolve",
                                use_presolve ? kHighsOnString : kHighsOffString);
  if (status != HighsStatus::kOk) {
    error_message =
        use_presolve ? "Failed to set presolve=on" : "Failed to set presolve=off";
    return HighsStatus::kError;
  }
  status = highs.setOptionValue("kkt_tolerance", kExperimentalKktTolerance);
  if (status != HighsStatus::kOk) {
    error_message = "Failed to set kkt_tolerance";
    return HighsStatus::kError;
  }
  status = highs.setOptionValue("pdlp_restart_strategy",
                                kExperimentalPdlpRestartStrategy);
  if (status != HighsStatus::kOk) {
    error_message = "Failed to set pdlp_restart_strategy";
    return HighsStatus::kError;
  }
  status = highs.setOptionValue("pdlp_scaling_mode",
                                kExperimentalPdlpScalingMode);
  if (status != HighsStatus::kOk) {
    error_message = "Failed to set pdlp_scaling_mode";
    return HighsStatus::kError;
  }
  status = highs.setOptionValue("pdlp_step_size_strategy",
                                kExperimentalPdlpStepSizeStrategy);
  if (status != HighsStatus::kOk) {
    error_message = "Failed to set pdlp_step_size_strategy";
    return HighsStatus::kError;
  }
  status = highs.setOptionValue("pdlp_iteration_limit",
                                kExperimentalPdlpIterationLimit);
  if (status != HighsStatus::kOk) {
    error_message = "Failed to set pdlp_iteration_limit";
    return HighsStatus::kError;
  }
  status = highs.setOptionValue("log_dev_level", kExperimentalLogDevLevel);
  if (status != HighsStatus::kOk) {
    error_message = "Failed to set log_dev_level";
    return HighsStatus::kError;
  }
  return HighsStatus::kOk;
}

HighsStatus loadModelAsLp(const std::string& model_file, HighsLp& lp,
                          std::string& error_message) {
  Highs highs;
  HighsStatus status = highs.setOptionValue("output_flag", kExperimentalOutputFlag);
  if (status != HighsStatus::kOk) {
    error_message = "Failed to set output_flag while reading model";
    return HighsStatus::kError;
  }
  status = highs.readModel(model_file);
  if (status == HighsStatus::kError) {
    error_message = "Failed to read model file: " + model_file;
    return HighsStatus::kError;
  }
  lp = highs.getLp();
  lp.ensureColwise();
  return HighsStatus::kOk;
}

HighsStatus solveWithHiPdlp(const HighsLp& lp, HighsModelStatus& model_status,
                            HighsInfo& info, HighsSolution& solution,
                            double& solve_time, std::string& error_message,
                            const bool use_presolve = false) {
  Highs highs;
  HighsStatus status = setCommonHiPdlpOptions(highs, use_presolve, error_message);
  if (status != HighsStatus::kOk) return status;

  status = highs.passModel(lp);
  if (status != HighsStatus::kOk) {
    error_message = "Failed to pass LP model to HiGHS";
    return HighsStatus::kError;
  }

  const auto solve_start = std::chrono::steady_clock::now();
  status = highs.run();
  solve_time =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - solve_start)
          .count();
  if (status == HighsStatus::kError) {
    error_message = "HiPDLP solve returned HighsStatus::kError";
    return HighsStatus::kError;
  }

  model_status = highs.getModelStatus();
  info = highs.getInfo();
  solution = highs.getSolution();
  return HighsStatus::kOk;
}

HighsStatus solvePresolvedWithHiPdlp(
    Highs& original_highs, HighsModelStatus& model_status, HighsInfo& info,
    HighsSolution& original_solution, double& presolve_time, double& solve_time,
    HighsInt& reduced_rows, HighsInt& reduced_cols, HighsInt& reduced_nnz,
    std::string& error_message) {
  const auto presolve_start = std::chrono::steady_clock::now();
  HighsStatus status = original_highs.presolve();
  presolve_time =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - presolve_start)
          .count();
  if (status != HighsStatus::kOk) {
    error_message = "HiGHS presolve returned an error";
    return HighsStatus::kError;
  }

  const HighsPresolveStatus presolve_status =
      original_highs.getModelPresolveStatus();
  if (presolve_status == HighsPresolveStatus::kInfeasible) {
    model_status = HighsModelStatus::kInfeasible;
    return HighsStatus::kOk;
  }
  if (presolve_status == HighsPresolveStatus::kUnboundedOrInfeasible) {
    model_status = HighsModelStatus::kUnboundedOrInfeasible;
    return HighsStatus::kOk;
  }
  if (presolve_status == HighsPresolveStatus::kTimeout) {
    model_status = original_highs.getModelStatus();
    error_message = "HiGHS presolve timed out";
    return HighsStatus::kError;
  }

  const HighsLp& presolved_lp = original_highs.getPresolvedLp();
  reduced_rows = presolved_lp.num_row_;
  reduced_cols = presolved_lp.num_col_;
  reduced_nnz = presolved_lp.a_matrix_.numNz();

  if (presolve_status == HighsPresolveStatus::kReducedToEmpty) {
    HighsSolution reduced_solution;
    reduced_solution.value_valid = true;
    reduced_solution.dual_valid = true;
    status = original_highs.postsolve(reduced_solution);
    if (status != HighsStatus::kOk) {
      error_message = "HiGHS postsolve failed for empty presolved LP";
      return HighsStatus::kError;
    }
    model_status = original_highs.getModelStatus();
    info = original_highs.getInfo();
    info.pdlp_iteration_count = 0;
    original_solution = original_highs.getSolution();
    solve_time = 0.0;
    return HighsStatus::kOk;
  }

  HighsSolution reduced_solution;
  status = solveWithHiPdlp(presolved_lp, model_status, info, reduced_solution,
                           solve_time, error_message, false);
  if (status != HighsStatus::kOk) return status;

  if (modelStatusHasReliableObjective(model_status)) {
    const HighsInt reduced_iteration_count = info.pdlp_iteration_count;
    status = original_highs.postsolve(reduced_solution);
    if (status != HighsStatus::kOk) {
      error_message = "HiGHS postsolve failed for presolved HiPDLP solution";
      return HighsStatus::kError;
    }
    model_status = original_highs.getModelStatus();
    original_solution = original_highs.getSolution();
    info.pdlp_iteration_count = reduced_iteration_count;
  } else {
    original_solution = reduced_solution;
  }

  return HighsStatus::kOk;
}

bool zeroColumnReducedProblemFeasible(const HighsLp& reduced_lp) {
  for (HighsInt row = 0; row < reduced_lp.num_row_; ++row) {
    if (reduced_lp.row_lower_[row] > 0.0) return false;
    if (reduced_lp.row_upper_[row] < 0.0) return false;
  }
  return true;
}

void initialiseZeroColumnReducedSolution(const HighsLp& reduced_lp,
                                         HighsSolution& reduced_solution) {
  reduced_solution.clear();
  reduced_solution.col_value.assign(reduced_lp.num_col_, 0.0);
  reduced_solution.col_dual.assign(reduced_lp.num_col_, 0.0);
  reduced_solution.row_value.assign(reduced_lp.num_row_, 0.0);
  reduced_solution.row_dual.assign(reduced_lp.num_row_, 0.0);
  reduced_solution.value_valid = true;
  reduced_solution.dual_valid = true;
}

void ensurePostsolveVectorsSized(const HighsLp& reduced_lp,
                                 HighsSolution& reduced_solution) {
  if (reduced_solution.col_value.size() !=
      static_cast<size_t>(reduced_lp.num_col_)) {
    reduced_solution.col_value.resize(reduced_lp.num_col_, 0.0);
  }
  if (reduced_solution.col_dual.size() !=
      static_cast<size_t>(reduced_lp.num_col_)) {
    reduced_solution.col_dual.resize(reduced_lp.num_col_, 0.0);
  }
  if (reduced_solution.row_dual.size() !=
      static_cast<size_t>(reduced_lp.num_row_)) {
    reduced_solution.row_dual.resize(reduced_lp.num_row_, 0.0);
  }
}

std::string modelStatusToString(const HighsModelStatus model_status) {
  Highs highs;
  return highs.modelStatusToString(model_status);
}

}  // namespace

HighsStatus experimentalConvertHighsLpToPslpCsr(
    const HighsLp& lp, ExperimentalPslpCsrMatrix& csr,
    std::string& error_message) {
  HighsSparseMatrix matrix = lp.a_matrix_;
  matrix.num_col_ = lp.num_col_;
  matrix.num_row_ = lp.num_row_;
  matrix.ensureRowwise();

  if (matrix.start_.size() != static_cast<size_t>(lp.num_row_ + 1)) {
    error_message = "Rowwise matrix start vector has invalid length";
    return HighsStatus::kError;
  }
  if (matrix.index_.size() != matrix.value_.size()) {
    error_message = "Matrix index and value arrays have different lengths";
    return HighsStatus::kError;
  }

  csr.Ax = matrix.value_;
  csr.Ai.resize(matrix.index_.size());
  csr.Ap.resize(matrix.start_.size());

  for (size_t i = 0; i < matrix.index_.size(); ++i) {
    int converted = 0;
    if (!checkedHighsIntToInt(matrix.index_[i], converted, error_message)) {
      return HighsStatus::kError;
    }
    csr.Ai[i] = converted;
  }
  for (size_t i = 0; i < matrix.start_.size(); ++i) {
    int converted = 0;
    if (!checkedHighsIntToInt(matrix.start_[i], converted, error_message)) {
      return HighsStatus::kError;
    }
    csr.Ap[i] = converted;
  }
  return HighsStatus::kOk;
}

HighsStatus experimentalConvertPslpReducedProblemToHighsLp(
    const PresolvedProblem& reduced_problem, const HighsLp& original_lp,
    HighsLp& reduced_lp, std::string& error_message) {
  HighsInt num_row = 0;
  HighsInt num_col = 0;
  HighsInt num_nz = 0;
  if (!checkedSizeTToHighsInt(reduced_problem.m, num_row, error_message) ||
      !checkedSizeTToHighsInt(reduced_problem.n, num_col, error_message) ||
      !checkedSizeTToHighsInt(reduced_problem.nnz, num_nz, error_message)) {
    return HighsStatus::kError;
  }

  reduced_lp.clear();
  reduced_lp.num_row_ = num_row;
  reduced_lp.num_col_ = num_col;
  reduced_lp.sense_ = original_lp.sense_;
  reduced_lp.offset_ = original_lp.offset_ + reduced_problem.obj_offset;
  reduced_lp.model_name_ = original_lp.model_name_;
  reduced_lp.origin_name_ = original_lp.origin_name_;
  reduced_lp.objective_name_ = original_lp.objective_name_;

  reduced_lp.col_cost_.assign(reduced_problem.c,
                              reduced_problem.c + reduced_problem.n);
  reduced_lp.col_lower_.assign(reduced_problem.lbs,
                               reduced_problem.lbs + reduced_problem.n);
  reduced_lp.col_upper_.assign(reduced_problem.ubs,
                               reduced_problem.ubs + reduced_problem.n);
  reduced_lp.row_lower_.assign(reduced_problem.lhs,
                               reduced_problem.lhs + reduced_problem.m);
  reduced_lp.row_upper_.assign(reduced_problem.rhs,
                               reduced_problem.rhs + reduced_problem.m);

  reduced_lp.a_matrix_.format_ = MatrixFormat::kColwise;
  reduced_lp.a_matrix_.num_row_ = reduced_lp.num_row_;
  reduced_lp.a_matrix_.num_col_ = reduced_lp.num_col_;
  reduced_lp.a_matrix_.start_.assign(reduced_lp.num_col_ + 1, 0);
  reduced_lp.a_matrix_.index_.assign(num_nz, 0);
  reduced_lp.a_matrix_.value_.assign(num_nz, 0.0);

  for (size_t row = 0; row < reduced_problem.m; ++row) {
    const int row_start = reduced_problem.Ap[row];
    const int row_end = reduced_problem.Ap[row + 1];
    if (row_start > row_end || row_start < 0 || row_end < 0) {
      error_message = "Invalid PSLP CSR row pointer range";
      return HighsStatus::kError;
    }
    for (int position = row_start; position < row_end; ++position) {
      const int col = reduced_problem.Ai[position];
      if (col < 0 || col >= static_cast<int>(reduced_problem.n)) {
        error_message = "Invalid PSLP CSR column index";
        return HighsStatus::kError;
      }
      reduced_lp.a_matrix_.start_[col + 1]++;
    }
  }

  for (HighsInt col = 0; col < reduced_lp.num_col_; ++col) {
    reduced_lp.a_matrix_.start_[col + 1] += reduced_lp.a_matrix_.start_[col];
  }

  std::vector<HighsInt> next = reduced_lp.a_matrix_.start_;
  for (size_t row = 0; row < reduced_problem.m; ++row) {
    for (int position = reduced_problem.Ap[row];
         position < reduced_problem.Ap[row + 1]; ++position) {
      const int col = reduced_problem.Ai[position];
      const HighsInt destination = next[col]++;
      reduced_lp.a_matrix_.index_[destination] = static_cast<HighsInt>(row);
      reduced_lp.a_matrix_.value_[destination] = reduced_problem.Ax[position];
    }
  }
  return HighsStatus::kOk;
}

HighsStatus runExperimentalHiPdlpPslpBenchmark(
    const std::string& model_file, const std::string& mode,
    ExperimentalHiPdlpPslpRun& result, std::string& error_message) {
  result = ExperimentalHiPdlpPslpRun();
  result.instance = basenameOf(model_file);
  result.mode = canonicalMode(mode);

  if (!isModeValid(mode)) {
    error_message = "Mode must be none, highs, or pslp";
    return HighsStatus::kError;
  }

  const auto total_start = std::chrono::steady_clock::now();

  HighsLp original_lp;
  HighsStatus status = loadModelAsLp(model_file, original_lp, error_message);
  if (status != HighsStatus::kOk) return status;

  result.orig_rows = original_lp.num_row_;
  result.orig_cols = original_lp.num_col_;
  result.orig_nnz = original_lp.a_matrix_.numNz();
  result.reduced_rows = result.orig_rows;
  result.reduced_cols = result.orig_cols;
  result.reduced_nnz = result.orig_nnz;

  if (result.mode == "none") {
    HighsInfo info;
    HighsSolution solution;
    status = solveWithHiPdlp(original_lp, result.model_status, info, solution,
                             result.solve_time, error_message);
    result.iterations = info.pdlp_iteration_count;
    if (status != HighsStatus::kOk) return status;
    if (modelStatusHasReliableObjective(result.model_status) &&
        solution.col_value.size() != static_cast<size_t>(original_lp.num_col_)) {
      error_message = "Plain HiPDLP solution has invalid column dimension";
      return HighsStatus::kError;
    }
    if (modelStatusHasReliableObjective(result.model_status)) {
      result.objective = original_lp.objectiveValue(solution.col_value);
    }
    result.total_time =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - total_start)
            .count();
    return HighsStatus::kOk;
  }

  if (result.mode == "highs") {
    Highs highs;
    status = setCommonHiPdlpOptions(highs, false, error_message);
    if (status != HighsStatus::kOk) return status;
    status = highs.passModel(original_lp);
    if (status != HighsStatus::kOk) {
      error_message = "Failed to pass LP model to HiGHS internal presolver";
      return HighsStatus::kError;
    }

    HighsInfo info;
    HighsSolution solution;
    status = solvePresolvedWithHiPdlp(
        highs, result.model_status, info, solution, result.presolve_time,
        result.solve_time, result.reduced_rows, result.reduced_cols,
        result.reduced_nnz, error_message);
    result.iterations = info.pdlp_iteration_count;
    if (status != HighsStatus::kOk) return status;
    if (modelStatusHasReliableObjective(result.model_status)) {
      if (solution.col_value.size() != static_cast<size_t>(original_lp.num_col_)) {
        error_message =
            "HiGHS postsolve solution has invalid column dimension";
        return HighsStatus::kError;
      }
      result.objective = original_lp.objectiveValue(solution.col_value);
    }
    result.total_time =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - total_start)
            .count();
    return HighsStatus::kOk;
  }

  ExperimentalPslpCsrMatrix csr;
  status = experimentalConvertHighsLpToPslpCsr(original_lp, csr, error_message);
  if (status != HighsStatus::kOk) return status;

  std::vector<double> lhs = original_lp.row_lower_;
  std::vector<double> rhs = original_lp.row_upper_;
  std::vector<double> lbs = original_lp.col_lower_;
  std::vector<double> ubs = original_lp.col_upper_;
  std::vector<double> c = original_lp.col_cost_;

  Settings* settings = default_settings();
  if (!settings) {
    error_message = "PSLP default_settings returned nullptr";
    return HighsStatus::kError;
  }
  settings->verbose = false;

  Presolver* presolver = nullptr;
  const auto pslp_start = std::chrono::steady_clock::now();
  presolver = new_presolver(
      dataOrNull(csr.Ax), dataOrNull(csr.Ai), dataOrNull(csr.Ap),
      static_cast<size_t>(original_lp.num_row_),
      static_cast<size_t>(original_lp.num_col_), csr.Ax.size(), dataOrNull(lhs),
      dataOrNull(rhs), dataOrNull(lbs), dataOrNull(ubs), dataOrNull(c),
      settings);
  if (!presolver) {
    free_settings(settings);
    error_message = "PSLP new_presolver returned nullptr";
    return HighsStatus::kError;
  }

  const PresolveStatus presolve_status = run_presolver(presolver);
  result.presolve_time =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - pslp_start)
          .count();

  if (!presolver->reduced_prob) {
    free_presolver(presolver);
    free_settings(settings);
    error_message = "PSLP did not populate reduced_prob";
    return HighsStatus::kError;
  }

  result.reduced_rows = static_cast<HighsInt>(presolver->reduced_prob->m);
  result.reduced_cols = static_cast<HighsInt>(presolver->reduced_prob->n);
  result.reduced_nnz = static_cast<HighsInt>(presolver->reduced_prob->nnz);

  if (presolve_status & INFEASIBLE) {
    result.model_status = HighsModelStatus::kInfeasible;
    result.total_time =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - total_start)
         .count();
  free_presolver(presolver);
  free_settings(settings);
  return HighsStatus::kOk;
}
  if (presolve_status & UNBNDORINFEAS) {
    result.model_status = HighsModelStatus::kUnboundedOrInfeasible;
    result.total_time =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - total_start)
            .count();
    free_presolver(presolver);
    free_settings(settings);
    return HighsStatus::kOk;
  }

  HighsLp reduced_lp;
  status = experimentalConvertPslpReducedProblemToHighsLp(
      *presolver->reduced_prob, original_lp, reduced_lp, error_message);
  if (status != HighsStatus::kOk) {
    free_presolver(presolver);
    free_settings(settings);
    return status;
  }

  HighsSolution reduced_solution;
  if (reduced_lp.num_col_ == 0) {
    if (!zeroColumnReducedProblemFeasible(reduced_lp)) {
      result.model_status = HighsModelStatus::kInfeasible;
      result.total_time =
          std::chrono::duration<double>(std::chrono::steady_clock::now() - total_start)
              .count();
      free_presolver(presolver);
      free_settings(settings);
      return HighsStatus::kOk;
    }
    initialiseZeroColumnReducedSolution(reduced_lp, reduced_solution);
    result.iterations = 0;
    result.model_status = HighsModelStatus::kOptimal;
  } else {
    HighsInfo reduced_info;
    status = solveWithHiPdlp(reduced_lp, result.model_status, reduced_info,
                             reduced_solution, result.solve_time, error_message);
    result.iterations = reduced_info.pdlp_iteration_count;
    if (status != HighsStatus::kOk) {
      free_presolver(presolver);
      free_settings(settings);
      return status;
    }
  }

  if (modelStatusHasReliableObjective(result.model_status)) {
    ensurePostsolveVectorsSized(reduced_lp, reduced_solution);
    postsolve(presolver, dataOrNull(reduced_solution.col_value),
              dataOrNull(reduced_solution.row_dual),
              dataOrNull(reduced_solution.col_dual));
    if (!presolver->sol) {
      free_presolver(presolver);
      free_settings(settings);
      error_message = "PSLP postsolve did not populate presolver->sol";
      return HighsStatus::kError;
    }
    if (presolver->sol->dim_x != static_cast<size_t>(original_lp.num_col_)) {
      free_presolver(presolver);
      free_settings(settings);
      error_message =
          "PSLP postsolve primal solution dimension does not match original LP";
      return HighsStatus::kError;
    }
    std::vector<double> original_col_value(
        presolver->sol->x, presolver->sol->x + presolver->sol->dim_x);
    result.objective = original_lp.objectiveValue(original_col_value);
  }

  result.total_time =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - total_start)
          .count();
  free_presolver(presolver);
  free_settings(settings);
  return HighsStatus::kOk;
}

bool experimentalObjectivesClose(double lhs, double rhs,
                                 double relative_tolerance,
                                 double absolute_tolerance) {
  const double difference = std::fabs(lhs - rhs);
  if (difference <= absolute_tolerance) return true;
  const double scale = std::max(1.0, std::max(std::fabs(lhs), std::fabs(rhs)));
  return difference <= relative_tolerance * scale;
}

std::string experimentalHiPdlpPslpCsvLine(
    const ExperimentalHiPdlpPslpRun& result) {
  std::ostringstream stream;
  stream << result.instance << "," << result.mode << "," << result.orig_rows
         << "," << result.orig_cols << "," << result.orig_nnz << ","
         << result.reduced_rows << "," << result.reduced_cols << ","
         << result.reduced_nnz << "," << result.presolve_time << ","
         << result.iterations << "," << result.solve_time << ","
         << result.total_time << "," << result.objective << ","
         << modelStatusToString(result.model_status);
  return stream.str();
}
