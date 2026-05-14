#include "ExperimentalHighsCupdlpx.h"

#include <chrono>
#include <cmath>
#include <limits>
#include <sstream>
#include <vector>

#include "cupdlpx.h"
#include "lp_data/HConst.h"

namespace {

const bool kExperimentalOutputFlag = false;
const double kExperimentalKktTolerance = 1e-4;
const HighsInt kExperimentalPdlpIterationLimit = 859400;
const double kCupdlpxInfiniteBoundThreshold = 1e20;

template <typename T>
const T* dataOrNull(const std::vector<T>& values) {
  return values.empty() ? nullptr : values.data();
}

bool checkedHighsIntToInt(const HighsInt value, int& converted,
                          std::string& error_message) {
  if (value < 0 || value > std::numeric_limits<int>::max()) {
    std::ostringstream stream;
    stream << "Value " << value << " does not fit in int";
    error_message = stream.str();
    return false;
  }
  converted = static_cast<int>(value);
  return true;
}

std::string basenameOf(const std::string& path) {
  const size_t slash = path.find_last_of("/\\");
  return slash == std::string::npos ? path : path.substr(slash + 1);
}

std::string modelStatusToString(const HighsModelStatus model_status) {
  Highs highs;
  return highs.modelStatusToString(model_status);
}

bool zeroColumnReducedProblemFeasible(const HighsLp& reduced_lp) {
  for (HighsInt row = 0; row < reduced_lp.num_row_; ++row) {
    if (reduced_lp.row_lower_[row] > 0.0) return false;
    if (reduced_lp.row_upper_[row] < 0.0) return false;
  }
  return true;
}

std::vector<double> normaliseBoundsForCupdlpx(const std::vector<double>& bounds,
                                              const bool is_lower) {
  std::vector<double> normalised = bounds;
  const double infinite_value = is_lower ? -INFINITY : INFINITY;
  for (double& value : normalised) {
    if (is_lower) {
      if (value < -kCupdlpxInfiniteBoundThreshold) value = infinite_value;
    } else {
      if (value > kCupdlpxInfiniteBoundThreshold) value = infinite_value;
    }
  }
  return normalised;
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

HighsModelStatus cupdlpxTerminationReasonToModelStatus(
    const termination_reason_t termination_reason) {
  switch (termination_reason) {
    case TERMINATION_REASON_OPTIMAL:
      return HighsModelStatus::kOptimal;
    case TERMINATION_REASON_PRIMAL_INFEASIBLE:
      return HighsModelStatus::kInfeasible;
    case TERMINATION_REASON_DUAL_INFEASIBLE:
      return HighsModelStatus::kUnbounded;
    case TERMINATION_REASON_INFEASIBLE_OR_UNBOUNDED:
      return HighsModelStatus::kUnboundedOrInfeasible;
    case TERMINATION_REASON_TIME_LIMIT:
      return HighsModelStatus::kTimeLimit;
    case TERMINATION_REASON_ITERATION_LIMIT:
      return HighsModelStatus::kIterationLimit;
    case TERMINATION_REASON_FEAS_POLISH_SUCCESS:
    case TERMINATION_REASON_UNSPECIFIED:
    default:
      return HighsModelStatus::kUnknown;
  }
}

const char* cupdlpxTerminationReasonToString(
    const termination_reason_t termination_reason) {
  switch (termination_reason) {
    case TERMINATION_REASON_OPTIMAL:
      return "OPTIMAL";
    case TERMINATION_REASON_PRIMAL_INFEASIBLE:
      return "PRIMAL_INFEASIBLE";
    case TERMINATION_REASON_DUAL_INFEASIBLE:
      return "DUAL_INFEASIBLE";
    case TERMINATION_REASON_INFEASIBLE_OR_UNBOUNDED:
      return "INFEASIBLE_OR_UNBOUNDED";
    case TERMINATION_REASON_TIME_LIMIT:
      return "TIME_LIMIT";
    case TERMINATION_REASON_ITERATION_LIMIT:
      return "ITERATION_LIMIT";
    case TERMINATION_REASON_FEAS_POLISH_SUCCESS:
      return "FEAS_POLISH_SUCCESS";
    case TERMINATION_REASON_UNSPECIFIED:
    default:
      return "UNSPECIFIED";
  }
}

HighsStatus loadModelAsLp(const std::string& model_file, HighsLp& lp,
                          std::string& error_message) {
  Highs highs;
  HighsStatus status =
      highs.setOptionValue("output_flag", kExperimentalOutputFlag);
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

HighsStatus convertHighsLpToCsr(const HighsLp& lp,
                                ExperimentalHighsCupdlpxCsrMatrix& csr,
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

  csr.row_ptr.resize(matrix.start_.size());
  csr.col_ind.resize(matrix.index_.size());
  csr.value = matrix.value_;

  for (size_t i = 0; i < matrix.start_.size(); ++i) {
    int converted = 0;
    if (!checkedHighsIntToInt(matrix.start_[i], converted, error_message)) {
      return HighsStatus::kError;
    }
    csr.row_ptr[i] = converted;
  }
  for (size_t i = 0; i < matrix.index_.size(); ++i) {
    int converted = 0;
    if (!checkedHighsIntToInt(matrix.index_[i], converted, error_message)) {
      return HighsStatus::kError;
    }
    csr.col_ind[i] = converted;
  }
  return HighsStatus::kOk;
}

}  // namespace

HighsStatus runExperimentalHighsCupdlpxBenchmark(
    const std::string& model_file, ExperimentalHighsCupdlpxRun& result,
    std::string& error_message) {
  result = ExperimentalHighsCupdlpxRun();
  result.instance = basenameOf(model_file);

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

  Highs highs;
  status = highs.setOptionValue("output_flag", kExperimentalOutputFlag);
  if (status != HighsStatus::kOk) {
    error_message = "Failed to set output_flag";
    return HighsStatus::kError;
  }
  status = highs.passModel(original_lp);
  if (status != HighsStatus::kOk) {
    error_message = "Failed to pass LP model to HiGHS presolver";
    return HighsStatus::kError;
  }

  const auto presolve_start = std::chrono::steady_clock::now();
  status = highs.presolve();
  result.presolve_time =
      std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                    presolve_start)
          .count();
  if (status != HighsStatus::kOk) {
    error_message = "HiGHS presolve returned an error";
    return HighsStatus::kError;
  }

  const HighsPresolveStatus presolve_status = highs.getModelPresolveStatus();
  if (presolve_status == HighsPresolveStatus::kInfeasible) {
    result.model_status = HighsModelStatus::kInfeasible;
    result.total_time =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      total_start)
            .count();
    return HighsStatus::kOk;
  }
  if (presolve_status == HighsPresolveStatus::kUnboundedOrInfeasible) {
    result.model_status = HighsModelStatus::kUnboundedOrInfeasible;
    result.total_time =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      total_start)
            .count();
    return HighsStatus::kOk;
  }
  if (presolve_status == HighsPresolveStatus::kTimeout) {
    error_message = "HiGHS presolve timed out";
    return HighsStatus::kError;
  }

  const HighsLp& presolved_lp = highs.getPresolvedLp();
  result.reduced_rows = presolved_lp.num_row_;
  result.reduced_cols = presolved_lp.num_col_;
  result.reduced_nnz = presolved_lp.a_matrix_.numNz();

  if (presolve_status == HighsPresolveStatus::kReducedToEmpty) {
    HighsSolution reduced_solution;
    reduced_solution.value_valid = true;
    reduced_solution.dual_valid = true;
    status = highs.postsolve(reduced_solution);
    if (status != HighsStatus::kOk) {
      error_message = "HiGHS postsolve failed for empty presolved LP";
      return HighsStatus::kError;
    }
    result.iterations = 0;
    result.solve_time = 0.0;
    result.model_status = highs.getModelStatus();
    result.objective = highs.getInfo().objective_function_value;
    result.total_time =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      total_start)
            .count();
    return HighsStatus::kOk;
  }

  if (presolved_lp.num_col_ == 0) {
    if (!zeroColumnReducedProblemFeasible(presolved_lp)) {
      result.model_status = HighsModelStatus::kInfeasible;
      result.total_time =
          std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                        total_start)
              .count();
      return HighsStatus::kOk;
    }
    HighsSolution reduced_solution;
    initialiseZeroColumnReducedSolution(presolved_lp, reduced_solution);
    const HighsStatus postsolve_status = highs.postsolve(reduced_solution);
    if (postsolve_status == HighsStatus::kError) {
      error_message = "HiGHS postsolve failed for zero-column presolved LP";
      return HighsStatus::kError;
    }
    result.iterations = 0;
    result.solve_time = 0.0;
    result.model_status = HighsModelStatus::kOptimal;
    result.objective = highs.getInfo().objective_function_value;
    result.total_time =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      total_start)
            .count();
    return HighsStatus::kOk;
  }

  ExperimentalHighsCupdlpxCsrMatrix csr;
  status = convertHighsLpToCsr(presolved_lp, csr, error_message);
  if (status != HighsStatus::kOk) return status;

  const std::vector<double> row_lower =
      normaliseBoundsForCupdlpx(presolved_lp.row_lower_, true);
  const std::vector<double> row_upper =
      normaliseBoundsForCupdlpx(presolved_lp.row_upper_, false);
  const std::vector<double> col_lower =
      normaliseBoundsForCupdlpx(presolved_lp.col_lower_, true);
  const std::vector<double> col_upper =
      normaliseBoundsForCupdlpx(presolved_lp.col_upper_, false);

  std::vector<double> objective = presolved_lp.col_cost_;
  double objective_constant = presolved_lp.offset_;
  if (presolved_lp.sense_ == ObjSense::kMaximize) {
    for (double& value : objective) value = -value;
    objective_constant = -objective_constant;
  }

  matrix_desc_t matrix_desc;
  matrix_desc.m = static_cast<int>(presolved_lp.num_row_);
  matrix_desc.n = static_cast<int>(presolved_lp.num_col_);
  matrix_desc.fmt = matrix_csr;
  matrix_desc.data.csr.nnz = static_cast<int>(csr.value.size());
  matrix_desc.data.csr.row_ptr = dataOrNull(csr.row_ptr);
  matrix_desc.data.csr.col_ind = dataOrNull(csr.col_ind);
  matrix_desc.data.csr.vals = dataOrNull(csr.value);

  lp_problem_t* problem = create_lp_problem(
      dataOrNull(objective), &matrix_desc, dataOrNull(row_lower),
      dataOrNull(row_upper), dataOrNull(col_lower), dataOrNull(col_upper),
      &objective_constant);
  if (!problem) {
    error_message = "cuPDLPx create_lp_problem returned nullptr";
    return HighsStatus::kError;
  }

  pdhg_parameters_t params;
  set_default_parameters(&params);
  params.verbose = false;
  params.presolve = false;
  params.feasibility_polishing = false;
  params.matrix_zero_tol = 0.0;
  params.termination_criteria.eps_optimal_relative =
      kExperimentalKktTolerance;
  params.termination_criteria.eps_feasible_relative =
      kExperimentalKktTolerance;
  params.termination_criteria.iteration_limit =
      static_cast<int>(kExperimentalPdlpIterationLimit);

  cupdlpx_result_t* cupdlpx_result = solve_lp_problem(problem, &params);
  lp_problem_free(problem);
  if (!cupdlpx_result) {
    error_message = "cuPDLPx solve_lp_problem returned nullptr";
    return HighsStatus::kError;
  }

  result.rescaling_time = cupdlpx_result->rescaling_time_sec;
  result.iterations = cupdlpx_result->total_count;
  result.solve_time = cupdlpx_result->rescaling_time_sec +
                      cupdlpx_result->cumulative_time_sec +
                      cupdlpx_result->feasibility_polishing_time;
  result.model_status = cupdlpxTerminationReasonToModelStatus(
      cupdlpx_result->termination_reason);
  result.termination_reason =
      cupdlpxTerminationReasonToString(cupdlpx_result->termination_reason);

  HighsSolution reduced_solution;
  if (cupdlpx_result->primal_solution) {
    reduced_solution.value_valid = true;
    reduced_solution.col_value.assign(
        cupdlpx_result->primal_solution,
        cupdlpx_result->primal_solution + presolved_lp.num_col_);
  }
  if (cupdlpx_result->dual_solution) {
    reduced_solution.row_dual.assign(
        cupdlpx_result->dual_solution,
        cupdlpx_result->dual_solution + presolved_lp.num_row_);
  }

  if (reduced_solution.value_valid) {
    const HighsStatus postsolve_status = highs.postsolve(reduced_solution);
    if (postsolve_status == HighsStatus::kError) {
      cupdlpx_result_free(cupdlpx_result);
      error_message =
          "HiGHS postsolve failed for cuPDLPx solution on presolved LP";
      return HighsStatus::kError;
    }
    const HighsSolution original_solution = highs.getSolution();
    if (original_solution.col_value.size() ==
        static_cast<size_t>(original_lp.num_col_)) {
      result.objective = original_lp.objectiveValue(original_solution.col_value);
    }
  }

  if (result.objective == kHighsUndefined &&
      cupdlpx_result->primal_solution &&
      static_cast<HighsInt>(cupdlpx_result->num_variables) ==
          presolved_lp.num_col_) {
    std::vector<double> presolved_solution(cupdlpx_result->primal_solution,
                                           cupdlpx_result->primal_solution +
                                               cupdlpx_result->num_variables);
    result.objective = presolved_lp.objectiveValue(presolved_solution);
  }

  result.total_time =
      std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                    total_start)
          .count();
  cupdlpx_result_free(cupdlpx_result);
  return HighsStatus::kOk;
}

std::string experimentalHighsCupdlpxCsvLine(
    const ExperimentalHighsCupdlpxRun& result) {
  std::ostringstream stream;
  stream << result.instance << "," << result.orig_rows << "," << result.orig_cols
         << "," << result.orig_nnz << "," << result.reduced_rows << ","
         << result.reduced_cols << "," << result.reduced_nnz << ","
         << result.presolve_time << "," << result.rescaling_time << ","
         << result.iterations << "," << result.solve_time << ","
         << result.total_time << "," << result.objective << ","
         << modelStatusToString(result.model_status) << ","
         << result.termination_reason;
  return stream.str();
}
