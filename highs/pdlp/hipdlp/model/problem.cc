#include "pdlp/hipdlp/model/problem.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <string>
#include <utility>

#include "lp_data/HConst.h"

namespace hipdlp {
namespace model {

namespace {

struct MatrixEntry {
  HighsInt row;
  HighsInt col;
  double value;
};

bool entryLess(const MatrixEntry& left, const MatrixEntry& right) {
  if (left.row != right.row) return left.row < right.row;
  if (left.col != right.col) return left.col < right.col;
  if (left.value != right.value) return left.value < right.value;
  return std::signbit(left.value) && !std::signbit(right.value);
}

api::Error invalidModel(const std::string& message) {
  return api::Error(api::ErrorCode::kInvalidModel, message);
}

bool isNegativeInfinity(const double value) { return value <= -kHighsInf; }

bool isPositiveInfinity(const double value) { return value >= kHighsInf; }

bool normalizeLowerBound(const double value, double& normalized) {
  if (std::isnan(value) || isPositiveInfinity(value)) return false;
  if (isNegativeInfinity(value)) {
    normalized = -kHighsInf;
    return true;
  }
  if (!std::isfinite(value)) return false;
  normalized = value;
  return true;
}

bool normalizeUpperBound(const double value, double& normalized) {
  if (std::isnan(value) || isNegativeInfinity(value)) return false;
  if (isPositiveInfinity(value)) {
    normalized = kHighsInf;
    return true;
  }
  if (!std::isfinite(value)) return false;
  normalized = value;
  return true;
}

bool finiteVector(const std::vector<double>& values) {
  for (std::size_t index = 0; index < values.size(); ++index) {
    if (!std::isfinite(values[index])) return false;
  }
  return true;
}

double stableNorm(const std::vector<double>& values, bool& valid) {
  double norm = 0.0;
  for (std::size_t index = 0; index < values.size(); ++index) {
    norm = std::hypot(norm, values[index]);
    if (!std::isfinite(norm)) {
      valid = false;
      return 0.0;
    }
  }
  return norm;
}

api::Error canonicalizeEntries(const HighsInt num_rows, const HighsInt num_cols,
                               std::vector<MatrixEntry>& entries) {
  for (std::size_t entry = 0; entry < entries.size(); ++entry) {
    if (entries[entry].row < 0 || entries[entry].row >= num_rows ||
        entries[entry].col < 0 || entries[entry].col >= num_cols)
      return invalidModel("HiPDLP matrix index is out of range");
    if (!std::isfinite(entries[entry].value))
      return invalidModel("HiPDLP matrix contains a non-finite value");
  }
  std::sort(entries.begin(), entries.end(), entryLess);
  std::vector<MatrixEntry> canonical;
  canonical.reserve(entries.size());
  std::size_t entry = 0;
  while (entry < entries.size()) {
    const HighsInt row = entries[entry].row;
    const HighsInt col = entries[entry].col;
    long double sum = 0.0;
    do {
      sum += static_cast<long double>(entries[entry].value);
      if (!std::isfinite(sum))
        return invalidModel("HiPDLP duplicate matrix sum is non-finite");
      ++entry;
    } while (entry < entries.size() && entries[entry].row == row &&
             entries[entry].col == col);

    const double value = static_cast<double>(sum);
    if (!std::isfinite(value))
      return invalidModel("HiPDLP duplicate matrix sum exceeds FP64 range");
    if (value != 0.0) canonical.push_back({row, col, value});
  }
  entries.swap(canonical);
  return api::Error();
}

SparseMatrix makeCompressedMatrix(const HighsInt num_rows,
                                  const HighsInt num_cols,
                                  const std::vector<MatrixEntry>& entries) {
  SparseMatrix matrix;
  matrix.num_rows = num_rows;
  matrix.num_cols = num_cols;
  matrix.starts.assign(static_cast<std::size_t>(num_rows) + 1, 0);
  matrix.indices.reserve(entries.size());
  matrix.values.reserve(entries.size());
  HighsInt current_row = 0;
  for (std::size_t entry = 0; entry < entries.size(); ++entry) {
    while (current_row < entries[entry].row) {
      matrix.starts[static_cast<std::size_t>(++current_row)] =
          static_cast<HighsInt>(matrix.indices.size());
    }
    matrix.indices.push_back(entries[entry].col);
    matrix.values.push_back(entries[entry].value);
  }
  while (current_row < num_rows) {
    matrix.starts[static_cast<std::size_t>(++current_row)] =
        static_cast<HighsInt>(matrix.indices.size());
  }
  return matrix;
}

api::Error validateLp(const HighsLp& lp) {
  if (lp.num_row_ < 0 || lp.num_col_ < 0)
    return invalidModel("HiPDLP model dimensions must be nonnegative");
  if (lp.sense_ != ObjSense::kMinimize && lp.sense_ != ObjSense::kMaximize)
    return invalidModel("HiPDLP objective sense is invalid");
  if (!std::isfinite(lp.offset_))
    return invalidModel("HiPDLP objective offset must be finite");
  if (lp.col_cost_.size() != static_cast<std::size_t>(lp.num_col_) ||
      lp.col_lower_.size() != static_cast<std::size_t>(lp.num_col_) ||
      lp.col_upper_.size() != static_cast<std::size_t>(lp.num_col_) ||
      lp.row_lower_.size() != static_cast<std::size_t>(lp.num_row_) ||
      lp.row_upper_.size() != static_cast<std::size_t>(lp.num_row_))
    return invalidModel("HiPDLP model vector dimensions are inconsistent");
  if (!finiteVector(lp.col_cost_))
    return invalidModel("HiPDLP objective contains a non-finite value");

  const HighsSparseMatrix& matrix = lp.a_matrix_;
  if (matrix.num_row_ != lp.num_row_ || matrix.num_col_ != lp.num_col_)
    return invalidModel("HiPDLP matrix dimensions are inconsistent");
  if (matrix.format_ != MatrixFormat::kColwise &&
      matrix.format_ != MatrixFormat::kRowwise)
    return invalidModel("HiPDLP matrix format must be rowwise or colwise");
  const HighsInt outer_size =
      matrix.format_ == MatrixFormat::kColwise ? lp.num_col_ : lp.num_row_;
  if (matrix.start_.size() != static_cast<std::size_t>(outer_size) + 1 ||
      matrix.start_.empty() || matrix.start_[0] != 0)
    return invalidModel("HiPDLP sparse starts are malformed");
  for (HighsInt outer = 0; outer < outer_size; ++outer) {
    if (matrix.start_[outer] > matrix.start_[outer + 1] ||
        matrix.start_[outer] < 0)
      return invalidModel("HiPDLP sparse starts must be monotone");
  }
  const HighsInt num_nonzeros = matrix.start_[outer_size];
  if (num_nonzeros < 0 ||
      matrix.index_.size() != static_cast<std::size_t>(num_nonzeros) ||
      matrix.value_.size() != static_cast<std::size_t>(num_nonzeros))
    return invalidModel("HiPDLP sparse arrays have inconsistent lengths");
  return api::Error();
}

api::Error copyAndNormalizeBounds(const std::vector<double>& lower,
                                  const std::vector<double>& upper,
                                  std::vector<double>& normalized_lower,
                                  std::vector<double>& normalized_upper,
                                  const char* kind) {
  normalized_lower.resize(lower.size());
  normalized_upper.resize(upper.size());
  for (std::size_t index = 0; index < lower.size(); ++index) {
    if (!normalizeLowerBound(lower[index], normalized_lower[index]) ||
        !normalizeUpperBound(upper[index], normalized_upper[index]))
      return invalidModel(std::string("HiPDLP ") + kind +
                          " bound has a non-finite or wrong-direction "
                          "infinity");
    if (normalized_lower[index] > normalized_upper[index])
      return invalidModel(std::string("HiPDLP ") + kind +
                          " lower bound exceeds its upper bound");
  }
  return api::Error();
}

api::Error validateFiniteInput(const std::vector<double>& values,
                               const HighsInt expected_size, const char* name) {
  if (values.size() != static_cast<std::size_t>(expected_size))
    return invalidModel(std::string("HiPDLP ") + name +
                        " has the wrong dimension");
  if (!finiteVector(values))
    return invalidModel(std::string("HiPDLP ") + name +
                        " contains a non-finite value");
  return api::Error();
}

bool scaleBound(const double value, const double factor, const bool multiply,
                double& result) {
  if (isNegativeInfinity(value)) {
    result = -kHighsInf;
    return true;
  }
  if (isPositiveInfinity(value)) {
    result = kHighsInf;
    return true;
  }
  const long double scaled = multiply
                                 ? static_cast<long double>(value) * factor
                                 : static_cast<long double>(value) / factor;
  if (!std::isfinite(scaled) ||
      std::abs(scaled) > std::numeric_limits<double>::max())
    return false;
  result = static_cast<double>(scaled);
  return std::isfinite(result) && (scaled == 0.0 || result != 0.0);
}

bool scaleProduct(const double value, const double first, const double second,
                  double& result) {
  const long double scaled = static_cast<long double>(value) * first * second;
  if (!std::isfinite(scaled) ||
      std::abs(scaled) > std::numeric_limits<double>::max())
    return false;
  result = static_cast<double>(scaled);
  return std::isfinite(result) && (scaled == 0.0 || result != 0.0);
}

bool divideValue(const double value, const double divisor, double& result) {
  if (!std::isfinite(value) || !std::isfinite(divisor) || divisor <= 0.0)
    return false;
  const long double quotient =
      static_cast<long double>(value) / static_cast<long double>(divisor);
  if (!std::isfinite(quotient) ||
      std::abs(quotient) > std::numeric_limits<double>::max())
    return false;
  result = static_cast<double>(quotient);
  return std::isfinite(result) && (quotient == 0.0 || result != 0.0);
}

bool addProduct(const double left, const double right, double& sum) {
  double product = 0.0;
  if (!scaleProduct(left, right, 1.0, product)) return false;
  const double updated = sum + product;
  if (!std::isfinite(updated)) return false;
  sum = updated;
  return true;
}

}  // namespace

SparseMatrix::SparseMatrix() : num_rows(0), num_cols(0), starts(1, 0) {}

HighsInt SparseMatrix::numNonzeros() const {
  return static_cast<HighsInt>(values.size());
}

bool SparseMatrix::samePattern(const SparseMatrix& other) const {
  return num_rows == other.num_rows && num_cols == other.num_cols &&
         starts == other.starts && indices == other.indices;
}

ScalingMetadata::ScalingMetadata()
    : applied(false), curtis_reid_sweeps(0), ruiz_iterations(0) {}

TransformationMetadata::TransformationMetadata()
    : original_num_rows(0),
      original_num_cols(0),
      original_sense(ObjSense::kMinimize),
      original_offset(0.0),
      objective_sign(1.0) {}

UnscaledNorms::UnscaledNorms()
    : objective_l2(0.0), row_bounds_l2(0.0), matrix_frobenius(0.0) {}

OriginalModel::OriginalModel(const HighsLp& lp) : lp_(lp) {}

const HighsLp& OriginalModel::lp() const { return lp_; }

struct Problem::Data {
  HighsInt num_rows;
  HighsInt num_cols;
  ObjSense objective_sense;
  double objective_offset;
  IndexWidth index_width;
  std::vector<double> costs;
  std::vector<double> col_lower;
  std::vector<double> col_upper;
  std::vector<double> row_lower;
  std::vector<double> row_upper;
  std::vector<HighsBool> equality_rows;
  SparseMatrix matrix;
  SparseMatrix transpose;
  ScalingMetadata scaling;
  TransformationMetadata transformation;
  UnscaledNorms unscaled_norms;

  Data()
      : num_rows(0),
        num_cols(0),
        objective_sense(ObjSense::kMinimize),
        objective_offset(0.0),
        index_width(sizeof(HighsInt) > 4 ? IndexWidth::kInt64
                                         : IndexWidth::kInt32) {}
};

Problem::Problem() {}

Problem::~Problem() {}

Problem::Problem(Problem&& other) : data_(std::move(other.data_)) {}

Problem& Problem::operator=(Problem&& other) {
  data_ = std::move(other.data_);
  return *this;
}

bool Problem::initialized() const { return static_cast<bool>(data_); }

HighsInt Problem::numRows() const {
  assert(initialized());
  return data_->num_rows;
}

HighsInt Problem::numCols() const {
  assert(initialized());
  return data_->num_cols;
}

ObjSense Problem::objectiveSense() const {
  assert(initialized());
  return data_->objective_sense;
}

double Problem::objectiveOffset() const {
  assert(initialized());
  return data_->objective_offset;
}

IndexWidth Problem::indexWidth() const {
  assert(initialized());
  return data_->index_width;
}

const std::vector<double>& Problem::costs() const {
  assert(initialized());
  return data_->costs;
}

const std::vector<double>& Problem::colLower() const {
  assert(initialized());
  return data_->col_lower;
}

const std::vector<double>& Problem::colUpper() const {
  assert(initialized());
  return data_->col_upper;
}

const std::vector<double>& Problem::rowLower() const {
  assert(initialized());
  return data_->row_lower;
}

const std::vector<double>& Problem::rowUpper() const {
  assert(initialized());
  return data_->row_upper;
}

const std::vector<HighsBool>& Problem::equalityRows() const {
  assert(initialized());
  return data_->equality_rows;
}

const SparseMatrix& Problem::matrix() const {
  assert(initialized());
  return data_->matrix;
}

const SparseMatrix& Problem::transpose() const {
  assert(initialized());
  return data_->transpose;
}

const ScalingMetadata& Problem::scaling() const {
  assert(initialized());
  return data_->scaling;
}

const TransformationMetadata& Problem::transformation() const {
  assert(initialized());
  return data_->transformation;
}

const UnscaledNorms& Problem::unscaledNorms() const {
  assert(initialized());
  return data_->unscaled_norms;
}

api::Error Problem::applyDiagonalScaling(const std::vector<double>& row_scale,
                                         const std::vector<double>& col_scale,
                                         const HighsInt curtis_reid_sweeps) {
  ScalingMetadata metadata;
  metadata.curtis_reid_sweeps = curtis_reid_sweeps;
  if (curtis_reid_sweeps > 0)
    metadata.stages.push_back(ScalingStage::kCurtisReid);
  return applyDiagonalScaling(row_scale, col_scale, metadata);
}

api::Error Problem::applyDiagonalScaling(const std::vector<double>& row_scale,
                                         const std::vector<double>& col_scale,
                                         const ScalingMetadata& metadata) {
  if (!initialized()) return invalidModel("HiPDLP problem is not initialized");
  if (data_->scaling.applied)
    return api::Error(api::ErrorCode::kInvalidOption,
                      "HiPDLP problem scaling is already applied");
  if (row_scale.size() != static_cast<std::size_t>(data_->num_rows) ||
      col_scale.size() != static_cast<std::size_t>(data_->num_cols))
    return invalidModel("HiPDLP scale vectors have the wrong dimensions");
  for (std::size_t row = 0; row < row_scale.size(); ++row) {
    if (!std::isfinite(row_scale[row]) || row_scale[row] <= 0.0)
      return invalidModel("HiPDLP row scale must be finite and positive");
  }
  for (std::size_t col = 0; col < col_scale.size(); ++col) {
    if (!std::isfinite(col_scale[col]) || col_scale[col] <= 0.0)
      return invalidModel("HiPDLP column scale must be finite and positive");
  }

  std::vector<double> scaled_costs = data_->costs;
  std::vector<double> scaled_col_lower = data_->col_lower;
  std::vector<double> scaled_col_upper = data_->col_upper;
  std::vector<double> scaled_row_lower = data_->row_lower;
  std::vector<double> scaled_row_upper = data_->row_upper;
  std::vector<double> scaled_values = data_->matrix.values;
  std::vector<double> scaled_transpose_values = data_->transpose.values;

  for (HighsInt col = 0; col < data_->num_cols; ++col) {
    if (!scaleProduct(data_->costs[col], col_scale[col], 1.0,
                      scaled_costs[col]) ||
        !scaleBound(data_->col_lower[col], col_scale[col], false,
                    scaled_col_lower[col]) ||
        !scaleBound(data_->col_upper[col], col_scale[col], false,
                    scaled_col_upper[col]))
      return invalidModel("HiPDLP column data became invalid while scaling");
  }
  for (HighsInt row = 0; row < data_->num_rows; ++row) {
    if (!scaleBound(data_->row_lower[row], row_scale[row], true,
                    scaled_row_lower[row]) ||
        !scaleBound(data_->row_upper[row], row_scale[row], true,
                    scaled_row_upper[row]))
      return invalidModel("HiPDLP row data became invalid while scaling");
    for (HighsInt entry = data_->matrix.starts[row];
         entry < data_->matrix.starts[row + 1]; ++entry) {
      const HighsInt col = data_->matrix.indices[entry];
      if (!scaleProduct(data_->matrix.values[entry], row_scale[row],
                        col_scale[col], scaled_values[entry]))
        return invalidModel(
            "HiPDLP matrix value became non-finite while "
            "scaling");
    }
  }
  for (HighsInt col = 0; col < data_->num_cols; ++col) {
    for (HighsInt entry = data_->transpose.starts[col];
         entry < data_->transpose.starts[col + 1]; ++entry) {
      const HighsInt row = data_->transpose.indices[entry];
      if (!scaleProduct(data_->transpose.values[entry], row_scale[row],
                        col_scale[col], scaled_transpose_values[entry]))
        return invalidModel(
            "HiPDLP transpose value became non-finite while "
            "scaling");
    }
  }

  data_->costs.swap(scaled_costs);
  data_->col_lower.swap(scaled_col_lower);
  data_->col_upper.swap(scaled_col_upper);
  data_->row_lower.swap(scaled_row_lower);
  data_->row_upper.swap(scaled_row_upper);
  data_->matrix.values.swap(scaled_values);
  data_->transpose.values.swap(scaled_transpose_values);
  data_->scaling = metadata;
  data_->scaling.applied = !metadata.stages.empty();
  data_->scaling.row_scale = row_scale;
  data_->scaling.col_scale = col_scale;
  return api::Error();
}

api::Error Problem::recoverOriginalPrimal(
    const std::vector<double>& transformed_primal,
    std::vector<double>& original_primal) const {
  if (!initialized()) return invalidModel("HiPDLP problem is not initialized");
  api::Error error = validateFiniteInput(transformed_primal, data_->num_cols,
                                         "transformed primal vector");
  if (!error.ok()) return error;
  original_primal.assign(
      transformed_primal.begin(),
      transformed_primal.begin() + data_->transformation.original_num_cols);
  return api::Error();
}

api::Error Problem::recoverOriginalDual(
    const std::vector<double>& transformed_dual,
    std::vector<double>& original_dual) const {
  if (!initialized()) return invalidModel("HiPDLP problem is not initialized");
  api::Error error = validateFiniteInput(transformed_dual, data_->num_rows,
                                         "transformed dual vector");
  if (!error.ok()) return error;
  const TransformationMetadata& transform = data_->transformation;
  original_dual.resize(transform.original_num_rows);
  for (HighsInt row = 0; row < transform.original_num_rows; ++row) {
    original_dual[row] =
        transform.objective_sign * transform.original_row_sign[row] *
        transformed_dual[transform.original_to_transformed_row[row]];
  }
  return api::Error();
}

api::Error Problem::recoverOriginalReducedCosts(
    const std::vector<double>& transformed_reduced_costs,
    std::vector<double>& original_reduced_costs) const {
  if (!initialized()) return invalidModel("HiPDLP problem is not initialized");
  api::Error error =
      validateFiniteInput(transformed_reduced_costs, data_->num_cols,
                          "transformed reduced-cost vector");
  if (!error.ok()) return error;
  const TransformationMetadata& transform = data_->transformation;
  original_reduced_costs.resize(transform.original_num_cols);
  for (HighsInt col = 0; col < transform.original_num_cols; ++col) {
    original_reduced_costs[col] =
        transform.objective_sign * transformed_reduced_costs[col];
  }
  return api::Error();
}

api::Error Problem::recoverPrimalDirection(
    const std::vector<double>& scaled_direction,
    std::vector<double>& original_direction) const {
  if (!initialized()) return invalidModel("HiPDLP problem is not initialized");
  api::Error error = validateFiniteInput(scaled_direction, data_->num_cols,
                                         "scaled primal direction");
  if (!error.ok()) return error;
  if (data_->scaling.col_scale.size() !=
      static_cast<std::size_t>(data_->num_cols))
    return invalidModel("HiPDLP direction recovery scale is invalid");
  const HighsInt original_num_cols =
      data_->transformation.original_num_cols;
  original_direction.assign(original_num_cols, 0.0);
  for (HighsInt col = 0; col < original_num_cols; ++col) {
    if (!scaleProduct(scaled_direction[col], data_->scaling.col_scale[col], 1.0,
                      original_direction[col]))
      return invalidModel("HiPDLP primal direction unscaling is unsafe");
  }
  return api::Error();
}

api::Error Problem::recoverDualDirection(
    const std::vector<double>& scaled_direction,
    std::vector<double>& original_direction) const {
  if (!initialized()) return invalidModel("HiPDLP problem is not initialized");
  api::Error error = validateFiniteInput(scaled_direction, data_->num_rows,
                                         "scaled dual direction");
  if (!error.ok()) return error;
  if (data_->scaling.row_scale.size() !=
      static_cast<std::size_t>(data_->num_rows))
    return invalidModel("HiPDLP direction recovery scale is invalid");
  const TransformationMetadata& transform = data_->transformation;
  original_direction.assign(transform.original_num_rows, 0.0);
  for (HighsInt row = 0; row < transform.original_num_rows; ++row) {
    const HighsInt transformed = transform.original_to_transformed_row[row];
    double unscaled = 0.0;
    if (!scaleProduct(scaled_direction[transformed],
                      data_->scaling.row_scale[transformed], 1.0, unscaled))
      return invalidModel("HiPDLP dual direction unscaling is unsafe");
    original_direction[row] = transform.original_row_sign[row] * unscaled;
    if (!std::isfinite(original_direction[row]))
      return invalidModel("HiPDLP recovered dual direction is non-finite");
  }
  return api::Error();
}

api::Error Problem::recoverReducedCostDirection(
    const std::vector<double>& scaled_direction,
    std::vector<double>& original_direction) const {
  if (!initialized()) return invalidModel("HiPDLP problem is not initialized");
  api::Error error = validateFiniteInput(scaled_direction, data_->num_cols,
                                         "scaled reduced-cost direction");
  if (!error.ok()) return error;
  if (data_->scaling.col_scale.size() !=
      static_cast<std::size_t>(data_->num_cols))
    return invalidModel("HiPDLP direction recovery scale is invalid");
  const HighsInt original_num_cols =
      data_->transformation.original_num_cols;
  original_direction.assign(original_num_cols, 0.0);
  for (HighsInt col = 0; col < original_num_cols; ++col) {
    if (!divideValue(scaled_direction[col], data_->scaling.col_scale[col],
                     original_direction[col]))
      return invalidModel("HiPDLP reduced-cost direction unscaling is unsafe");
  }
  return api::Error();
}

ProblemBuildResult::ProblemBuildResult() {}
ProblemBuildResult::ProblemBuildResult(ProblemBuildResult&& other)
    : error(std::move(other.error)), problem(std::move(other.problem)) {}

ProblemBuildResult& ProblemBuildResult::operator=(ProblemBuildResult&& other) {
  error = std::move(other.error);
  problem = std::move(other.problem);
  return *this;
}

bool ProblemBuildResult::ok() const {
  return error.ok() && problem.initialized();
}

PostsolveResult::PostsolveResult() : objective_value(0.0) {}

bool PostsolveResult::ok() const {
  return error.ok() && solution.value_valid && solution.dual_valid;
}

PostsolveResult postsolve(const Problem& problem,
                          const std::vector<double>& scaled_primal,
                          const std::vector<double>& scaled_dual,
                          const std::vector<double>& scaled_reduced_costs) {
  PostsolveResult result;
  if (!problem.initialized()) {
    result.error = invalidModel("HiPDLP problem is not initialized");
    return result;
  }
  api::Error error =
      validateFiniteInput(scaled_primal, problem.numCols(), "scaled primal");
  if (!error.ok()) {
    result.error = error;
    return result;
  }
  error = validateFiniteInput(scaled_dual, problem.numRows(), "scaled dual");
  if (!error.ok()) {
    result.error = error;
    return result;
  }
  error = validateFiniteInput(scaled_reduced_costs, problem.numCols(),
                              "scaled reduced-cost vector");
  if (!error.ok()) {
    result.error = error;
    return result;
  }

  const ScalingMetadata& scaling = problem.scaling();
  if (scaling.row_scale.size() != static_cast<std::size_t>(problem.numRows()) ||
      scaling.col_scale.size() != static_cast<std::size_t>(problem.numCols())) {
    result.error =
        invalidModel("HiPDLP postsolve scale dimensions are invalid");
    return result;
  }
  std::vector<double> transformed_primal(problem.numCols());
  std::vector<double> transformed_dual(problem.numRows());
  std::vector<double> transformed_reduced_costs(problem.numCols());
  for (HighsInt col = 0; col < problem.numCols(); ++col) {
    if (!scaleProduct(scaled_primal[col], scaling.col_scale[col], 1.0,
                      transformed_primal[col]) ||
        !divideValue(scaled_reduced_costs[col], scaling.col_scale[col],
                     transformed_reduced_costs[col])) {
      result.error =
          invalidModel("HiPDLP column postsolve unscaling is unsafe");
      return result;
    }
  }
  for (HighsInt row = 0; row < problem.numRows(); ++row) {
    if (!scaleProduct(scaled_dual[row], scaling.row_scale[row], 1.0,
                      transformed_dual[row])) {
      result.error = invalidModel("HiPDLP row postsolve unscaling is unsafe");
      return result;
    }
  }

  HighsSolution solution;
  error = problem.recoverOriginalPrimal(transformed_primal, solution.col_value);
  if (!error.ok()) {
    result.error = error;
    return result;
  }
  error = problem.recoverOriginalDual(transformed_dual, solution.row_dual);
  if (!error.ok()) {
    result.error = error;
    return result;
  }
  error = problem.recoverOriginalReducedCosts(transformed_reduced_costs,
                                              solution.col_dual);
  if (!error.ok()) {
    result.error = error;
    return result;
  }

  const TransformationMetadata& transform = problem.transformation();
  solution.row_value.assign(transform.original_num_rows, 0.0);
  const SparseMatrix& matrix = problem.matrix();
  for (HighsInt transformed_row = 0; transformed_row < problem.numRows();
       ++transformed_row) {
    double activity = 0.0;
    for (HighsInt entry = matrix.starts[transformed_row];
         entry < matrix.starts[transformed_row + 1]; ++entry) {
      const HighsInt col = matrix.indices[entry];
      if (col >= transform.original_num_cols) continue;
      double scaled_product = 0.0;
      double contribution = 0.0;
      if (!scaleProduct(matrix.values[entry], scaled_primal[col], 1.0,
                        scaled_product) ||
          !divideValue(scaled_product, scaling.row_scale[transformed_row],
                       contribution) ||
          !std::isfinite(activity + contribution)) {
        result.error = invalidModel("HiPDLP row activity postsolve is unsafe");
        return result;
      }
      activity += contribution;
    }
    const HighsInt original_row =
        transform.transformed_to_original_row[transformed_row];
    solution.row_value[original_row] =
        transform.original_row_sign[original_row] * activity;
    if (!std::isfinite(solution.row_value[original_row])) {
      result.error =
          invalidModel("HiPDLP recovered row activity is non-finite");
      return result;
    }
  }

  double objective = transform.original_offset;
  for (HighsInt col = 0; col < transform.original_num_cols; ++col) {
    double transformed_cost = 0.0;
    if (!divideValue(problem.costs()[col], scaling.col_scale[col],
                     transformed_cost) ||
        !addProduct(transform.objective_sign * transformed_cost,
                    solution.col_value[col], objective)) {
      result.error = invalidModel("HiPDLP recovered objective is unsafe");
      return result;
    }
  }

  solution.value_valid = true;
  solution.dual_valid = true;
  result.solution = solution;
  result.objective_value = objective;
  return result;
}

ProblemBuildResult buildProblem(const HighsLp& lp) {
  ProblemBuildResult result;
  result.error = validateLp(lp);
  if (!result.error.ok()) return result;

  std::vector<double> original_col_lower;
  std::vector<double> original_col_upper;
  result.error =
      copyAndNormalizeBounds(lp.col_lower_, lp.col_upper_, original_col_lower,
                             original_col_upper, "column");
  if (!result.error.ok()) return result;
  std::vector<double> original_row_lower;
  std::vector<double> original_row_upper;
  result.error =
      copyAndNormalizeBounds(lp.row_lower_, lp.row_upper_, original_row_lower,
                             original_row_upper, "row");
  if (!result.error.ok()) return result;

  std::vector<MatrixEntry> original_entries;
  original_entries.reserve(lp.a_matrix_.value_.size());
  if (lp.a_matrix_.format_ == MatrixFormat::kColwise) {
    for (HighsInt col = 0; col < lp.num_col_; ++col) {
      for (HighsInt entry = lp.a_matrix_.start_[col];
           entry < lp.a_matrix_.start_[col + 1]; ++entry) {
        original_entries.push_back(
            {lp.a_matrix_.index_[entry], col, lp.a_matrix_.value_[entry]});
      }
    }
  } else {
    for (HighsInt row = 0; row < lp.num_row_; ++row) {
      for (HighsInt entry = lp.a_matrix_.start_[row];
           entry < lp.a_matrix_.start_[row + 1]; ++entry) {
        original_entries.push_back(
            {row, lp.a_matrix_.index_[entry], lp.a_matrix_.value_[entry]});
      }
    }
  }
  result.error =
      canonicalizeEntries(lp.num_row_, lp.num_col_, original_entries);
  if (!result.error.ok()) return result;

  std::unique_ptr<Problem::Data> data(new Problem::Data());
  TransformationMetadata& transform = data->transformation;
  transform.original_num_rows = lp.num_row_;
  transform.original_num_cols = lp.num_col_;
  transform.original_sense = lp.sense_;
  transform.original_offset = lp.offset_;
  transform.objective_sign = lp.sense_ == ObjSense::kMinimize ? 1.0 : -1.0;
  transform.original_row_type.resize(lp.num_row_);
  transform.original_to_transformed_row.resize(lp.num_row_);
  transform.transformed_to_original_row.resize(lp.num_row_);
  transform.original_row_sign.assign(lp.num_row_, 1.0);
  transform.slack_col_for_original_row.assign(lp.num_row_, -1);

  HighsInt num_equalities = 0;
  HighsInt num_slacks = 0;
  for (HighsInt row = 0; row < lp.num_row_; ++row) {
    const bool has_lower = !isNegativeInfinity(original_row_lower[row]);
    const bool has_upper = !isPositiveInfinity(original_row_upper[row]);
    if (has_lower && has_upper &&
        original_row_lower[row] == original_row_upper[row]) {
      transform.original_row_type[row] = RowType::kEquality;
      ++num_equalities;
    } else if (has_lower && has_upper) {
      transform.original_row_type[row] = RowType::kBoxed;
      ++num_equalities;
      ++num_slacks;
    } else if (has_lower) {
      transform.original_row_type[row] = RowType::kLower;
    } else if (has_upper) {
      transform.original_row_type[row] = RowType::kUpper;
      transform.original_row_sign[row] = -1.0;
    } else {
      transform.original_row_type[row] = RowType::kFree;
      ++num_equalities;
      ++num_slacks;
    }
  }

  HighsInt equality_position = 0;
  HighsInt inequality_position = num_equalities;
  for (HighsInt row = 0; row < lp.num_row_; ++row) {
    const RowType type = transform.original_row_type[row];
    const bool equality = type == RowType::kEquality ||
                          type == RowType::kBoxed || type == RowType::kFree;
    const HighsInt transformed_row =
        equality ? equality_position++ : inequality_position++;
    transform.original_to_transformed_row[row] = transformed_row;
    transform.transformed_to_original_row[transformed_row] = row;
  }

  data->num_rows = lp.num_row_;
  data->num_cols = lp.num_col_ + num_slacks;
  data->objective_sense = ObjSense::kMinimize;
  data->objective_offset = transform.objective_sign * lp.offset_;
  data->costs.assign(data->num_cols, 0.0);
  data->col_lower.assign(data->num_cols, -kHighsInf);
  data->col_upper.assign(data->num_cols, kHighsInf);
  for (HighsInt col = 0; col < lp.num_col_; ++col) {
    data->costs[col] = transform.objective_sign * lp.col_cost_[col];
    data->col_lower[col] = original_col_lower[col];
    data->col_upper[col] = original_col_upper[col];
  }
  data->row_lower.resize(data->num_rows);
  data->row_upper.resize(data->num_rows);
  data->equality_rows.assign(data->num_rows, false);

  HighsInt next_slack = lp.num_col_;
  for (HighsInt original_row = 0; original_row < lp.num_row_; ++original_row) {
    const HighsInt row = transform.original_to_transformed_row[original_row];
    const RowType type = transform.original_row_type[original_row];
    switch (type) {
      case RowType::kEquality:
        data->row_lower[row] = original_row_lower[original_row];
        data->row_upper[row] = original_row_upper[original_row];
        data->equality_rows[row] = true;
        break;
      case RowType::kBoxed:
      case RowType::kFree:
        data->row_lower[row] = 0.0;
        data->row_upper[row] = 0.0;
        data->equality_rows[row] = true;
        transform.slack_col_for_original_row[original_row] = next_slack;
        data->col_lower[next_slack] = original_row_lower[original_row];
        data->col_upper[next_slack] = original_row_upper[original_row];
        ++next_slack;
        break;
      case RowType::kLower:
        data->row_lower[row] = original_row_lower[original_row];
        data->row_upper[row] = kHighsInf;
        break;
      case RowType::kUpper:
        data->row_lower[row] = -original_row_upper[original_row];
        data->row_upper[row] = kHighsInf;
        break;
    }
  }

  std::vector<MatrixEntry> transformed_entries;
  transformed_entries.reserve(original_entries.size() + num_slacks);
  for (std::size_t entry = 0; entry < original_entries.size(); ++entry) {
    const HighsInt original_row = original_entries[entry].row;
    transformed_entries.push_back(
        {transform.original_to_transformed_row[original_row],
         original_entries[entry].col,
         transform.original_row_sign[original_row] *
             original_entries[entry].value});
  }
  for (HighsInt original_row = 0; original_row < lp.num_row_; ++original_row) {
    const HighsInt slack = transform.slack_col_for_original_row[original_row];
    if (slack >= 0) {
      transformed_entries.push_back(
          {transform.original_to_transformed_row[original_row], slack, -1.0});
    }
  }
  result.error =
      canonicalizeEntries(data->num_rows, data->num_cols, transformed_entries);
  if (!result.error.ok()) return result;
  data->matrix =
      makeCompressedMatrix(data->num_rows, data->num_cols, transformed_entries);

  std::vector<MatrixEntry> transpose_entries;
  transpose_entries.reserve(transformed_entries.size());
  for (std::size_t entry = 0; entry < transformed_entries.size(); ++entry) {
    transpose_entries.push_back({transformed_entries[entry].col,
                                 transformed_entries[entry].row,
                                 transformed_entries[entry].value});
  }
  result.error =
      canonicalizeEntries(data->num_cols, data->num_rows, transpose_entries);
  if (!result.error.ok()) return result;
  data->transpose =
      makeCompressedMatrix(data->num_cols, data->num_rows, transpose_entries);

  data->scaling.row_scale.assign(data->num_rows, 1.0);
  data->scaling.col_scale.assign(data->num_cols, 1.0);
  bool norms_valid = true;
  data->unscaled_norms.objective_l2 = stableNorm(data->costs, norms_valid);
  data->unscaled_norms.row_bounds_l2 = stableNorm(data->row_lower, norms_valid);
  data->unscaled_norms.matrix_frobenius =
      stableNorm(data->matrix.values, norms_valid);
  if (!norms_valid) {
    result.error = invalidModel("HiPDLP unscaled norm is non-finite");
    return result;
  }

  result.problem.data_ = std::move(data);
  return result;
}

}  // namespace model
}  // namespace hipdlp
