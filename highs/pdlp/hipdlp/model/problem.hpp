#ifndef PDLP_HIPDLP_MODEL_PROBLEM_HPP_
#define PDLP_HIPDLP_MODEL_PROBLEM_HPP_

#include <memory>
#include <vector>

#include "lp_data/HighsLp.h"
#include "lp_data/HighsSolution.h"
#include "pdlp/hipdlp/api/status.hpp"
#include "pdlp/hipdlp/model/preprocessing.hpp"

namespace hipdlp {
namespace model {

enum class RowType { kEquality, kLower, kUpper, kBoxed, kFree };

enum class IndexWidth { kInt32, kInt64 };

enum class ScalingStage { kCurtisReid, kRuiz, kPockChambolle, kL2 };

struct SparseMatrix {
  HighsInt num_rows;
  HighsInt num_cols;
  std::vector<HighsInt> starts;
  std::vector<HighsInt> indices;
  std::vector<double> values;

  SparseMatrix();

  HighsInt numNonzeros() const;
  bool samePattern(const SparseMatrix& other) const;
};

struct ScalingMetadata {
  bool applied;
  // True when auto rejected CR's raw log-scale spread before applying it.
  bool curtis_reid_rejected;
  double curtis_reid_row_log_range;
  double curtis_reid_col_log_range;
  HighsInt curtis_reid_sweeps;
  HighsInt ruiz_iterations;
  std::vector<ScalingStage> stages;
  std::vector<double> row_scale;
  std::vector<double> col_scale;

  ScalingMetadata();
};

struct TransformationMetadata {
  HighsInt original_num_rows;
  HighsInt original_num_cols;
  ObjSense original_sense;
  double original_offset;
  double objective_sign;
  std::vector<RowType> original_row_type;
  std::vector<HighsInt> original_to_transformed_row;
  std::vector<HighsInt> transformed_to_original_row;
  std::vector<double> original_row_sign;
  std::vector<HighsInt> slack_col_for_original_row;

  TransformationMetadata();
};

struct UnscaledNorms {
  double objective_l2;
  double row_bounds_l2;
  double matrix_frobenius;

  UnscaledNorms();
};

struct ProblemBuildResult;

struct PostsolveResult {
  api::Error error;
  HighsSolution solution;
  double objective_value;

  PostsolveResult();
  bool ok() const;
};

class OriginalModel {
 public:
  explicit OriginalModel(const HighsLp& lp);

  const HighsLp& lp() const;

 private:
  const HighsLp& lp_;
};

class Problem {
 public:
  Problem();
  ~Problem();

  Problem(Problem&& other);
  Problem& operator=(Problem&& other);

  Problem(const Problem&) = delete;
  Problem& operator=(const Problem&) = delete;

  bool initialized() const;

  HighsInt numRows() const;
  HighsInt numCols() const;
  ObjSense objectiveSense() const;
  double objectiveOffset() const;
  IndexWidth indexWidth() const;

  const std::vector<double>& costs() const;
  const std::vector<double>& colLower() const;
  const std::vector<double>& colUpper() const;
  const std::vector<double>& rowLower() const;
  const std::vector<double>& rowUpper() const;
  const std::vector<HighsBool>& equalityRows() const;
  const SparseMatrix& matrix() const;
  const SparseMatrix& transpose() const;
  const ScalingMetadata& scaling() const;
  const TransformationMetadata& transformation() const;
  const UnscaledNorms& unscaledNorms() const;

  api::Error applyDiagonalScaling(const std::vector<double>& row_scale,
                                  const std::vector<double>& col_scale,
                                  HighsInt curtis_reid_sweeps);
  api::Error applyDiagonalScaling(const std::vector<double>& row_scale,
                                  const std::vector<double>& col_scale,
                                  const ScalingMetadata& metadata);

  api::Error recoverOriginalPrimal(
      const std::vector<double>& transformed_primal,
      std::vector<double>& original_primal) const;
  api::Error recoverOriginalDual(const std::vector<double>& transformed_dual,
                                 std::vector<double>& original_dual) const;
  api::Error recoverOriginalReducedCosts(
      const std::vector<double>& transformed_reduced_costs,
      std::vector<double>& original_reduced_costs) const;

  // Direction recovery. Unlike the solution recovery above these map a scaled
  // *direction* back to the original space: they undo the diagonal scaling and
  // the row permutation/negation, drop synthetic slack columns, and never
  // apply the objective offset or the maximization sign flip, because a
  // recession or Farkas direction has no affine part and no objective sense.
  api::Error recoverPrimalDirection(
      const std::vector<double>& scaled_direction,
      std::vector<double>& original_direction) const;
  api::Error recoverDualDirection(
      const std::vector<double>& scaled_direction,
      std::vector<double>& original_direction) const;
  api::Error recoverReducedCostDirection(
      const std::vector<double>& scaled_direction,
      std::vector<double>& original_direction) const;

 private:
  struct Data;
  std::unique_ptr<Data> data_;

  friend ProblemBuildResult buildProblem(const HighsLp& lp,
                                         PreprocessingDiagnostics* diagnostics);
};

struct ProblemBuildResult {
  api::Error error;
  Problem problem;

  ProblemBuildResult();
  ProblemBuildResult(ProblemBuildResult&& other);
  ProblemBuildResult& operator=(ProblemBuildResult&& other);

  ProblemBuildResult(const ProblemBuildResult&) = delete;
  ProblemBuildResult& operator=(const ProblemBuildResult&) = delete;

  bool ok() const;
};

ProblemBuildResult buildProblem(const HighsLp& lp,
                                PreprocessingDiagnostics* diagnostics = nullptr);
PostsolveResult postsolve(const Problem& problem,
                          const std::vector<double>& scaled_primal,
                          const std::vector<double>& scaled_dual,
                          const std::vector<double>& scaled_reduced_costs);

}  // namespace model
}  // namespace hipdlp

#endif
