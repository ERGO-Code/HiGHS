#include <cmath>
#include <limits>
#include <type_traits>
#include <vector>

#include "catch.hpp"
#include "lp_data/HConst.h"
#include "pdlp/hipdlp/model/problem.hpp"
#include "pdlp/hipdlp/model/scaling.hpp"

namespace {

HighsLp makeTransformationLp() {
  HighsLp lp;
  lp.num_col_ = 3;
  lp.num_row_ = 5;
  lp.col_cost_ = {1.0, -2.0, 3.0};
  lp.col_lower_ = {2.0, -kHighsInf, 0.0};
  lp.col_upper_ = {2.0, kHighsInf, 5.0};
  lp.row_lower_ = {1.0, 2.0, -kHighsInf, -1.0, -kHighsInf};
  lp.row_upper_ = {1.0, kHighsInf, 3.0, 4.0, kHighsInf};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = lp.num_col_;
  lp.a_matrix_.num_row_ = lp.num_row_;
  lp.a_matrix_.start_ = {0, 4, 6, 9};
  lp.a_matrix_.index_ = {0, 1, 1, 2, 3, 4, 2, 2, 4};
  lp.a_matrix_.value_ = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, -7.0, -2.0};
  lp.sense_ = ObjSense::kMaximize;
  lp.offset_ = 7.0;
  return lp;
}

hipdlp::model::SparseMatrix makeCurtisReidFixture() {
  hipdlp::model::SparseMatrix matrix;
  matrix.num_rows = 2;
  matrix.num_cols = 2;
  matrix.starts = {0, 2, 4};
  matrix.indices = {0, 1, 0, 1};
  matrix.values = {1.0, 2.0, 4.0, 8.0};
  return matrix;
}

hipdlp::model::SparseMatrix makeDiagonalScalingFixture() {
  hipdlp::model::SparseMatrix matrix;
  matrix.num_rows = 3;
  matrix.num_cols = 3;
  matrix.starts = {0, 1, 2, 2};
  matrix.indices = {0, 1};
  matrix.values = {4.0, 9.0};
  return matrix;
}

HighsLp makeUnsafeScalingLp() {
  HighsLp lp;
  lp.num_col_ = 1;
  lp.num_row_ = 1;
  lp.col_cost_ = {std::numeric_limits<double>::max()};
  lp.col_lower_ = {-kHighsInf};
  lp.col_upper_ = {kHighsInf};
  lp.row_lower_ = {0.0};
  lp.row_upper_ = {0.0};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 1;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1};
  lp.a_matrix_.index_ = {0};
  lp.a_matrix_.value_ = {1e-300};
  lp.sense_ = ObjSense::kMinimize;
  lp.offset_ = 0.0;
  return lp;
}

void checkVectorApprox(const std::vector<double>& actual,
                       const std::vector<double>& expected) {
  REQUIRE(actual.size() == expected.size());
  for (std::size_t index = 0; index < actual.size(); ++index)
    CHECK(actual[index] == Approx(expected[index]).epsilon(1e-14));
}

}  // namespace

TEST_CASE("hi-pdlp-owned-problem-transformation",
          "[pdlp][hipdlp][phase3][model]") {
  HighsLp lp = makeTransformationLp();
  const HighsLp original = lp;

  hipdlp::model::ProblemBuildResult built = hipdlp::model::buildProblem(lp);
  INFO(built.error.message);
  REQUIRE(built.ok());
  const hipdlp::model::Problem& problem = built.problem;

  CHECK(lp.num_col_ == original.num_col_);
  CHECK(lp.num_row_ == original.num_row_);
  CHECK(lp.col_cost_ == original.col_cost_);
  CHECK(lp.col_lower_ == original.col_lower_);
  CHECK(lp.col_upper_ == original.col_upper_);
  CHECK(lp.row_lower_ == original.row_lower_);
  CHECK(lp.row_upper_ == original.row_upper_);
  CHECK(lp.a_matrix_ == original.a_matrix_);

  CHECK(problem.numRows() == 5);
  CHECK(problem.numCols() == 5);
  CHECK(problem.objectiveSense() == ObjSense::kMinimize);
  CHECK(problem.objectiveOffset() == -7.0);
  CHECK(problem.costs() == std::vector<double>({-1.0, 2.0, -3.0, 0.0, 0.0}));
  CHECK(problem.colLower() ==
        std::vector<double>({2.0, -kHighsInf, 0.0, -1.0, -kHighsInf}));
  CHECK(problem.colUpper() ==
        std::vector<double>({2.0, kHighsInf, 5.0, 4.0, kHighsInf}));
  CHECK(problem.rowLower() == std::vector<double>({1.0, 0.0, 0.0, 2.0, -3.0}));
  CHECK(problem.rowUpper() ==
        std::vector<double>({1.0, 0.0, 0.0, kHighsInf, kHighsInf}));
  CHECK(problem.equalityRows() ==
        std::vector<HighsBool>({true, true, true, false, false}));

  const hipdlp::model::SparseMatrix& matrix = problem.matrix();
  CHECK(matrix.starts == std::vector<HighsInt>({0, 1, 3, 6, 7, 8}));
  CHECK(matrix.indices == std::vector<HighsInt>({0, 1, 3, 1, 2, 4, 0, 0}));
  CHECK(matrix.values ==
        std::vector<double>({1.0, 5.0, -1.0, 6.0, -2.0, -1.0, 5.0, -4.0}));
  CHECK(problem.transpose().num_rows == problem.numCols());
  CHECK(problem.transpose().num_cols == problem.numRows());
  CHECK(problem.transpose().numNonzeros() == matrix.numNonzeros());

  HighsLp rowwise_lp = lp;
  rowwise_lp.a_matrix_.ensureRowwise();
  hipdlp::model::ProblemBuildResult rowwise_built =
      hipdlp::model::buildProblem(rowwise_lp);
  REQUIRE(rowwise_built.ok());
  CHECK(rowwise_built.problem.matrix().starts == matrix.starts);
  CHECK(rowwise_built.problem.matrix().indices == matrix.indices);
  CHECK(rowwise_built.problem.matrix().values == matrix.values);

  const hipdlp::model::TransformationMetadata& transform =
      problem.transformation();
  CHECK(transform.original_sense == ObjSense::kMaximize);
  CHECK(transform.original_offset == 7.0);
  CHECK(transform.objective_sign == -1.0);
  CHECK(transform.original_to_transformed_row ==
        std::vector<HighsInt>({0, 3, 4, 1, 2}));
  CHECK(transform.transformed_to_original_row ==
        std::vector<HighsInt>({0, 3, 4, 1, 2}));
  CHECK(transform.original_row_sign ==
        std::vector<double>({1.0, 1.0, -1.0, 1.0, 1.0}));
  CHECK(transform.slack_col_for_original_row ==
        std::vector<HighsInt>({-1, -1, -1, 3, 4}));
  CHECK(transform.original_row_type[0] == hipdlp::model::RowType::kEquality);
  CHECK(transform.original_row_type[1] == hipdlp::model::RowType::kLower);
  CHECK(transform.original_row_type[2] == hipdlp::model::RowType::kUpper);
  CHECK(transform.original_row_type[3] == hipdlp::model::RowType::kBoxed);
  CHECK(transform.original_row_type[4] == hipdlp::model::RowType::kFree);
  CHECK(problem.indexWidth() == (sizeof(HighsInt) > 4
                                     ? hipdlp::model::IndexWidth::kInt64
                                     : hipdlp::model::IndexWidth::kInt32));

  std::vector<double> original_primal;
  REQUIRE(problem
              .recoverOriginalPrimal({10.0, 20.0, 30.0, 40.0, 50.0},
                                     original_primal)
              .ok());
  CHECK(original_primal == std::vector<double>({10.0, 20.0, 30.0}));
  std::vector<double> original_dual;
  REQUIRE(problem.recoverOriginalDual({1.0, 2.0, 3.0, 4.0, 5.0}, original_dual)
              .ok());
  CHECK(original_dual == std::vector<double>({-1.0, -4.0, 5.0, -2.0, -3.0}));
  std::vector<double> original_reduced_costs;
  REQUIRE(problem
              .recoverOriginalReducedCosts({1.0, 2.0, 3.0, 4.0, 5.0},
                                           original_reduced_costs)
              .ok());
  CHECK(original_reduced_costs == std::vector<double>({-1.0, -2.0, -3.0}));
}

TEST_CASE("hi-pdlp-curtis-reid-fixed-scales",
          "[pdlp][hipdlp][phase3][model][scaling]") {
  const hipdlp::model::SparseMatrix matrix = makeCurtisReidFixture();
  const hipdlp::model::CurtisReidResult first =
      hipdlp::model::computeCurtisReidScaling(matrix, 20);
  INFO(first.error.message);
  REQUIRE(first.ok());
  checkVectorApprox(first.row_scale, {0.7071067811865476, 0.1767766952966369});
  checkVectorApprox(first.col_scale, {1.4142135623730951, 0.7071067811865476});

  const hipdlp::model::CurtisReidResult second =
      hipdlp::model::computeCurtisReidScaling(matrix, 20);
  REQUIRE(second.ok());
  CHECK(second.row_scale == first.row_scale);
  CHECK(second.col_scale == first.col_scale);
}

TEST_CASE("hi-pdlp-scaling-pipeline-fixed-vectors",
          "[pdlp][hipdlp][phase3][phase3b][model][scaling]") {
  const hipdlp::model::SparseMatrix matrix = makeCurtisReidFixture();

  const hipdlp::model::ScalingResult ruiz =
      hipdlp::model::computeRuizScaling(matrix, 1);
  REQUIRE(ruiz.ok());
  checkVectorApprox(ruiz.scaling.row_scale,
                    {0.70710678118654746, 0.35355339059327373});
  checkVectorApprox(ruiz.scaling.col_scale, {0.5, 0.35355339059327373});

  const hipdlp::model::ScalingResult pock_chambolle =
      hipdlp::model::computePockChambolleScaling(matrix);
  REQUIRE(pock_chambolle.ok());
  checkVectorApprox(pock_chambolle.scaling.row_scale,
                    {0.57735026918962584, 0.28867513459481292});
  checkVectorApprox(pock_chambolle.scaling.col_scale,
                    {0.44721359549995793, 0.31622776601683794});

  const hipdlp::model::ScalingResult l2 =
      hipdlp::model::computeL2Scaling(matrix);
  REQUIRE(l2.ok());
  checkVectorApprox(l2.scaling.row_scale,
                    {0.66874030497642201, 0.334370152488211});
  checkVectorApprox(l2.scaling.col_scale,
                    {0.49247906050545237, 0.3482352832757854});

  hipdlp::model::ScalingOptions legacy_options;
  legacy_options.mode = hipdlp::model::ScalingMode::kLegacy;
  legacy_options.legacy_scaling_mode = kPdlpScalingMax;
  legacy_options.ruiz_iterations = 1;
  const hipdlp::model::ScalingResult legacy =
      hipdlp::model::computeScalingPipeline(matrix, legacy_options);
  REQUIRE(legacy.ok());
  checkVectorApprox(legacy.scaling.row_scale,
                    {1.0072794410463148, 0.29946596952206039});
  checkVectorApprox(legacy.scaling.col_scale,
                    {0.60516777658976795, 0.32997006996549161});
  CHECK(legacy.scaling.stages ==
        std::vector<hipdlp::model::ScalingStage>(
            {hipdlp::model::ScalingStage::kRuiz,
             hipdlp::model::ScalingStage::kPockChambolle,
             hipdlp::model::ScalingStage::kL2}));

  hipdlp::model::ScalingOptions auto_options = legacy_options;
  auto_options.mode = hipdlp::model::ScalingMode::kAuto;
  const hipdlp::model::ScalingResult automatic =
      hipdlp::model::computeScalingPipeline(matrix, auto_options);
  REQUIRE(automatic.ok());
  checkVectorApprox(automatic.scaling.row_scale,
                    {0.59460355750136062, 0.14865088937534016});
  checkVectorApprox(automatic.scaling.col_scale,
                    {1.1892071150027208, 0.5946035575013604});
  CHECK(automatic.scaling.curtis_reid_sweeps == 20);
  CHECK(automatic.scaling.stages ==
        std::vector<hipdlp::model::ScalingStage>(
            {hipdlp::model::ScalingStage::kCurtisReid,
             hipdlp::model::ScalingStage::kRuiz,
             hipdlp::model::ScalingStage::kPockChambolle,
             hipdlp::model::ScalingStage::kL2}));

  const hipdlp::model::SparseMatrix diagonal = makeDiagonalScalingFixture();
  const std::vector<double> split_scale = {0.5, 1.0 / 3.0, 1.0};
  const hipdlp::model::ScalingResult legacy_first =
      hipdlp::model::computeScalingPipeline(diagonal, legacy_options);
  REQUIRE(legacy_first.ok());
  checkVectorApprox(legacy_first.scaling.row_scale, split_scale);
  checkVectorApprox(legacy_first.scaling.col_scale, split_scale);
  hipdlp::model::SparseMatrix legacy_scaled = diagonal;
  for (HighsInt row = 0; row < legacy_scaled.num_rows; ++row) {
    for (HighsInt entry = legacy_scaled.starts[row];
         entry < legacy_scaled.starts[row + 1]; ++entry) {
      const HighsInt col = legacy_scaled.indices[entry];
      legacy_scaled.values[entry] *= legacy_first.scaling.row_scale[row] *
                                     legacy_first.scaling.col_scale[col];
    }
  }
  const hipdlp::model::CurtisReidResult curtis_reid_last =
      hipdlp::model::computeCurtisReidScaling(legacy_scaled, 20);
  REQUIRE(curtis_reid_last.ok());
  std::vector<double> reversed_rows(3);
  std::vector<double> reversed_cols(3);
  for (HighsInt index = 0; index < 3; ++index) {
    reversed_rows[index] = legacy_first.scaling.row_scale[index] *
                           curtis_reid_last.row_scale[index];
    reversed_cols[index] = legacy_first.scaling.col_scale[index] *
                           curtis_reid_last.col_scale[index];
  }
  const hipdlp::model::ScalingResult auto_diagonal =
      hipdlp::model::computeScalingPipeline(diagonal, auto_options);
  REQUIRE(auto_diagonal.ok());
  checkVectorApprox(auto_diagonal.scaling.row_scale, {0.25, 1.0 / 9.0, 1.0});
  checkVectorApprox(auto_diagonal.scaling.col_scale, {1.0, 1.0, 1.0});
  CHECK(auto_diagonal.scaling.row_scale != reversed_rows);
  CHECK(auto_diagonal.scaling.col_scale != reversed_cols);
}

TEST_CASE("hi-pdlp-scaling-modes-determinism-and-pattern",
          "[pdlp][hipdlp][phase3][phase3b][model][scaling]") {
  struct ModeCase {
    hipdlp::model::ScalingMode mode;
    HighsInt features_off;
    HighsInt legacy_scaling_mode;
    std::vector<hipdlp::model::ScalingStage> stages;
  };
  const std::vector<ModeCase> cases = {
      {hipdlp::model::ScalingMode::kAuto,
       kPdlpAllFeaturesOn,
       kPdlpScalingMax,
       {hipdlp::model::ScalingStage::kCurtisReid,
        hipdlp::model::ScalingStage::kRuiz,
        hipdlp::model::ScalingStage::kPockChambolle,
        hipdlp::model::ScalingStage::kL2}},
      {hipdlp::model::ScalingMode::kOff,
       kPdlpAllFeaturesOn,
       kPdlpScalingMax,
       {}},
      {hipdlp::model::ScalingMode::kCurtisReid,
       kPdlpAllFeaturesOn,
       kPdlpScalingMax,
       {hipdlp::model::ScalingStage::kCurtisReid}},
      {hipdlp::model::ScalingMode::kLegacy,
       kPdlpAllFeaturesOn,
       kPdlpScalingRuiz,
       {hipdlp::model::ScalingStage::kRuiz}},
      {hipdlp::model::ScalingMode::kLegacy,
       kPdlpAllFeaturesOn,
       kPdlpScalingPC,
       {hipdlp::model::ScalingStage::kPockChambolle}},
      {hipdlp::model::ScalingMode::kLegacy,
       kPdlpAllFeaturesOn,
       kPdlpScalingL2,
       {hipdlp::model::ScalingStage::kL2}},
      {hipdlp::model::ScalingMode::kLegacy,
       kPdlpAllFeaturesOn,
       kPdlpScalingMax,
       {hipdlp::model::ScalingStage::kRuiz,
        hipdlp::model::ScalingStage::kPockChambolle,
        hipdlp::model::ScalingStage::kL2}},
      {hipdlp::model::ScalingMode::kAuto,
       kPdlpScalingOff,
       kPdlpScalingMax,
       {}}};

  for (const ModeCase& test_case : cases) {
    CAPTURE(static_cast<int>(test_case.mode), test_case.features_off);
    hipdlp::model::ProblemBuildResult built =
        hipdlp::model::buildProblem(makeTransformationLp());
    REQUIRE(built.ok());
    const hipdlp::model::SparseMatrix original_matrix = built.problem.matrix();
    const hipdlp::model::SparseMatrix original_transpose =
        built.problem.transpose();
    hipdlp::model::ScalingOptions options;
    options.mode = test_case.mode;
    options.features_off = test_case.features_off;
    options.legacy_scaling_mode = test_case.legacy_scaling_mode;
    options.ruiz_iterations = 2;
    const hipdlp::model::ScalingResult first =
        hipdlp::model::computeScalingPipeline(built.problem.matrix(), options);
    const hipdlp::model::ScalingResult second =
        hipdlp::model::computeScalingPipeline(built.problem.matrix(), options);
    REQUIRE(first.ok());
    REQUIRE(second.ok());
    CHECK(first.scaling.row_scale == second.scaling.row_scale);
    CHECK(first.scaling.col_scale == second.scaling.col_scale);
    CHECK(first.scaling.stages == test_case.stages);
    CHECK(
        first.scaling.curtis_reid_sweeps ==
        (test_case.stages.empty() ||
                 test_case.stages[0] != hipdlp::model::ScalingStage::kCurtisReid
             ? 0
             : 20));

    REQUIRE(hipdlp::model::applyScalingPipeline(built.problem, options).ok());
    CHECK(built.problem.matrix().samePattern(original_matrix));
    CHECK(built.problem.transpose().samePattern(original_transpose));
    CHECK(built.problem.matrix().starts == original_matrix.starts);
    CHECK(built.problem.matrix().indices == original_matrix.indices);
    CHECK(built.problem.scaling().stages == test_case.stages);
    CHECK(built.problem.scaling().applied == !test_case.stages.empty());
    if (test_case.stages.empty()) {
      CHECK(built.problem.matrix().values == original_matrix.values);
      CHECK(built.problem.transpose().values == original_transpose.values);
    }
    CHECK(built.problem.colLower()[1] == -kHighsInf);
    CHECK(built.problem.colUpper()[1] == kHighsInf);
    CHECK(built.problem.rowUpper()[3] == kHighsInf);
  }
}

TEST_CASE("hi-pdlp-scaling-validation-and-transactionality",
          "[pdlp][hipdlp][phase3][phase3b][model][scaling]") {
  const hipdlp::model::SparseMatrix matrix = makeCurtisReidFixture();
  hipdlp::model::ScalingOptions invalid;
  invalid.mode = static_cast<hipdlp::model::ScalingMode>(-1);
  CHECK_FALSE(hipdlp::model::computeScalingPipeline(matrix, invalid).ok());
  invalid = hipdlp::model::ScalingOptions();
  invalid.curtis_reid_sweeps = 0;
  CHECK_FALSE(hipdlp::model::computeScalingPipeline(matrix, invalid).ok());
  invalid = hipdlp::model::ScalingOptions();
  invalid.legacy_scaling_mode = kPdlpScalingMax + 1;
  CHECK_FALSE(hipdlp::model::computeScalingPipeline(matrix, invalid).ok());
  invalid = hipdlp::model::ScalingOptions();
  invalid.ruiz_iterations = -1;
  CHECK_FALSE(hipdlp::model::computeScalingPipeline(matrix, invalid).ok());
  invalid = hipdlp::model::ScalingOptions();
  invalid.features_off = kPdlpAllFeaturesOff + 1;
  CHECK_FALSE(hipdlp::model::computeScalingPipeline(matrix, invalid).ok());

  hipdlp::model::SparseMatrix malformed = matrix;
  malformed.starts = {0, 4};
  CHECK_FALSE(hipdlp::model::computeRuizScaling(malformed, 1).ok());
  CHECK_FALSE(hipdlp::model::computePockChambolleScaling(malformed).ok());
  CHECK_FALSE(hipdlp::model::computeL2Scaling(malformed).ok());
  hipdlp::model::SparseMatrix non_finite = matrix;
  non_finite.values[0] = std::numeric_limits<double>::infinity();
  CHECK_FALSE(hipdlp::model::computeRuizScaling(non_finite, 1).ok());
  CHECK_FALSE(hipdlp::model::computePockChambolleScaling(non_finite).ok());
  CHECK_FALSE(hipdlp::model::computeL2Scaling(non_finite).ok());
  hipdlp::model::SparseMatrix overflowing_sum;
  overflowing_sum.num_rows = 1;
  overflowing_sum.num_cols = 2;
  overflowing_sum.starts = {0, 2};
  overflowing_sum.indices = {0, 1};
  overflowing_sum.values = {1e308, 1e308};
  CHECK_FALSE(hipdlp::model::computePockChambolleScaling(overflowing_sum).ok());

  hipdlp::model::SparseMatrix extreme;
  extreme.num_rows = 3;
  extreme.num_cols = 3;
  extreme.starts = {0, 1, 2, 2};
  extreme.indices = {0, 1};
  extreme.values = {1e-300, 1e300};
  const hipdlp::model::SparseMatrix original_extreme = extreme;
  hipdlp::model::ScalingOptions extreme_options;
  extreme_options.legacy_scaling_mode = kPdlpScalingMax;
  extreme_options.ruiz_iterations = 2;
  const hipdlp::model::ScalingResult extreme_result =
      hipdlp::model::computeScalingPipeline(extreme, extreme_options);
  REQUIRE(extreme_result.ok());
  CHECK(extreme.samePattern(original_extreme));
  CHECK(extreme.values == original_extreme.values);
  for (double scale : extreme_result.scaling.row_scale) {
    CHECK(std::isfinite(scale));
    CHECK(scale > 0.0);
  }
  for (double scale : extreme_result.scaling.col_scale) {
    CHECK(std::isfinite(scale));
    CHECK(scale > 0.0);
  }

  HighsLp lp = makeUnsafeScalingLp();
  const HighsLp original_lp = lp;
  hipdlp::model::ProblemBuildResult built = hipdlp::model::buildProblem(lp);
  REQUIRE(built.ok());
  const hipdlp::model::SparseMatrix original_matrix = built.problem.matrix();
  const std::vector<double> original_costs = built.problem.costs();
  const std::vector<double> original_col_lower = built.problem.colLower();
  hipdlp::model::ScalingOptions unsafe;
  unsafe.mode = hipdlp::model::ScalingMode::kLegacy;
  unsafe.legacy_scaling_mode = kPdlpScalingRuiz;
  unsafe.ruiz_iterations = 1;
  CHECK_FALSE(hipdlp::model::applyScalingPipeline(built.problem, unsafe).ok());
  CHECK_FALSE(built.problem.scaling().applied);
  CHECK(built.problem.scaling().stages.empty());
  CHECK(built.problem.matrix().starts == original_matrix.starts);
  CHECK(built.problem.matrix().indices == original_matrix.indices);
  CHECK(built.problem.matrix().values == original_matrix.values);
  CHECK(built.problem.costs() == original_costs);
  CHECK(built.problem.colLower() == original_col_lower);
  CHECK(lp.col_cost_ == original_lp.col_cost_);
  CHECK(lp.a_matrix_ == original_lp.a_matrix_);

  hipdlp::model::ProblemBuildResult malformed_options_problem =
      hipdlp::model::buildProblem(makeTransformationLp());
  REQUIRE(malformed_options_problem.ok());
  const hipdlp::model::SparseMatrix malformed_options_matrix =
      malformed_options_problem.problem.matrix();
  hipdlp::model::ScalingOptions malformed_options;
  malformed_options.curtis_reid_sweeps = 0;
  CHECK_FALSE(hipdlp::model::applyScalingPipeline(
                  malformed_options_problem.problem, malformed_options)
                  .ok());
  CHECK_FALSE(malformed_options_problem.problem.scaling().applied);
  CHECK(malformed_options_problem.problem.scaling().stages.empty());
  CHECK(malformed_options_problem.problem.matrix().starts ==
        malformed_options_matrix.starts);
  CHECK(malformed_options_problem.problem.matrix().indices ==
        malformed_options_matrix.indices);
  CHECK(malformed_options_problem.problem.matrix().values ==
        malformed_options_matrix.values);
}

TEST_CASE("hi-pdlp-curtis-reid-safeguards",
          "[pdlp][hipdlp][phase3][model][scaling]") {
  SECTION("empty rows and columns have unit scales") {
    hipdlp::model::SparseMatrix matrix;
    matrix.num_rows = 3;
    matrix.num_cols = 3;
    matrix.starts = {0, 1, 1, 1};
    matrix.indices = {1};
    matrix.values = {2.0};
    const hipdlp::model::CurtisReidResult result =
        hipdlp::model::computeCurtisReidScaling(matrix, 20);
    REQUIRE(result.ok());
    CHECK(result.row_scale[1] == 1.0);
    CHECK(result.row_scale[2] == 1.0);
    CHECK(result.col_scale[0] == 1.0);
    CHECK(result.col_scale[2] == 1.0);
  }

  SECTION("extreme coefficients remain finite and clamped") {
    hipdlp::model::SparseMatrix matrix;
    matrix.num_rows = 2;
    matrix.num_cols = 2;
    matrix.starts = {0, 1, 2};
    matrix.indices = {0, 1};
    matrix.values = {1e-320, 1e300};
    const hipdlp::model::CurtisReidResult result =
        hipdlp::model::computeCurtisReidScaling(matrix, 20);
    REQUIRE(result.ok());
    for (double scale : result.row_scale) {
      CHECK(std::isfinite(scale));
      CHECK(scale >= 1e-30);
      CHECK(scale <= 1e30);
    }
    for (double scale : result.col_scale) CHECK(std::isfinite(scale));
  }

  SECTION("malformed and non-finite scaling inputs are rejected") {
    hipdlp::model::SparseMatrix malformed = makeCurtisReidFixture();
    malformed.starts = {0, 4};
    CHECK_FALSE(hipdlp::model::computeCurtisReidScaling(malformed, 20).ok());

    hipdlp::model::SparseMatrix non_finite = makeCurtisReidFixture();
    non_finite.values[0] = std::numeric_limits<double>::infinity();
    CHECK_FALSE(hipdlp::model::computeCurtisReidScaling(non_finite, 20).ok());
    CHECK_FALSE(
        hipdlp::model::computeCurtisReidScaling(makeCurtisReidFixture(), 0)
            .ok());
  }

  SECTION("non-finite scaling intermediates do not partially mutate") {
    HighsLp lp;
    lp.num_row_ = 1;
    lp.num_col_ = 1;
    lp.col_cost_ = {0.0};
    lp.col_lower_ = {-kHighsInf};
    lp.col_upper_ = {kHighsInf};
    lp.row_lower_ = {0.0};
    lp.row_upper_ = {0.0};
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    lp.a_matrix_.num_row_ = 1;
    lp.a_matrix_.num_col_ = 1;
    lp.a_matrix_.start_ = {0, 1};
    lp.a_matrix_.index_ = {0};
    lp.a_matrix_.value_ = {std::numeric_limits<double>::max()};
    hipdlp::model::ProblemBuildResult built = hipdlp::model::buildProblem(lp);
    REQUIRE(built.ok());
    const std::vector<double> original_values = built.problem.matrix().values;
    CHECK_FALSE(built.problem.applyDiagonalScaling({1e30}, {1e30}, 20).ok());
    CHECK_FALSE(built.problem.scaling().applied);
    CHECK(built.problem.matrix().values == original_values);
  }
}

TEST_CASE("hi-pdlp-scaling-pattern-and-round-trips",
          "[pdlp][hipdlp][phase3][model][scaling]") {
  hipdlp::model::ProblemBuildResult built =
      hipdlp::model::buildProblem(makeTransformationLp());
  REQUIRE(built.ok());
  const hipdlp::model::SparseMatrix original_matrix = built.problem.matrix();
  const hipdlp::model::SparseMatrix original_transpose =
      built.problem.transpose();
  REQUIRE(hipdlp::model::applyCurtisReidScaling(built.problem, 20).ok());
  CHECK(built.problem.scaling().applied);
  CHECK(built.problem.scaling().curtis_reid_sweeps == 20);
  CHECK(built.problem.matrix().samePattern(original_matrix));
  CHECK(built.problem.transpose().samePattern(original_transpose));

  const std::vector<double> primal = {2.0, -3.0, 4.0, 1.5, -2.5};
  const std::vector<double> dual = {1.0, -2.0, 3.0, 4.0, -5.0};
  const std::vector<double> reduced_costs = {6.0, -7.0, 8.0, 9.0, -10.0};
  std::vector<double> scaled;
  std::vector<double> recovered;
  REQUIRE(hipdlp::model::scalePrimal(built.problem, primal, scaled).ok());
  REQUIRE(hipdlp::model::unscalePrimal(built.problem, scaled, recovered).ok());
  checkVectorApprox(recovered, primal);
  REQUIRE(hipdlp::model::scaleDual(built.problem, dual, scaled).ok());
  REQUIRE(hipdlp::model::unscaleDual(built.problem, scaled, recovered).ok());
  checkVectorApprox(recovered, dual);
  REQUIRE(hipdlp::model::scaleReducedCosts(built.problem, reduced_costs, scaled)
              .ok());
  REQUIRE(hipdlp::model::unscaleReducedCosts(built.problem, scaled, recovered)
              .ok());
  checkVectorApprox(recovered, reduced_costs);
}

TEST_CASE("hi-pdlp-full-scaled-postsolve-round-trip",
          "[pdlp][hipdlp][phase3][phase3b][model][postsolve]") {
  HighsLp lp = makeTransformationLp();
  const HighsLp original = lp;
  hipdlp::model::ProblemBuildResult built = hipdlp::model::buildProblem(lp);
  REQUIRE(built.ok());

  hipdlp::model::ScalingOptions options;
  options.legacy_scaling_mode = kPdlpScalingMax;
  options.ruiz_iterations = 2;
  REQUIRE(hipdlp::model::applyScalingPipeline(built.problem, options).ok());
  REQUIRE(built.problem.scaling().curtis_reid_sweeps == 20);

  const std::vector<double> transformed_primal = {2.0, 0.5, 4.0, 2.5, -5.0};
  const std::vector<double> transformed_dual = {-1.0, 4.0, -5.0, 2.0, 3.0};
  const std::vector<double> transformed_reduced_costs = {-6.0, 7.0, -8.0, 0.0,
                                                         0.0};
  std::vector<double> scaled_primal;
  std::vector<double> scaled_dual;
  std::vector<double> scaled_reduced_costs;
  REQUIRE(hipdlp::model::scalePrimal(built.problem, transformed_primal,
                                     scaled_primal)
              .ok());
  REQUIRE(hipdlp::model::scaleDual(built.problem, transformed_dual, scaled_dual)
              .ok());
  REQUIRE(hipdlp::model::scaleReducedCosts(
              built.problem, transformed_reduced_costs, scaled_reduced_costs)
              .ok());

  const hipdlp::model::PostsolveResult recovered = hipdlp::model::postsolve(
      built.problem, scaled_primal, scaled_dual, scaled_reduced_costs);
  INFO(recovered.error.message);
  REQUIRE(recovered.ok());
  checkVectorApprox(recovered.solution.col_value, {2.0, 0.5, 4.0});
  checkVectorApprox(recovered.solution.row_value, {2.0, 10.0, 8.0, 2.5, -5.0});
  checkVectorApprox(recovered.solution.row_dual, {1.0, -2.0, 3.0, -4.0, 5.0});
  checkVectorApprox(recovered.solution.col_dual, {6.0, -7.0, 8.0});
  CHECK(recovered.objective_value == Approx(20.0).epsilon(1e-14));
  CHECK(recovered.solution.value_valid);
  CHECK(recovered.solution.dual_valid);
  CHECK(recovered.solution.col_value.size() == 3);
  CHECK(recovered.solution.col_dual.size() == 3);
  CHECK(recovered.solution.row_value.size() == 5);
  CHECK(recovered.solution.row_dual.size() == 5);

  const hipdlp::model::PostsolveResult repeated = hipdlp::model::postsolve(
      built.problem, scaled_primal, scaled_dual, scaled_reduced_costs);
  REQUIRE(repeated.ok());
  CHECK(repeated.solution.col_value == recovered.solution.col_value);
  CHECK(repeated.solution.row_value == recovered.solution.row_value);
  CHECK(repeated.solution.col_dual == recovered.solution.col_dual);
  CHECK(repeated.solution.row_dual == recovered.solution.row_dual);
  CHECK(repeated.objective_value == recovered.objective_value);

  CHECK(lp.num_col_ == original.num_col_);
  CHECK(lp.num_row_ == original.num_row_);
  CHECK(lp.col_cost_ == original.col_cost_);
  CHECK(lp.col_lower_ == original.col_lower_);
  CHECK(lp.col_upper_ == original.col_upper_);
  CHECK(lp.row_lower_ == original.row_lower_);
  CHECK(lp.row_upper_ == original.row_upper_);
  CHECK(lp.a_matrix_ == original.a_matrix_);
  CHECK(lp.sense_ == original.sense_);
  CHECK(lp.offset_ == original.offset_);
}

TEST_CASE("hi-pdlp-all-scaling-mode-round-trips",
          "[pdlp][hipdlp][phase3][phase3b][model][scaling][postsolve]") {
  struct Selection {
    hipdlp::model::ScalingMode mode;
    HighsInt features_off;
  };
  const std::vector<Selection> selections = {
      {hipdlp::model::ScalingMode::kAuto, kPdlpAllFeaturesOn},
      {hipdlp::model::ScalingMode::kOff, kPdlpAllFeaturesOn},
      {hipdlp::model::ScalingMode::kCurtisReid, kPdlpAllFeaturesOn},
      {hipdlp::model::ScalingMode::kLegacy, kPdlpAllFeaturesOn},
      {hipdlp::model::ScalingMode::kAuto, kPdlpScalingOff}};
  const std::vector<double> transformed_primal = {2.0, 0.5, 4.0, 2.5, -5.0};
  const std::vector<double> transformed_dual = {-1.0, 4.0, -5.0, 2.0, 3.0};
  const std::vector<double> transformed_reduced_costs = {-6.0, 7.0, -8.0, 0.0,
                                                         0.0};

  for (const Selection& selection : selections) {
    CAPTURE(static_cast<int>(selection.mode), selection.features_off);
    hipdlp::model::ProblemBuildResult built =
        hipdlp::model::buildProblem(makeTransformationLp());
    REQUIRE(built.ok());
    hipdlp::model::ScalingOptions options;
    options.mode = selection.mode;
    options.features_off = selection.features_off;
    options.legacy_scaling_mode = kPdlpScalingMax;
    options.ruiz_iterations = 2;
    REQUIRE(hipdlp::model::applyScalingPipeline(built.problem, options).ok());

    std::vector<double> scaled_primal;
    std::vector<double> scaled_dual;
    std::vector<double> scaled_reduced_costs;
    std::vector<double> recovered_vector;
    REQUIRE(hipdlp::model::scalePrimal(built.problem, transformed_primal,
                                       scaled_primal)
                .ok());
    REQUIRE(hipdlp::model::unscalePrimal(built.problem, scaled_primal,
                                         recovered_vector)
                .ok());
    checkVectorApprox(recovered_vector, transformed_primal);
    REQUIRE(
        hipdlp::model::scaleDual(built.problem, transformed_dual, scaled_dual)
            .ok());
    REQUIRE(
        hipdlp::model::unscaleDual(built.problem, scaled_dual, recovered_vector)
            .ok());
    checkVectorApprox(recovered_vector, transformed_dual);
    REQUIRE(hipdlp::model::scaleReducedCosts(
                built.problem, transformed_reduced_costs, scaled_reduced_costs)
                .ok());
    REQUIRE(hipdlp::model::unscaleReducedCosts(
                built.problem, scaled_reduced_costs, recovered_vector)
                .ok());
    checkVectorApprox(recovered_vector, transformed_reduced_costs);

    const hipdlp::model::PostsolveResult recovered = hipdlp::model::postsolve(
        built.problem, scaled_primal, scaled_dual, scaled_reduced_costs);
    INFO(recovered.error.message);
    REQUIRE(recovered.ok());
    checkVectorApprox(recovered.solution.col_value, {2.0, 0.5, 4.0});
    checkVectorApprox(recovered.solution.row_value,
                      {2.0, 10.0, 8.0, 2.5, -5.0});
    checkVectorApprox(recovered.solution.row_dual, {1.0, -2.0, 3.0, -4.0, 5.0});
    checkVectorApprox(recovered.solution.col_dual, {6.0, -7.0, 8.0});
    CHECK(recovered.objective_value == Approx(20.0).epsilon(1e-14));
  }
}

TEST_CASE("hi-pdlp-postsolve-transactional-rejection",
          "[pdlp][hipdlp][phase3][phase3b][model][postsolve]") {
  hipdlp::model::ProblemBuildResult built =
      hipdlp::model::buildProblem(makeTransformationLp());
  REQUIRE(built.ok());
  const std::vector<double> primal(built.problem.numCols(), 0.0);
  const std::vector<double> dual(built.problem.numRows(), 0.0);
  const std::vector<double> reduced_costs(built.problem.numCols(), 0.0);

  hipdlp::model::PostsolveResult rejected = hipdlp::model::postsolve(
      built.problem, std::vector<double>(), dual, reduced_costs);
  CHECK_FALSE(rejected.ok());
  CHECK_FALSE(rejected.solution.value_valid);
  CHECK_FALSE(rejected.solution.dual_valid);
  CHECK(rejected.solution.col_value.empty());
  std::vector<double> non_finite = primal;
  non_finite[0] = std::numeric_limits<double>::quiet_NaN();
  rejected =
      hipdlp::model::postsolve(built.problem, non_finite, dual, reduced_costs);
  CHECK_FALSE(rejected.ok());
  CHECK_FALSE(rejected.solution.value_valid);
  CHECK(rejected.solution.row_value.empty());

  HighsLp unsafe_lp = makeUnsafeScalingLp();
  unsafe_lp.col_cost_[0] = 0.0;
  unsafe_lp.a_matrix_.value_[0] = std::numeric_limits<double>::max();
  hipdlp::model::ProblemBuildResult unsafe =
      hipdlp::model::buildProblem(unsafe_lp);
  REQUIRE(unsafe.ok());
  rejected = hipdlp::model::postsolve(
      unsafe.problem, {std::numeric_limits<double>::max()}, {0.0}, {0.0});
  CHECK_FALSE(rejected.ok());
  CHECK_FALSE(rejected.solution.value_valid);
  CHECK(rejected.solution.col_value.empty());
}

TEST_CASE("hi-pdlp-minimization-postsolve-with-offset",
          "[pdlp][hipdlp][phase3][phase3b][model][postsolve]") {
  HighsLp lp = makeUnsafeScalingLp();
  lp.col_cost_[0] = 2.0;
  lp.col_lower_[0] = 0.0;
  lp.a_matrix_.value_[0] = 1.0;
  lp.row_lower_[0] = 1.0;
  lp.row_upper_[0] = 1.0;
  lp.offset_ = 3.0;
  const HighsLp original = lp;
  hipdlp::model::ProblemBuildResult built = hipdlp::model::buildProblem(lp);
  REQUIRE(built.ok());

  hipdlp::model::ScalingOptions options;
  options.mode = hipdlp::model::ScalingMode::kCurtisReid;
  REQUIRE(hipdlp::model::applyScalingPipeline(built.problem, options).ok());
  const hipdlp::model::PostsolveResult recovered =
      hipdlp::model::postsolve(built.problem, {1.0}, {2.0}, {0.0});
  INFO(recovered.error.message);
  REQUIRE(recovered.ok());
  CHECK(recovered.solution.col_value == std::vector<double>({1.0}));
  CHECK(recovered.solution.row_value == std::vector<double>({1.0}));
  CHECK(recovered.solution.col_dual == std::vector<double>({0.0}));
  CHECK(recovered.solution.row_dual == std::vector<double>({2.0}));
  CHECK(recovered.objective_value == 5.0);
  CHECK(lp.col_cost_ == original.col_cost_);
  CHECK(lp.row_lower_ == original.row_lower_);
  CHECK(lp.a_matrix_ == original.a_matrix_);
  CHECK(lp.offset_ == original.offset_);
}

TEST_CASE("hi-pdlp-owned-problem-rejects-invalid-models",
          "[pdlp][hipdlp][phase3][model]") {
  SECTION("malformed sparse storage") {
    HighsLp lp = makeTransformationLp();
    lp.a_matrix_.start_.pop_back();
    CHECK_FALSE(hipdlp::model::buildProblem(lp).ok());

    lp = makeTransformationLp();
    lp.a_matrix_.index_[0] = lp.num_row_;
    CHECK_FALSE(hipdlp::model::buildProblem(lp).ok());
  }

  SECTION("non-finite model data") {
    HighsLp lp = makeTransformationLp();
    lp.col_cost_[0] = std::numeric_limits<double>::quiet_NaN();
    CHECK_FALSE(hipdlp::model::buildProblem(lp).ok());

    lp = makeTransformationLp();
    lp.a_matrix_.value_[0] = std::numeric_limits<double>::infinity();
    CHECK_FALSE(hipdlp::model::buildProblem(lp).ok());

    lp = makeTransformationLp();
    lp.a_matrix_.value_[0] = std::numeric_limits<double>::quiet_NaN();
    CHECK_FALSE(hipdlp::model::buildProblem(lp).ok());

    lp = makeTransformationLp();
    lp.offset_ = std::numeric_limits<double>::infinity();
    CHECK_FALSE(hipdlp::model::buildProblem(lp).ok());
  }

  SECTION("invalid bounds and actual infinities") {
    HighsLp lp = makeTransformationLp();
    lp.col_lower_[0] = 3.0;
    lp.col_upper_[0] = 2.0;
    CHECK_FALSE(hipdlp::model::buildProblem(lp).ok());

    lp = makeTransformationLp();
    lp.row_lower_[0] = std::numeric_limits<double>::infinity();
    CHECK_FALSE(hipdlp::model::buildProblem(lp).ok());
  }
}

static_assert(!std::is_copy_constructible<hipdlp::model::Problem>::value,
              "Problem must remain move-only");
static_assert(std::is_move_constructible<hipdlp::model::Problem>::value,
              "Problem must remain movable");
