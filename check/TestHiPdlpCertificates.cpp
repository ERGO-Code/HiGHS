#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "catch.hpp"
#include "lp_data/HConst.h"
#include "lp_data/HighsLp.h"
#include "pdlp/hipdlp/algorithm/solver.hpp"
#include "pdlp/hipdlp/api/options.hpp"
#include "pdlp/hipdlp/api/status.hpp"
#include "pdlp/hipdlp/backend/backend_factory.hpp"
#include "pdlp/hipdlp/model/problem.hpp"
#include "pdlp/hipdlp/verification/certificate_verifier.hpp"

namespace {

const double kTolerance = 1e-6;

// x >= 0, x1 + x2 <= -1. Farkas: y = -1 gives z = -A^T y = (1, 1) >= 0 and
// bound objective y * u_row = 1 > 0.
HighsLp makeInfeasibleLp() {
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 1;
  lp.col_cost_ = {-1.0, -2.0};
  lp.col_lower_ = {0.0, 0.0};
  lp.col_upper_ = {kHighsInf, kHighsInf};
  lp.row_lower_ = {-kHighsInf};
  lp.row_upper_ = {-1.0};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 2;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1, 2};
  lp.a_matrix_.index_ = {0, 0};
  lp.a_matrix_.value_ = {1.0, 1.0};
  lp.sense_ = ObjSense::kMinimize;
  return lp;
}

// x >= 0, x1 + x2 >= 1, minimize -x1 - 2 x2. Ray d = (0, 1).
HighsLp makeUnboundedLp() {
  HighsLp lp = makeInfeasibleLp();
  lp.row_lower_ = {1.0};
  lp.row_upper_ = {kHighsInf};
  return lp;
}

HighsLp makeMaximizationUnboundedLp() {
  HighsLp lp = makeUnboundedLp();
  lp.col_cost_ = {1.0, 2.0};
  lp.sense_ = ObjSense::kMaximize;
  lp.offset_ = 7.0;
  return lp;
}

HighsLp makeLargeFixedFeasibleLp() {
  HighsLp lp;
  lp.num_col_ = 1;
  lp.num_row_ = 1;
  lp.col_cost_ = {0.0};
  lp.col_lower_ = {10000.0};
  lp.col_upper_ = {10000.0};
  lp.row_lower_ = {0.0};
  lp.row_upper_ = {0.0};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 1;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 0};
  lp.sense_ = ObjSense::kMinimize;
  return lp;
}

// One row of each type over four columns of each type. Row 5 is an extra
// equality row that touches only the free column, so a multiplier can move the
// free column's bound multiplier without moving the other three.
HighsLp makeMixedTypeLp() {
  HighsLp lp;
  lp.num_col_ = 4;
  lp.num_row_ = 6;
  lp.col_cost_ = {1.0, 1.0, 1.0, 1.0};
  //                    lower-only  upper-only  boxed  free
  lp.col_lower_ = {0.0, -kHighsInf, 1.0, -kHighsInf};
  lp.col_upper_ = {kHighsInf, 0.0, 3.0, kHighsInf};
  //               equality  lower  upper  boxed  free  equality
  lp.row_lower_ = {2.0, 1.0, -kHighsInf, -1.0, -kHighsInf, 0.0};
  lp.row_upper_ = {2.0, kHighsInf, 4.0, 1.0, kHighsInf, 0.0};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 4;
  lp.a_matrix_.num_row_ = 6;
  lp.a_matrix_.start_ = {0, 5, 10, 15, 21};
  lp.a_matrix_.index_ = {0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0,
                         1, 2, 3, 4, 0, 1, 2, 3, 4, 5};
  lp.a_matrix_.value_.assign(21, 1.0);
  lp.sense_ = ObjSense::kMinimize;
  return lp;
}

std::vector<double> transposeProduct(const HighsLp& lp,
                                     const std::vector<double>& row_values) {
  std::vector<double> result(static_cast<std::size_t>(lp.num_col_), 0.0);
  for (HighsInt col = 0; col < lp.num_col_; ++col) {
    for (HighsInt entry = lp.a_matrix_.start_[col];
         entry < lp.a_matrix_.start_[col + 1]; ++entry) {
      result[col] +=
          lp.a_matrix_.value_[entry] * row_values[lp.a_matrix_.index_[entry]];
    }
  }
  return result;
}

std::vector<double> negatedTransposeProduct(
    const HighsLp& lp, const std::vector<double>& row_values) {
  std::vector<double> result = transposeProduct(lp, row_values);
  for (std::size_t col = 0; col < result.size(); ++col)
    result[col] = -result[col];
  return result;
}

// Row 0 equality, row 1 upper-only (sign flipped), row 2 boxed (synthetic
// slack), row 3 lower-only. Maximization, so the solution dual mapping carries
// an objective sign that a direction must not carry.
HighsLp makePermutedMaximizationLp() {
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 4;
  lp.col_cost_ = {1.0, 2.0};
  lp.col_lower_ = {0.0, 0.0};
  lp.col_upper_ = {kHighsInf, kHighsInf};
  lp.row_lower_ = {3.0, -kHighsInf, -1.0, 1.0};
  lp.row_upper_ = {3.0, 5.0, 2.0, kHighsInf};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 2;
  lp.a_matrix_.num_row_ = 4;
  lp.a_matrix_.start_ = {0, 3, 6};
  lp.a_matrix_.index_ = {0, 2, 3, 1, 2, 3};
  lp.a_matrix_.value_ = {1.0, 1.0, 2.0, 1.0, 1.0, 3.0};
  lp.sense_ = ObjSense::kMaximize;
  lp.offset_ = 11.0;
  return lp;
}

bool cudaBackendUsable() {
  return hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda).ok();
}

hipdlp::api::SolverOptions makeSolveOptions(const hipdlp::backend::Kind kind) {
  hipdlp::api::SolverOptions options;
  options.implementation = hipdlp::api::SolverImplementation::kRewrite;
  options.backend = kind;
  options.time_limit = kHighsInf;
  options.kkt_tolerance = 1e-4;
  options.features_off = kPdlpAllFeaturesOn;
  options.iteration_limit = 20000;
  options.scaling.mode = hipdlp::model::ScalingMode::kOff;
  options.policy.step_size = hipdlp::api::StepSizePolicy::kFixed;
  options.policy.restart = hipdlp::api::RestartPolicy::kHalpern;
  options.optimality_tolerance = 1e-4;
  return options;
}

// A bounded feasible LP: no certificate may ever be produced for it.
HighsLp makeBoundedFeasibleLp() {
  HighsLp lp = makeInfeasibleLp();
  lp.col_upper_ = {4.0, 4.0};
  lp.row_lower_ = {1.0};
  lp.row_upper_ = {2.0};
  return lp;
}

}  // namespace

TEST_CASE("hi-pdlp-primal-ray-accepts-hand-derived-unbounded-direction",
          "[pdlp][hipdlp][phase8][certificate][ray]") {
  const HighsLp lp = makeUnboundedLp();
  const std::vector<double> direction = {0.0, 1.0};
  const hipdlp::verification::PrimalRayVerificationResult result =
      hipdlp::verification::verifyPrimalRay(lp, direction, kTolerance);
  REQUIRE(result.input_valid);
  REQUIRE(result.metrics_finite);
  REQUIRE(result.direction_nonzero);
  REQUIRE(result.column_cone_ok);
  REQUIRE(result.row_cone_ok);
  REQUIRE(result.objective_improving);
  REQUIRE(result.pass);
  REQUIRE(result.metrics.direction_norm == Approx(1.0));
  REQUIRE(result.metrics.raw_objective_direction == Approx(-2.0));
  REQUIRE(result.metrics.normalized_objective_direction == Approx(-2.0));
  REQUIRE(result.metrics.objective_margin == Approx(2.0));
  REQUIRE(result.metrics.max_raw_column_violation == Approx(0.0));
  REQUIRE(result.metrics.max_raw_row_violation == Approx(0.0));
}

TEST_CASE("hi-pdlp-primal-ray-normalizes-a-scaled-direction",
          "[pdlp][hipdlp][phase8][certificate][ray]") {
  const HighsLp lp = makeUnboundedLp();
  const std::vector<double> direction = {0.0, 10.0};
  const hipdlp::verification::PrimalRayVerificationResult result =
      hipdlp::verification::verifyPrimalRay(lp, direction, kTolerance);
  REQUIRE(result.pass);
  REQUIRE(result.metrics.direction_norm == Approx(10.0));
  REQUIRE(result.metrics.raw_objective_direction == Approx(-20.0));
  REQUIRE(result.metrics.normalized_objective_direction == Approx(-2.0));
  REQUIRE(result.metrics.objective_margin == Approx(2.0));
}

TEST_CASE("hi-pdlp-primal-ray-uses-the-original-objective-sense",
          "[pdlp][hipdlp][phase8][certificate][ray]") {
  const HighsLp lp = makeMaximizationUnboundedLp();
  const std::vector<double> direction = {0.0, 1.0};
  const hipdlp::verification::PrimalRayVerificationResult result =
      hipdlp::verification::verifyPrimalRay(lp, direction, kTolerance);
  REQUIRE(result.pass);
  REQUIRE(result.metrics.raw_objective_direction == Approx(2.0));
  REQUIRE(result.metrics.objective_margin == Approx(2.0));

  HighsLp minimized = lp;
  minimized.sense_ = ObjSense::kMinimize;
  const hipdlp::verification::PrimalRayVerificationResult flipped =
      hipdlp::verification::verifyPrimalRay(minimized, direction, kTolerance);
  REQUIRE_FALSE(flipped.objective_improving);
  REQUIRE_FALSE(flipped.pass);
  REQUIRE(flipped.metrics.objective_margin == Approx(-2.0));
}

TEST_CASE("hi-pdlp-primal-ray-enforces-every-column-and-row-recession-cone",
          "[pdlp][hipdlp][phase8][certificate][ray]") {
  const HighsLp lp = makeMixedTypeLp();
  SECTION("lower-bounded column rejects a negative component") {
    std::vector<double> direction = {-1.0, 0.0, 0.0, 0.0};
    const hipdlp::verification::PrimalRayVerificationResult result =
        hipdlp::verification::verifyPrimalRay(lp, direction, kTolerance);
    REQUIRE_FALSE(result.column_cone_ok);
    REQUIRE(result.metrics.worst_col == 0);
    REQUIRE(result.metrics.max_raw_column_violation == Approx(1.0));
    REQUIRE_FALSE(result.pass);
  }
  SECTION("upper-bounded column rejects a positive component") {
    std::vector<double> direction = {0.0, 1.0, 0.0, 0.0};
    const hipdlp::verification::PrimalRayVerificationResult result =
        hipdlp::verification::verifyPrimalRay(lp, direction, kTolerance);
    REQUIRE_FALSE(result.column_cone_ok);
    REQUIRE(result.metrics.worst_col == 1);
    REQUIRE(result.metrics.max_raw_column_violation == Approx(1.0));
  }
  SECTION("boxed column rejects any nonzero component") {
    std::vector<double> direction = {0.0, 0.0, -1.0, 0.0};
    const hipdlp::verification::PrimalRayVerificationResult result =
        hipdlp::verification::verifyPrimalRay(lp, direction, kTolerance);
    REQUIRE_FALSE(result.column_cone_ok);
    REQUIRE(result.metrics.worst_col == 2);
    REQUIRE(result.metrics.max_raw_column_violation == Approx(1.0));
  }
  SECTION("free column allows any component but the rows still bind") {
    std::vector<double> direction = {0.0, 0.0, 0.0, 1.0};
    const hipdlp::verification::PrimalRayVerificationResult result =
        hipdlp::verification::verifyPrimalRay(lp, direction, kTolerance);
    REQUIRE(result.column_cone_ok);
    // Equality, upper and boxed rows all see activity 1.
    REQUIRE_FALSE(result.row_cone_ok);
    REQUIRE(result.metrics.max_raw_row_violation == Approx(1.0));
    REQUIRE(result.metrics.worst_row == 0);
  }
}

TEST_CASE("hi-pdlp-primal-ray-rejects-degenerate-input",
          "[pdlp][hipdlp][phase8][certificate][ray]") {
  const HighsLp lp = makeUnboundedLp();
  SECTION("zero ray") {
    const std::vector<double> direction = {0.0, 0.0};
    const hipdlp::verification::PrimalRayVerificationResult result =
        hipdlp::verification::verifyPrimalRay(lp, direction, kTolerance);
    REQUIRE(result.input_valid);
    REQUIRE_FALSE(result.direction_nonzero);
    REQUIRE_FALSE(result.pass);
  }
  SECTION("wrong dimension") {
    const std::vector<double> direction = {1.0};
    const hipdlp::verification::PrimalRayVerificationResult result =
        hipdlp::verification::verifyPrimalRay(lp, direction, kTolerance);
    REQUIRE_FALSE(result.input_valid);
    REQUIRE_FALSE(result.pass);
    REQUIRE_FALSE(result.error.empty());
  }
  SECTION("NaN component") {
    const std::vector<double> direction = {
        std::numeric_limits<double>::quiet_NaN(), 1.0};
    const hipdlp::verification::PrimalRayVerificationResult result =
        hipdlp::verification::verifyPrimalRay(lp, direction, kTolerance);
    REQUIRE_FALSE(result.input_valid);
    REQUIRE_FALSE(result.pass);
  }
  SECTION("infinite component") {
    const std::vector<double> direction = {
        std::numeric_limits<double>::infinity(), 1.0};
    const hipdlp::verification::PrimalRayVerificationResult result =
        hipdlp::verification::verifyPrimalRay(lp, direction, kTolerance);
    REQUIRE_FALSE(result.input_valid);
    REQUIRE_FALSE(result.pass);
  }
}

TEST_CASE("hi-pdlp-primal-ray-rejects-a-near-threshold-direction",
          "[pdlp][hipdlp][phase8][certificate][ray]") {
  const HighsLp lp = makeUnboundedLp();
  // Unit direction so raw and normalized violations coincide.
  SECTION("just inside the cone tolerance") {
    const std::vector<double> direction = {-0.5e-6, 1.0};
    const hipdlp::verification::PrimalRayVerificationResult result =
        hipdlp::verification::verifyPrimalRay(lp, direction, kTolerance);
    REQUIRE(result.column_cone_ok);
    REQUIRE(result.pass);
  }
  SECTION("just outside the cone tolerance") {
    const std::vector<double> direction = {-2.0e-6, 1.0};
    const hipdlp::verification::PrimalRayVerificationResult result =
        hipdlp::verification::verifyPrimalRay(lp, direction, kTolerance);
    REQUIRE_FALSE(result.column_cone_ok);
    REQUIRE_FALSE(result.pass);
  }
  SECTION("objective improvement inside the tolerance is not strict") {
    HighsLp flat = lp;
    flat.col_cost_ = {0.0, -1.0e-9};
    const std::vector<double> direction = {0.0, 1.0};
    const hipdlp::verification::PrimalRayVerificationResult result =
        hipdlp::verification::verifyPrimalRay(flat, direction, kTolerance);
    REQUIRE(result.column_cone_ok);
    REQUIRE(result.row_cone_ok);
    REQUIRE_FALSE(result.objective_improving);
    REQUIRE_FALSE(result.pass);
  }
}

TEST_CASE("hi-pdlp-farkas-accepts-a-hand-derived-infeasibility-certificate",
          "[pdlp][hipdlp][phase8][certificate][farkas]") {
  const HighsLp lp = makeInfeasibleLp();
  const std::vector<double> row_direction = {-1.0};
  const std::vector<double> column_direction =
      negatedTransposeProduct(lp, row_direction);
  REQUIRE(column_direction[0] == Approx(1.0));
  REQUIRE(column_direction[1] == Approx(1.0));
  const hipdlp::verification::FarkasVerificationResult result =
      hipdlp::verification::verifyFarkasCertificate(
          lp, row_direction, column_direction, kTolerance);
  REQUIRE(result.input_valid);
  REQUIRE(result.metrics_finite);
  REQUIRE(result.direction_nonzero);
  REQUIRE(result.row_signs_ok);
  REQUIRE(result.column_signs_ok);
  REQUIRE(result.stationarity_ok);
  REQUIRE(result.bound_objective_positive);
  REQUIRE(result.pass);
  REQUIRE(result.metrics.raw_bound_objective == Approx(1.0));
  REQUIRE(result.metrics.direction_norm == Approx(std::sqrt(3.0)));
  REQUIRE(result.metrics.certificate_margin == Approx(1.0 / std::sqrt(3.0)));
  REQUIRE(result.metrics.max_raw_stationarity_residual == Approx(0.0));
}

TEST_CASE("hi-pdlp-farkas-rejects-the-wrong-row-multiplier-sign",
          "[pdlp][hipdlp][phase8][certificate][farkas]") {
  const HighsLp lp = makeInfeasibleLp();
  const std::vector<double> row_direction = {1.0};
  const std::vector<double> column_direction =
      negatedTransposeProduct(lp, row_direction);
  const hipdlp::verification::FarkasVerificationResult result =
      hipdlp::verification::verifyFarkasCertificate(
          lp, row_direction, column_direction, kTolerance);
  REQUIRE_FALSE(result.row_signs_ok);
  REQUIRE(result.metrics.worst_row == 0);
  REQUIRE(result.metrics.max_raw_row_sign_violation == Approx(1.0));
  REQUIRE_FALSE(result.pass);
}

TEST_CASE("hi-pdlp-farkas-enforces-every-row-and-column-sign-condition",
          "[pdlp][hipdlp][phase8][certificate][farkas]") {
  const HighsLp lp = makeMixedTypeLp();
  SECTION("free row forces a zero multiplier") {
    std::vector<double> row_direction = {0.0, 0.0, 0.0, 0.0, 1.0, 0.0};
    const std::vector<double> column_direction =
        negatedTransposeProduct(lp, row_direction);
    const hipdlp::verification::FarkasVerificationResult result =
        hipdlp::verification::verifyFarkasCertificate(
            lp, row_direction, column_direction, kTolerance);
    REQUIRE_FALSE(result.row_signs_ok);
    REQUIRE(result.metrics.worst_row == 4);
    REQUIRE(result.metrics.max_raw_row_sign_violation == Approx(1.0));
  }
  SECTION("lower-bounded row rejects a negative multiplier") {
    std::vector<double> row_direction = {0.0, -1.0, 0.0, 0.0, 0.0, 0.0};
    const std::vector<double> column_direction =
        negatedTransposeProduct(lp, row_direction);
    const hipdlp::verification::FarkasVerificationResult result =
        hipdlp::verification::verifyFarkasCertificate(
            lp, row_direction, column_direction, kTolerance);
    REQUIRE_FALSE(result.row_signs_ok);
    REQUIRE(result.metrics.worst_row == 1);
  }
  SECTION("upper-bounded row rejects a positive multiplier") {
    std::vector<double> row_direction = {0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
    const std::vector<double> column_direction =
        negatedTransposeProduct(lp, row_direction);
    const hipdlp::verification::FarkasVerificationResult result =
        hipdlp::verification::verifyFarkasCertificate(
            lp, row_direction, column_direction, kTolerance);
    REQUIRE_FALSE(result.row_signs_ok);
    REQUIRE(result.metrics.worst_row == 2);
  }
  SECTION("equality and boxed rows accept either multiplier sign") {
    std::vector<double> row_direction = {-1.0, 0.0, 0.0, 1.0, 0.0, 0.0};
    const std::vector<double> column_direction =
        negatedTransposeProduct(lp, row_direction);
    const hipdlp::verification::FarkasVerificationResult result =
        hipdlp::verification::verifyFarkasCertificate(
            lp, row_direction, column_direction, kTolerance);
    REQUIRE(result.row_signs_ok);
    // -A^T y is exactly zero here, so the free column sign condition holds.
    REQUIRE(result.column_signs_ok);
    REQUIRE(result.stationarity_ok);
    // Bound objective is (-1)(2) + (1)(-1) = -3, so the certificate fails.
    REQUIRE(result.metrics.raw_bound_objective == Approx(-3.0));
    REQUIRE_FALSE(result.bound_objective_positive);
    REQUIRE_FALSE(result.pass);
  }
  SECTION("free column forces a zero bound multiplier") {
    std::vector<double> row_direction = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
    const std::vector<double> column_direction =
        negatedTransposeProduct(lp, row_direction);
    const hipdlp::verification::FarkasVerificationResult result =
        hipdlp::verification::verifyFarkasCertificate(
            lp, row_direction, column_direction, kTolerance);
    REQUIRE(result.row_signs_ok);
    REQUIRE_FALSE(result.column_signs_ok);
    REQUIRE(result.metrics.worst_col == 3);
    REQUIRE(result.metrics.max_raw_column_sign_violation == Approx(1.0));
    REQUIRE_FALSE(result.pass);
  }
}

TEST_CASE("hi-pdlp-farkas-checks-stationarity-independently",
          "[pdlp][hipdlp][phase8][certificate][farkas]") {
  const HighsLp lp = makeInfeasibleLp();
  const std::vector<double> row_direction = {-1.0};
  std::vector<double> column_direction =
      negatedTransposeProduct(lp, row_direction);
  SECTION("a perturbation inside the tolerance still passes") {
    column_direction[1] += 0.5e-6;
    const hipdlp::verification::FarkasVerificationResult result =
        hipdlp::verification::verifyFarkasCertificate(
            lp, row_direction, column_direction, kTolerance);
    REQUIRE(result.stationarity_ok);
    REQUIRE(result.pass);
  }
  SECTION("a perturbation outside the tolerance fails") {
    column_direction[1] += 1.0e-3;
    const hipdlp::verification::FarkasVerificationResult result =
        hipdlp::verification::verifyFarkasCertificate(
            lp, row_direction, column_direction, kTolerance);
    REQUIRE_FALSE(result.stationarity_ok);
    REQUIRE(result.metrics.worst_col == 1);
    REQUIRE(result.metrics.max_raw_stationarity_residual == Approx(1.0e-3));
    REQUIRE_FALSE(result.pass);
  }
  SECTION("an unrelated column direction fails stationarity") {
    const std::vector<double> unrelated = {-1.0, -1.0};
    const hipdlp::verification::FarkasVerificationResult result =
        hipdlp::verification::verifyFarkasCertificate(lp, row_direction,
                                                      unrelated, kTolerance);
    REQUIRE_FALSE(result.stationarity_ok);
    REQUIRE_FALSE(result.pass);
  }
}

TEST_CASE("hi-pdlp-farkas-accounts-for-stationarity-error-at-bound-scale",
          "[pdlp][hipdlp][phase7][certificate][farkas][regression]") {
  const HighsLp lp = makeLargeFixedFeasibleLp();
  const std::vector<double> row_direction = {1.0};
  const std::vector<double> column_direction = {5.0e-5};
  const hipdlp::verification::FarkasVerificationResult result =
      hipdlp::verification::verifyFarkasCertificate(
          lp, row_direction, column_direction, 1.0e-4);

  REQUIRE(result.input_valid);
  REQUIRE(result.metrics_finite);
  REQUIRE(result.row_signs_ok);
  REQUIRE(result.column_signs_ok);
  REQUIRE(result.stationarity_ok);
  REQUIRE_FALSE(result.bound_objective_positive);
  REQUIRE_FALSE(result.pass);
}

TEST_CASE("hi-pdlp-farkas-rejects-degenerate-input",
          "[pdlp][hipdlp][phase8][certificate][farkas]") {
  const HighsLp lp = makeInfeasibleLp();
  SECTION("zero certificate") {
    const std::vector<double> row_direction = {0.0};
    const std::vector<double> column_direction = {0.0, 0.0};
    const hipdlp::verification::FarkasVerificationResult result =
        hipdlp::verification::verifyFarkasCertificate(
            lp, row_direction, column_direction, kTolerance);
    REQUIRE(result.input_valid);
    REQUIRE_FALSE(result.direction_nonzero);
    REQUIRE_FALSE(result.pass);
  }
  SECTION("wrong row dimension") {
    const std::vector<double> row_direction = {-1.0, -1.0};
    const std::vector<double> column_direction = {1.0, 1.0};
    const hipdlp::verification::FarkasVerificationResult result =
        hipdlp::verification::verifyFarkasCertificate(
            lp, row_direction, column_direction, kTolerance);
    REQUIRE_FALSE(result.input_valid);
    REQUIRE_FALSE(result.pass);
    REQUIRE_FALSE(result.error.empty());
  }
  SECTION("wrong column dimension") {
    const std::vector<double> row_direction = {-1.0};
    const std::vector<double> column_direction = {1.0};
    const hipdlp::verification::FarkasVerificationResult result =
        hipdlp::verification::verifyFarkasCertificate(
            lp, row_direction, column_direction, kTolerance);
    REQUIRE_FALSE(result.input_valid);
    REQUIRE_FALSE(result.pass);
  }
  SECTION("non-finite multiplier") {
    const std::vector<double> row_direction = {
        -std::numeric_limits<double>::infinity()};
    const std::vector<double> column_direction = {1.0, 1.0};
    const hipdlp::verification::FarkasVerificationResult result =
        hipdlp::verification::verifyFarkasCertificate(
            lp, row_direction, column_direction, kTolerance);
    REQUIRE_FALSE(result.input_valid);
    REQUIRE_FALSE(result.pass);
  }
  SECTION("NaN bound multiplier") {
    const std::vector<double> row_direction = {-1.0};
    const std::vector<double> column_direction = {
        std::numeric_limits<double>::quiet_NaN(), 1.0};
    const hipdlp::verification::FarkasVerificationResult result =
        hipdlp::verification::verifyFarkasCertificate(
            lp, row_direction, column_direction, kTolerance);
    REQUIRE_FALSE(result.input_valid);
    REQUIRE_FALSE(result.pass);
  }
}

TEST_CASE("hi-pdlp-farkas-rejects-a-near-threshold-certificate",
          "[pdlp][hipdlp][phase8][certificate][farkas]") {
  HighsLp lp = makeInfeasibleLp();
  // Moving the row bound to zero makes the bound objective exactly zero.
  lp.row_upper_ = {0.0};
  const std::vector<double> row_direction = {-1.0};
  const std::vector<double> column_direction =
      negatedTransposeProduct(lp, row_direction);
  const hipdlp::verification::FarkasVerificationResult result =
      hipdlp::verification::verifyFarkasCertificate(
          lp, row_direction, column_direction, kTolerance);
  REQUIRE(result.row_signs_ok);
  REQUIRE(result.column_signs_ok);
  REQUIRE(result.stationarity_ok);
  REQUIRE(result.metrics.raw_bound_objective == Approx(0.0));
  REQUIRE_FALSE(result.bound_objective_positive);
  REQUIRE_FALSE(result.pass);
}

TEST_CASE("hi-pdlp-certificates-reject-each-other",
          "[pdlp][hipdlp][phase8][certificate]") {
  const HighsLp infeasible = makeInfeasibleLp();
  const HighsLp unbounded = makeUnboundedLp();
  // The unbounded ray is not a recession direction of the infeasible LP row.
  const std::vector<double> ray = {0.0, 1.0};
  const hipdlp::verification::PrimalRayVerificationResult wrong_ray =
      hipdlp::verification::verifyPrimalRay(infeasible, ray, kTolerance);
  REQUIRE_FALSE(wrong_ray.row_cone_ok);
  REQUIRE_FALSE(wrong_ray.pass);

  // The Farkas multiplier of the infeasible LP has the wrong sign once the row
  // becomes a lower-bounded row.
  const std::vector<double> row_direction = {-1.0};
  const std::vector<double> column_direction =
      negatedTransposeProduct(unbounded, row_direction);
  const hipdlp::verification::FarkasVerificationResult wrong_farkas =
      hipdlp::verification::verifyFarkasCertificate(
          unbounded, row_direction, column_direction, kTolerance);
  REQUIRE_FALSE(wrong_farkas.row_signs_ok);
  REQUIRE_FALSE(wrong_farkas.pass);
}

TEST_CASE("hi-pdlp-certificate-verifiers-accept-a-row-wise-matrix",
          "[pdlp][hipdlp][phase8][certificate]") {
  HighsLp lp = makeInfeasibleLp();
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, 2};
  lp.a_matrix_.index_ = {0, 1};
  lp.a_matrix_.value_ = {1.0, 1.0};
  const std::vector<double> row_direction = {-1.0};
  const std::vector<double> column_direction = {1.0, 1.0};
  const hipdlp::verification::FarkasVerificationResult farkas =
      hipdlp::verification::verifyFarkasCertificate(
          lp, row_direction, column_direction, kTolerance);
  REQUIRE(farkas.pass);

  HighsLp unbounded = lp;
  unbounded.row_lower_ = {1.0};
  unbounded.row_upper_ = {kHighsInf};
  const std::vector<double> ray = {0.0, 1.0};
  const hipdlp::verification::PrimalRayVerificationResult primal =
      hipdlp::verification::verifyPrimalRay(unbounded, ray, kTolerance);
  REQUIRE(primal.pass);
  REQUIRE(primal.metrics.max_raw_row_violation == Approx(0.0));
}

TEST_CASE("hi-pdlp-direction-recovery-undoes-scaling-and-row-permutation",
          "[pdlp][hipdlp][phase8][certificate][recovery]") {
  const HighsLp lp = makePermutedMaximizationLp();
  hipdlp::model::ProblemBuildResult built = hipdlp::model::buildProblem(lp);
  REQUIRE(built.ok());
  const hipdlp::model::TransformationMetadata& transform =
      built.problem.transformation();
  // Equality-like rows come first, so the boxed row 2 moves ahead of rows 1
  // and 3 and the upper row 1 keeps a negative sign.
  REQUIRE(built.problem.numRows() == 4);
  REQUIRE(built.problem.numCols() == 3);
  REQUIRE(transform.original_to_transformed_row[0] == 0);
  REQUIRE(transform.original_to_transformed_row[1] == 2);
  REQUIRE(transform.original_to_transformed_row[2] == 1);
  REQUIRE(transform.original_to_transformed_row[3] == 3);
  REQUIRE(transform.original_row_sign[1] == Approx(-1.0));
  REQUIRE(transform.slack_col_for_original_row[2] == 2);
  REQUIRE(transform.objective_sign == Approx(-1.0));

  const std::vector<double> row_scale = {2.0, 4.0, 8.0, 16.0};
  const std::vector<double> col_scale = {0.5, 0.25, 2.0};
  REQUIRE(built.problem.applyDiagonalScaling(row_scale, col_scale, 1).ok());

  SECTION("primal direction drops slacks and rescales columns") {
    const std::vector<double> scaled = {3.0, 5.0, 7.0};
    std::vector<double> original;
    REQUIRE(built.problem.recoverPrimalDirection(scaled, original).ok());
    REQUIRE(original.size() == 2);
    REQUIRE(original[0] == Approx(1.5));
    REQUIRE(original[1] == Approx(1.25));
  }

  SECTION("dual direction permutes rows and applies only the row sign") {
    const std::vector<double> scaled = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> original;
    REQUIRE(built.problem.recoverDualDirection(scaled, original).ok());
    REQUIRE(original.size() == 4);
    REQUIRE(original[0] == Approx(2.0));
    REQUIRE(original[1] == Approx(-24.0));
    REQUIRE(original[2] == Approx(8.0));
    REQUIRE(original[3] == Approx(64.0));

    // The solution mapping additionally flips every sign for a maximization
    // model, so postsolve is not a valid direction recovery here.
    std::vector<double> transformed(4, 0.0);
    for (HighsInt row = 0; row < 4; ++row)
      transformed[row] = scaled[row] * row_scale[row];
    std::vector<double> solution_dual;
    REQUIRE(built.problem.recoverOriginalDual(transformed, solution_dual).ok());
    for (HighsInt row = 0; row < 4; ++row)
      REQUIRE(solution_dual[row] == Approx(-original[row]));
  }

  SECTION("reduced-cost direction divides by the column scale") {
    const std::vector<double> scaled = {3.0, 5.0, 7.0};
    std::vector<double> original;
    REQUIRE(built.problem.recoverReducedCostDirection(scaled, original).ok());
    REQUIRE(original.size() == 2);
    REQUIRE(original[0] == Approx(6.0));
    REQUIRE(original[1] == Approx(20.0));
  }

  SECTION("recovery rejects malformed directions") {
    std::vector<double> original;
    const std::vector<double> short_col = {1.0, 2.0};
    const std::vector<double> short_row = {1.0, 2.0, 3.0};
    const std::vector<double> nan_col = {
        1.0, std::numeric_limits<double>::quiet_NaN(), 3.0};
    CHECK_FALSE(built.problem.recoverPrimalDirection(short_col, original).ok());
    CHECK_FALSE(built.problem.recoverDualDirection(short_row, original).ok());
    CHECK_FALSE(
        built.problem.recoverReducedCostDirection(nan_col, original).ok());
    const std::vector<double> inf_row = {
        1.0, std::numeric_limits<double>::infinity(), 3.0, 4.0};
    CHECK_FALSE(built.problem.recoverDualDirection(inf_row, original).ok());
  }
}

TEST_CASE("hi-pdlp-recovered-farkas-certificate-verifies-in-original-space",
          "[pdlp][hipdlp][phase8][certificate][recovery]") {
  // Row 0 fixes x1 = 3 while boxed row 2 caps x1 + x2 at 2 with x2 >= 0.
  const HighsLp lp = makePermutedMaximizationLp();
  hipdlp::model::ProblemBuildResult built = hipdlp::model::buildProblem(lp);
  REQUIRE(built.ok());
  const std::vector<double> row_scale = {2.0, 4.0, 8.0, 16.0};
  const std::vector<double> col_scale = {0.5, 0.25, 2.0};
  REQUIRE(built.problem.applyDiagonalScaling(row_scale, col_scale, 1).ok());

  // Transformed row order is [eq r0, boxed r2, negated upper r1, lower r3].
  // The transformed multiplier (1, -1, 0, 0) recovers the original certificate
  // y = (1, 0, -1, 0) with z = -A^T y = (0, 1) and bound objective
  // 3 * 1 + 2 * (-1) = 1 > 0.
  const std::vector<double> transformed_dual = {1.0, -1.0, 0.0, 0.0};
  std::vector<double> scaled_dual(4, 0.0);
  for (HighsInt row = 0; row < 4; ++row)
    scaled_dual[row] = transformed_dual[row] / row_scale[row];

  std::vector<double> original_dual;
  REQUIRE(built.problem.recoverDualDirection(scaled_dual, original_dual).ok());
  REQUIRE(original_dual[0] == Approx(1.0));
  REQUIRE(original_dual[1] == Approx(0.0));
  REQUIRE(original_dual[2] == Approx(-1.0));
  REQUIRE(original_dual[3] == Approx(0.0));

  const std::vector<double> original_column =
      negatedTransposeProduct(lp, original_dual);
  REQUIRE(original_column[0] == Approx(0.0));
  REQUIRE(original_column[1] == Approx(1.0));
  const hipdlp::verification::FarkasVerificationResult result =
      hipdlp::verification::verifyFarkasCertificate(
          lp, original_dual, original_column, kTolerance);
  REQUIRE(result.row_signs_ok);
  REQUIRE(result.column_signs_ok);
  REQUIRE(result.stationarity_ok);
  REQUIRE(result.bound_objective_positive);
  REQUIRE(result.pass);
  REQUIRE(result.metrics.raw_bound_objective == Approx(1.0));
}

TEST_CASE("hi-pdlp-certificate-result-carries-kind-and-metrics",
          "[pdlp][hipdlp][phase8][certificate][status]") {
  hipdlp::api::Certificate certificate;
  REQUIRE(certificate.kind == hipdlp::api::CertificateKind::kNone);
  REQUIRE_FALSE(certificate.valid);
  REQUIRE(certificate.column_direction.empty());
  REQUIRE(certificate.row_direction.empty());
  REQUIRE(certificate.metrics.direction_norm == Approx(0.0));
  REQUIRE(certificate.metrics.worst_row == -1);
  REQUIRE(certificate.metrics.worst_col == -1);

  hipdlp::api::SolveResult result;
  REQUIRE(result.certificate.kind == hipdlp::api::CertificateKind::kNone);
}

TEST_CASE("hi-pdlp-certificate-status-maps-to-highs-model-status",
          "[pdlp][hipdlp][phase8][certificate][status]") {
  HighsSolution solution;
  solution.col_value = {0.0};
  solution.col_dual = {0.0};
  solution.row_value = {0.0};
  solution.row_dual = {0.0};
  solution.value_valid = true;
  solution.dual_valid = true;

  const hipdlp::api::SolverStatus statuses[] = {
      hipdlp::api::SolverStatus::kInfeasible,
      hipdlp::api::SolverStatus::kUnbounded,
      hipdlp::api::SolverStatus::kUnboundedOrInfeasible};
  const HighsModelStatus expected[] = {
      HighsModelStatus::kInfeasible, HighsModelStatus::kUnbounded,
      HighsModelStatus::kUnboundedOrInfeasible};
  for (int index = 0; index < 3; ++index) {
    hipdlp::api::SolveResult result = hipdlp::api::SolveResult::success(
        statuses[index], solution, 7);
    REQUIRE(result.ok());
    HighsBasis basis;
    HighsSolution applied;
    HighsModelStatus model_status = HighsModelStatus::kNotset;
    HighsInfo info;
    REQUIRE(hipdlp::api::applySolveResult(result, basis, applied, model_status,
                                          info) == HighsStatus::kOk);
    REQUIRE(model_status == expected[index]);
    REQUIRE(info.pdlp_iteration_count == 7);
    REQUIRE(applied.value_valid);
  }
}

TEST_CASE("hi-pdlp-solve-certifies-an-infeasible-lp",
          "[pdlp][hipdlp][phase8][certificate][solver]") {
  const HighsLp lp = makeInfeasibleLp();
  const hipdlp::model::OriginalModel model(lp);
  const hipdlp::api::SolveResult result =
      hipdlp::algorithm::solve(makeSolveOptions(hipdlp::backend::Kind::kCpu),
                               model);
  REQUIRE(result.ok());
  REQUIRE(result.status == hipdlp::api::SolverStatus::kInfeasible);
  REQUIRE(result.certificate.kind ==
          hipdlp::api::CertificateKind::kPrimalInfeasibility);
  REQUIRE(result.certificate.valid);
  REQUIRE(result.certificate.row_direction.size() == 1);
  REQUIRE(result.certificate.column_direction.size() == 2);
  REQUIRE(result.certificate.metrics.margin > 1e-4);
  REQUIRE(result.certificate.metrics.direction_norm > 0.0);
  // The row multiplier must respect the upper-bounded row sign condition.
  REQUIRE(result.certificate.row_direction[0] < 0.0);
  // A complete postsolved iterate is still returned.
  REQUIRE(result.has_solution);
  REQUIRE(result.solution.value_valid);
  REQUIRE(result.solution.col_value.size() == 2);
}

TEST_CASE("hi-pdlp-solve-certifies-an-unbounded-lp",
          "[pdlp][hipdlp][phase8][certificate][solver]") {
  const HighsLp lp = makeUnboundedLp();
  const hipdlp::model::OriginalModel model(lp);
  const hipdlp::api::SolveResult result =
      hipdlp::algorithm::solve(makeSolveOptions(hipdlp::backend::Kind::kCpu),
                               model);
  REQUIRE(result.ok());
  REQUIRE(result.status == hipdlp::api::SolverStatus::kUnbounded);
  REQUIRE(result.certificate.kind ==
          hipdlp::api::CertificateKind::kPrimalUnboundedness);
  REQUIRE(result.certificate.valid);
  REQUIRE(result.certificate.column_direction.size() == 2);
  REQUIRE(result.certificate.metrics.margin > 1e-4);
  // Both columns are nonnegative, so the ray must be nonnegative.
  REQUIRE(result.certificate.column_direction[0] >= -1e-9);
  REQUIRE(result.certificate.column_direction[1] >= -1e-9);
  REQUIRE(result.has_solution);
  REQUIRE(result.solution.value_valid);
}

TEST_CASE("hi-pdlp-solve-certifies-a-maximization-unbounded-lp",
          "[pdlp][hipdlp][phase8][certificate][solver]") {
  const HighsLp lp = makeMaximizationUnboundedLp();
  const hipdlp::model::OriginalModel model(lp);
  const hipdlp::api::SolveResult result =
      hipdlp::algorithm::solve(makeSolveOptions(hipdlp::backend::Kind::kCpu),
                               model);
  REQUIRE(result.status == hipdlp::api::SolverStatus::kUnbounded);
  REQUIRE(result.certificate.kind ==
          hipdlp::api::CertificateKind::kPrimalUnboundedness);
  const hipdlp::verification::PrimalRayVerificationResult replay =
      hipdlp::verification::verifyPrimalRay(
          lp, result.certificate.column_direction, 1e-4);
  REQUIRE(replay.pass);
  REQUIRE(replay.metrics.raw_objective_direction > 0.0);
}

TEST_CASE("hi-pdlp-solve-certificates-survive-every-scaling-mode",
          "[pdlp][hipdlp][phase8][certificate][solver]") {
  const hipdlp::model::ScalingMode modes[] = {
      hipdlp::model::ScalingMode::kOff,
      hipdlp::model::ScalingMode::kAuto,
      hipdlp::model::ScalingMode::kCurtisReid,
      hipdlp::model::ScalingMode::kLegacy};
  for (int index = 0; index < 4; ++index) {
    hipdlp::api::SolverOptions options =
        makeSolveOptions(hipdlp::backend::Kind::kCpu);
    options.scaling.mode = modes[index];
    {
      const HighsLp lp = makeInfeasibleLp();
      const hipdlp::model::OriginalModel model(lp);
      const hipdlp::api::SolveResult result =
          hipdlp::algorithm::solve(options, model);
      CHECK(result.status == hipdlp::api::SolverStatus::kInfeasible);
      CHECK(result.certificate.kind ==
            hipdlp::api::CertificateKind::kPrimalInfeasibility);
      const hipdlp::verification::FarkasVerificationResult replay =
          hipdlp::verification::verifyFarkasCertificate(
              lp, result.certificate.row_direction,
              result.certificate.column_direction, 1e-4);
      CHECK(replay.pass);
    }
    {
      const HighsLp lp = makeUnboundedLp();
      const hipdlp::model::OriginalModel model(lp);
      const hipdlp::api::SolveResult result =
          hipdlp::algorithm::solve(options, model);
      CHECK(result.status == hipdlp::api::SolverStatus::kUnbounded);
      const hipdlp::verification::PrimalRayVerificationResult replay =
          hipdlp::verification::verifyPrimalRay(
              lp, result.certificate.column_direction, 1e-4);
      CHECK(replay.pass);
    }
  }
}

TEST_CASE("hi-pdlp-solve-certificates-use-the-auto-backend",
          "[pdlp][hipdlp][phase8][certificate][solver]") {
  const HighsLp lp = makeInfeasibleLp();
  const hipdlp::model::OriginalModel model(lp);
  const hipdlp::api::SolveResult result =
      hipdlp::algorithm::solve(makeSolveOptions(hipdlp::backend::Kind::kAuto),
                               model);
  REQUIRE(result.status == hipdlp::api::SolverStatus::kInfeasible);
  REQUIRE(result.certificate.valid);
}

TEST_CASE("hi-pdlp-solve-never-certifies-a-bounded-feasible-lp",
          "[pdlp][hipdlp][phase8][certificate][solver]") {
  const HighsLp lp = makeBoundedFeasibleLp();
  const hipdlp::model::OriginalModel model(lp);
  hipdlp::api::SolverOptions options =
      makeSolveOptions(hipdlp::backend::Kind::kCpu);
  // Stop long before convergence so the loop keeps rejecting candidates.
  options.iteration_limit = 40;
  const hipdlp::api::SolveResult limited =
      hipdlp::algorithm::solve(options, model);
  REQUIRE(limited.ok());
  REQUIRE(limited.status == hipdlp::api::SolverStatus::kIterationLimit);
  REQUIRE(limited.certificate.kind == hipdlp::api::CertificateKind::kNone);
  REQUIRE_FALSE(limited.certificate.valid);
  REQUIRE(limited.solution.value_valid);

  options.iteration_limit = 20000;
  const hipdlp::api::SolveResult solved =
      hipdlp::algorithm::solve(options, model);
  REQUIRE(solved.status == hipdlp::api::SolverStatus::kOptimal);
  REQUIRE(solved.certificate.kind == hipdlp::api::CertificateKind::kNone);
}

TEST_CASE("hi-pdlp-solve-prefers-a-verified-certificate-over-a-limit",
          "[pdlp][hipdlp][phase8][certificate][solver]") {
  const HighsLp lp = makeInfeasibleLp();
  const hipdlp::model::OriginalModel model(lp);
  hipdlp::api::SolverOptions options =
      makeSolveOptions(hipdlp::backend::Kind::kCpu);
  const hipdlp::api::SolveResult reference =
      hipdlp::algorithm::solve(options, model);
  REQUIRE(reference.status == hipdlp::api::SolverStatus::kInfeasible);
  REQUIRE(reference.iteration_count > 0);

  // The iteration limit now lands exactly on the certifying boundary.
  options.iteration_limit = reference.iteration_count;
  const hipdlp::api::SolveResult tied = hipdlp::algorithm::solve(options, model);
  CHECK(tied.status == hipdlp::api::SolverStatus::kInfeasible);
  CHECK(tied.iteration_count == reference.iteration_count);

  // One boundary earlier there is no certificate yet, so the limit wins.
  options.iteration_limit = reference.iteration_count - 40;
  if (options.iteration_limit > 0) {
    const hipdlp::api::SolveResult early =
        hipdlp::algorithm::solve(options, model);
    CHECK(early.status == hipdlp::api::SolverStatus::kIterationLimit);
    CHECK(early.certificate.kind == hipdlp::api::CertificateKind::kNone);
  }
}

TEST_CASE("hi-pdlp-solve-certificates-are-deterministic",
          "[pdlp][hipdlp][phase8][certificate][solver]") {
  const HighsLp lp = makeInfeasibleLp();
  const hipdlp::model::OriginalModel model(lp);
  const hipdlp::api::SolverOptions options =
      makeSolveOptions(hipdlp::backend::Kind::kCpu);
  const hipdlp::api::SolveResult first =
      hipdlp::algorithm::solve(options, model);
  const hipdlp::api::SolveResult second =
      hipdlp::algorithm::solve(options, model);
  REQUIRE(first.status == second.status);
  REQUIRE(first.iteration_count == second.iteration_count);
  REQUIRE(first.certificate.row_direction == second.certificate.row_direction);
  REQUIRE(first.certificate.column_direction ==
          second.certificate.column_direction);
}

TEST_CASE("hi-pdlp-solve-certificates-agree-on-cpu-and-cuda",
          "[pdlp][hipdlp][phase8][certificate][solver][cuda]") {
  if (!cudaBackendUsable()) {
    WARN("CUDA backend unavailable; skipping the CPU/CUDA certificate check");
    return;
  }
  const HighsLp models[] = {makeInfeasibleLp(), makeUnboundedLp()};
  for (int index = 0; index < 2; ++index) {
    const hipdlp::model::OriginalModel model(models[index]);
    const hipdlp::api::SolveResult cpu = hipdlp::algorithm::solve(
        makeSolveOptions(hipdlp::backend::Kind::kCpu), model);
    const hipdlp::api::SolveResult cuda = hipdlp::algorithm::solve(
        makeSolveOptions(hipdlp::backend::Kind::kCuda), model);
    CHECK(cpu.status == cuda.status);
    CHECK(cpu.certificate.kind == cuda.certificate.kind);
    CHECK(cpu.certificate.valid == cuda.certificate.valid);
    CHECK(cuda.certificate.valid);
    REQUIRE(cpu.certificate.column_direction.size() ==
            cuda.certificate.column_direction.size());
    for (std::size_t col = 0; col < cpu.certificate.column_direction.size();
         ++col) {
      CHECK(cuda.certificate.column_direction[col] ==
            Approx(cpu.certificate.column_direction[col]).margin(1e-9));
    }
    REQUIRE(cpu.certificate.row_direction.size() ==
            cuda.certificate.row_direction.size());
    for (std::size_t row = 0; row < cpu.certificate.row_direction.size();
         ++row) {
      CHECK(cuda.certificate.row_direction[row] ==
            Approx(cpu.certificate.row_direction[row]).margin(1e-9));
    }
  }
}

TEST_CASE("hi-pdlp-solve-certificate-logging-reports-the-certificate",
          "[pdlp][hipdlp][phase8][certificate][solver]") {
  class Sink : public hipdlp::algorithm::SolveLogSink {
   public:
    void info(const std::string& message) override {
      messages.push_back(message);
    }
    void warning(const std::string& message) override {
      messages.push_back(message);
    }
    std::vector<std::string> messages;
  };
  Sink sink;
  hipdlp::algorithm::SolveContext context;
  context.log = &sink;
  const HighsLp lp = makeInfeasibleLp();
  const hipdlp::model::OriginalModel model(lp);
  const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
      makeSolveOptions(hipdlp::backend::Kind::kCpu), model, context);
  REQUIRE(result.status == hipdlp::api::SolverStatus::kInfeasible);
  bool reported = false;
  bool claimed_optimal_kkt = false;
  for (std::size_t index = 0; index < sink.messages.size(); ++index) {
    if (sink.messages[index].find("certificate") != std::string::npos &&
        sink.messages[index].find("margin") != std::string::npos) {
      reported = true;
    }
    if (sink.messages[index].find("finished") != std::string::npos &&
        sink.messages[index].find("relative KKT") != std::string::npos) {
      claimed_optimal_kkt = true;
    }
  }
  CHECK(reported);
  CHECK_FALSE(claimed_optimal_kkt);
}
