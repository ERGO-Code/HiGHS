#ifndef PDLP_HIPDLP_VERIFICATION_CERTIFICATE_VERIFIER_HPP_
#define PDLP_HIPDLP_VERIFICATION_CERTIFICATE_VERIFIER_HPP_

#include <string>
#include <vector>

#include "lp_data/HConst.h"
#include "lp_data/HStruct.h"

class HighsLp;

// Independent certificate verification for the HiGHS LP form
//
//   l_x <= x <= u_x,   l_r <= A x <= u_r,   optimize c^T x + offset.
//
// The two verifiers below take a candidate direction in the *original* LP
// space and decide, without any solver state, whether that direction proves
// unboundedness or infeasibility. They are deliberately separate from
// kkt_verifier, which only ever validates a complete primal-dual iterate.
//
// Primal unboundedness. A column direction d proves unboundedness when it lies
// in the recession cone of both bound systems and strictly improves the
// objective:
//   d_j = 0            if l_x[j] and u_x[j] are both finite,
//   d_j >= 0           if only l_x[j] is finite,
//   d_j <= 0           if only u_x[j] is finite,
//   (A d)_i likewise for l_r[i] and u_r[i],
//   sigma * c^T d < 0  with sigma = +1 for minimize and -1 for maximize.
// Violations are reported raw and divided by ||d||_2; the objective margin is
// -sigma * c^T d / ||d||_2 and must exceed the tolerance.
//
// Primal infeasibility (Farkas). A row multiplier direction y together with a
// column-bound multiplier direction z proves infeasibility when
//   z = -A^T y,
//   y_i = 0 for a free row, y_i >= 0 if only l_r[i] is finite,
//   y_i <= 0 if only u_r[i] is finite, y_i free if both are finite,
//   z_j = 0 for a free column, z_j >= 0 if only l_x[j] is finite,
//   z_j <= 0 if only u_x[j] is finite, z_j free if both are finite,
//   B = sum_i y_i b_i(y_i) + sum_j z_j g_j(z_j) > 0,
// where b_i(t) is l_r[i] for t > 0 and u_r[i] for t < 0, and g_j(t) is l_x[j]
// for t > 0 and u_x[j] for t < 0. The verifier recomputes -A^T y itself and
// reports ||z + A^T y||_inf as the stationarity residual, so a caller cannot
// certify infeasibility by supplying a consistent-looking pair. The
// certificate margin is B / ||(y, z)||_2 and must exceed the tolerance.
//
// Neither verifier applies the objective offset: a direction has no affine
// part.

namespace hipdlp {
namespace verification {

struct PrimalRayMetrics {
  double direction_norm = 0.0;
  double max_raw_column_violation = 0.0;
  double max_normalized_column_violation = 0.0;
  double max_raw_row_violation = 0.0;
  double max_normalized_row_violation = 0.0;
  double raw_objective_direction = 0.0;
  double normalized_objective_direction = 0.0;
  double objective_margin = 0.0;
  HighsInt worst_col = -1;
  HighsInt worst_row = -1;
};

struct PrimalRayVerificationResult {
  bool input_valid = false;
  bool metrics_finite = false;
  bool direction_nonzero = false;
  bool column_cone_ok = false;
  bool row_cone_ok = false;
  bool objective_improving = false;
  bool pass = false;
  std::string error;
  PrimalRayMetrics metrics;
};

struct FarkasMetrics {
  double direction_norm = 0.0;
  double max_raw_row_sign_violation = 0.0;
  double max_normalized_row_sign_violation = 0.0;
  double max_raw_column_sign_violation = 0.0;
  double max_normalized_column_sign_violation = 0.0;
  double max_raw_stationarity_residual = 0.0;
  double max_normalized_stationarity_residual = 0.0;
  double raw_bound_objective = 0.0;
  double certificate_margin = 0.0;
  HighsInt worst_row = -1;
  HighsInt worst_col = -1;
};

struct FarkasVerificationResult {
  bool input_valid = false;
  bool metrics_finite = false;
  bool direction_nonzero = false;
  bool row_signs_ok = false;
  bool column_signs_ok = false;
  bool stationarity_ok = false;
  bool bound_objective_finite = false;
  bool bound_objective_positive = false;
  bool pass = false;
  std::string error;
  FarkasMetrics metrics;
};

PrimalRayVerificationResult verifyPrimalRay(
    const HighsLp& lp, const std::vector<double>& column_direction,
    double tolerance = 1e-6);

FarkasVerificationResult verifyFarkasCertificate(
    const HighsLp& lp, const std::vector<double>& row_direction,
    const std::vector<double>& column_direction, double tolerance = 1e-6);

}  // namespace verification
}  // namespace hipdlp

#endif
