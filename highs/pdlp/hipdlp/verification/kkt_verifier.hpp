#ifndef PDLP_HIPDLP_VERIFICATION_KKT_VERIFIER_HPP_
#define PDLP_HIPDLP_VERIFICATION_KKT_VERIFIER_HPP_

#include <string>

#include "lp_data/HConst.h"
#include "lp_data/HStruct.h"

class HighsLp;

namespace hipdlp {
namespace verification {

struct KktMetrics {
  double max_primal_violation = 0.0;
  double max_dual_violation = 0.0;
  double max_stationarity_residual = 0.0;
  double max_row_activity_residual = 0.0;
  double max_complementarity = 0.0;
  double primal_objective = 0.0;
  double dual_objective = 0.0;
  double absolute_objective_gap = 0.0;
  double relative_primal_residual = 0.0;
  double relative_dual_residual = 0.0;
  double relative_objective_gap = 0.0;
};

struct KktVerificationResult {
  bool input_valid = false;
  bool metrics_finite = false;
  bool primal_feasible = false;
  bool dual_feasible = false;
  bool objective_gap_ok = false;
  bool pass = false;
  std::string error;
  KktMetrics metrics;
  HighsInt worst_primal_col = -1;
  HighsInt worst_primal_row = -1;
  HighsInt worst_dual_col = -1;
  HighsInt worst_dual_row = -1;
  HighsInt worst_complementarity_col = -1;
  HighsInt worst_complementarity_row = -1;
};

KktVerificationResult verifyKkt(const HighsLp& lp,
                                const HighsSolution& solution,
                                double tolerance = 1e-4);

}  // namespace verification
}  // namespace hipdlp

#endif  // PDLP_HIPDLP_VERIFICATION_KKT_VERIFIER_HPP_