#ifndef HIPO_ITERATION_HISTORY_H
#define HIPO_ITERATION_HISTORY_H

#include <string>
#include <vector>

#include "IterationReport.h"
#include "Recovery.h"
#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

// Record of each ipm iteration, used for the log and to detect numerical
// trouble

struct IterationRecord {
  double sigma_aff = 0.0;
  double sigma = 0.0;
  Int correctors = 0;
  Int num_solves = 0;
  double min_theta = 0.0;
  double max_theta = 0.0;
  double min_prod = 0.0;
  double max_prod = 0.0;
  Int num_small_prod = 0;
  Int num_large_prod = 0;
  Int factorisation_used = 0;

  double bw_error = 0.0;
  double cw_error = 0.0;

  // pivots statistics
  Int num_dynamic_reg = 0;    // perturbed beyond the static regularisation
  Int num_wrong_sign = 0;     // wrong sign after regularisation
  Int num_reg_dominated = 0;  // smaller than their regularisation
  double max_reg = 0.0;
  double min_pivot = 0.0;  // smallest |regularised pivot|
  double max_pivot = 0.0;  // largest |regularised pivot|
  Int num_2x2 = 0;
  Int wrong_inertia = 0;  // distance of the inertia from the expected one

  double alpha_primal = 0.0;
  double alpha_dual = 0.0;
  double res_norm_before[4]{};
  double res_norm_after[4]{};

  Int num_near_free = 0;
  Int num_absorbed_rows = 0;

  // issues detected at the end of the iteration, and the detectors' messages,
  // one per line
  IterationReport report;
  std::string messages;

  RecoveryAction action = RecoveryAction::kNone;
  bool reset_to_best = false;
};

struct IterationHistory {
  std::vector<IterationRecord> records;
  void append();
  IterationRecord& back();
};

}  // namespace hipo

#endif