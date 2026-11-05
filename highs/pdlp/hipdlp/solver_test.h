#pragma once

#include <cstdint>
#include <vector>

#include "solver.h"

namespace hipdlp {
constexpr int kDecisionCapacity = 4096;
struct Decision {
  int64_t iteration = 0, epoch = 0;
  Status status = Status::kRunning;
  bool restart = false;
  double omega = 1, fpe = 0, reference_fpe = 0;
  int64_t balance_epoch_limit = 0;
  double primal = 0, dual = 0, gap = 0;
};

// Test-only observation: every buffer is read after the entire graph completes. Record all
// checks through iteration 20,000, then restarts, every 100,000 steps and the terminal check,
// up to kDecisionCapacity entries. Each decision precedes its restart/weight update.
struct Snapshot {
  std::vector<double> row_multiplier, column_multiplier, scaled_value;
  std::vector<double> x, y, xp, yp, xa, ya;
  std::vector<double> epoch_dx, epoch_dy;
  int64_t epoch_iteration = 0;
  double omega = 1, fpe = 0;
  std::vector<Decision> decisions;
};

SolveResult SolveForTesting(const HostLp& lp, const SolveOptions& options, Snapshot* snapshot, bool fixed_weight = true);

SolveResult SolveCpuForTesting(const HostLp& lp, const SolveOptions& options, Snapshot* snapshot, bool fixed_weight = true);

}  // namespace hipdlp
