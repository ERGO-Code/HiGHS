#pragma once
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

}  // namespace hipdlp