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

SolveResult SolveForTesting(const HostLp& lp, const SolveOptions& options, Snapshot* snapshot,
                            bool fixed_weight = true);
// The reduced CSR (arrays only) that the device builds from presolve masks (nonzero = kept).
// Rows and columns keep their original relative order; explicit zeros are dropped.
HostLp PackOnDeviceForTesting(const HostLp& original, const std::vector<unsigned char>& row_active,
                              const std::vector<unsigned char>& column_active);
// The model SolveHybrid hands to the reduced GPU solve (CPU-packed matrix) and the original index
// of every reduced row/column. `uncertain`: presolve gave up and the original is solved instead.
struct PresolvedForTesting {
  HostLp reduced;
  std::vector<int> row_map, column_map;
  bool uncertain = false;
  size_t reductions = 0;
};
PresolvedForTesting PresolveForTesting(const HostLp& lp, double budget_seconds);
// PresolveForTesting, then the hybrid recovery of a reduced solution (primal, dual and
// reduced_cost in reduced coordinates) into original primal and dual vectors.
void RecoverForTesting(const HostLp& lp, double budget_seconds, const SolveResult& reduced,
                       std::vector<double>* primal, std::vector<double>* dual);
// The CSC of lp's matrix computed on the device, returned as the CSR of A^T (`transposed` has
// columns rows, rows columns), and the CSR slot each transposed position came from.
void TransposeOnDeviceForTesting(const HostLp& lp, HostLp* transposed,
                                 std::vector<int>* source_slot);
// Faults injected into SolveHybrid by tests. `corrupt_recovery` perturbs the recovered primal
// solution after postsolve, so the final original-model check fails and the (single,
// remaining-budget) original-model fallback must run. `fail_pack` makes the device construction of
// the reduced matrix fail with that status: kOutOfMemory takes the CPU construction, any other
// status must be reported as the error it is.
struct HybridFault {
  bool corrupt_recovery = false;
  Status fail_pack = Status::kRunning;  // kRunning: no injected failure
};
SolveResult SolveHybridWithFaultForTesting(const HostLp& lp, const SolveOptions& options,
                                           const HybridFault& fault);
// ValidateSolution with the original matrix borrowed from the device (uploaded and kept by the
// caller of the real pipeline); the result must equal ValidateSolution's.
SolveResult ValidateResidentForTesting(const HostLp& lp, const std::vector<double>& primal,
                                       const std::vector<double>& dual, double tolerance);

}  // namespace hipdlp