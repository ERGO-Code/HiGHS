#ifndef PDLP_HIPDLP_ALGORITHM_STEP_SIZE_POLICY_HPP_
#define PDLP_HIPDLP_ALGORITHM_STEP_SIZE_POLICY_HPP_

#include "pdlp/hipdlp/api/status.hpp"

namespace hipdlp {
namespace algorithm {

struct StepSizeInitializationParameters {
  double safety_factor;

  StepSizeInitializationParameters();
};

// Consumes an already-computed squared operator norm; it never computes one.
struct StepSizeInitializationInput {
  double operator_norm_squared;
  double unscaled_rhs_norm;
  double unscaled_cost_norm;
  StepSizeInitializationParameters parameters;

  StepSizeInitializationInput();
};

struct InitialStepSizes {
  double primal_weight;
  double base_step;
  double primal_step;
  double dual_step;

  InitialStepSizes();
};

struct StepSizeInitializationResult {
  api::Error error;
  InitialStepSizes step_sizes;

  StepSizeInitializationResult();

  bool ok() const;
};

StepSizeInitializationResult initializeStepSizes(
    const StepSizeInitializationInput& input);

enum class PrimalWeightOutcome { kNone, kPidUpdate, kFallback };

enum class PrimalWeightUpdateKind { kNone, kControlled, kRevertedToBest };

struct PrimalWeightParameters {
  double integral_smoothing;
  double k_p;
  double k_i;
  double k_d;

  PrimalWeightParameters();
};

struct PrimalWeightState {
  double primal_weight;
  double error_sum;
  double previous_error;
  double best_primal_weight;
  double best_residual_gap;
  bool best_gap_valid;
  bool history_valid;

  PrimalWeightState();
};

// Consumes only precomputed scalars; it owns no iterate, backend, or reduction.
struct PrimalWeightInput {
  PrimalWeightState previous;
  double primal_step;
  double dual_step;
  double primal_distance;
  double dual_distance;
  double relative_primal_residual;
  double relative_dual_residual;
  PrimalWeightParameters parameters;

  PrimalWeightInput();
};

struct PrimalWeightUpdateInput {
  PrimalWeightState previous;
  double primal_anchor_distance;
  double dual_anchor_distance;
  double relative_primal_residual;
  double relative_dual_residual;
  double primal_step_size;
  double dual_step_size;
  PrimalWeightParameters parameters;

  PrimalWeightUpdateInput();
};

struct PrimalWeightUpdateResult {
  api::Error error;
  PrimalWeightState next;
  PrimalWeightUpdateKind kind;
  InitialStepSizes step_sizes;

  // Compatibility fields for the unintegrated Tier A6 entry point.
  PrimalWeightOutcome outcome;
  double base_step;
  double primal_step;
  double dual_step;

  PrimalWeightUpdateResult();

  bool ok() const;
};

PrimalWeightUpdateResult updatePrimalWeight(const PrimalWeightInput& input);
PrimalWeightUpdateResult updatePrimalWeightAtRestart(
    const PrimalWeightUpdateInput& input);
PrimalWeightState resetPrimalWeight(double primal_weight);

}  // namespace algorithm
}  // namespace hipdlp

#endif
