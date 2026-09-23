#include "pdlp/hipdlp/algorithm/halpern_pdhg_step.hpp"

#include <cmath>
#include <limits>
#include <utility>

namespace hipdlp {
namespace algorithm {

HalpernPdhgState::HalpernPdhgState()
    : problem(nullptr),
      workspace(nullptr),
      current_primal(nullptr),
      current_dual(nullptr),
      next_primal(nullptr),
      next_dual(nullptr),
      primal_anchor(nullptr),
      dual_anchor(nullptr),
      primal_candidate(nullptr),
      dual_candidate(nullptr),
      reflected_primal(nullptr),
      reflected_dual(nullptr),
      current_aty(nullptr),
      reflected_ax(nullptr),
      primal_step_size(std::numeric_limits<double>::quiet_NaN()),
      dual_step_size(std::numeric_limits<double>::quiet_NaN()),
      iteration(0) {}

api::Error enqueueHalpernPdhgStep(backend::Backend& backend,
                                  HalpernPdhgState& state,
                                  const backend::VectorHandle* scalar_parameters,
                                  const HighsInt scalar_offset) {
  if (state.problem == nullptr || state.workspace == nullptr ||
      state.current_primal == nullptr || state.current_dual == nullptr ||
      state.next_primal == nullptr || state.next_dual == nullptr ||
      state.primal_anchor == nullptr || state.dual_anchor == nullptr ||
      state.primal_candidate == nullptr || state.dual_candidate == nullptr ||
      state.reflected_primal == nullptr || state.reflected_dual == nullptr ||
      state.current_aty == nullptr ||
      state.reflected_ax == nullptr || !std::isfinite(state.primal_step_size) ||
      state.primal_step_size <= 0.0 || !std::isfinite(state.dual_step_size) ||
      state.dual_step_size <= 0.0 || state.iteration < 0 ||
    state.iteration == std::numeric_limits<HighsInt>::max() ||
    scalar_offset < 0 || scalar_offset > std::numeric_limits<HighsInt>::max() - 2) {
    return api::Error(api::ErrorCode::kInvalidArgument,
                      "Halpern-PDHG step state is invalid");
  }

  const HighsInt next_iteration = state.iteration + 1;
  const backend::ScalarParameter primal_step(state.primal_step_size,
                                             scalar_parameters, scalar_offset);
  api::Error error = scalar_parameters != nullptr
      ? backend.fusedProjectedPrimalUpdateAndReflect(
          *state.problem, *state.workspace, *state.current_primal,
          *state.current_aty, primal_step, *state.primal_candidate,
          *state.reflected_primal)
      : backend.fusedProjectedPrimalUpdateAndReflect(
      *state.problem, *state.workspace, *state.current_primal,
      *state.current_aty, state.primal_step_size, *state.primal_candidate,
      *state.reflected_primal);
  if (!error.ok()) return error;
  const backend::ScalarParameter dual_step(state.dual_step_size,
                                           scalar_parameters, scalar_offset + 1);
  error = backend.multiplyAndProjectedDualUpdateAndReflect(
      *state.problem, *state.workspace, *state.reflected_primal,
      *state.current_dual, dual_step, *state.reflected_ax, *state.dual_candidate,
      *state.reflected_dual);
  if (!error.ok()) return error;
  const double anchor_weight =
      1.0 / (static_cast<double>(next_iteration) + 1.0);
  const backend::ScalarParameter blend_weight(anchor_weight, scalar_parameters,
                                              scalar_offset + 2);
  error = scalar_parameters != nullptr
      ? backend.fusedVectorBlend(*state.workspace, *state.primal_anchor,
          *state.dual_anchor, *state.reflected_primal, *state.reflected_dual,
          blend_weight, *state.next_primal, *state.next_dual)
      : backend.fusedVectorBlend(*state.workspace, *state.primal_anchor,
                                   *state.dual_anchor, *state.reflected_primal,
                                   *state.reflected_dual, anchor_weight,
                                   *state.next_primal, *state.next_dual);
  if (!error.ok()) return error;
  error = backend.multiplyTranspose(*state.problem, *state.next_dual,
                                    *state.current_aty);
  if (!error.ok()) return error;

  std::swap(state.current_primal, state.next_primal);
  std::swap(state.current_dual, state.next_dual);
  state.iteration = next_iteration;
  return api::Error();
}

api::Error restartHalpernPdhgStateFromCandidate(backend::Backend& backend,
                                                                                                HalpernPdhgState& state) {
    if (state.problem == nullptr || state.current_primal == nullptr ||
            state.current_dual == nullptr || state.primal_anchor == nullptr ||
            state.dual_anchor == nullptr || state.primal_candidate == nullptr ||
            state.dual_candidate == nullptr ||
            state.current_aty == nullptr) {
        return api::Error(api::ErrorCode::kInvalidArgument,
                                            "Halpern-PDHG restart state is invalid");
    }

    api::Error error =
            backend.copy(*state.primal_candidate, *state.current_primal);
    if (!error.ok()) return error;
    error = backend.copy(*state.dual_candidate, *state.current_dual);
    if (!error.ok()) return error;
    error = backend.copy(*state.primal_candidate, *state.primal_anchor);
    if (!error.ok()) return error;
    error = backend.copy(*state.dual_candidate, *state.dual_anchor);
    if (!error.ok()) return error;
    error = backend.multiplyTranspose(*state.problem, *state.current_dual,
                                                                        *state.current_aty);
    if (!error.ok()) return error;
    state.iteration = 0;
    return api::Error();
}

}  // namespace algorithm
}  // namespace hipdlp