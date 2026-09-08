#ifndef PDLP_HIPDLP_BACKEND_CHECK_HPP_
#define PDLP_HIPDLP_BACKEND_CHECK_HPP_

#include <cmath>
#include <cstddef>

#include "pdlp/hipdlp/algorithm/convergence.hpp"
#include "pdlp/hipdlp/backend/types.hpp"

namespace hipdlp {
namespace backend {

struct CheckInput {
  const VectorHandle* primal;
  const VectorHandle* dual;
  const VectorHandle* ax;
  const VectorHandle* aty;
  const VectorHandle* delta_primal;
  const VectorHandle* delta_dual;
  const VectorHandle* aty_delta;
  const VectorHandle* primal_anchor;
  const VectorHandle* dual_anchor;
  bool fixed_point;
  bool anchor_distances;
  double certificate_tolerance;
  double primal_tolerance;
  double dual_tolerance;
  double gap_tolerance;

  CheckInput()
      : primal(nullptr),
        dual(nullptr),
        ax(nullptr),
        aty(nullptr),
        delta_primal(nullptr),
        delta_dual(nullptr),
        aty_delta(nullptr),
        primal_anchor(nullptr),
        dual_anchor(nullptr),
        fixed_point(false),
        anchor_distances(false),
        certificate_tolerance(1e-4),
        primal_tolerance(0.0),
        dual_tolerance(0.0),
        gap_tolerance(0.0) {}

  bool validTolerances() const {
    return std::isfinite(primal_tolerance) && primal_tolerance >= 0.0 &&
           std::isfinite(dual_tolerance) && dual_tolerance >= 0.0 &&
           std::isfinite(gap_tolerance) && gap_tolerance >= 0.0;
  }
};

struct CheckResult {
  bool valid;
  algorithm::RawConvergenceMetrics raw;
  double primal_delta_norm;
  double dual_delta_norm;
  double delta_dot;
  double primal_anchor_distance;
  double dual_anchor_distance;
  bool requires_host_oracle;
  bool primal_ray_possible;
  bool farkas_possible;
  std::size_t scalar_count;
  std::size_t device_to_host_transfers;
  std::size_t device_to_host_bytes;
  std::size_t synchronizations;

  CheckResult()
      : valid(false),
        primal_delta_norm(0.0),
        dual_delta_norm(0.0),
        delta_dot(0.0),
        primal_anchor_distance(0.0),
        dual_anchor_distance(0.0),
        requires_host_oracle(false),
        primal_ray_possible(true),
        farkas_possible(true),
        scalar_count(0),
        device_to_host_transfers(0),
        device_to_host_bytes(0),
        synchronizations(0) {}
};

}  // namespace backend
}  // namespace hipdlp

#endif