#include "Solver.h"

// Iterative refinement on the 6x6 system

namespace hipo {

void Solver::refine(NewtonDir& delta) {
  Clock clock;

  NewtonDir correction(m_, n_);

  double inf_norm_rhs = infNorm(it_->res.r1);
  inf_norm_rhs = std::max(inf_norm_rhs, infNorm(it_->res.r2));
  inf_norm_rhs = std::max(inf_norm_rhs, infNorm(it_->res.r3));
  inf_norm_rhs = std::max(inf_norm_rhs, infNorm(it_->res.r4));
  inf_norm_rhs = std::max(inf_norm_rhs, infNorm(it_->res.r5));
  inf_norm_rhs = std::max(inf_norm_rhs, infNorm(it_->res.r6));

  // compute the residuals of the linear system in it_->ires
  clock.start();
  double inf_norm_res = it_->residuals6x6(delta);
  info_.residual_time += clock.stop();

  double old_norm_res{};

  for (Int iter = 0; iter < kMaxIterRefine; ++iter) {
    if (inf_norm_res < kRefinTol + kRefineMult * inf_norm_rhs) break;

    correction.clear();
    solve6x6(correction, it_->ires);

    correction.add(delta);

    if (correction.isNan() || correction.isInf()) break;

    clock.start();
    old_norm_res = inf_norm_res;
    inf_norm_res = it_->residuals6x6(correction);
    info_.residual_time += clock.stop();

    if (old_norm_res < inf_norm_res * kRefineRatio) {
      if (inf_norm_res < old_norm_res) {
        delta = correction;
      } else {
        inf_norm_res = old_norm_res;
      }
      break;
    }
    delta = correction;
  }

  // save worst residual seen in this ipm iteration
  if (inf_norm_res > it_->data.back().worst_res) {
    it_->data.back().worst_res = inf_norm_res;
    it_->data.back().worst_res_rhs = kRefinTol + kRefineMult * inf_norm_rhs;
  }
}

}  // namespace hipo