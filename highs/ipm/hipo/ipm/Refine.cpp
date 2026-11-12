#include "Solver.h"

// Iterative refinement on the 6x6 system

namespace hipo {

void Solver::refine(NewtonDir& delta) {
  Clock clock;

  NewtonDir correction(m_, n_);
  NewtonDir temp(m_, n_);

  // compute the residuals of the linear system in it_->ires
  clock.start();
  it_->residuals6x6(delta);
  info_.times[kRefinementResTime] += clock.stop();

  clock.start();
  double omega = computeOmega(delta);
  info_.times[kRefinementOmegaTime] += clock.stop();

  double old_omega{};

  for (Int iter = 0; iter < kRefineMaxIter; ++iter) {
    if (omega < kRefineTol) break;

    correction.clear();
    solve6x6(correction, it_->ires);

    temp = delta;
    temp.add(correction);

    clock.start();
    it_->residuals6x6(temp);
    info_.times[kRefinementResTime] += clock.stop();

    old_omega = omega;

    clock.start();
    omega = computeOmega(temp);
    info_.times[kRefinementOmegaTime] += clock.stop();

    if (omega < old_omega) {
      delta = temp;
    } else {
      omega = old_omega;
      break;
    }
  }

  // save worst residual seen in this ipm iteration
  it_->data.back().omega = std::max(it_->data.back().omega, omega);
}

static void updateOmega(double tau_const, double& omega1, double& omega2,
                        double num, double den1, double den2, double den3) {
  // den1 = (|A| * |x|)_i
  // den2 = |b_i|
  // den3 = (|A| * e)_i * ||x||_inf
  // tau_i = tau_const * ((|A| * e)_i * ||x||_inf + |b_i|)
  const double tau = tau_const * (den3 + den2);
  if (den1 + den2 > tau) {
    double omega = num / (den1 + den2);
    omega1 = std::max(omega1, omega);
  } else {
    double omega = num / (den1 + den3);
    omega2 = std::max(omega2, omega);
  }
}

double Solver::computeOmega(const NewtonDir& delta) const {
  // Evaluate iterative refinement progress. Based on "Solving sparse linear
  // systems with sparse backward error", Arioli, Demmel, Duff.

  // |6x6 matrix| * |delta| is in it_->iabs, computed together with the
  // residuals in it_->residuals6x6(delta).

  double inf_norm_delta{};
  inf_norm_delta = std::max(inf_norm_delta, infNorm(delta.x));
  inf_norm_delta = std::max(inf_norm_delta, infNorm(delta.xl));
  inf_norm_delta = std::max(inf_norm_delta, infNorm(delta.xu));
  inf_norm_delta = std::max(inf_norm_delta, infNorm(delta.y));
  inf_norm_delta = std::max(inf_norm_delta, infNorm(delta.zl));
  inf_norm_delta = std::max(inf_norm_delta, infNorm(delta.zu));

  // rhs is in it_->res
  // residual of linear system is in it_->ires

  // tau_i =
  // tau_const * (one_norm_rows(big 6x6 matrix)_i * inf_norm_delta + |rhs_i|)

  const double tau_const = 1000.0 * (5 * n_ + m_) * 1e-16;
  double omega_1{}, omega_2{};

  assert(it_->Rd || m_ == 0);

  // First block
  for (Int i = 0; i < m_; ++i) {
    updateOmega(
        tau_const, omega_1, omega_2,
        // numerator
        std::abs(it_->ires.r1[i]),
        // denominators
        it_->iabs.r1[i], std::abs(it_->res.r1[i]),
        (model_.oneNormRows(i) + std::abs(it_->Rd[i])) * inf_norm_delta);
  }

  // Second block
  for (Int i = 0; i < n_; ++i) {
    if (model_.hasLb(i)) {
      updateOmega(tau_const, omega_1, omega_2,
                  // numerator
                  std::abs(it_->ires.r2[i]),
                  // denominators
                  it_->iabs.r2[i], std::abs(it_->res.r2[i]),
                  2 * inf_norm_delta);
    }
  }

  // Third block
  for (Int i = 0; i < n_; ++i) {
    if (model_.hasUb(i)) {
      updateOmega(tau_const, omega_1, omega_2,
                  // numerator
                  std::abs(it_->ires.r3[i]),
                  // denominators
                  it_->iabs.r3[i], std::abs(it_->res.r3[i]),
                  2 * inf_norm_delta);
    }
  }

  // Fourth block
  for (Int i = 0; i < n_; ++i) {
    const double reg_p = it_->Rp ? it_->Rp[i] : regul_.primal;

    double one_norm_row = model_.oneNormCols(i) + 2 + std::abs(reg_p);
    if (model_.qp()) one_norm_row += model_.oneNormQ(i);

    updateOmega(tau_const, omega_1, omega_2,
                // numerator
                std::abs(it_->ires.r4[i]),
                // denominators
                it_->iabs.r4[i], std::abs(it_->res.r4[i]),
                one_norm_row * inf_norm_delta);
  }

  // Fifth block
  for (Int i = 0; i < n_; ++i) {
    if (model_.hasLb(i)) {
      updateOmega(
          tau_const, omega_1, omega_2,
          // numerator
          std::abs(it_->ires.r5[i]),
          // denominators
          it_->iabs.r5[i], std::abs(it_->res.r5[i]),
          (std::abs(it_->zl[i]) + std::abs(it_->xl[i])) * inf_norm_delta);
    }
  }

  // Sixth block
  for (Int i = 0; i < n_; ++i) {
    if (model_.hasUb(i)) {
      updateOmega(
          tau_const, omega_1, omega_2,
          // numerator
          std::abs(it_->ires.r6[i]),
          // denominators
          it_->iabs.r6[i], std::abs(it_->res.r6[i]),
          (std::abs(it_->zu[i]) + std::abs(it_->xu[i])) * inf_norm_delta);
    }
  }

  return omega_1 + omega_2;
}

}  // namespace hipo