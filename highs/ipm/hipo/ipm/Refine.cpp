#include "Solver.h"

// Iterative refinement on the 6x6 system

namespace hipo {

void Solver::refine(NewtonDir& delta) {
  Clock clock;

  // measure of the error of the solve, selected by it_->refine_method
  const bool componentwise =
      it_->refine_method == RefineMethod::kComponentwise;
  const double tol =
      componentwise ? kRefineTolComponentwise : kRefineTolBlockwise;
  auto computeError = [&](const NewtonDir& d) {
    return componentwise ? computeComponentwiseError(d)
                         : computeBlockwiseError(d);
  };

  NewtonDir correction(m_, n_);
  NewtonDir temp(m_, n_);

  // compute the residuals of the linear system in it_->ires
  clock.start();
  it_->residuals6x6(delta, componentwise);
  info_.times[kRefinementResTime] += clock.stop();

  clock.start();
  double error = computeError(delta);
  info_.times[kRefinementErrorTime] += clock.stop();

  double old_error{};

  for (Int iter = 0; iter < kRefineMaxIter; ++iter) {
    if (error < tol) break;

    correction.clear();
    solve6x6(correction, it_->ires);

    temp = delta;
    temp.add(correction);

    clock.start();
    it_->residuals6x6(temp, componentwise);
    info_.times[kRefinementResTime] += clock.stop();

    old_error = error;

    clock.start();
    error = computeError(temp);
    info_.times[kRefinementErrorTime] += clock.stop();

    // accept the correction if the error decreased
    if (error < old_error) {
      delta = temp;
    } else {
      error = old_error;
      break;
    }

    // stop if the error did not decrease enough
    if (error > kRefineMinDecrease * old_error) break;
  }

  // save worst error seen in this ipm iteration
  it_->data.back().refinement_error =
      std::max(it_->data.back().refinement_error, error);
}

static double safeRatio(double num, double den) {
  // 0/0 is 0, num/0 is inf
  if (num == 0.0) return 0.0;
  if (den == 0.0) return kHighsInf;
  return num / den;
}

double Solver::computeBlockwiseError(const NewtonDir& delta) const {
  // Normwise backward error of each block k of the 6x6 system:
  //  eta_k = ||r_k|| / (||M_k|| * ||delta|| + ||rhs_k||),
  // where r_k, M_k and rhs_k are the residual, matrix and rhs of block k, and
  // all norms are infinity norms. ||M_k|| * ||delta|| is replaced by
  // sum_j ||M_kj|| * ||delta_j||, where M_kj is the sub-block that multiplies
  // component delta_j, so that components of different size are not mixed.
  // Returns max_k eta_k.
  //
  // It does not require |A| * |dx| and similar products, so it costs O(m + n)
  // once the residual is known.

  // rhs is in it_->res
  // residual of linear system is in it_->ires

  const double norm_dx = infNorm(delta.x);
  const double norm_dy = infNorm(delta.y);
  const double norm_dxl = infNorm(delta.xl);
  const double norm_dxu = infNorm(delta.xu);
  const double norm_dzl = infNorm(delta.zl);
  const double norm_dzu = infNorm(delta.zu);

  // norms of diagonal blocks
  double norm_Rd{}, norm_Rp{}, norm_xl{}, norm_xu{}, norm_zl{}, norm_zu{};
  assert(it_->Rd || m_ == 0);
  for (Int i = 0; i < m_; ++i)
    norm_Rd = std::max(norm_Rd, std::abs(it_->Rd[i]));
  for (Int i = 0; i < n_; ++i) {
    const double reg_p = it_->Rp ? it_->Rp[i] : regul_.primal;
    norm_Rp = std::max(norm_Rp, std::abs(reg_p));
    if (model_.hasLb(i)) {
      norm_xl = std::max(norm_xl, std::abs(it_->xl[i]));
      norm_zl = std::max(norm_zl, std::abs(it_->zl[i]));
    }
    if (model_.hasUb(i)) {
      norm_xu = std::max(norm_xu, std::abs(it_->xu[i]));
      norm_zu = std::max(norm_zu, std::abs(it_->zu[i]));
    }
  }

  double eta{};

  // First block: A * dx + Rd * dy
  eta = std::max(eta, safeRatio(infNorm(it_->ires.r1),
                                model_.normA() * norm_dx + norm_Rd * norm_dy +
                                    infNorm(it_->res.r1)));

  // Second block: dx - dxl
  eta = std::max(eta, safeRatio(infNorm(it_->ires.r2),
                                norm_dx + norm_dxl + infNorm(it_->res.r2)));

  // Third block: dx + dxu
  eta = std::max(eta, safeRatio(infNorm(it_->ires.r3),
                                norm_dx + norm_dxu + infNorm(it_->res.r3)));

  // Fourth block: A^T * dy + dzl - dzu - (Q + Rp) * dx
  eta =
      std::max(eta, safeRatio(infNorm(it_->ires.r4),
                              model_.normAt() * norm_dy + norm_dzl + norm_dzu +
                                  (model_.normQ() + norm_Rp) * norm_dx +
                                  infNorm(it_->res.r4)));

  // Fifth block: zl * dxl + xl * dzl
  eta = std::max(eta, safeRatio(infNorm(it_->ires.r5),
                                norm_zl * norm_dxl + norm_xl * norm_dzl +
                                    infNorm(it_->res.r5)));

  // Sixth block: zu * dxu + xu * dzu
  eta = std::max(eta, safeRatio(infNorm(it_->ires.r6),
                                norm_zu * norm_dxu + norm_xu * norm_dzu +
                                    infNorm(it_->res.r6)));

  return eta;
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

double Solver::computeComponentwiseError(const NewtonDir& delta) const {
  // Evaluate iterative refinement progress. Based on "Solving sparse linear
  // systems with sparse backward error", Arioli, Demmel, Duff.

  // |6x6 matrix| * |delta| is in it_->iabs, computed together with the
  // residuals in it_->residuals6x6(delta, true).

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