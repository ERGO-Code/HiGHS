#include "Iterate.h"

#include "Parameters.h"
#include "ipm/hipo/factorhighs/FactorHiGHSSettings.h"

namespace hipo {

NewtonDir::NewtonDir(Int m, Int n)
    : x(n, 0.0), y(m, 0.0), xl(n, 0.0), xu(n, 0.0), zl(n, 0.0), zu(n, 0.0) {}

Iterate::Iterate(const Model& model_input, Regularisation& r)
    : model{&model_input}, delta(model->m(), model->n()), regul{r} {
  clearIter();
  clearRes();
  best_mu = 0;
}

bool Iterate::isNan() const {
  if (isNanVector(x) || isNanVector(xl) || isNanVector(xu) || isNanVector(y) ||
      isNanVector(zl) || isNanVector(zu))
    return true;
  return false;
}
bool Iterate::isInf() const {
  if (isInfVector(x) || isInfVector(xl) || isInfVector(xu) || isInfVector(y) ||
      isInfVector(zl) || isInfVector(zu))
    return true;
  return false;
}
bool Iterate::isResNan() const {
  if (isNanVector(res1) || isNanVector(res2) || isNanVector(res3) ||
      isNanVector(res4) || isNanVector(res5) || isNanVector(res6))
    return true;
  return false;
}
bool Iterate::isResInf() const {
  if (isInfVector(res1) || isInfVector(res2) || isInfVector(res3) ||
      isInfVector(res4) || isInfVector(res5) || isInfVector(res6))
    return true;
  return false;
}
bool Iterate::isDirNan() const {
  if (isNanVector(delta.x) || isNanVector(delta.xl) || isNanVector(delta.xu) ||
      isNanVector(delta.y) || isNanVector(delta.zl) || isNanVector(delta.zu))
    return true;
  return false;
}
bool Iterate::isDirInf() const {
  if (isInfVector(delta.x) || isInfVector(delta.xl) || isInfVector(delta.xu) ||
      isInfVector(delta.y) || isInfVector(delta.zl) || isInfVector(delta.zu))
    return true;
  return false;
}

void Iterate::computeMu() {
  mu = 0.0;
  Int number_finite_bounds{};
  for (Int i = 0; i < model->n(); ++i) {
    if (model->hasLb(i)) {
      mu += xl[i] * zl[i];
      ++number_finite_bounds;
    }
    if (model->hasUb(i)) {
      mu += xu[i] * zu[i];
      ++number_finite_bounds;
    }
  }
  mu /= number_finite_bounds;

  if (best_mu > 0.0)
    best_mu = std::min(best_mu, mu);
  else
    best_mu = mu;
}
void Iterate::computeScaling() {
  scaling.assign(model->n(), 0.0);

  for (Int i = 0; i < model->n(); ++i) {
    if (model->hasLb(i)) scaling[i] += zl[i] / xl[i];
    if (model->hasUb(i)) scaling[i] += zu[i] / xu[i];

    // slow down the growth of theta
    if (scaling[i] < 1e-12) scaling[i] = sqrt(1e-12 * scaling[i]);
  }

  // compute extremes of theta
  data.back().min_theta = std::numeric_limits<double>::infinity();
  data.back().max_theta = 0.0;
  for (double d : scaling) {
    if (d != 0.0) {
      data.back().min_theta = std::min(data.back().min_theta, 1.0 / d);
      data.back().max_theta = std::max(data.back().max_theta, 1.0 / d);
    }
  }
}
void Iterate::products() {
  double& min_prod = data.back().min_prod;
  double& max_prod = data.back().max_prod;
  Int& num_small = data.back().num_small_prod;
  Int& num_large = data.back().num_large_prod;

  min_prod = std::numeric_limits<double>::max();
  max_prod = 0.0;
  num_small = 0;
  num_large = 0;

  for (Int i = 0; i < model->n(); ++i) {
    if (model->hasLb(i)) {
      double prod = xl[i] * zl[i] / mu;
      min_prod = std::min(min_prod, prod);
      max_prod = std::max(max_prod, prod);
      if (prod < kSmallProduct) ++num_small;
      if (prod > kLargeProduct) ++num_large;
    }
    if (model->hasUb(i)) {
      double prod = xu[i] * zu[i] / mu;
      min_prod = std::min(min_prod, prod);
      max_prod = std::max(max_prod, prod);
      if (prod < kSmallProduct) ++num_small;
      if (prod > kLargeProduct) ++num_large;
    }
  }
}

void Iterate::indicators() {
  primalObj();
  dualObj();
  primalInfeasUnscaled();
  dualInfeasUnscaled();
  pdGap();
  products();
}

void Iterate::primalObj() { pobj = model->offset() + dotProd(x, model->c()); }
void Iterate::dualObj() {
  dobj = model->offset() + dotProd(y, model->b());
  for (Int i = 0; i < model->n(); ++i) {
    if (model->hasLb(i)) dobj += model->lb(i) * zl[i];
    if (model->hasUb(i)) dobj -= model->ub(i) * zu[i];
  }
}
void Iterate::pdGap() {
  // relative primal-dual gap
  pdgap = std::abs(pobj - dobj) / (1.0 + 0.5 * std::abs(pobj + dobj));
}

void Iterate::primalInfeas() {
  // relative infinity norm of scaled primal residuals
  pinf = infNorm(res1);
  pinf = std::max(pinf, infNorm(res2));
  pinf = std::max(pinf, infNorm(res3));
  pinf /= (1.0 + model->normScaledRhs());
}
void Iterate::dualInfeas() {
  // relative infinity norm of scaled dual residual
  dinf = infNorm(res4) / (1.0 + model->normScaledObj());
}
void Iterate::primalInfeasUnscaled() {
  // relative infinity norm of unscaled primal residuals
  pinf = 0.0;
  for (Int i = 0; i < model->m(); ++i) {
    double val = std::abs(res1[i]);
    if (model->scaled()) val /= model->rowScale(i);
    pinf = std::max(pinf, val);
  }
  for (Int i = 0; i < model->n(); ++i) {
    double val = std::abs(res2[i]);
    if (model->scaled()) val *= model->colScale(i);
    pinf = std::max(pinf, val);

    val = std::abs(res3[i]);
    if (model->scaled()) val *= model->colScale(i);
    pinf = std::max(pinf, val);
  }
  pinf /= (1.0 + model->normUnscaledRhs());
}
void Iterate::dualInfeasUnscaled() {
  // relative infinity norm of unscaled dual residual
  dinf = 0.0;
  for (Int i = 0; i < model->n(); ++i) {
    double val = std::abs(res4[i]);
    if (model->scaled()) val /= model->colScale(i);
    dinf = std::max(dinf, val);
  }
  dinf /= (1.0 + model->normUnscaledObj());
}

void Iterate::residual1234() {
  // res1
  res1 = model->b();
  model->A().alphaProductPlusY(-1.0, x, res1);

  // res2
  for (Int i = 0; i < model->n(); ++i) {
    if (model->hasLb(i))
      res2[i] = model->lb(i) - x[i] + xl[i];
    else
      res2[i] = 0.0;
  }

  // res3
  for (Int i = 0; i < model->n(); ++i) {
    if (model->hasUb(i))
      res3[i] = model->ub(i) - x[i] - xu[i];
    else
      res3[i] = 0.0;
  }

  // res4
  res4 = model->c();
  model->A().alphaProductPlusY(-1.0, y, res4, true);
  for (Int i = 0; i < model->n(); ++i) {
    if (model->hasLb(i)) res4[i] -= zl[i];
    if (model->hasUb(i)) res4[i] += zu[i];
  }
}
void Iterate::residual56(double sigma) {
  for (Int i = 0; i < model->n(); ++i) {
    // res5
    if (model->hasLb(i))
      res5[i] = sigma * mu - xl[i] * zl[i];
    else
      res5[i] = 0.0;

    // res6
    if (model->hasUb(i))
      res6[i] = sigma * mu - xu[i] * zu[i];
    else
      res6[i] = 0.0;
  }
}

std::vector<double> Iterate::residual7() const {
  std::vector<double> res7(res4);
  for (Int i = 0; i < model->n(); ++i) {
    if (model->hasLb(i)) res7[i] -= ((res5[i] + zl[i] * res2[i]) / xl[i]);
    if (model->hasUb(i)) res7[i] += ((res6[i] - zu[i] * res3[i]) / xu[i]);
  }
  return res7;
}
std::vector<double> Iterate::residual8(const std::vector<double>& res7) const {
  std::vector<double> res8(res1);
  std::vector<double> temp(res7);

  // temp = (Theta^-1+Rp)^-1 * res7
  for (Int i = 0; i < model->n(); ++i) temp[i] /= scaling[i] + regul.primal;

  // res8 += A * temp
  model->A().alphaProductPlusY(1.0, temp, res8);

  return res8;
}

void Iterate::clearIter() {
  x.assign(model->n(), 0.0);
  xl.assign(model->n(), 0.0);
  xu.assign(model->n(), 0.0);
  y.assign(model->m(), 0.0);
  zl.assign(model->n(), 0.0);
  zu.assign(model->n(), 0.0);
}
void Iterate::clearRes() {
  res1.assign(model->m(), 0.0);
  res2.assign(model->n(), 0.0);
  res3.assign(model->n(), 0.0);
  res4.assign(model->n(), 0.0);
  res5.assign(model->n(), 0.0);
  res6.assign(model->n(), 0.0);
}
void Iterate::clearDir() {
  delta.x.assign(model->n(), 0.0);
  delta.xl.assign(model->n(), 0.0);
  delta.xu.assign(model->n(), 0.0);
  delta.y.assign(model->m(), 0.0);
  delta.zl.assign(model->n(), 0.0);
  delta.zu.assign(model->n(), 0.0);
}

void Iterate::extract(std::vector<double>& x_user, std::vector<double>& xl_user,
                      std::vector<double>& xu_user,
                      std::vector<double>& slack_user,
                      std::vector<double>& y_user, std::vector<double>& zl_user,
                      std::vector<double>& zu_user) const {
  // Extract solution with internal format

  // Copy x, xl, xu, zl, zu without slacks
  x_user = std::vector<double>(x.begin(), x.begin() + model->n_orig());
  xl_user = std::vector<double>(xl.begin(), xl.begin() + model->n_orig());
  xu_user = std::vector<double>(xu.begin(), xu.begin() + model->n_orig());
  zl_user = std::vector<double>(zl.begin(), zl.begin() + model->n_orig());
  zu_user = std::vector<double>(zu.begin(), zu.begin() + model->n_orig());

  // For the Lagrange multipliers, use slacks from zl and zu, to get correct
  // sign. NB: there is no explicit slack stored for equality constraints.
  y_user.resize(model->m());
  Int slack_pos = 0;
  for (Int i = 0; i < model->m(); ++i) {
    switch (model->constraint(i)) {
      case '=':
        y_user[i] = y[i];
        break;
      case '>':
        y_user[i] = zu[model->n_orig() + slack_pos];
        ++slack_pos;
        break;
      case '<':
        y_user[i] = -zl[model->n_orig() + slack_pos];
        ++slack_pos;
        break;
    }
  }

  // For x-slacks, use slacks from xl and xu, to get correct sign.
  // NB: there is no explicit slack stored for equality constraints.
  slack_user.resize(model->m());
  slack_pos = 0;
  for (Int i = 0; i < model->m(); ++i) {
    switch (model->constraint(i)) {
      case '=':
        slack_user[i] = 0.0;
        break;
      case '>':
        slack_user[i] = -xu[model->n_orig() + slack_pos];
        ++slack_pos;
        break;
      case '<':
        slack_user[i] = xl[model->n_orig() + slack_pos];
        ++slack_pos;
        break;
    }
  }
}

void Iterate::extract(std::vector<double>& x_user,
                      std::vector<double>& slack_user,
                      std::vector<double>& y_user,
                      std::vector<double>& z_user) const {
  // Extract solution with format for crossover

  // Construct complementary point (x_temp, y_temp, z_temp)
  std::vector<double> x_temp, y_temp, z_temp;
  dropToComplementarity(x_temp, y_temp, z_temp);

  // Both x_temp and z_temp include slacks.
  // They are removed from x and z, but they are used to compute slack and y.

  // Remove slacks from x and z
  x_user =
      std::vector<double>(x_temp.begin(), x_temp.begin() + model->n_orig());
  z_user =
      std::vector<double>(z_temp.begin(), z_temp.begin() + model->n_orig());

  // For inequality constraints, the corresponding z-slack may have been dropped
  // to zero, so build y from z-slacks.
  // NB: there is no explicit slack stored for equality constraints.
  y_user.resize(model->m());
  Int slack_pos = 0;
  for (Int i = 0; i < model->m(); ++i) {
    switch (model->constraint(i)) {
      case '=':
        y_user[i] = y_temp[i];
        break;
      case '>':
      case '<':
        y_user[i] = -z_temp[model->n_orig() + slack_pos];
        ++slack_pos;
        break;
    }
  }

  // Use slacks from x_temp and add slack for equality constraints.
  // NB: there is no explicit slack stored for equality constraints.
  slack_user.resize(model->m());
  slack_pos = 0;
  for (Int i = 0; i < model->m(); ++i) {
    switch (model->constraint(i)) {
      case '=':
        slack_user[i] = 0.0;
        break;
      case '>':
      case '<':
        slack_user[i] = x_temp[model->n_orig() + slack_pos];
        ++slack_pos;
        break;
    }
  }
}

void Iterate::dropToComplementarity(std::vector<double>& x_cmp,
                                    std::vector<double>& y_cmp,
                                    std::vector<double>& z_cmp) const {
  x_cmp.assign(model->n(), 0.0);
  z_cmp.assign(model->n(), 0.0);
  y_cmp = y;

  for (Int j = 0; j < model->n(); ++j) {
    // value of x_[j] within bounds
    double xj = std::max(x[j], model->lb(j));
    xj = std::min(xj, model->ub(j));

    // FIXED VARIABLE
    if (model->lb(j) == model->ub(j)) {
      x_cmp[j] = model->lb(j);
      z_cmp[j] = zl[j] - zu[j];
    }

    // BOTH BOUNDS FINITE
    else if (model->hasLb(j) && model->hasUb(j)) {
      if (zl[j] * xu[j] >= zu[j] * xl[j]) {
        // xlj/zlj <= xuj/zuj
        // Primal lower is smaller than primal upper, wrt respective duals
        if (zl[j] >= xl[j]) {
          // drop x to lower bound, set z positive
          x_cmp[j] = model->lb(j);
          z_cmp[j] = std::max(0.0, zl[j] - zu[j]);
        } else {
          // drop z to zero, set x within bounds
          x_cmp[j] = xj;
          z_cmp[j] = 0.0;
        }
      } else {
        // xuj/zuj < xlj/zlj
        // Primal upper is smaller than primal lower, wrt respective duals
        if (zu[j] >= xu[j]) {
          // drop x to upper bound, set z negative
          x_cmp[j] = model->ub(j);
          z_cmp[j] = std::min(0.0, zl[j] - zu[j]);
        } else {
          // drop z to zero, set x within bounds
          x_cmp[j] = xj;
          z_cmp[j] = 0.0;
        }
      }
    }

    // LOWER BOUND FINITE
    else if (model->hasLb(j)) {
      if (zl[j] >= xl[j]) {
        // drop x to lower bound, set z positive
        x_cmp[j] = model->lb(j);
        z_cmp[j] = std::max(0.0, zl[j] - zu[j]);
      } else {
        // drop z to zero, set x within bounds
        x_cmp[j] = xj;
        z_cmp[j] = 0.0;
      }
    }

    // UPPER BOUND FINITE
    else if (model->hasUb(j)) {
      if (zu[j] >= xu[j]) {
        // drop x to upper bound, set z negative
        x_cmp[j] = model->ub(j);
        z_cmp[j] = std::min(0.0, zl[j] - zu[j]);
      } else {
        // drop z to zero, set x within bounds
        x_cmp[j] = xj;
        z_cmp[j] = 0.0;
      }
    }

    // NO BOUNDS
    else {
      x_cmp[j] = xj;
      z_cmp[j] = 0.0;
    }
  }
}

double Iterate::infeasAfterDropping() const {
  // Compute estimate of residuals after dropping to complementarity (taken from
  // ipx).

  double pinf_max = 0.0;
  double dinf_max = 0.0;

  for (Int j = 0; j < model->n(); ++j) {
    double xdrop = 0.0;
    double zdrop = 0.0;

    if (model->hasLb(j) && model->hasUb(j)) {
      // BOTH BOUNDS FINITE
      if (zl[j] * xu[j] >= zu[j] * xl[j]) {
        if (zl[j] >= xl[j])
          xdrop = x[j] - model->lb(j);
        else
          zdrop = zl[j] - zu[j];
      } else {
        if (zu[j] >= xu[j])
          xdrop = x[j] - model->ub(j);
        else
          zdrop = zl[j] - zu[j];
      }
    }

    // LOWER BOUND FINITE
    else if (model->hasLb(j)) {
      if (zl[j] >= xl[j])
        xdrop = x[j] - model->lb(j);
      else {
        // zu[j] is zero
        zdrop = zl[j];
      }
    }

    // UPPER BOUND FINITE
    else if (model->hasUb(j)) {
      if (zu[j] >= xu[j])
        xdrop = x[j] - model->ub(j);
      else {
        // zl[j] is zero
        zdrop = -zu[j];
      }
    }

    // largest entry in column j of A
    double Amax = 0.0;
    for (Int el = model->A().start_[j]; el < model->A().start_[j + 1]; ++el)
      Amax = std::max(Amax, std::abs(model->A().value_[el]));

    pinf_max = std::max(pinf_max, std::abs(xdrop) * Amax);
    dinf_max = std::max(dinf_max, std::abs(zdrop));
  }

  pinf_max /= (1.0 + model->normScaledRhs());
  dinf_max /= (1.0 + model->normScaledObj());

  return std::max(pinf_max, dinf_max);
}

void Iterate::finalResiduals(Info& info) const {
  // If ipx has been used, the information is already available, otherwise,
  // compute it.

  if (!info.ipx_used) {
    std::vector<double> x_local, xl_local, xu_local, y_local, zl_local,
        zu_local, slack_local;

    extract(x_local, xl_local, xu_local, slack_local, y_local, zl_local,
            zu_local);

    const Int m = model->m();
    const Int n_orig = model->n_orig();

    // res1 = b - slack - A*x
    std::vector<double> res1(m);
    model->multWithoutSlack(-1.0, x_local, res1);
    for (Int i = 0; i < m; ++i) {
      res1[i] = res1[i] - slack_local[i] + model->b()[i];
      if (model->scaled()) res1[i] /= model->rowScale(i);
    }

    // res2 = lower - x + xl
    std::vector<double> res2(n_orig);
    for (Int i = 0; i < n_orig; ++i) {
      if (model->hasLb(i)) res2[i] = model->lb(i) - x_local[i] + xl_local[i];
      if (model->scaled()) res2[i] *= model->colScale(i);
    }

    // res3 = upper - x - xu
    std::vector<double> res3(n_orig);
    for (Int i = 0; i < n_orig; ++i) {
      if (model->hasUb(i)) res3[i] = model->ub(i) - x_local[i] - xu_local[i];
      if (model->scaled()) res3[i] *= model->colScale(i);
    }

    // res4 = c - A^T * y - zl + zu
    std::vector<double> res4(n_orig);
    model->multWithoutSlack(-1.0, y_local, res4, true);
    for (Int i = 0; i < n_orig; ++i) {
      res4[i] = res4[i] - zl_local[i] + zu_local[i] + model->c()[i];
      if (model->scaled()) res4[i] /= model->colScale(i);
    }

    info.p_res_abs = infNorm(res1);
    info.p_res_abs = std::max(info.p_res_abs, infNorm(res2));
    info.p_res_abs = std::max(info.p_res_abs, infNorm(res3));
    info.p_res_rel = info.p_res_abs / (1.0 + model->normUnscaledRhs());

    info.d_res_abs = infNorm(res4);
    info.d_res_rel = info.d_res_abs / (1.0 + model->normUnscaledObj());

    double pobj = model->offset();
    for (Int i = 0; i < n_orig; ++i) pobj += model->c()[i] * x_local[i];

    double dobj = model->offset();
    dobj += dotProd(y_local, model->b());
    for (Int i = 0; i < n_orig; ++i) {
      if (model->hasLb(i)) dobj += model->lb(i) * zl_local[i];
      if (model->hasUb(i)) dobj -= model->ub(i) * zu_local[i];
    }

    info.p_obj = pobj;
    info.d_obj = dobj;
    info.pd_gap = std::abs(pobj - dobj) / (1.0 + 0.5 * std::abs(pobj + dobj));
  } else {
    info.p_res_abs = info.ipx_info.abs_presidual;
    info.p_res_rel = info.ipx_info.rel_presidual;
    info.d_res_abs = info.ipx_info.abs_dresidual;
    info.d_res_rel = info.ipx_info.rel_dresidual;
    info.p_obj = info.ipx_info.pobjval;
    info.d_obj = info.ipx_info.dobjval;
    info.pd_gap = info.ipx_info.rel_objgap;
  }
}

}  // namespace hipo