#include "mip/HighsSeparation.h"

#include <algorithm>
#include <cassert>
#include <queue>

#include "mip/HighsCliqueTable.h"
#include "mip/HighsDomain.h"
#include "mip/HighsImplications.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolverData.h"

enum class RowType : int8_t {
  Unusuable = -2,
  Geq = -1,
  Eq = 0,
  Leq = 1,
};

static double complementWithUpper(double coef, double ub, HighsCDouble& rhs) {
  rhs -= coef * ub;
  return -coef;
}

static double complementWithLower(double coef, double lb, HighsCDouble& rhs) {
  rhs -= coef * lb;
  return coef;
}

static bool separatePureBinaryKnapsackCover(
    const std::vector<double>& solvals, HighsCDouble coverweight,
    HighsCDouble lambda, std::vector<int>& cover, std::vector<int>& inds,
    std::vector<double>& vals, HighsCDouble& rhs, double feastol) {
  int coversize = cover.size();
  int rowlen = inds.size();
  std::vector<double> S;
  S.resize(coversize);
  std::vector<int8_t> coverflag;
  coverflag.resize(rowlen);
  std::sort(cover.begin(), cover.end(),
            [&](int a, int b) { return vals[a] > vals[b]; });

  HighsCDouble abartmp = vals[cover[0]];
  HighsCDouble sigma = lambda;
  for (int i = 1; i != coversize; ++i) {
    HighsCDouble delta = abartmp - vals[cover[i]];
    HighsCDouble kdelta = double(i) * delta;
    if (double(kdelta) < double(sigma)) {
      abartmp = vals[cover[i]];
      sigma -= kdelta;
    } else {
      abartmp -= sigma * (1.0 / i);
      sigma = 0.0;
      break;
    }
  }

  if (double(sigma) > 0) abartmp = rhs / double(coversize);

  double abar = double(abartmp);

  HighsCDouble sum = 0.0;
  int cplussize = 0;
  for (int i = 0; i != coversize; ++i) {
    sum += std::min(abar, vals[cover[i]]);
    S[i] = double(sum);

    if (vals[cover[i]] > abar + feastol) {
      ++cplussize;
      coverflag[cover[i]] = 1;
    } else
      coverflag[cover[i]] = -1;
  }
  assert(std::abs(double(sum - rhs) / double(rhs)) <= 1e-14);
  bool halfintegral = false;

  /* define the lifting function */
  auto g = [&](double z) {
    double hfrac = z / abar;
    double coef = 0.0;

    int h = std::floor(hfrac + 0.5);
    if (std::abs(hfrac - h) <= 1e-10 && h <= cplussize - 1) {
      halfintegral = true;
      coef = 0.5;
    }

    h = std::max(h - 1, 0);
    for (; h < coversize; ++h) {
      if (z <= S[h] + feastol) break;
    }

    return coef + h;
  };

  rhs = coversize - 1;

  for (int i = 0; i != rowlen; ++i) {
    if (coverflag[i] == -1) {
      vals[i] = 1;
    } else {
      vals[i] = g(vals[i]);
    }
  }

  if (halfintegral) {
    rhs *= 2;
    for (int i = 0; i != rowlen; ++i) vals[i] *= 2;
  }

  HighsCDouble viol = -rhs;
  for (int i = 0; i != rowlen; ++i) {
    viol += vals[i] * solvals[i];
  }

  if (double(viol) > 1e-5) {
    // printf("found pure 0-1 cover cut with violation %g\n", double(viol));
    return true;
  }

  return false;
}

static bool separateMixedBinaryKnapsackCover(
    const HighsMipSolver& mip, const std::vector<double>& solvals,
    HighsCDouble coverweight, HighsCDouble lambda, std::vector<int>& cover,
    std::vector<int>& inds, std::vector<double>& vals, HighsCDouble& rhs) {
  int coversize = cover.size();
  int rowlen = inds.size();
  std::vector<double> S;
  S.resize(coversize);
  std::vector<uint8_t> coverflag;
  coverflag.resize(rowlen);

  if (coversize == 0) return false;

  for (int i = 0; i != coversize; ++i) coverflag[cover[i]] = 1;

  std::sort(cover.begin(), cover.end(),
            [&](int a, int b) { return vals[a] > vals[b]; });
  HighsCDouble sum = 0;

  int p = coversize;
  for (int i = 0; i != coversize; ++i) {
    if (vals[cover[i]] - lambda < 0) {
      p = i;
      break;
    }
    sum += vals[cover[i]];
    S[i] = double(sum);
  }
  if (p == 0) return false;
  /* define the lifting function */
  auto phi = [&](double a) {
    for (int i = 0; i < p; ++i) {
      if (a <= S[i] - lambda) return double(i * lambda);

      if (a <= S[i]) return double((i + 1) * lambda + (HighsCDouble(a) - S[i]));
    }

    return double(p * lambda + (HighsCDouble(a) - S[p - 1]));
  };

  rhs = -lambda;

  for (int i = 0; i != rowlen; ++i) {
    if (inds[i] >= mip.numCol() ||
        mip.variableType(inds[i]) == HighsVarType::CONTINUOUS)
      continue;

    if (coverflag[i]) {
      vals[i] = std::min(vals[i], double(lambda));
      rhs += vals[i];
    } else {
      vals[i] = phi(vals[i]);
    }
  }

  HighsCDouble viol = -rhs;
  for (int i = 0; i != rowlen; ++i) {
    viol += vals[i] * solvals[i];
  }

  if (double(viol) > 1e-5) {
    // printf("found mixed 0-1 cover cut with violation %g\n", double(viol));
    return true;
  }

  return false;
}

static bool separateMixedIntegerKnapsackCover(
    const HighsMipSolver& mip, const std::vector<double>& solvals,
    const std::vector<double>& upper, HighsCDouble coverweight,
    HighsCDouble lambda, std::vector<int>& cover, std::vector<int>& inds,
    std::vector<double>& vals, HighsCDouble& rhs) {
  int coversize = cover.size();
  int rowlen = inds.size();

  int l = -1;

  HighsCDouble mu;
  for (int i = coversize - 1; i >= 0; --i) {
    int j = cover[i];
    // we divide by the coefficient for computing
    // r and eta before rounding. Function values of the lifting function
    // depend on those values and using too small values can lead to numerical
    // blow up.
    if (vals[j] < 10 * mip.mipdata_->feastol) continue;
    mu = upper[j] * vals[j] - lambda;

    if (mu > 10 * mip.mipdata_->feastol) {
      l = j;
      break;
    }
  }

  if (l == -1) return false;

  assert(mu > 10 * mip.mipdata_->feastol);

  double al = vals[l];
  double mudival = double(mu / al);
  double eta = ceil(mudival);
  HighsCDouble r = mu - floor(mudival) * al;
  // we multiply with r and it is important that it does not flip the sign
  // so we safe guard against tiny numerical errors here
  if (r < 0) r = 0;

  int kmin = floor(eta - upper[l] - 0.5);

  auto phi_l = [&](double a) {
    assert(a < 0);

    for (int k = -1; k >= kmin; --k) {
      if (a >= k * al + r) {
        assert(a < (k + 1) * al);
        return double(a - (k + 1) * r);
      }

      if (a >= k * al) {
        assert(a < k * al + r);
        return double(k * (al - r));
      }
    }

    assert(a <= -lambda + mip.mipdata_->epsilon);
    return double(kmin * (al - r));
  };

  std::sort(cover.begin(), cover.end(),
            [&](int a, int b) { return vals[a] > vals[b]; });

  HighsCDouble ulminusetaplusone = HighsCDouble(upper[l]) - eta + 1.0;
  HighsCDouble cplusthreshold = ulminusetaplusone * al;

  std::vector<HighsCDouble> a;
  std::vector<HighsCDouble> u;
  std::vector<HighsCDouble> m;

  a.resize(coversize);
  u.resize(coversize + 1);
  m.resize(coversize + 1);

  HighsCDouble usum = 0.0;
  HighsCDouble msum = 0.0;

  int cplussize;
  for (cplussize = 0; cplussize != coversize; ++cplussize) {
    int i = cover[cplussize];
    if (vals[i] < cplusthreshold) break;

    u[cplussize] = usum;
    m[cplussize] = msum;
    a[cplussize] = vals[i];
    usum += upper[i];
    msum += upper[i] * a[cplussize];
  }

  u[cplussize] = usum;
  m[cplussize] = msum;
  int kmax = floor(upper[l] - eta + 0.5);

  auto gamma_l = [&](double z) {
    assert(z > 0);
    for (int i = 0; i < cplussize; ++i) {
      int upperi = upper[cover[i]];

      for (int h = 0; h <= upperi; ++h) {
        HighsCDouble mih = m[i] + h * a[i];
        HighsCDouble uih = u[i] + h;
        HighsCDouble mihplusdeltai = mih + a[i] - cplusthreshold;
        if (z <= mihplusdeltai) {
          assert(mih <= z);
          return double(uih * ulminusetaplusone * (al - r));
        }

        for (int k = 0; k <= kmax; ++k) {
          if (z <= mihplusdeltai + k * al + r) {
            assert(mihplusdeltai + k * al <= z);
            return double((uih * ulminusetaplusone + k) * (al - r));
          }

          if (z <= mihplusdeltai + (k + 1) * al) {
            assert(mihplusdeltai + k * al + r <= z);
            return double((uih * ulminusetaplusone) * (al - r) + z - mih -
                          a[i] + cplusthreshold - (k + 1) * r);
          }
        }
      }
    }

    for (int p = 0;; ++p) {
      if (z <= m[cplussize] + p * al + r) {
        assert(m[cplussize] + p * al <= z);
        return double((u[cplussize] * ulminusetaplusone + p) * (al - r));
      }

      if (z <= m[cplussize] + (p + 1) * al) {
        assert(m[cplussize] + p * al + r <= z);
        return double((u[cplussize] * ulminusetaplusone) * (al - r) + z -
                      m[cplussize] - (p + 1) * r);
      }
    }
  };

  std::vector<uint8_t> coverflag;
  coverflag.resize(rowlen);
  for (int i : cover) coverflag[i] = 1;

  // for computing the right hand side value we recompute eta and r but use
  // tolerances for floor/ceil such that the right hand side is possibly weaker
  // and numerically safer. The largest(weakest) value for the right hand side
  // is attained when eta is as small as possible and r as large as possible.
  double safe_eta = ceil(mudival - mip.mipdata_->epsilon);
  HighsCDouble safe_r = mu - floor(mudival - mip.mipdata_->epsilon) * al;
  if (safe_r < 0) safe_r = 0;

  rhs = (HighsCDouble(upper[l]) - safe_eta) * safe_r - lambda;

  for (int i = 0; i != rowlen; ++i) {
    int col = inds[i];
    if (col >= mip.numCol() ||
        mip.variableType(col) == HighsVarType::CONTINUOUS)
      continue;

    if (coverflag[i]) {
      vals[i] = -phi_l(-vals[i]);
      rhs += vals[i] * upper[i];
    } else {
      vals[i] = gamma_l(vals[i]);
    }
  }

  HighsCDouble viol = -rhs;
  for (int i = 0; i != rowlen; ++i) {
    viol += vals[i] * solvals[i];
  }

  if (double(viol) > 1e-5) {
    // printf("found mixed integer cover cut with violation %g\n",
    // double(viol));
    // printf("al: %g  r: %g  eta: %g  mu: %g  lambda: %g\n", al, double(r),
    // eta,
    //        double(mu), double(lambda));
    return true;
  }

  return false;
}

static bool transformBaseEquation(
    const HighsMipSolver& mip, const HighsDomain& domain,
    const HighsSolution& lpsol, const int* Rindex, const double* Rvalue,
    int Rlen, double scale, HighsCDouble& rhs,
    std::vector<int8_t>& complementation, std::vector<int>& inds,
    std::vector<double>& vals, std::vector<double>& upper,
    std::vector<double>& solvals, int& nbin, int& nint, int& ncont,
    int& nunbndint) {
  complementation.clear();
  inds.clear();
  vals.clear();
  upper.clear();
  solvals.clear();

  nbin = 0;
  nint = 0;
  ncont = 0;
  nunbndint = 0;

  rhs *= scale;

  double minabsval = 1.0;
  for (int j = 0; j != Rlen; ++j) {
    double absval = std::abs(Rvalue[j]);
    if (absval > minabsval) minabsval = absval;
  }

  minabsval *= mip.mipdata_->feastol;

  for (int j = 0; j != Rlen; ++j) {
    int col = Rindex[j];
    if (col >= mip.numCol()) {
      int row = col - mip.numCol();
      bool rowintegral;
      if (row < mip.numRow())
        rowintegral = mip.mipdata_->rowintegral[row];
      else
        rowintegral = mip.mipdata_->cutpool.cutIsIntegral(
            mip.mipdata_->lp.getCutIndex(row));
      double val = Rvalue[j] * scale;
      if (false && rowintegral) {
        // todo
      } else {
        if (val > 0.0) continue;

        inds.push_back(col);
        upper.push_back(HIGHS_CONST_INF);
        vals.push_back(val);
        complementation.push_back(0);
        solvals.push_back(0);

        ++ncont;
      }

      continue;
    } else if (domain.colLower_[col] == -HIGHS_CONST_INF &&
               domain.colUpper_[col] == HIGHS_CONST_INF) {
      return false;
    }

    double val = Rvalue[j] * scale;
    if (domain.colLower_[col] == domain.colUpper_[col]) {
      rhs -= val * domain.colLower_[col];
      continue;
    }

    if (std::abs(val) < minabsval) {
      if (val > 0) {
        if (domain.colLower_[col] == -HIGHS_CONST_INF) return false;

        rhs -= val * domain.colLower_[col];
      } else {
        if (domain.colUpper_[col] == HIGHS_CONST_INF) return false;

        rhs -= val * domain.colUpper_[col];
      }
      continue;
    }

    int8_t c = 0;
    double ub = HIGHS_CONST_INF;
    double solval = lpsol.col_value[col];

    if (mip.variableType(col) != HighsVarType::CONTINUOUS) {
      if (domain.colLower_[col] != -HIGHS_CONST_INF &&
          domain.colUpper_[col] != HIGHS_CONST_INF)
        ub = floor(domain.colUpper_[col] - domain.colLower_[col] + 0.5);
      // complement integer variables to have a positive coefficient
      if (domain.colUpper_[col] != HIGHS_CONST_INF && val < 0) {
        rhs -= val * domain.colUpper_[col];
        c = -1;
        val = -val;
        solval = domain.colUpper_[col] - solval;
      } else if (domain.colLower_[col] != 0.0 &&
                 domain.colLower_[col] != -HIGHS_CONST_INF) {
        c = 1;
        rhs -= val * domain.colLower_[col];
        solval = solval - domain.colLower_[col];
      }

      if (ub == 1.0)
        ++nbin;
      else if (ub == HIGHS_CONST_INF)
        ++nunbndint;
      else
        ++nint;
    } else {
      if (domain.colLower_[col] == -HIGHS_CONST_INF) {
        rhs -= val * domain.colUpper_[col];
        c = -1;
        val = -val;
        solval = domain.colUpper_[col] - solval;
      } else if (domain.colUpper_[col] == HIGHS_CONST_INF) {
        if (domain.colLower_[col] != 0.0) {
          c = 1;
          rhs -= val * domain.colLower_[col];
          solval = solval - domain.colLower_[col];
        }
      } else {
        ub = domain.colUpper_[col] - domain.colLower_[col];
        double bounddist = (lpsol.col_value[col] - domain.colLower_[col]) / ub;
        if (bounddist < 0.5) {
          if (domain.colLower_[col] != 0.0) {
            c = 1;
            rhs -= val * domain.colLower_[col];
            solval = solval - domain.colLower_[col];
          }
        } else {
          rhs -= val * domain.colUpper_[col];
          c = -1;
          val = -val;
          solval = domain.colUpper_[col] - solval;
        }
      }

      if (val > 0.0) continue;

      ++ncont;
    }

    inds.push_back(col);
    upper.push_back(ub);
    vals.push_back(val);
    complementation.push_back(c);
    solvals.push_back(solval);
  }

  // on the inequality relaxed from the base equality we perform coefficient
  // tightening before the cut generation
  int len = inds.size();
  HighsCDouble maxact = 0.0;
  bool unbnd = false;
  for (int i = 0; i != len; ++i) {
    if (upper[i] == HIGHS_CONST_INF) {
      unbnd = true;
      break;
    }
    if (vals[i] > 0) maxact += vals[i] * upper[i];
  }

  HighsCDouble maxabscoef = maxact - rhs;
  if (!unbnd && maxabscoef > mip.mipdata_->feastol) {
    int ntightened = 0;
    for (int i = 0; i != len; ++i) {
      if (inds[i] >= mip.numCol() ||
          mip.variableType(inds[i]) == HighsVarType::CONTINUOUS)
        continue;

      if (vals[i] > maxabscoef) {
        HighsCDouble delta = vals[i] - maxabscoef;
        rhs -= delta * upper[i];
        vals[i] = double(maxabscoef);
        ++ntightened;
      }
    }
  }

  /* relax the right hand side a little to ensure numerical safety */
  rhs += HIGHS_CONST_TINY;
  rhs *= (1.0 + std::copysign(HIGHS_CONST_TINY, double(rhs)));

  rhs.renormalize();
  return true;
}

static bool determineCover(const HighsMipSolver& mip, std::vector<int>& inds,
                           const std::vector<double>& vals,
                           const std::vector<double>& upper,
                           const std::vector<double>& solvals,
                           const HighsCDouble& rhs, std::vector<int>& cover,
                           HighsCDouble& coverweight, HighsCDouble& lambda) {
  cover.clear();
  int len = inds.size();
  for (int j = 0; j != len; ++j) {
    if (inds[j] >= mip.numCol()) continue;

    if (solvals[j] > mip.mipdata_->feastol &&
        mip.variableType(inds[j]) != HighsVarType::CONTINUOUS) {
      cover.push_back(j);
    }
  }

  if (rhs <= 1e-5) return false;

  std::sort(cover.begin(), cover.end(), [&](int a, int b) {
    if (solvals[a] * upper[b] > solvals[b] * upper[a]) return true;

    if (solvals[a] * upper[b] < solvals[b] * upper[a]) return false;

    if (vals[a] * upper[a] > vals[b] * upper[b]) return true;

    return false;
  });

  coverweight = 0.0;
  int coversize = 0;
  double maxcontribution = 0;
  for (int j : cover) {
    double lambda = double(coverweight - rhs);
    if (lambda > 1e-5) break;

    double contribution = vals[j] * upper[j];
    coverweight += contribution;
    ++coversize;
    maxcontribution = std::max(contribution, maxcontribution);
  }

  if (coversize == 0) return false;

  cover.resize(coversize);

  coverweight.renormalize();
  lambda = coverweight - rhs;

  if (lambda <= 1e-5) return false;

  assert(lambda > 1e-5);

  if (maxcontribution <= lambda + 1e-5) {
    for (int j = coversize - 2; j >= 0; --j) {
      int i = cover[j];
      double contribution = vals[i] * upper[i];
      if (contribution < lambda - 1e-5) {
        coverweight -= contribution;
        lambda -= contribution;

        assert(lambda > 1e-5);
        assert(coverweight > rhs);

        cover.erase(cover.begin() + j);
        return true;
      }
    }
  }

  return true;
}

bool cmirCutGenerationHeuristic(const HighsMipSolver& mip,
                                const std::vector<double>& upper,
                                HighsCDouble& efficacy,
                                std::vector<double>& solvals,
                                std::vector<int8_t>& complementation,
                                std::vector<int>& inds,
                                std::vector<double>& vals, HighsCDouble& rhs) {
  int len = inds.size();

  /* first complement every variable to its closest bound */
  for (int i = 0; i != len; ++i) {
    if (upper[i] == HIGHS_CONST_INF) continue;
    if (upper[i] - solvals[i] < solvals[i]) {
      complementation[i] = -complementation[i];
      rhs -= upper[i] * vals[i];
      vals[i] = -vals[i];
      solvals[i] = upper[i] - solvals[i];
    }
  }

  std::vector<double> deltas;

  HighsCDouble continuouscontribution = 0.0;
  HighsCDouble continuoussqrnorm = 0.0;
  std::vector<int> integerinds;
  integerinds.reserve(inds.size());
  double maxabsdelta = 0.0;

  for (int i = 0; i != len; ++i) {
    if (inds[i] < mip.numCol() &&
        mip.variableType(inds[i]) != HighsVarType::CONTINUOUS) {
      integerinds.push_back(i);
      if (solvals[i] > mip.mipdata_->feastol) {
        double delta = std::abs(vals[i]);
        if (delta <= 1e-4 || delta >= 1e4) continue;
        maxabsdelta = std::max(maxabsdelta, delta);
        deltas.push_back(delta);
      }
    } else {
      continuouscontribution += vals[i] * solvals[i];
      continuoussqrnorm += vals[i] * vals[i];
    }
  }

  if (maxabsdelta + 1.0 > 1e-4 && maxabsdelta + 1.0 < 1e4)
    deltas.push_back(maxabsdelta + 1.0);
  deltas.push_back(1.0);

  if (deltas.empty()) return false;

  std::sort(deltas.begin(), deltas.end());
  double curdelta = deltas[0];
  for (size_t i = 1; i < deltas.size(); ++i) {
    if (deltas[i] - curdelta <= mip.mipdata_->feastol)
      deltas[i] = 0.0;
    else
      curdelta = deltas[i];
  }

  deltas.erase(std::remove(deltas.begin(), deltas.end(), 0.0), deltas.end());
  double bestdelta = -1;
  double bestefficacy = 0.0;

  for (double delta : deltas) {
    HighsCDouble scale = 1.0 / HighsCDouble(delta);
    HighsCDouble scalrhs = rhs * scale;
    double downrhs = std::floor(double(scalrhs));

    HighsCDouble f0 = scalrhs - downrhs;
    if (f0 < 0.01 || f0 > 0.99) continue;
    HighsCDouble oneoveroneminusf0 = 1.0 / (1.0 - f0);
    if (double(oneoveroneminusf0) * double(scale) > 1e4) continue;

    HighsCDouble sqrnorm = scale * scale * continuoussqrnorm;
    HighsCDouble viol = continuouscontribution * oneoveroneminusf0 - scalrhs;

    for (int j : integerinds) {
      HighsCDouble scalaj = vals[j] * scale;
      double downaj = std::floor(double(scalaj));
      HighsCDouble fj = scalaj - downaj;
      double aj;
      if (fj > f0)
        aj = double(downaj + fj - f0);
      else
        aj = downaj;

      viol += aj * solvals[j];
      sqrnorm += aj * aj;
    }

    double efficacy = double(viol / sqrt(sqrnorm));
    if (efficacy > bestefficacy) {
      bestdelta = delta;
      bestefficacy = efficacy;
    }
  }

  if (bestdelta == -1) return false;

  /* try if multiplying best delta by 2 4 or 8 gives a better efficacy */
  for (int k = 1; k <= 3; ++k) {
    double delta = bestdelta * (1 << k);
    if (delta <= 1e-4 || delta >= 1e4) continue;
    HighsCDouble scale = 1.0 / HighsCDouble(delta);
    HighsCDouble scalrhs = rhs * scale;
    double downrhs = std::floor(double(scalrhs));
    HighsCDouble f0 = scalrhs - downrhs;
    if (f0 < 0.01 || f0 > 0.99) continue;

    HighsCDouble oneoveroneminusf0 = 1.0 / (1.0 - f0);
    if (double(oneoveroneminusf0) * double(scale) > 1e4) continue;

    HighsCDouble sqrnorm = scale * scale * continuoussqrnorm;
    HighsCDouble viol = continuouscontribution * oneoveroneminusf0 - scalrhs;

    for (int j : integerinds) {
      HighsCDouble scalaj = vals[j] * scale;
      double downaj = std::floor(double(scalaj));
      HighsCDouble fj = scalaj - downaj;
      double aj;
      if (fj > f0)
        aj = double(downaj + fj - f0);
      else
        aj = downaj;

      viol += aj * solvals[j];
      sqrnorm += aj * aj;
    }

    double efficacy = double(viol / sqrt(sqrnorm));
    if (efficacy > bestefficacy) {
      bestdelta = delta;
      bestefficacy = efficacy;
    }
  }

  if (bestdelta == -1) return false;

  // try to flip complementation of integers to increase efficacy

  for (int k : integerinds) {
    if (upper[k] == HIGHS_CONST_INF) continue;

    complementation[k] = -complementation[k];
    rhs -= upper[k] * vals[k];
    vals[k] = -vals[k];
    solvals[k] = upper[k] - solvals[k];

    double delta = bestdelta;
    HighsCDouble scale = 1.0 / HighsCDouble(delta);
    HighsCDouble scalrhs = rhs * scale;
    double downrhs = std::floor(double(scalrhs));

    HighsCDouble f0 = scalrhs - downrhs;
    if (f0 < 0.01 || f0 > 0.99) {
      complementation[k] = -complementation[k];
      rhs -= upper[k] * vals[k];
      vals[k] = -vals[k];
      solvals[k] = upper[k] - solvals[k];
      continue;
    }

    HighsCDouble oneoveroneminusf0 = 1.0 / (1.0 - f0);
    if (double(oneoveroneminusf0) * double(scale) > 1e4) {
      complementation[k] = -complementation[k];
      rhs -= upper[k] * vals[k];
      vals[k] = -vals[k];
      solvals[k] = upper[k] - solvals[k];
      continue;
    }

    HighsCDouble sqrnorm = scale * scale * continuoussqrnorm;
    HighsCDouble viol = continuouscontribution * oneoveroneminusf0 - scalrhs;

    for (int j : integerinds) {
      HighsCDouble scalaj = vals[j] * scale;
      double downaj = std::floor(double(scalaj));
      HighsCDouble fj = scalaj - downaj;
      double aj;
      if (fj > f0)
        aj = double(downaj + fj - f0);
      else
        aj = downaj;

      viol += aj * solvals[j];
      sqrnorm += aj * aj;
    }

    double efficacy = double(viol / sqrt(sqrnorm));
    if (efficacy > bestefficacy) {
      bestefficacy = efficacy;
    } else {
      complementation[k] = -complementation[k];
      rhs -= upper[k] * vals[k];
      vals[k] = -vals[k];
      solvals[k] = upper[k] - solvals[k];
    }
  }

  if (bestefficacy <= efficacy) return false;
  efficacy = bestefficacy;
  HighsCDouble scale = 1.0 / HighsCDouble(bestdelta);
  HighsCDouble scalrhs = rhs * scale;
  double downrhs = std::floor(double(scalrhs));

  HighsCDouble f0 = scalrhs - downrhs;
  HighsCDouble oneoveroneminusf0 = 1.0 / (1.0 - f0);

  rhs = downrhs * bestdelta;

  for (int j = 0; j != len; ++j) {
    if (inds[j] >= mip.numCol() ||
        mip.variableType(inds[j]) == HighsVarType::CONTINUOUS) {
      vals[j] = double(vals[j] * oneoveroneminusf0);
    } else {
      HighsCDouble scalaj = scale * vals[j];
      double downaj = std::floor(double(scalaj));
      HighsCDouble fj = scalaj - downaj;
      HighsCDouble aj;
      if (fj > f0)
        aj = downaj + fj - f0;
      else
        aj = downaj;
      vals[j] = double(aj * bestdelta);
    }
  }

  return true;
}

bool generateCut(const HighsMipSolver& mip, const std::vector<double>& upper,
                 int nbin, int nint, int ncont, int nunbndint,
                 std::vector<double>& solvals,
                 std::vector<int8_t>& complementation, std::vector<int>& inds,
                 std::vector<double>& vals, HighsCDouble& rhs, bool& integral) {
  std::vector<int> cover;

  bool success;
  integral = false;
#if 0
  auto tmpinds = inds;
  auto tmpvals = vals;
  auto tmpcompl = complementation;
  HighsCDouble tmprhs = rhs;
  HighsCDouble efficacy = feastol;
#endif
  if (nunbndint != 0) {
    HighsCDouble efficacy = 0;
    success = cmirCutGenerationHeuristic(mip, upper, efficacy, solvals,
                                         complementation, inds, vals, rhs);
  } else {
    HighsCDouble coverweight;
    HighsCDouble lambda;

    success = determineCover(mip, inds, vals, upper, solvals, rhs, cover,
                             coverweight, lambda);

    if (success) {
      if (nint == 0 && ncont == 0) {
        integral = true;
        success = separatePureBinaryKnapsackCover(solvals, coverweight, lambda,
                                                  cover, inds, vals, rhs,
                                                  mip.mipdata_->feastol);
      } else if (nint == 0)
        success = separateMixedBinaryKnapsackCover(
            mip, solvals, coverweight, lambda, cover, inds, vals, rhs);
      else
        success = separateMixedIntegerKnapsackCover(
            mip, solvals, upper, coverweight, lambda, cover, inds, vals, rhs);
    }
  }

#if 0
  if (success) {
    efficacy = -rhs;
    HighsCDouble sqrnorm = 0.0;
    for (size_t i = 0; i != inds.size(); ++i) {
      sqrnorm += vals[i] * vals[i];
      efficacy += vals[i] * solvals[i];
    }
    efficacy = efficacy / sqrt(sqrnorm);
  }

  if (cmirCutGenerationHeuristic(mip, upper, efficacy, solvals, tmpcompl,
                                 tmpinds, tmpvals, tmprhs)) {
    inds = tmpinds;
    vals = tmpvals;
    complementation = tmpcompl;
    rhs = tmprhs;
    return true;
  }
#endif
  return success;
}

class ImpliedBounds {
  HighsImplications& implications;

 public:
  ImpliedBounds(HighsImplications& implications) : implications(implications) {}

  void separateImplBounds(const HighsLpRelaxation& lp, HighsCutPool& cutpool,
                          HighsDomain& propdomain) {
    HighsDomain& globaldomain = lp.getMipSolver().mipdata_->domain;

    const double feastol = lp.getMipSolver().mipdata_->feastol;
    int inds[2];
    double vals[2];
    double rhs;
    const auto& sol = lp.getLpSolver().getSolution().col_value;
    lp.getMipSolver().mipdata_->cliquetable.cleanupFixed(globaldomain);
    if (globaldomain.infeasible()) return;
    int numboundchgs = 0;

    // first do probing on all candidates that have not been probed yet
    for (std::pair<int, double> fracint : lp.getFractionalIntegers()) {
      int col = fracint.first;
      if (globaldomain.colLower_[col] != 0.0 ||
          globaldomain.colUpper_[col] != 1.0)
        continue;

      if (implications.runProbing(col, numboundchgs)) {
        ++numboundchgs;
        if (globaldomain.infeasible()) return;
      }
    }

    for (std::pair<int, double> fracint : lp.getFractionalIntegers()) {
      int col = fracint.first;
      // skip non binary variables
      if (globaldomain.colLower_[col] != 0.0 ||
          globaldomain.colUpper_[col] != 1.0)
        continue;

      bool infeas;
      const HighsDomainChange* implics = nullptr;
      int nimplics = implications.getImplications(col, 1, implics, infeas);
      if (globaldomain.infeasible()) return;
      if (infeas) {
        vals[0] = 1.0;
        inds[0] = col;
        cutpool.addCut(lp.getMipSolver(), inds, vals, 1, 0.0);
        continue;
      }

      for (int i = 0; i != nimplics; ++i) {
        if (implics[i].boundtype == HighsBoundType::Upper) {
          if (implics[i].boundval + feastol >=
              globaldomain.colUpper_[implics[i].column])
            continue;

          vals[0] = 1.0;
          inds[0] = implics[i].column;
          vals[1] =
              globaldomain.colUpper_[implics[i].column] - implics[i].boundval;
          inds[1] = col;
          rhs = globaldomain.colUpper_[implics[i].column];

        } else {
          if (implics[i].boundval - feastol <=
              globaldomain.colLower_[implics[i].column])
            continue;

          vals[0] = -1.0;
          inds[0] = implics[i].column;
          vals[1] =
              globaldomain.colLower_[implics[i].column] - implics[i].boundval;
          inds[1] = col;
          rhs = -globaldomain.colLower_[implics[i].column];
        }

        double viol = sol[inds[0]] * vals[0] + sol[inds[1]] * vals[1] - rhs;

        if (viol > feastol) {
          // printf("added implied bound cut to pool\n");
          cutpool.addCut(lp.getMipSolver(), inds, vals, 2, rhs);
        }
      }

      nimplics = implications.getImplications(col, 0, implics, infeas);
      if (globaldomain.infeasible()) return;
      if (infeas) {
        vals[0] = -1.0;
        inds[0] = col;
        cutpool.addCut(lp.getMipSolver(), inds, vals, 1, -1.0);
        continue;
      }

      for (int i = 0; i != nimplics; ++i) {
        if (implics[i].boundtype == HighsBoundType::Upper) {
          if (implics[i].boundval + feastol >=
              globaldomain.colUpper_[implics[i].column])
            continue;

          vals[0] = 1.0;
          inds[0] = implics[i].column;
          vals[1] =
              implics[i].boundval - globaldomain.colUpper_[implics[i].column];
          inds[1] = col;
          rhs = implics[i].boundval;
        } else {
          if (implics[i].boundval - feastol <=
              globaldomain.colLower_[implics[i].column])
            continue;

          vals[0] = -1.0;
          inds[0] = implics[i].column;
          vals[1] =
              globaldomain.colLower_[implics[i].column] - implics[i].boundval;
          inds[1] = col;
          rhs = -implics[i].boundval;
        }

        double viol = sol[inds[0]] * vals[0] + sol[inds[1]] * vals[1] - rhs;

        if (viol > feastol) {
          // printf("added implied bound cut to pool\n");
          cutpool.addCut(lp.getMipSolver(), inds, vals, 2, rhs);
        }
      }
    }
  }
};

class AggregationHeuristic {
  const HighsLpRelaxation& lprelaxation;
  const HighsDomain& domain;
  HighsCutPool& cutpool;
  HighsDomain& propdomain;

  const HighsMipSolver& mip;
  const HighsSolution& lpsol;
  const HighsLp& lp;

  std::vector<RowType> rowtype;
  std::vector<uint8_t> rowusable;
  std::vector<int> numcontinuous;

  std::vector<int> colubtype;
  std::vector<double> colubscale;
  std::vector<int> collbtype;
  std::vector<double> collbscale;
  std::vector<double> bounddistance;

  // vector sum object for bound substitution and retransformation
  HighsSparseVectorSum vectorsum;

  // current indices and values and rhs for the aggregated row we are looking
  // at
  std::vector<int> baseinds;
  std::vector<double> basevals;
  HighsCDouble baserhs;

  // position of continuous variables with positive bound distance
  std::vector<int> bounddistpos;

  // temporary indices and values for the row aggregating the next row
  std::vector<int> tmpinds;
  std::vector<double> tmpvals;

  // current transformed row we are looking at
  HighsCDouble rhs;
  std::vector<int> inds;
  std::vector<double> vals;
  std::vector<double> upper;
  std::vector<double> solvals;
  std::vector<int8_t> complementation;
  int ncont;
  int nbin;
  int nint;
  int nunbndint;
  bool freevar;

  // maximal number of aggregations we allow
  static constexpr int maxAggr = 6;
  // array for the current path and its length
  int currpath[maxAggr];
  double currpathweights[maxAggr];
  int currpathlen;

  int numcuts;

 public:
  AggregationHeuristic(const HighsLpRelaxation& lprelaxation,
                       const HighsDomain& domain, HighsCutPool& cutpool,
                       HighsDomain& propdomain)
      : lprelaxation(lprelaxation),
        domain(domain),
        cutpool(cutpool),
        propdomain(propdomain),
        mip(lprelaxation.getMipSolver()),
        lpsol(lprelaxation.getLpSolver().getSolution()),
        lp(lprelaxation.getLpSolver().getLp()) {
    numcuts = 0;
  }

  void determineRowTypes() {
    rowtype.resize(lp.numRow_);
    rowusable.resize(lp.numRow_);
    for (int i = 0; i != lp.numRow_; ++i) {
      if (lp.rowLower_[i] == lp.rowUpper_[i]) {
        rowtype[i] = RowType::Eq;
        rowusable[i] = true;
        continue;
      }

      double lowerslack = HIGHS_CONST_INF;
      double upperslack = HIGHS_CONST_INF;

      if (lp.rowLower_[i] != -HIGHS_CONST_INF)
        lowerslack = lpsol.row_value[i] - lp.rowLower_[i];

      if (lp.rowUpper_[i] != HIGHS_CONST_INF)
        upperslack = lp.rowUpper_[i] - lpsol.row_value[i];

      if (lowerslack > mip.mipdata_->feastol &&
          upperslack > mip.mipdata_->feastol)
        rowtype[i] = RowType::Unusuable;
      else if (lowerslack < upperslack)
        rowtype[i] = RowType::Geq;
      else
        rowtype[i] = RowType::Leq;

      rowusable[i] = rowtype[i] != RowType::Unusuable;
    }
  }

  void detectContinuousBounds() {
    // count number of continuous variables
    numcontinuous.assign(lp.numRow_, 0);
    for (int i = 0; i != lp.numCol_; ++i) {
      if (mip.variableType(i) != HighsVarType::CONTINUOUS) continue;

      int start = lp.Astart_[i];
      int end = lp.Astart_[i + 1];

      for (int j = start; j != end; ++j) ++numcontinuous[lp.Aindex_[j]];
    }

    // now detect variable bound constraints with one or multiple integral
    // variables and store the sparsest one for each constinuous column Also
    // compute the bound distance to the best simple or variable bound
    colubtype.assign(lp.numCol_, -1);
    colubscale.resize(lp.numCol_);
    collbtype.assign(lp.numCol_, -1);
    collbscale.resize(lp.numCol_);
    bounddistance.assign(lp.numCol_, 0.0);

    for (int i = 0; i != lp.numCol_; ++i) {
      if (mip.variableType(i) != HighsVarType::CONTINUOUS) continue;

      int start = lp.Astart_[i];
      int end = lp.Astart_[i + 1];

      int lblen = HIGHS_CONST_I_INF;
      int ublen = HIGHS_CONST_I_INF;

      for (int j = start; j != end; ++j) {
        int row = lp.Aindex_[j];
        if (numcontinuous[row] != 1) continue;

        int rowlen =
            row < mip.model_->numRow_
                ? mip.mipdata_->ARstart_[row + 1] - mip.mipdata_->ARstart_[row]
                : cutpool.getRowLength(lprelaxation.getCutIndex(row));
        if (rowlen != 2) continue;
        switch (rowtype[row]) {
          case RowType::Unusuable:
            break;
          case RowType::Leq:
            if (lp.Avalue_[j] < 0 && rowlen < lblen) {
              lblen = rowlen;
              collbtype[i] = row;
              collbscale[i] = lp.Avalue_[j];
            } else if (lp.Avalue_[j] > 0 && rowlen < ublen) {
              ublen = rowlen;
              colubtype[i] = row;
              colubscale[i] = lp.Avalue_[j];
            }
            break;
          case RowType::Geq:
            if (lp.Avalue_[j] > 0 && rowlen < lblen) {
              lblen = rowlen;
              collbtype[i] = row;
              collbscale[i] = lp.Avalue_[j];
            } else if (lp.Avalue_[j] < 0 && rowlen < ublen) {
              ublen = rowlen;
              colubtype[i] = row;
              colubscale[i] = lp.Avalue_[j];
            }
            break;
          case RowType::Eq:
            if (rowlen < ublen) {
              ublen = rowlen;
              colubtype[i] = row;
              colubscale[i] = lp.Avalue_[j];
            }
            if (rowlen < lblen) {
              lblen = rowlen;
              collbtype[i] = row;
              collbscale[i] = lp.Avalue_[j];
            }
        }
      }

      double lbdist;
      double ubdist;

      assert(collbtype[i] == -1 || rowusable[collbtype[i]]);
      assert(colubtype[i] == -1 || rowusable[colubtype[i]]);
      if (collbtype[i] != -1) {
        lbdist = 0.0;
        rowusable[collbtype[i]] = false;
      } else
        lbdist = domain.colLower_[i] != -HIGHS_CONST_INF
                     ? lpsol.col_value[i] - domain.colLower_[i]
                     : HIGHS_CONST_INF;

      if (colubtype[i] != -1) {
        ubdist = 0.0;
        rowusable[colubtype[i]] = false;
      } else
        ubdist = domain.colUpper_[i] != HIGHS_CONST_INF
                     ? domain.colUpper_[i] - lpsol.col_value[i]
                     : HIGHS_CONST_INF;

      bounddistance[i] = std::min(lbdist, ubdist);
      if (bounddistance[i] < mip.mipdata_->feastol) bounddistance[i] = 0;
    }

    // todo :mark all continuous variables that have fractional integer
    // variables in their variable bound constraints then mark all rows that
    // contain fractional integer variables or marked continuous variables

    // now only count the number of continuous variables
    // not at their bound
    numcontinuous.assign(lp.numRow_, 0);
    for (int i = 0; i != lp.numCol_; ++i) {
      if (bounddistance[i] == 0.0) continue;

      int start = lp.Astart_[i];
      int end = lp.Astart_[i + 1];

      for (int j = start; j != end; ++j) ++numcontinuous[lp.Aindex_[j]];
    }

    // mark each of the variables that has an equation row where it is the
    // only continuous variable by flipping its bound distance such variables
    // will be substituted in the aggregation heuristic first
  }

  void setupStartRow(int row) {
    baseinds.clear();
    basevals.clear();

    currpath[0] = row;
    currpathlen = 1;

    const int* Rindex;
    const double* Rvalue;
    int Rlen;
    if (row < mip.model_->numRow_)
      mip.mipdata_->getRow(row, Rlen, Rindex, Rvalue);
    else
      cutpool.getCut(lprelaxation.getCutIndex(row), Rlen, Rindex, Rvalue);

    baseinds.insert(baseinds.end(), Rindex, Rindex + Rlen);

    if (rowtype[row] == RowType::Leq) {
      baserhs = lp.rowUpper_[row];
      basevals.insert(basevals.end(), Rvalue, Rvalue + Rlen);
      currpathweights[0] = 1.0;
    } else {
      baserhs = -lp.rowLower_[row];
      basevals.resize(Rlen);
      currpathweights[0] = -1.0;
      for (int i = 0; i != Rlen; ++i) basevals[i] = -Rvalue[i];
    }
  }

  void computeTransformedRow(bool negate) {
    vectorsum.clear();
    bounddistpos.clear();
    inds.clear();
    vals.clear();
    upper.clear();
    solvals.clear();
    complementation.clear();

    freevar = false;

    vectorsum.setDimension(lp.numCol_);

    double basescale;

    if (negate) {
      for (int i = 0; i < currpathlen; ++i) {
        if (rowtype[currpath[i]] == RowType::Eq) continue;
        // add a slack variable for the rows in the path if they are not
        // equations
        inds.push_back(lp.numCol_ + currpath[i]);
        vals.push_back(-std::abs(currpathweights[i]));
        upper.push_back(HIGHS_CONST_INF);
        complementation.push_back(0);
        solvals.push_back(0);
      }
      basescale = -1.0;
    } else
      basescale = 1.0;

    rhs = basescale * baserhs;

    for (size_t j = 0; j != baseinds.size(); ++j) {
      int col = baseinds[j];
      double baseval = basescale * basevals[j];

      // add things to the vectorsum
      // store solution values and upper bounds
      // the coefficient values area added up with the sparse vector sum
      // as variable bound usage can alter things

      // we first need to complement the continuous variables
      if (mip.variableType(col) != HighsVarType::CONTINUOUS) {
        vectorsum.add(col, baseval);
        continue;
      }

      // complement positive if possible without slack using a simple or
      // variable bound otherwise use the tightest bound and prefer simple
      // bounds adjust the right hand side accordingly
      if (bounddistance[col] > mip.mipdata_->feastol) bounddistpos.push_back(j);

      if (bounddistance[col] == HIGHS_CONST_INF) freevar = true;

      if (freevar) continue;

      double simplelbdist = domain.colLower_[col] != -HIGHS_CONST_INF
                                ? lpsol.col_value[col] - domain.colLower_[col]
                                : HIGHS_CONST_INF;
      double simpleubdist = domain.colUpper_[col] != HIGHS_CONST_INF
                                ? domain.colUpper_[col] - lpsol.col_value[col]
                                : HIGHS_CONST_INF;

      if (baseval < 0) {
        if (colubtype[col] != -1) {
          // use variable upper bound
          int UBrow = colubtype[col];
          double UBscale = colubscale[col];
          int UBlen;
          const int* UBindex;
          const double* UBvalue;
          double UBconst;

          if (UBscale < 0)
            UBconst = lp.rowLower_[UBrow];
          else
            UBconst = lp.rowUpper_[UBrow];

          if (UBrow < mip.numRow()) {
            mip.mipdata_->getRow(UBrow, UBlen, UBindex, UBvalue);
          } else {
            int cut = lprelaxation.getCutIndex(UBrow);
            cutpool.getCut(cut, UBlen, UBindex, UBvalue);
          }

          HighsCDouble scale = HighsCDouble(-baseval) / UBscale;

          // the upper bound constraint adds a conmtinuous slack variable
          // with positive coefficient that we relax out

          rhs += scale * UBconst;

          for (int k = 0; k != UBlen; ++k) {
            if (UBindex[k] == col)  // the column cancels out
            {
              assert(std::abs(double(baseval + scale * UBvalue[k])) < 1e-10);
              continue;
            }
            assert(mip.variableType(UBindex[k]) != HighsVarType::CONTINUOUS);
            vectorsum.add(UBindex[k], scale * UBvalue[k]);
          }
        } else if (simpleubdist - mip.mipdata_->feastol <= bounddistance[col]) {
          // use  the simple upper bound for complementation and then relax
          // the positive continuous variable as it has a positive
          // coefficient
          complementWithUpper(baseval, domain.colUpper_[col], rhs);
        } else if (simplelbdist - mip.mipdata_->feastol <= bounddistance[col]) {
          inds.push_back(col);
          complementation.push_back(1);
          if (domain.colUpper_[col] == HIGHS_CONST_INF)
            upper.push_back(HIGHS_CONST_INF);
          else
            upper.push_back(domain.colUpper_[col] - domain.colLower_[col]);
          vals.push_back(
              complementWithLower(baseval, domain.colLower_[col], rhs));
          solvals.push_back(lpsol.col_value[col] - domain.colLower_[col]);
        } else {
          assert(collbtype[col] != -1);

          // use variable lower bound
          int LBrow = collbtype[col];
          double LBscale = collbscale[col];
          int LBlen;
          const int* LBindex;
          const double* LBvalue;
          double LBconst;

          if (LBscale < 0)
            LBconst = lp.rowUpper_[LBrow];
          else
            LBconst = lp.rowLower_[LBrow];

          if (LBrow < mip.numRow()) {
            mip.mipdata_->getRow(LBrow, LBlen, LBindex, LBvalue);
          } else {
            int cut = lprelaxation.getCutIndex(LBrow);
            cutpool.getCut(cut, LBlen, LBindex, LBvalue);
          }

          HighsCDouble scale = HighsCDouble(-baseval) / LBscale;

          rhs += scale * LBconst;

          // the lb constraint does add a slack variable with negative sign
          // as the inequality orientations do not match
          if (LBscale > 0)
            vals.push_back(-double(scale));
          else
            vals.push_back(double(scale));
          assert(vals.back() < 0);

          inds.push_back(lp.numCol_ + LBrow);
          complementation.push_back(0);
          upper.push_back(HIGHS_CONST_INF);
          solvals.push_back(0);

          for (int k = 0; k != LBlen; ++k) {
            if (LBindex[k] == col)  // the column cancels out
            {
              assert(std::abs(double(baseval + scale * LBvalue[k])) < 1e-10);
              continue;
            }

            assert(mip.variableType(LBindex[k]) != HighsVarType::CONTINUOUS);
            vectorsum.add(LBindex[k], scale * LBvalue[k]);
          }
        }
      } else {
        assert(baseval > 0);
        if (collbtype[col] != -1) {
          int LBrow = collbtype[col];
          double LBscale = collbscale[col];
          int LBlen;
          const int* LBindex;
          const double* LBvalue;
          double LBconst;

          if (LBscale < 0)
            LBconst = lp.rowUpper_[LBrow];
          else
            LBconst = lp.rowLower_[LBrow];

          if (LBrow < mip.numRow()) {
            mip.mipdata_->getRow(LBrow, LBlen, LBindex, LBvalue);
          } else {
            int cut = lprelaxation.getCutIndex(LBrow);
            cutpool.getCut(cut, LBlen, LBindex, LBvalue);
          }

          HighsCDouble scale = HighsCDouble(-baseval) / LBscale;

          rhs += scale * LBconst;

          assert(int(rowtype[LBrow]) * double(scale) >= 0);

          // the lower bound constraint does not add a slack variable we
          // need to remember as the inequalities match in orientation and
          // the positive slack variable is relaxed out

          for (int k = 0; k != LBlen; ++k) {
            if (LBindex[k] == col)  // the column cancels out
            {
              assert(std::abs(double(baseval + scale * LBvalue[k])) < 1e-10);
              continue;
            }
            assert(mip.variableType(LBindex[k]) != HighsVarType::CONTINUOUS);
            vectorsum.add(LBindex[k], scale * LBvalue[k]);
          }
        } else if (simplelbdist - mip.mipdata_->feastol <= bounddistance[col]) {
          // use  the simple lower bound for complementation and then relax
          // the positive continuous variable as it has a positive
          // coefficient
          complementWithLower(baseval, domain.colLower_[col], rhs);
        } else if (simpleubdist - mip.mipdata_->feastol <= bounddistance[col]) {
          // use the simple upper bound and add the negative complemented
          // variable
          inds.push_back(col);
          complementation.push_back(-1);
          if (domain.colLower_[col] == -HIGHS_CONST_INF)
            upper.push_back(HIGHS_CONST_INF);
          else
            upper.push_back(domain.colUpper_[col] - domain.colLower_[col]);
          solvals.push_back(domain.colUpper_[col] - lpsol.col_value[col]);

          vals.push_back(
              complementWithUpper(baseval, domain.colUpper_[col], rhs));
        } else {
          assert(colubtype[col] != -1);
          // use variable upper bound
          int UBrow = colubtype[col];
          double UBscale = colubscale[col];
          int UBlen;
          const int* UBindex;
          const double* UBvalue;
          double UBconst;

          if (UBscale < 0)
            UBconst = lp.rowLower_[UBrow];
          else
            UBconst = lp.rowUpper_[UBrow];

          if (UBrow < mip.numRow()) {
            mip.mipdata_->getRow(UBrow, UBlen, UBindex, UBvalue);
          } else {
            int cut = lprelaxation.getCutIndex(UBrow);
            cutpool.getCut(cut, UBlen, UBindex, UBvalue);
          }

          HighsCDouble scale = HighsCDouble(-baseval) / UBscale;
          // the ub constraint does add a slack variable as the inequality
          // orientations do not match

          if (UBscale > 0)
            vals.push_back(double(scale));
          else
            vals.push_back(-double(scale));
          assert(vals.back() < 0);

          inds.push_back(lp.numCol_ + UBrow);
          complementation.push_back(0);
          upper.push_back(HIGHS_CONST_INF);
          solvals.push_back(0);
          rhs += scale * UBconst;

          for (int k = 0; k != UBlen; ++k) {
            if (UBindex[k] == col)  // the column cancels out
            {
              assert(std::abs(double(baseval + scale * UBvalue[k])) < 1e-10);
              continue;
            }
            assert(mip.variableType(UBindex[k]) != HighsVarType::CONTINUOUS);
            vectorsum.add(UBindex[k], scale * UBvalue[k]);
          }
        }
      }
    }

    // so far we only added continuous variables as the values of integer
    // variables could still change due to variable bound usage
    ncont = inds.size();
    nbin = 0;
    nint = 0;
    nunbndint = 0;

    // now we add the integer variables and complement them to be positive
    for (int col : vectorsum.getNonzeros()) {
      double val = vectorsum.getValue(col);

      if (std::abs(val) <= 1e-10) continue;

      if (std::abs(val) <= mip.mipdata_->feastol) {
        if (val < 0)
          rhs -= domain.colUpper_[col] * val;
        else
          rhs -= domain.colLower_[col] * val;
      }

      assert(mip.variableType(col) != HighsVarType::CONTINUOUS);
      assert(domain.colLower_[col] != -HIGHS_CONST_INF ||
             domain.colUpper_[col] != HIGHS_CONST_INF);

      if (domain.colLower_[col] == -HIGHS_CONST_INF ||
          domain.colUpper_[col] == HIGHS_CONST_INF)
        upper.push_back(HIGHS_CONST_INF);
      else
        upper.push_back(
            std::floor(domain.colUpper_[col] - domain.colLower_[col] + 0.5));

      inds.push_back(col);

      if (upper.back() == 1.0)
        ++nbin;
      else if (upper.back() == HIGHS_CONST_INF)
        ++nunbndint;
      else
        ++nint;
      if (val < 0 && domain.colUpper_[col] != HIGHS_CONST_INF) {
        solvals.push_back(domain.colUpper_[col] - lpsol.col_value[col]);
        complementation.push_back(-1);
        vals.push_back(complementWithUpper(val, domain.colUpper_[col], rhs));
      } else {
        solvals.push_back(lpsol.col_value[col] - domain.colLower_[col]);
        complementation.push_back(1);
        vals.push_back(complementWithLower(val, domain.colLower_[col], rhs));
      }
    }

    // perform coefficient tightening on the transformed row
    int len = inds.size();
    HighsCDouble maxact = 0.0;
    bool unbnd = false;
    for (int i = 0; i != len; ++i) {
      if (upper[i] == HIGHS_CONST_INF) {
        unbnd = true;
        break;
      }
      if (vals[i] > 0) maxact += vals[i] * upper[i];
    }

    HighsCDouble maxabscoef = maxact - rhs;
    if (!unbnd && maxabscoef > mip.mipdata_->feastol) {
      int ntightened = 0;
      for (int i = 0; i != len; ++i) {
        if (inds[i] >= mip.numCol() ||
            mip.variableType(inds[i]) == HighsVarType::CONTINUOUS)
          continue;

        if (vals[i] > maxabscoef) {
          HighsCDouble delta = vals[i] - maxabscoef;
          rhs -= delta * upper[i];
          vals[i] = double(maxabscoef);
          ++ntightened;
        }
      }
    }
  }

  void computeCut() {
    bool cutintegral;
    bool foundcut =
        generateCut(mip, upper, nbin, nint, ncont, nunbndint, solvals,
                    complementation, inds, vals, rhs, cutintegral);

    if (foundcut) {
      vectorsum.clear();
      vectorsum.setDimension(mip.numCol());

      int cutlen = inds.size();

      for (int i = 0; i != cutlen; ++i) {
        if (vals[i] == 0.0) continue;

        if (inds[i] >= mip.numCol()) {
          int row = inds[i] - mip.numCol();

          const int* Rindex;
          const double* Rvalue;
          int Rlen;
          double Rside;

          assert(rowtype[row] == RowType::Leq || rowtype[row] == RowType::Geq);

          if (row < mip.numRow()) {
            if (rowtype[row] == RowType::Leq)
              Rside = lp.rowUpper_[row];
            else
              Rside = -lp.rowLower_[row];
            mip.mipdata_->getRow(row, Rlen, Rindex, Rvalue);
          } else {
            assert(rowtype[row] == RowType::Leq);
            int cut = lprelaxation.getCutIndex(row);
            cutpool.getCut(cut, Rlen, Rindex, Rvalue);
            Rside = lp.rowUpper_[row];
          }

          rhs -= vals[i] * Rside;

          double slackval = -int(rowtype[row]) * vals[i];

          for (int j = 0; j != Rlen; ++j)
            vectorsum.add(Rindex[j], slackval * Rvalue[j]);

          continue;
        }

        if (complementation[i] == 1)
          rhs += vals[i] * domain.colLower_[inds[i]];
        else if (complementation[i] == -1) {
          vals[i] = -vals[i];
          rhs += vals[i] * domain.colUpper_[inds[i]];
        }

        vectorsum.add(inds[i], vals[i]);
      }

      inds.clear();
      vals.clear();

      vectorsum.sort();

      for (int col : vectorsum.getNonzeros()) {
        double val = vectorsum.getValue(col);

        if (std::abs(val) > 1e-10) {
          vals.push_back(val);
          inds.push_back(col);
        }
      }

      double upper = double(rhs);

      cutpool.addCut(mip, inds.data(), vals.data(), inds.size(), upper,
                     cutintegral);

      ++numcuts;
    }
  }

  bool aggregateNextRow() {
    if (currpathlen == maxAggr || bounddistpos.empty()) return false;

    // sort by decreasing bound distance
    std::sort(bounddistpos.begin(), bounddistpos.end(), [&](int a, int b) {
      int col1 = baseinds[a];
      int col2 = baseinds[b];
      return bounddistance[col1] > bounddistance[col2];
    });

    int nextaggrow = -1;
    HighsCDouble nextaggscale;
    int naggbounddist = HIGHS_CONST_I_INF;

    for (int pos : bounddistpos) {
      int col = baseinds[pos];

      int start = lp.Astart_[col];
      int end = lp.Astart_[col + 1];

      for (int j = start; j != end; ++j) {
        int row = lp.Aindex_[j];
        // if the row is marked as unusable or the inequality orientation does
        // not match we skip it
        if (!rowusable[row] ||
            basevals[pos] * lp.Avalue_[j] * int(rowtype[row]) > 0.0)
          continue;

        bool incurrentpath = false;
        for (int k = 0; k < currpathlen; ++k) {
          if (currpath[k] == row) {
            incurrentpath = true;
            break;
          }
        }

        if (incurrentpath) continue;

        bool nextaggusedbefore = currpath[0] > nextaggrow;
        bool rowusedbefore = currpath[0] > row;

        if ((nextaggusedbefore && !rowusedbefore) ||
            (nextaggusedbefore == rowusedbefore &&
             numcontinuous[row] < naggbounddist)) {
          nextaggscale = HighsCDouble(-basevals[pos]) / lp.Avalue_[j];
          nextaggrow = row;
          naggbounddist = numcontinuous[row];
          if (numcontinuous[row] == 1 && rowtype[row] == RowType::Eq) break;
        }
      }

      if (nextaggrow != -1) break;
    }

    if (nextaggrow == -1) return false;

    const int* nextRindex;
    const double* nextRvalue;
    int nextRlen;

    if (rowtype[nextaggrow] == RowType::Geq) {
      baserhs += nextaggscale * lp.rowLower_[nextaggrow];
      assert(nextaggscale < 0);
    } else {
      assert(rowtype[nextaggrow] == RowType::Eq || nextaggscale > 0);
      baserhs += nextaggscale * lp.rowUpper_[nextaggrow];
    }
    if (nextaggrow < mip.numRow())
      mip.mipdata_->getRow(nextaggrow, nextRlen, nextRindex, nextRvalue);
    else
      cutpool.getCut(lprelaxation.getCutIndex(nextaggrow), nextRlen, nextRindex,
                     nextRvalue);

    tmpinds.clear();
    tmpvals.clear();

    int a = 0;
    int b = 0;

    int baselen = baseinds.size();

    while (a != nextRlen && b != baselen) {
      if (nextRindex[a] < baseinds[b]) {
        tmpinds.push_back(nextRindex[a]);
        tmpvals.push_back(double(nextaggscale * nextRvalue[a]));
        ++a;

      } else if (baseinds[b] < nextRindex[a]) {
        tmpinds.push_back(baseinds[b]);
        tmpvals.push_back(basevals[b]);
        ++b;
      } else {
        assert(nextRindex[a] == baseinds[b]);
        double val = double(basevals[b] + nextaggscale * nextRvalue[a]);
        if (std::abs(val) > 1e-10) {
          tmpinds.push_back(baseinds[b]);
          tmpvals.push_back(val);
        }
        ++a;
        ++b;
      }
    }

    if (a != nextRlen) {
      tmpinds.insert(tmpinds.end(), nextRindex + a, nextRindex + nextRlen);
      do {
        tmpvals.push_back(double(nextaggscale * nextRvalue[a]));
        ++a;
      } while (a != nextRlen);
    }

    if (b != baselen) {
      tmpinds.insert(tmpinds.end(), baseinds.begin() + b, baseinds.end());
      tmpvals.insert(tmpvals.end(), basevals.begin() + b, basevals.end());
    }

    tmpinds.swap(baseinds);
    tmpvals.swap(basevals);
    currpath[currpathlen] = nextaggrow;
    currpathweights[currpathlen] = double(nextaggscale);
    ++currpathlen;
    return true;
  }

  void run() {
    determineRowTypes();
    detectContinuousBounds();

    for (int row = 0; row != lp.numRow_; ++row) {
      if (!rowusable[row]) continue;

      int currnumcuts = numcuts;
      setupStartRow(row);

      do {
        computeTransformedRow(false);

        if (!freevar) {
          computeCut();
          computeTransformedRow(true);
          computeCut();

          if (currnumcuts < numcuts) {
            break;
          }
        }

        if (!aggregateNextRow()) break;
      } while (true);
    }
  }
};

static void doSeparate(const HighsDomain& domain, const HighsLpRelaxation& lp,
                       const HighsSolution& lpsol,
                       const HighsSeparation::BaseRows& baserows,
                       HighsCutPool& cutpool, HighsDomain& propdomain) {
  std::vector<int8_t> complementation;
  std::vector<int> inds;
  std::vector<double> vals;
  std::vector<double> upper;
  std::vector<double> solvals;
  HighsCDouble rhs;

  const HighsMipSolver& mip = lp.getMipSolver();
  int numbaserows = baserows.rhs_.size();

  for (int i = 0; i != numbaserows; ++i) {
    bool success;
    bool cutintegral;

    int start = baserows.ARstart_[i];
    int end = baserows.ARstart_[i + 1];

    int nbin;
    int nint;
    int ncont;
    int nunbndint;

    rhs = baserows.rhs_[i];
    success = transformBaseEquation(
        mip, domain, lpsol, &baserows.ARindex_[start],
        &baserows.ARvalue_[start], end - start, 1.0, rhs, complementation, inds,
        vals, upper, solvals, nbin, nint, ncont, nunbndint);
    if (success)
      success = generateCut(mip, upper, nbin, nint, ncont, nunbndint, solvals,
                            complementation, inds, vals, rhs, cutintegral);

    if (success) {
      baserows.retransformAndAddCut(domain, lp, inds, vals, complementation,
                                    rhs, cutpool, propdomain, cutintegral);
    }

    rhs = baserows.rhs_[i];
    success = transformBaseEquation(
        mip, domain, lpsol, &baserows.ARindex_[start],
        &baserows.ARvalue_[start], end - start, -1.0, rhs, complementation,
        inds, vals, upper, solvals, nbin, nint, ncont, nunbndint);

    if (success)
      success = generateCut(mip, upper, nbin, nint, ncont, nunbndint, solvals,
                            complementation, inds, vals, rhs, cutintegral);

    if (success) {
      baserows.retransformAndAddCut(domain, lp, inds, vals, complementation,
                                    rhs, cutpool, propdomain, cutintegral);
    }
  }
}

static void tableauaggregator(HighsLpRelaxation& lp,
                              const HighsCutPool& cutpool,
                              HighsSeparation::BaseRows& baserows) {
  std::vector<int> basisinds;
  int numrow = lp.getNumLpRows();
  basisinds.resize(numrow);
  lp.getLpSolver().getBasicVariables(basisinds.data());

  std::vector<int> aggrinds;
  std::vector<double> aggrvals;
  int naggrinds;

  const HighsMipSolver& mip = lp.getMipSolver();

  aggrinds.resize(numrow);
  aggrvals.resize(numrow);
  baserows.slacktype_.resize(numrow);

  for (int i = 0; i != numrow; ++i) {
    if (lp.rowLower(i) == lp.rowUpper(i)) {
      baserows.slacktype_[i] = 0;
      continue;
    }

    switch (lp.getLpSolver().getBasis().row_status[i]) {
      case HighsBasisStatus::LOWER:
        baserows.slacktype_[i] = -1;
        break;
      case HighsBasisStatus::UPPER:
        baserows.slacktype_[i] = 1;
        break;
      default:
        continue;
    }
  }

  baserows.ARstart_.push_back(0);

  auto& lpsol = lp.getLpSolver().getSolution();

  for (int i = 0; i != int(basisinds.size()); ++i) {
    if (basisinds[i] < 0) {
      continue;
      // todo, cuts off solutions
      int row = -basisinds[i] - 1;

      bool rowintegral;
      if (row < mip.numRow())
        rowintegral = mip.mipdata_->rowintegral[row];
      else
        rowintegral = cutpool.cutIsIntegral(lp.getCutIndex(row));

      if (!rowintegral) continue;

      double solval = lpsol.row_value[row];
      double frac = std::abs(floor(solval + 0.5) - solval);
      if (frac < 1e-1) continue;
    } else {
      int col = basisinds[i];
      if (mip.variableType(col) == HighsVarType::CONTINUOUS) continue;
      double solval = lpsol.col_value[col];

      double frac = std::abs(std::floor(solval + 0.5) - solval);
      if (frac < 1e-4) continue;
    }

    if (lp.getLpSolver().getBasisInverseRow(i, aggrvals.data(), &naggrinds,
                                            aggrinds.data()) != HighsStatus::OK)
      continue;

    baserows.addAggregation(lp, cutpool, aggrvals.data(), aggrinds.data(),
                            naggrinds);
  }
}

void HighsSeparation::BaseRows::addAggregation(const HighsLpRelaxation& lp,
                                               const HighsCutPool& cutpool,
                                               const double* aggrvals,
                                               const int* aggrinds,
                                               int naggrinds) {
  HighsCDouble rhs = 0.0;

  const HighsMipSolver& mip = lp.getMipSolver();

  vectorsum.setDimension(mip.numCol() + lp.getNumLpRows());

  assert(std::all_of(vectorsum.nonzeroflag.begin(), vectorsum.nonzeroflag.end(),
                     [](uint8_t f) { return f == 0; }));

  double sum = 0.0;
  for (int k = 0; k != naggrinds; ++k) {
    double val = std::abs(aggrvals[aggrinds[k]]);
    sum += val;
  }

  int expscal;
  std::frexp(sum, &expscal);

  for (int k = 0; k != naggrinds; ++k) {
    int j = aggrinds[k];
    double aggval = std::ldexp(aggrvals[j], -expscal);
    if (std::abs(aggval) <= mip.mipdata_->feastol) continue;

    int rowlen;
    const int* rowinds;
    const double* rowvals;
    if (j < mip.numRow()) {
      mip.mipdata_->getRow(j, rowlen, rowinds, rowvals);
    } else {
      int cut = lp.getCutIndex(j);
      cutpool.getCut(cut, rowlen, rowinds, rowvals);
    }

    HighsCDouble scale = aggval;

    assert(slacktype_[j] != 0 || lp.rowLower(j) == lp.rowUpper(j));

    if (slacktype_[j] != 0) {
      vectorsum.set(mip.numCol() + j, slacktype_[j] * scale);
    }

    if (slacktype_[j] == -1)
      rhs += scale * lp.rowLower(j);
    else
      rhs += scale * lp.rowUpper(j);

    for (int k = 0; k != rowlen; ++k)
      vectorsum.add(rowinds[k], scale * rowvals[k]);
  }

  double minval = HIGHS_CONST_INF;
  double maxval = 0;
  int len = 0;
  for (int j : vectorsum.getNonzeros()) {
    double val = std::abs(vectorsum.getValue(j));

    if (val <= 1e-12) {
      vectorsum.chgValue(j, 0.0);
      continue;
    }

    ++len;
    if (val < minval) minval = val;

    if (val > maxval) maxval = val;
  }

  /* reject baserows that have a too large dynamism */
  if (maxval <= 1e6 * minval) {
    rhs_.push_back(double(rhs));
    for (int j : vectorsum.getNonzeros()) {
      double val = vectorsum.getValue(j);

      if (val == 0.0) continue;

      ARindex_.push_back(j);
      ARvalue_.push_back(val);
    }

    ARstart_.push_back(ARindex_.size());
  }

  vectorsum.clear();
}

void HighsSeparation::BaseRows::retransformAndAddCut(
    const HighsDomain& domain, const HighsLpRelaxation& lp,
    std::vector<int>& inds, std::vector<double>& vals,
    std::vector<int8_t>& complementation, HighsCDouble rhs,
    HighsCutPool& cutpool, HighsDomain& propdomain, bool cutintegral) const {
  const HighsMipSolver& mip = lp.getMipSolver();
  vectorsum.setDimension(mip.numCol() + lp.getNumLpRows());

  assert(std::all_of(vectorsum.nonzeroflag.begin(), vectorsum.nonzeroflag.end(),
                     [](uint8_t f) { return f == 0; }));
  int cutlen = inds.size();

  for (int i = 0; i != cutlen; ++i) {
    if (vals[i] == 0.0) continue;

    if (inds[i] >= mip.numCol()) {
      int row = inds[i] - mip.numCol();

      const int* Rindex;
      const double* Rvalue;
      int Rlen;
      double Rside;

      assert(slacktype_[row] == 1 || slacktype_[row] == -1);

      if (row < mip.numRow()) {
        if (slacktype_[row] == 1)
          Rside = lp.rowUpper(row);
        else
          Rside = -lp.rowLower(row);
        mip.mipdata_->getRow(row, Rlen, Rindex, Rvalue);
      } else {
        int cut = lp.getCutIndex(row);
        assert(slacktype_[row] == 1);
        cutpool.getCut(cut, Rlen, Rindex, Rvalue);
        Rside = lp.rowUpper(row);
      }

      rhs -= vals[i] * Rside;

      double slackval = -slacktype_[row] * vals[i];

      for (int j = 0; j != Rlen; ++j)
        vectorsum.add(Rindex[j], slackval * Rvalue[j]);

      continue;
    }

    if (complementation[i] == 1)
      rhs += vals[i] * domain.colLower_[inds[i]];
    else if (complementation[i] == -1) {
      vals[i] = -vals[i];
      rhs += vals[i] * domain.colUpper_[inds[i]];
    }

    vectorsum.add(inds[i], vals[i]);
  }

  inds.clear();
  vals.clear();

  vectorsum.sort();

  for (int col : vectorsum.getNonzeros()) {
    double val = vectorsum.getValue(col);

    if (std::abs(val) > 1e-10) {
      vals.push_back(val);
      inds.push_back(col);
    }
  }

  vectorsum.clear();

  double upper = double(rhs);

  cutpool.addCut(mip, inds.data(), vals.data(), inds.size(), upper,
                 cutintegral);
}

int HighsSeparation::separationRound(HighsDomain& propdomain,
                                     HighsLpRelaxation::Status& status) {
  const HighsSolution& sol = lp->getLpSolver().getSolution();

  HighsMipSolverData& mipdata = *lp->getMipSolver().mipdata_;
  ImpliedBounds implbound(mipdata.implications);
  implbound.separateImplBounds(*lp, mipdata.cutpool, propdomain);
  propdomain.propagate();

  if (propdomain.infeasible() || mipdata.domain.infeasible()) {
    status = HighsLpRelaxation::Status::Infeasible;
    propdomain.clearChangedCols();
    return 0;
  }

  if (!propdomain.getChangedCols().empty()) {
    lp->flushDomain(propdomain);
    status = lp->resolveLp();

    if (!lp->scaledOptimal(status)) return 0;
  }

  mipdata.cliquetable.separateCliques(lp->getMipSolver(), sol.col_value,
                                      mipdata.cutpool, mipdata.feastol);

  AggregationHeuristic aggheur(*lp, mipdata.domain, mipdata.cutpool,
                               propdomain);
  aggheur.run();

  baserows.clear();
  tableauaggregator(*lp, mipdata.cutpool, baserows);
  doSeparate(mipdata.domain, *lp, sol, baserows, mipdata.cutpool, propdomain);
  baserows.clear();

  propdomain.propagate();

  if (propdomain.infeasible()) {
    status = HighsLpRelaxation::Status::Infeasible;
    propdomain.clearChangedCols();
    return 0;
  }

  if (!propdomain.getChangedCols().empty()) {
    lp->flushDomain(propdomain);
    status = lp->resolveLp();
  }

  if (lp->scaledOptimal(status)) {
    mipdata.cutpool.separate(sol.col_value, propdomain, cutset,
                             mipdata.feastol);

    int ncuts = cutset.numCuts();

    if (ncuts > 0) {
      lp->addCuts(cutset);
      status = lp->resolveLp();
      mipdata.cutpool.ageLPRows(*lp);
    }
    return ncuts;
  }

  return 0;
}

void HighsSeparation::computeAndAddConflictCut(HighsMipSolver& mipsolver,
                                               HighsDomain& localdomain,
                                               std::vector<int>& inds,
                                               std::vector<double>& vals,
                                               double rowupper) {
  int len = inds.size();
  std::vector<double> solvals(len);
  std::vector<int8_t> complementation(len);
  std::vector<double> upper(len);
  HighsCDouble rhs = rowupper;

  mipsolver.mipdata_->debugSolution.checkCut(inds.data(), vals.data(),
                                             inds.size(), rowupper);

  int nbin = 0;
  int nint = 0;
  int ncont = 0;
  int nunbndint = 0;

  const HighsDomain& globaldomain = mipsolver.mipdata_->domain;

  double minact = 0.0;

  for (int i = 0; i != len; ++i) {
    int col = inds[i];
    assert(globaldomain.colUpper_[col] != HIGHS_CONST_INF ||
           globaldomain.colLower_[col] != -HIGHS_CONST_INF);

    if (globaldomain.colUpper_[col] == HIGHS_CONST_INF ||
        globaldomain.colLower_[col] == -HIGHS_CONST_INF) {
      upper[i] = HIGHS_CONST_INF;
    } else {
      upper[i] = globaldomain.colUpper_[col] - globaldomain.colLower_[col];
    }

    if (mipsolver.variableType(col) != HighsVarType::CONTINUOUS) {
      if (upper[i] < 1.5) {
        upper[i] = 1.0;
        ++nbin;
      } else if (upper[i] != HIGHS_CONST_INF) {
        upper[i] = std::floor(upper[i] + 0.5);
        ++nint;
      } else {
        ++nunbndint;
      }
    } else {
      ++ncont;
    }

    if (vals[i] < 0 && globaldomain.colUpper_[col] != HIGHS_CONST_INF) {
      vals[i] = complementWithUpper(vals[i], globaldomain.colUpper_[col], rhs);
      complementation[i] = -1;
      solvals[i] = globaldomain.colUpper_[col] - localdomain.colUpper_[col];
    } else {
      vals[i] = complementWithLower(vals[i], globaldomain.colLower_[col], rhs);
      complementation[i] = 1;
      solvals[i] = localdomain.colLower_[col] - globaldomain.colLower_[col];
    }

    minact += solvals[i] * vals[i];
  }

  bool cutintegral;
  bool success =
      generateCut(mipsolver, upper, nbin, nint, ncont, nunbndint, solvals,
                  complementation, inds, vals, rhs, cutintegral);

  if (success) {
    int offset = 0;

    for (int i = 0; i != len; ++i) {
      // skip zeros
      if (vals[i] == 0) continue;

      // undo complementation
      if (complementation[i] == 1)
        rhs += vals[i] * globaldomain.colLower_[inds[i]];
      else {
        assert(complementation[i] == -1);
        vals[i] = -vals[i];
        rhs += vals[i] * globaldomain.colUpper_[inds[i]];
      }

      // store back
      if (offset < i) {
        vals[offset] = vals[i];
        inds[offset] = inds[i];
      }
      ++offset;
    }

    mipsolver.mipdata_->cutpool.addCut(mipsolver, inds.data(), vals.data(),
                                       offset, double(rhs), cutintegral);
  }
}

void HighsSeparation::separate(HighsDomain& propdomain) {
  HighsLpRelaxation::Status status = lp->getStatus();
  const HighsMipSolver& mipsolver = lp->getMipSolver();

  if (lp->scaledOptimal(status) && !lp->getFractionalIntegers().empty()) {
    double firstobj = mipsolver.mipdata_->rootlpsolobj;

    while (lp->getObjective() < mipsolver.mipdata_->upper_limit) {
      double lastobj = lp->getObjective();

      size_t nlpiters = -lp->getNumLpIterations();
      int ncuts = separationRound(propdomain, status);
      nlpiters += lp->getNumLpIterations();
      mipsolver.mipdata_->sepa_lp_iterations += nlpiters;
      mipsolver.mipdata_->total_lp_iterations += nlpiters;
      // printf("separated %d cuts\n", ncuts);

      // printf(
      //     "separation round %d at node %d added %d cuts objective changed "
      //     "from %g to %g, first obj is %g\n",
      //     nrounds, (int)nnodes, ncuts, lastobj, lp->getObjective(),
      //     firstobj);
      if (ncuts == 0 || !lp->scaledOptimal(status) ||
          lp->getFractionalIntegers().empty())
        break;

      // if the objective improved considerably we continue
      if ((lp->getObjective() - firstobj) <=
          std::max((lastobj - firstobj), 0.0) * 1.01)
        break;
    }

    // printf("done separating\n");
  } else {
    // printf("no separation, just aging. status: %d\n", (int)status);
    mipsolver.mipdata_->cutpool.ageLPRows(*lp);
    mipsolver.mipdata_->cutpool.ageNonLPRows();
  }
}