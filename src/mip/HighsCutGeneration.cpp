#include "mip/HighsCutGeneration.h"

#include "mip/HighsMipSolverData.h"
#include "mip/HighsTransformedLp.h"
#include "util/HighsIntegers.h"

HighsCutGeneration::HighsCutGeneration(const HighsLpRelaxation& lpRelaxation,
                                       HighsCutPool& cutpool)
    : lpRelaxation(lpRelaxation),
      cutpool(cutpool),
      feastol(lpRelaxation.getMipSolver().mipdata_->feastol),
      epsilon(lpRelaxation.getMipSolver().mipdata_->epsilon) {}

bool HighsCutGeneration::determineCover() {
  if (rhs <= 10 * feastol) return false;

  cover.clear();
  cover.reserve(rowlen);

  for (int j = 0; j != rowlen; ++j) {
    if (!lpRelaxation.isColIntegral(inds[j])) continue;

    if (solval[j] <= feastol) continue;

    cover.push_back(j);
  }

  // take all variables that sit at their upper bound always into the cover
  int coversize =
      std::partition(cover.begin(), cover.end(),
                     [&](int j) { return solval[j] >= upper[j] - feastol; }) -
      cover.begin();
  int maxCoverSize = cover.size();

  coverweight = 0.0;
  for (int i = 0; i != coversize; ++i) {
    int j = cover[i];

    assert(solval[j] >= upper[j] - feastol);

    coverweight += vals[j] * upper[j];
  }

  // sort the remaining variables by the contribution to the rows activity in
  // the current solution
  std::sort(&cover[coversize], &cover[maxCoverSize], [&](int i, int j) {
    double contributionA = solval[i] * vals[i];
    double contributionB = solval[j] * vals[j];

    // for equal contributions take the larger coefficients first because
    // this makes some of the lifting functions more likely to generate
    // a facet
    if (std::abs(contributionA - contributionB) <= feastol)
      return vals[i] > vals[j];

    return contributionA > contributionB;
  });

  const double minlambda = 10 * feastol;

  for (; coversize != maxCoverSize; ++coversize) {
    double lambda = double(coverweight - rhs);
    if (lambda > minlambda) break;

    int j = cover[coversize];
    coverweight += vals[j] * upper[j];
  }
  if (coversize == 0) return false;

  coverweight.renormalize();
  lambda = coverweight - rhs;

  if (lambda <= 10 * feastol) return false;

  cover.resize(coversize);
  assert(lambda > feastol);

  return true;
}

void HighsCutGeneration::separateLiftedKnapsackCover() {
  const double feastol = lpRelaxation.getMipSolver().mipdata_->feastol;

  const int coversize = cover.size();

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

  if (double(sigma) > 0) abartmp = HighsCDouble(rhs) / double(coversize);

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
    if (vals[i] == 0.0) continue;
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

  // resulting cut is always integral
  integralSupport = true;
  integralCoefficients = true;
}

bool HighsCutGeneration::separateLiftedMixedBinaryCover() {
  int coversize = cover.size();
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
    if (vals[cover[i]] - lambda <= epsilon) {
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

  integralCoefficients = false;
  integralSupport = true;
  for (int i = 0; i != rowlen; ++i) {
    if (!lpRelaxation.isColIntegral(inds[i])) {
      if (vals[i] < 0)
        integralSupport = false;
      else
        vals[i] = 0;
      continue;
    }

    if (coverflag[i]) {
      vals[i] = std::min(vals[i], double(lambda));
      rhs += vals[i];
    } else {
      vals[i] = phi(vals[i]);
    }
  }

  return true;
}

bool HighsCutGeneration::separateLiftedMixedIntegerCover() {
  int coversize = cover.size();

  int l = -1;

  std::vector<uint8_t> coverflag;
  coverflag.resize(rowlen);
  for (int i : cover) coverflag[i] = 1;

  auto comp = [&](int a, int b) { return vals[a] > vals[b]; };
  std::sort(cover.begin(), cover.end(), comp);

  std::vector<HighsCDouble> a;
  std::vector<HighsCDouble> u;
  std::vector<HighsCDouble> m;

  a.resize(coversize);
  u.resize(coversize + 1);
  m.resize(coversize + 1);

  HighsCDouble usum = 0.0;
  HighsCDouble msum = 0.0;
  // set up the partial sums of the upper bounds, and the contributions
  for (int c = 0; c != coversize; ++c) {
    int i = cover[c];

    u[c] = usum;
    m[c] = msum;
    a[c] = vals[i];
    double ub = upper[i];
    usum += ub;
    msum += ub * a[c];
  }

  u[coversize] = usum;
  m[coversize] = msum;

  // determine which variable in the cover we want to create the MIR inequality
  // from which we lift we try to select a variable to have the highest chance
  // of satisfying the facet conditions for the superadditive lifting function
  // gamma to be satisfied.
  int lpos = -1;
  int bestlCplusend = -1;
  double bestlVal = 0.0;

  for (int i = 0; i != coversize; ++i) {
    int j = cover[i];
    double ub = upper[j];

    if (solval[j] >= ub - feastol) continue;

    double mju = ub * vals[j];
    HighsCDouble mu = mju - lambda;

    if (mu <= 10 * feastol) continue;

    double mudival = double(mu / vals[j]);
    if (std::abs(std::round(mudival) - mudival) <= feastol) continue;
    double eta = ceil(mudival);

    HighsCDouble ulminusetaplusone = HighsCDouble(ub) - eta + 1.0;
    HighsCDouble cplusthreshold = ulminusetaplusone * vals[j];

    int cplusend =
        std::upper_bound(cover.begin(), cover.end(), double(cplusthreshold),
                         [&](double cplusthreshold, int i) {
                           return cplusthreshold > vals[i];
                         }) -
        cover.begin();

    HighsCDouble mcplus = m[cplusend];
    if (i < cplusend) mcplus -= mju;

    double jlVal = double(mcplus + eta * vals[j]);

    if (jlVal > bestlVal) {
      lpos = i;
      bestlCplusend = cplusend;
      bestlVal = jlVal;
    }
  }

  if (lpos == -1) return false;

  l = cover[lpos];
  HighsCDouble al = vals[l];
  double upperl = upper[l];
  HighsCDouble mlu = upperl * al;
  HighsCDouble mu = mlu - lambda;

  a.resize(bestlCplusend);
  cover.resize(bestlCplusend);
  u.resize(bestlCplusend + 1);
  m.resize(bestlCplusend + 1);

  if (lpos < bestlCplusend) {
    a.erase(a.begin() + lpos);
    cover.erase(cover.begin() + lpos);
    u.erase(u.begin() + lpos + 1);
    m.erase(m.begin() + lpos + 1);
    for (int i = lpos + 1; i < bestlCplusend; ++i) {
      u[i] -= upperl;
      m[i] -= mlu;
    }
  }

  int cplussize = a.size();

  assert(mu > 10 * feastol);

  double mudival = double(mu / al);
  double eta = ceil(mudival);
  HighsCDouble r = mu - floor(mudival) * HighsCDouble(al);
  // we multiply with r and it is important that it does not flip the sign
  // so we safe guard against tiny numerical errors here
  if (r < 0) r = 0;

  HighsCDouble ulminusetaplusone = HighsCDouble(upperl) - eta + 1.0;
  HighsCDouble cplusthreshold = ulminusetaplusone * al;

  int kmin = floor(eta - upperl - 0.5);

  auto phi_l = [&](double a) {
    assert(a < 0);

    int64_t k = std::min(int64_t(a / double(al)), int64_t(-1));

    for (; k >= kmin; --k) {
      if (a >= k * al + r) {
        assert(a < (k + 1) * al);
        return double(a - (k + 1) * r);
      }

      if (a >= k * al) {
        assert(a < k * al + r);
        return double(k * (al - r));
      }
    }

    assert(a <= -lambda + epsilon);
    return double(kmin * (al - r));
  };

  int64_t kmax = floor(upperl - eta + 0.5);

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

        int64_t k = ((int64_t)(double)((z - mihplusdeltai) / al)) - 1;
        for (; k <= kmax; ++k) {
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

    int64_t p = ((int64_t)(double)((z - m[cplussize]) / al)) - 1;
    for (;; ++p) {
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

  rhs = (HighsCDouble(upperl) - eta) * r - lambda;
  integralSupport = true;
  integralCoefficients = false;
  for (int i = 0; i != rowlen; ++i) {
    if (vals[i] == 0.0) continue;
    int col = inds[i];

    if (!lpRelaxation.isColIntegral(col)) {
      if (vals[i] < 0.0)
        integralSupport = false;
      else
        vals[i] = 0.0;
      continue;
    }

    if (coverflag[i]) {
      vals[i] = -phi_l(-vals[i]);
      rhs += vals[i] * upper[i];
    } else {
      vals[i] = gamma_l(vals[i]);
    }
  }

  return true;
}

bool HighsCutGeneration::cmirCutGenerationHeuristic() {
  std::vector<double> deltas;

  HighsCDouble continuouscontribution = 0.0;
  HighsCDouble continuoussqrnorm = 0.0;
  std::vector<int> integerinds;
  integerinds.reserve(rowlen);
  double maxabsdelta = 0.0;

  std::vector<uint8_t> complementation(rowlen);

  for (int i = 0; i != rowlen; ++i) {
    if (lpRelaxation.isColIntegral(inds[i])) {
      integerinds.push_back(i);

      if (upper[i] < 2 * solval[i]) {
        complementation[i] = 1;
        rhs -= upper[i] * vals[i];
        vals[i] = -vals[i];
      }

      if (solval[i] > feastol) {
        double delta = std::abs(vals[i]);
        if (delta <= 1e-4 || delta >= 1e4) continue;
        maxabsdelta = std::max(maxabsdelta, delta);
        deltas.push_back(delta);
      }
    } else {
      continuouscontribution += vals[i] * solval[i];
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
    if (deltas[i] - curdelta <= feastol)
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

      viol += aj * solval[j];
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

      viol += aj * solval[j];
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

    complementation[k] = 1 - complementation[k];
    solval[k] = upper[k] - solval[k];
    rhs -= upper[k] * vals[k];
    vals[k] = -vals[k];

    double delta = bestdelta;
    HighsCDouble scale = 1.0 / HighsCDouble(delta);
    HighsCDouble scalrhs = rhs * scale;
    double downrhs = std::floor(double(scalrhs));

    HighsCDouble f0 = scalrhs - downrhs;
    if (f0 < 0.01 || f0 > 0.99) {
      complementation[k] = 1 - complementation[k];
      solval[k] = upper[k] - solval[k];
      rhs -= upper[k] * vals[k];
      vals[k] = -vals[k];

      continue;
    }

    HighsCDouble oneoveroneminusf0 = 1.0 / (1.0 - f0);
    if (double(oneoveroneminusf0) * double(scale) > 1e4) {
      complementation[k] = 1 - complementation[k];
      solval[k] = upper[k] - solval[k];
      rhs -= upper[k] * vals[k];
      vals[k] = -vals[k];

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

      viol += aj * solval[j];
      sqrnorm += aj * aj;
    }

    double efficacy = double(viol / sqrt(sqrnorm));
    if (efficacy > bestefficacy) {
      bestefficacy = efficacy;
    } else {
      complementation[k] = 1 - complementation[k];
      solval[k] = upper[k] - solval[k];
      rhs -= upper[k] * vals[k];
      vals[k] = -vals[k];
    }
  }

  HighsCDouble scale = 1.0 / HighsCDouble(bestdelta);
  HighsCDouble scalrhs = rhs * scale;
  double downrhs = std::floor(double(scalrhs));

  HighsCDouble f0 = scalrhs - downrhs;
  HighsCDouble oneoveroneminusf0 = 1.0 / (1.0 - f0);

  rhs = downrhs * bestdelta;
  integralSupport = true;
  integralCoefficients = false;
  for (int j = 0; j != rowlen; ++j) {
    if (vals[j] == 0.0) continue;
    if (!lpRelaxation.isColIntegral(inds[j])) {
      if (vals[j] > 0.0)
        vals[j] = 0.0;
      else {
        vals[j] = double(vals[j] * oneoveroneminusf0);
        integralSupport = false;
      }
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

  for (int j = 0; j != rowlen; ++j) {
    if (complementation[j]) {
      rhs -= upper[j] * vals[j];
      vals[j] = -vals[j];
    }
  }

  return true;
}

bool HighsCutGeneration::postprocessCut() {
  double maxAbsValue;
  if (integralSupport) {
    if (integralCoefficients) return true;

    // if the support is integral, allow a maximal dynamism of 1e4
    maxAbsValue = 0.0;
    for (int i = 0; i != rowlen; ++i)
      maxAbsValue = std::max(std::abs(vals[i]), maxAbsValue);

    double minCoefficientValue = std::max(maxAbsValue * 100 * feastol, epsilon);

    for (int i = 0; i != rowlen; ++i) {
      if (vals[i] == 0) continue;
      if (std::abs(vals[i]) <= minCoefficientValue) {
        if (vals[i] < 0) {
          double ub = upper[i];
          if (ub == HIGHS_CONST_INF)
            return false;
          else
            rhs -= ub * vals[i];
        }

        vals[i] = 0.0;
      }
    }

    std::vector<double> nonzerovals;
    nonzerovals.reserve(rowlen);

    for (int i = 0; i != rowlen; ++i)
      if (vals[i] != 0) nonzerovals.push_back(vals[i]);

    double intscale =
        HighsIntegers::integralScale(nonzerovals, feastol, epsilon);

    bool scaleSmallestValToOne = true;

    if (intscale != 0.0 &&
        intscale * std::max(1.0, maxAbsValue) <= (double)(uint64_t{1} << 53)) {
      // A scale to make all value integral was found. The scale is only
      // rejected if it is in a range where not all integral values are
      // representable in double precision anymore. Otherwise we want to always
      // use the scale to adjust the coefficients and right hand side for
      // numerical safety reasons. If the resulting integral values are too
      // large, however, we scale the cut down by shifting the exponent.
      rhs.renormalize();
      rhs *= intscale;
      maxAbsValue = std::round(maxAbsValue * intscale);
      for (int i = 0; i != rowlen; ++i) {
        if (vals[i] == 0.0) continue;

        HighsCDouble scaleval = intscale * HighsCDouble(vals[i]);
        HighsCDouble intval = round(scaleval);
        double delta = double(scaleval - intval);

        vals[i] = (double)intval;

        // if the coefficient would be strengthened by rounding, we add the
        // upperbound constraint to make it exactly integral instead and
        // therefore weaken the right hand side
        if (delta < 0.0) {
          if (upper[i] == HIGHS_CONST_INF) return false;

          rhs -= delta * upper[i];
        }
      }

      // finally we can round down the right hand side. Therefore in most cases
      // small errors for which the upper bound constraints where used and the
      // right hand side was weakened, do not weaken the final cut.
      rhs = floor(rhs + epsilon);

      if (intscale * maxAbsValue * feastol <= 1.0) {
        scaleSmallestValToOne = false;
        integralCoefficients = true;
      }
    }

    if (scaleSmallestValToOne) {
      double minAbsValue = HIGHS_CONST_INF;
      for (int i = 0; i != rowlen; ++i) {
        if (vals[i] == 0.0) continue;
        minAbsValue = std::min(std::abs(vals[i]), minAbsValue);
      }

      int expshift;
      std::frexp(minAbsValue - epsilon, &expshift);
      expshift = -expshift;

      maxAbsValue = std::ldexp(maxAbsValue, expshift);

      rhs = std::ldexp((double)rhs, expshift);

      for (int i = 0; i != rowlen; ++i) {
        if (vals[i] == 0) continue;

        vals[i] = std::ldexp(vals[i], expshift);
      }
    }
  } else {
    maxAbsValue = 0.0;
    for (int i = 0; i != rowlen; ++i)
      maxAbsValue = std::max(std::abs(vals[i]), maxAbsValue);

    double minCoefficientValue = maxAbsValue * 100 * feastol;

    // now remove small coefficients and determine the smallest absolute
    // coefficient of an integral variable
    double minIntCoef = HIGHS_CONST_INF;
    for (int i = 0; i != rowlen; ++i) {
      if (vals[i] == 0.0) continue;

      double val = std::abs(vals[i]);

      if (val <= minCoefficientValue) {
        if (vals[i] < 0.0) {
          if (upper[i] == HIGHS_CONST_INF) return false;
          rhs -= vals[i] * upper[i];
        } else
          vals[i] = 0.0;
      } else if (lpRelaxation.isColIntegral(inds[i]))
        minIntCoef = std::min(minIntCoef, val);
    }

    // if the smallest coefficient of an integral variable is not integral
    // we want to scale the cut so that it becomes 1. To avoid rounding errors
    // we instead shift the exponent.
    if (std::abs(std::round(minIntCoef) - minIntCoef) > epsilon) {
      int expshift;
      std::frexp(minIntCoef, &expshift);
      expshift = -expshift;

      if (expshift != 0) {
        maxAbsValue = std::ldexp(maxAbsValue, expshift);
        rhs = std::ldexp((double)rhs, expshift);
        for (int i = 0; i != rowlen; ++i)
          vals[i] = std::ldexp(vals[i], expshift);
      }
    }
  }

  return true;
}

bool HighsCutGeneration::preprocessBaseInequality(bool& hasUnboundedInts,
                                                  bool& hasGeneralInts,
                                                  bool& hasContinuous) {
  // preprocess the inequality before cut generation
  // 1. Determine the maximal activity to check for trivial redundancy and
  // tighten coefficients
  // 2. Check for presence of continuous variables and unbounded integers as not
  // all methods for cut generation are applicable in that case
  // 3. Remove coefficients that are below the feasibility tolerance to avoid
  // numerical troubles, use bound constraints to cancel them and
  // reject base inequalities where that is not possible due to unbounded
  // variables
  hasUnboundedInts = false;
  hasContinuous = false;
  hasGeneralInts = false;
  int numZeros = 0;

  HighsCDouble maxact = 0.0;
  bool maxactinf = false;
  for (int i = 0; i != rowlen; ++i) {
    if (std::abs(vals[i]) <= feastol) {
      if (vals[i] < 0) {
        if (upper[i] == HIGHS_CONST_INF) return false;
        rhs -= vals[i] * upper[i];
      }

      ++numZeros;
      vals[i] = 0.0;
      continue;
    }

    if (!lpRelaxation.isColIntegral(inds[i])) {
      hasContinuous = true;

      if (vals[i] > 0) {
        if (upper[i] == HIGHS_CONST_INF)
          maxactinf = true;
        else
          maxact += vals[i] * upper[i];
      }
    } else {
      if (upper[i] == HIGHS_CONST_INF) {
        hasUnboundedInts = true;
        hasGeneralInts = true;
        if (vals[i] > 0.0) maxactinf = true;
        if (maxactinf) break;
      } else if (upper[i] != 1.0) {
        hasGeneralInts = true;
      }

      if (vals[i] > 0) maxact += vals[i] * upper[i];
    }
  }

  int maxLen = 100 + 0.15 * (lpRelaxation.numCols());

  if (rowlen - numZeros > maxLen) {
    int numCancel = rowlen - numZeros - maxLen;
    std::vector<int> cancelNzs;

    for (int i = 0; i != rowlen; ++i) {
      double cancelSlack = vals[i] > 0 ? solval[i] : upper[i] - solval[i];
      if (cancelSlack <= feastol) cancelNzs.push_back(i);
    }

    if ((int)cancelNzs.size() < numCancel) return false;
    if ((int)cancelNzs.size() > numCancel)
      std::partial_sort(
          cancelNzs.begin(), cancelNzs.begin() + numCancel, cancelNzs.end(),
          [&](int a, int b) { return std::abs(vals[a]) < std::abs(vals[b]); });

    for (int i = 0; i < numCancel; ++i) {
      int j = cancelNzs[i];

      if (vals[j] < 0) {
        rhs -= vals[j] * upper[j];
      } else
        maxact -= vals[j] * upper[j];

      vals[j] = 0.0;
    }

    numZeros += numCancel;
  }

  if (numZeros != 0) {
    // remove zeros in place
    for (int i = rowlen - 1; i >= 0; --i) {
      if (vals[i] == 0.0) {
        --rowlen;
        inds[i] = inds[rowlen];
        vals[i] = vals[rowlen];
        upper[i] = upper[rowlen];
        solval[i] = solval[rowlen];
        if (--numZeros == 0) break;
      }
    }
  }

  if (!maxactinf) {
    double maxabscoef = double(maxact - rhs);
    if (maxabscoef <= feastol) return false;

    int ntightened = 0;
    for (int i = 0; i != rowlen; ++i) {
      if (!lpRelaxation.isColIntegral(inds[i])) continue;

      if (vals[i] > maxabscoef) {
        HighsCDouble delta = vals[i] - maxabscoef;
        rhs -= delta * upper[i];
        vals[i] = maxabscoef;
        ++ntightened;
      } else if (vals[i] < -maxabscoef) {
        vals[i] = -maxabscoef;
      }
    }
  }

  return true;
}

bool HighsCutGeneration::generateCut(HighsTransformedLp& transLp,
                                     std::vector<int>& inds_,
                                     std::vector<double>& vals_, double& rhs_) {
  if (!transLp.transform(vals_, upper, solval, inds_, rhs_, true)) return false;

  rowlen = inds_.size();
  this->inds = inds_.data();
  this->vals = vals_.data();
  this->rhs = rhs_;

  bool hasUnboundedInts = false;
  bool hasGeneralInts = false;
  bool hasContinuous = false;
  if (!preprocessBaseInequality(hasUnboundedInts, hasGeneralInts,
                                hasContinuous))
    return false;

  if (hasUnboundedInts) {
    if (!cmirCutGenerationHeuristic()) return false;
  } else {
    // 1. Determine a cover, cover does not need to be minimal as neither of
    // the
    //    lifting functions have minimality of the cover as necessary facet
    //    condition
    if (!determineCover()) return false;

    // 2. use superadditive lifting function depending on structure of base
    //    inequality:
    //    We have 3 lifting functions available for pure binary knapsack sets,
    //    for mixed-binary knapsack sets and for mixed integer knapsack sets.
    if (!hasContinuous && !hasGeneralInts)
      separateLiftedKnapsackCover();
    else if (hasGeneralInts) {
      if (!separateLiftedMixedIntegerCover()) return false;
    } else {
      assert(hasContinuous);
      assert(!hasGeneralInts);
      if (!separateLiftedMixedBinaryCover()) return false;
    }
  }

  // apply cut postprocessing including scaling and removal of small
  // coeffiicents
  if (!postprocessCut()) return false;

  // transform the cut back into the original space, i.e. remove the bound
  // substitution and replace implicit slack variables
  rhs_ = (double)rhs;
  bool cutintegral = integralSupport && integralCoefficients;
  vals_.resize(rowlen);
  inds_.resize(rowlen);
  if (!transLp.untransform(vals_, inds_, rhs_, cutintegral)) return false;

  // finally check whether the cut is violated
  rowlen = inds_.size();
  inds = inds_.data();
  vals = vals_.data();
  lpRelaxation.getMipSolver().mipdata_->debugSolution.checkCut(inds, vals,
                                                               rowlen, rhs_);

  // finally determine the violation of the cut in the original space
  HighsCDouble violation = -rhs_;
  const auto& sol = lpRelaxation.getSolution().col_value;
  for (int i = 0; i != rowlen; ++i) violation += sol[inds[i]] * vals_[i];

  if (violation <= 10 * feastol) return false;

  lpRelaxation.getMipSolver().mipdata_->domain.tightenCoefficients(
      inds, vals, rowlen, rhs_);

  // if the cut is violated by a small factor above the feasibility
  // tolerance, add it to the cutpool
  int cutindex = cutpool.addCut(lpRelaxation.getMipSolver(), inds_.data(),
                                vals_.data(), inds_.size(), rhs_, cutintegral);

  // only return true if cut was accepted by the cutpool, i.e. not a duplicate
  // of a cut already in the pool
  return cutindex != -1;
}

bool HighsCutGeneration::generateConflict(HighsDomain& localdomain,
                                          std::vector<int>& proofinds,
                                          std::vector<double>& proofvals,
                                          double& proofrhs) {
  this->inds = proofinds.data();
  this->vals = proofvals.data();
  this->rhs = proofrhs;
  rowlen = proofinds.size();

  std::vector<uint8_t> complementation(rowlen);

  upper.resize(rowlen);
  solval.resize(rowlen);

  HighsDomain& globaldomain = lpRelaxation.getMipSolver().mipdata_->domain;
  HighsNodeQueue& nodequeue = lpRelaxation.getMipSolver().mipdata_->nodequeue;
  std::vector<std::pair<size_t, int>> nonzeroSolvals;
  nonzeroSolvals.reserve(rowlen);
  for (int i = 0; i != rowlen; ++i) {
    int col = inds[i];

    upper[i] = globaldomain.colUpper_[col] - globaldomain.colLower_[col];
    size_t numAffectedNodes;
    if (vals[i] < 0 && globaldomain.colUpper_[col] != HIGHS_CONST_INF) {
      rhs -= globaldomain.colUpper_[col] * vals[i];
      vals[i] = -vals[i];
      complementation[i] = 1;

      if (globaldomain.isBinary(col))
        numAffectedNodes = nodequeue.numNodesDown(col);
      else
        numAffectedNodes =
            nodequeue.numNodesDown(col, localdomain.colUpper_[col]);
      solval[i] = globaldomain.colUpper_[col] - localdomain.colUpper_[col];
    } else {
      rhs -= globaldomain.colLower_[col] * vals[i];
      complementation[i] = 0;
      solval[i] = localdomain.colLower_[col] - globaldomain.colLower_[col];

      if (globaldomain.isBinary(col))
        numAffectedNodes = nodequeue.numNodesUp(col);
      else
        numAffectedNodes =
            nodequeue.numNodesUp(col, localdomain.colLower_[col]);
    }

    if (solval[i] > feastol) nonzeroSolvals.emplace_back(numAffectedNodes, i);
  }

  std::sort(
      nonzeroSolvals.begin(), nonzeroSolvals.end(),
      [&](std::pair<size_t, int> const& a, std::pair<size_t, int> const& b) {
        if (a.first > b.first) return true;
        if (a.first < b.first) return false;
        return vals[a.second] * solval[a.second] >
               vals[b.second] * solval[b.second];
      });
  HighsCDouble capacity = rhs;
  int nzSolvalsLen = nonzeroSolvals.size();
  for (int i = 0; i != nzSolvalsLen; ++i) {
    int j = nonzeroSolvals[i].second;
    double contribution = solval[j] * vals[j];
    HighsCDouble reducedCapacity = capacity - contribution;
    if (reducedCapacity >= 0) {
      capacity = reducedCapacity;
    } else {
      solval[j] = double(capacity / vals[j]);
      capacity = 0;
    }
  }

  bool hasUnboundedInts = false;
  bool hasGeneralInts = false;
  bool hasContinuous = false;

  if (!preprocessBaseInequality(hasUnboundedInts, hasGeneralInts,
                                hasContinuous))
    return false;

  if (hasUnboundedInts) {
    if (!cmirCutGenerationHeuristic()) return false;
  } else {
    // 1. Determine a cover, cover does not need to be minimal as neither of the
    //    lifting functions have minimality of the cover as necessary facet
    //    condition
    if (!determineCover()) return false;

    // 2. use superadditive lifting function depending on structure of base
    //    inequality:
    //    We have 3 lifting functions available for pure binary knapsack sets,
    //    for mixed-binary knapsack sets and for mixed integer knapsack sets.
    if (!hasContinuous && !hasGeneralInts)
      separateLiftedKnapsackCover();
    else if (hasGeneralInts) {
      if (!separateLiftedMixedIntegerCover()) return false;
    } else {
      assert(hasContinuous);
      assert(!hasGeneralInts);
      if (!separateLiftedMixedBinaryCover()) return false;
    }
  }

  // apply cut postprocessing including scaling and removal of small
  // coefficients
  if (!postprocessCut()) return false;

  // remove the complementation
  for (int i = 0; i != rowlen; ++i) {
    if (complementation[i]) {
      rhs -= globaldomain.colUpper_[inds[i]] * vals[i];
      vals[i] = -vals[i];
    } else
      rhs += globaldomain.colLower_[inds[i]] * vals[i];
  }

  // remove zeros in place
  for (int i = rowlen - 1; i >= 0; --i) {
    if (vals[i] == 0.0) {
      --rowlen;
      proofinds[i] = proofinds[rowlen];
      proofvals[i] = proofvals[rowlen];
    }
  }

  proofvals.resize(rowlen);
  proofrhs = (double)rhs;

  bool cutintegral = integralSupport && integralCoefficients;

  // finally determine the violation of the local domain
  HighsCDouble violation = -proofrhs;

  for (int i = 0; i != rowlen; ++i) {
    if (vals[i] < 0)
      violation += localdomain.colUpper_[inds[i]] * proofvals[i];
    else
      violation += localdomain.colLower_[inds[i]] * proofvals[i];
  }

  // if the cut is violated above the feasibility tolerance we add it
  if (violation <= feastol) return false;

  lpRelaxation.getMipSolver().mipdata_->domain.tightenCoefficients(
      proofinds.data(), proofvals.data(), proofinds.size(), proofrhs);

  int cutindex =
      cutpool.addCut(lpRelaxation.getMipSolver(), proofinds.data(),
                     proofvals.data(), rowlen, proofrhs, cutintegral);

  // only return true if cut was accepted by the cutpool, i.e. not a duplicate
  // of a cut already in the pool
  return cutindex != -1;
}