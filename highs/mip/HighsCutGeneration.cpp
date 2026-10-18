/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsCutGeneration.h"

#include "../extern/pdqsort/pdqsort.h"
#include "mip/HighsDomain.h"
#include "mip/HighsMipSolverData.h"
#include "mip/HighsTransformedLp.h"
#include "util/HighsIntegers.h"

HighsCutGeneration::HighsCutGeneration(const HighsLpRelaxation& lpRelaxation,
                                       HighsCutPool& cutpool)
    : lpRelaxation(lpRelaxation),
      cutpool(cutpool),
      randgen(lpRelaxation.getMipSolver().options_mip_->random_seed +
              lpRelaxation.getNumLpIterations() + cutpool.getNumCuts()),
      feastol(lpRelaxation.getMipSolver().mipdata_->feastol),
      epsilon(lpRelaxation.getMipSolver().mipdata_->epsilon) {}

bool HighsCutGeneration::determineCover(bool lpSol) {
  if (cut.rhs <= 10 * feastol) return false;

  const HighsInt numVars = cut.vars.size();
  cover.clear();
  cover.reserve(numVars);

  for (HighsInt j = 0; j != numVars; ++j) {
    const TransformedVar& var = cut.vars[j];
    if (!var.integral) continue;

    if (lpSol && var.solval <= feastol) continue;

    cover.push_back(j);
  }

  HighsInt maxCoverSize = cover.size();
  HighsInt coversize = 0;
  HighsInt r = randgen.integer();
  coverweight = 0.0;
  if (lpSol) {
    // take all variables that sit at their upper bound always into the cover
    coversize = std::partition(cover.begin(), cover.end(),
                               [&](HighsInt j) {
                                 return cut.vars[j].flippedSolval <= feastol;
                               }) -
                cover.begin();

    for (HighsInt i = 0; i != coversize; ++i) {
      const TransformedVar& var = cut.vars[cover[i]];

      assert(var.flippedSolval <= feastol);

      coverweight += static_cast<HighsCDouble>(var.upper) * var.val;
    }

    // sort the remaining variables by the contribution to the rows activity in
    // the current solution
    pdqsort(
        cover.begin() + coversize, cover.begin() + maxCoverSize,
        [&](HighsInt i, HighsInt j) {
          if (cut.vars[i].upper < 1.5 && cut.vars[j].upper > 1.5) return true;
          if (cut.vars[i].upper > 1.5 && cut.vars[j].upper < 1.5) return false;

          double contributionA = cut.vars[i].solval * cut.vars[i].val;
          double contributionB = cut.vars[j].solval * cut.vars[j].val;

          if (contributionA > contributionB + feastol) return true;
          if (contributionA < contributionB - feastol) return false;
          // for equal contributions take the larger coefficients first
          // because this makes some of the lifting functions more likely
          // to generate a facet
          // if the value is equal too, choose a random tiebreaker based
          // on hashing the column index and the current number of pool
          // cuts
          if (std::abs(cut.vars[i].val - cut.vars[j].val) <= feastol)
            return HighsHashHelpers::hash(std::make_pair(cut.vars[i].col, r)) >
                   HighsHashHelpers::hash(std::make_pair(cut.vars[j].col, r));
          return cut.vars[i].val > cut.vars[j].val;
        });
  } else {
    // the current solution
    const auto& nodequeue = lpRelaxation.getMipSolver().mipdata_->nodequeue;

    pdqsort(
        cover.begin() + coversize, cover.begin() + maxCoverSize,
        [&](HighsInt i, HighsInt j) {
          if (cut.vars[i].solval > feastol && cut.vars[j].solval <= feastol)
            return true;
          if (cut.vars[i].solval <= feastol && cut.vars[j].solval > feastol)
            return false;

          int64_t numNodesA;
          int64_t numNodesB;

          numNodesA = cut.vars[i].complemented
                          ? nodequeue.numNodesDown(cut.vars[i].col)
                          : nodequeue.numNodesUp(cut.vars[i].col);

          numNodesB = cut.vars[j].complemented
                          ? nodequeue.numNodesDown(cut.vars[j].col)
                          : nodequeue.numNodesUp(cut.vars[j].col);

          if (numNodesA > numNodesB) return true;
          if (numNodesA < numNodesB) return false;

          return HighsHashHelpers::hash(std::make_pair(cut.vars[i].col, r)) >
                 HighsHashHelpers::hash(std::make_pair(cut.vars[j].col, r));
        });
  }

  const double minlambda =
      std::max(10 * feastol, feastol * std::abs(static_cast<double>(cut.rhs)));

  for (; coversize != maxCoverSize; ++coversize) {
    double lambda = static_cast<double>(coverweight - cut.rhs);
    if (lambda > minlambda) break;

    const TransformedVar& var = cut.vars[cover[coversize]];
    coverweight += static_cast<HighsCDouble>(var.upper) * var.val;
  }
  if (coversize == 0) return false;

  coverweight.renormalize();
  lambda = coverweight - cut.rhs;

  if (lambda <= minlambda) return false;

  cover.resize(coversize);
  assert(lambda > feastol);

  return true;
}

void HighsCutGeneration::separateLiftedKnapsackCover() {
  const HighsInt coversize = cover.size();
  const HighsInt numVars = cut.vars.size();

  std::vector<double> S;
  S.resize(coversize);
  std::vector<int8_t> coverflag;
  coverflag.resize(numVars);
  pdqsort_branchless(cover.begin(), cover.end(), [&](HighsInt a, HighsInt b) {
    return cut.vars[a].val > cut.vars[b].val;
  });

  HighsCDouble abartmp = cut.vars[cover[0]].val;
  HighsCDouble sigma = lambda;
  for (HighsInt i = 1; i != coversize; ++i) {
    HighsCDouble delta = abartmp - cut.vars[cover[i]].val;
    HighsCDouble kdelta = static_cast<double>(i) * delta;
    if (kdelta < sigma) {
      abartmp = cut.vars[cover[i]].val;
      sigma -= kdelta;
    } else {
      abartmp -= sigma * (1.0 / i);
      sigma = 0.0;
      break;
    }
  }

  if (sigma > 0)
    abartmp =
        static_cast<HighsCDouble>(cut.rhs) / static_cast<double>(coversize);

  double abar = static_cast<double>(abartmp);

  HighsCDouble sum = 0.0;
  HighsInt cplussize = 0;
  for (HighsInt i = 0; i != coversize; ++i) {
    sum += std::min(abar, cut.vars[cover[i]].val);
    S[i] = static_cast<double>(sum);

    if (cut.vars[cover[i]].val > abar + feastol) {
      ++cplussize;
      coverflag[cover[i]] = 1;
    } else
      coverflag[cover[i]] = -1;
  }
  assert(std::abs(static_cast<double>((sum - cut.rhs) / cut.rhs)) <= 1e-14);
  bool halfintegral = false;

  /* define the lifting function */
  auto g = [&](double z) {
    double hfrac = z / abar;
    double coef = 0.0;

    HighsInt h = std::floor(hfrac + 0.5);
    if (h != 0 && std::abs(hfrac - h) * std::max(1.0, abar) <= epsilon &&
        h <= cplussize - 1) {
      halfintegral = true;
      coef = 0.5;
    }

    h = std::max(h - 1, HighsInt{0});
    for (; h < coversize; ++h) {
      if (z <= S[h] + feastol) break;
    }

    return coef + h;
  };

  cut.rhs = coversize - 1;

  for (HighsInt i = 0; i != numVars; ++i) {
    TransformedVar& var = cut.vars[i];
    if (var.val == 0.0) continue;
    if (coverflag[i] == -1) {
      var.val = 1;
    } else {
      var.val = g(var.val);
    }
  }

  if (halfintegral) cut.scale(2);

  // resulting cut is always integral
  integralSupport = true;
  integralCoefficients = true;
}

bool HighsCutGeneration::separateLiftedMixedBinaryCover() {
  // initialize indicators
  integralSupport = false;
  integralCoefficients = false;

  HighsInt coversize = cover.size();
  const HighsInt numVars = cut.vars.size();
  std::vector<double> S;
  S.resize(coversize);
  std::vector<HighsBool> coverflag;
  coverflag.resize(numVars);

  if (coversize == 0) return false;

  for (HighsInt i : cover) coverflag[i] = true;

  pdqsort_branchless(cover.begin(), cover.end(), [&](HighsInt a, HighsInt b) {
    return cut.vars[a].val > cut.vars[b].val;
  });
  HighsCDouble sum = 0;

  HighsInt p = coversize;
  for (HighsInt i = 0; i != coversize; ++i) {
    if (cut.vars[cover[i]].val - lambda <= epsilon) {
      p = i;
      break;
    }
    sum += cut.vars[cover[i]].val;
    S[i] = static_cast<double>(sum);
  }
  if (p == 0) return false;
  /* define the lifting function */
  auto phi = [&](double a) {
    for (HighsInt i = 0; i < p; ++i) {
      if (a <= S[i] - lambda) return static_cast<double>(i * lambda);

      if (a <= S[i])
        return static_cast<double>((i + 1) * lambda +
                                   (static_cast<HighsCDouble>(a) - S[i]));
    }

    return static_cast<double>(p * lambda +
                               (static_cast<HighsCDouble>(a) - S[p - 1]));
  };

  cut.rhs = -lambda;

  integralCoefficients = false;
  integralSupport = true;
  for (HighsInt i = 0; i != numVars; ++i) {
    TransformedVar& var = cut.vars[i];
    if (!var.integral) {
      if (var.val < 0)
        integralSupport = false;
      else
        var.val = 0;
      continue;
    }

    if (coverflag[i]) {
      var.val = std::min(var.val, static_cast<double>(lambda));
      cut.rhs += var.val;
    } else {
      var.val = phi(var.val);
    }
  }

  return true;
}

bool HighsCutGeneration::separateLiftedMixedIntegerCover() {
  // initialize indicators
  integralSupport = false;
  integralCoefficients = false;

  HighsInt coversize = cover.size();
  const HighsInt numVars = cut.vars.size();

  HighsInt l = -1;

  std::vector<HighsBool> coverflag;
  coverflag.resize(numVars);
  for (HighsInt i : cover) coverflag[i] = true;

  auto comp = [&](HighsInt a, HighsInt b) {
    return cut.vars[a].val > cut.vars[b].val;
  };
  pdqsort_branchless(cover.begin(), cover.end(), comp);

  std::vector<HighsCDouble> a;
  std::vector<HighsCDouble> u;
  std::vector<HighsCDouble> m;

  a.resize(coversize);
  u.resize(coversize + 1);
  m.resize(coversize + 1);

  HighsCDouble usum = 0.0;
  HighsCDouble msum = 0.0;
  // set up the partial sums of the upper bounds, and the contributions
  for (HighsInt c = 0; c != coversize; ++c) {
    HighsInt i = cover[c];

    u[c] = usum;
    m[c] = msum;
    a[c] = cut.vars[i].val;
    double ub = cut.vars[i].upper;
    usum += ub;
    msum += ub * a[c];
  }

  u[coversize] = usum;
  m[coversize] = msum;

  // determine which variable in the cover we want to create the MIR inequality
  // from which we lift we try to select a variable to have the highest chance
  // of satisfying the facet conditions for the superadditive lifting function
  // gamma to be satisfied.
  HighsInt lpos = -1;
  HighsInt bestlCplusend = -1;
  double bestlVal = 0.0;
  bool bestlAtUpper = true;

  for (HighsInt i = 0; i != coversize; ++i) {
    HighsInt j = cover[i];
    double ub = cut.vars[j].upper;

    bool atUpper = cut.vars[j].flippedSolval <= feastol;
    if (atUpper && !bestlAtUpper) continue;

    double valj = cut.vars[j].val;
    double mju = ub * valj;
    HighsCDouble mu = mju - lambda;

    if (mu <= 10 * feastol) continue;
    if (std::abs(valj) < 1000 * feastol) continue;

    double mudival = static_cast<double>(mu / valj);
    if (HighsIntegers::isIntegral(mudival, feastol)) continue;
    double eta = ceil(mudival);

    HighsCDouble ulminusetaplusone = static_cast<HighsCDouble>(ub) - eta + 1.0;
    HighsCDouble cplusthreshold = ulminusetaplusone * valj;

    HighsInt cplusend =
        std::upper_bound(cover.begin(), cover.end(),
                         static_cast<double>(cplusthreshold),
                         [&](double cplusthreshold, HighsInt i) {
                           return cplusthreshold > cut.vars[i].val;
                         }) -
        cover.begin();

    HighsCDouble mcplus = m[cplusend];
    if (i < cplusend) mcplus -= mju;

    double jlVal = static_cast<double>(mcplus + eta * valj);

    if (jlVal > bestlVal || (!atUpper && bestlAtUpper)) {
      lpos = i;
      bestlCplusend = cplusend;
      bestlVal = jlVal;
      bestlAtUpper = atUpper;
    }
  }

  if (lpos == -1) return false;

  l = cover[lpos];
  HighsCDouble al = cut.vars[l].val;
  double upperl = cut.vars[l].upper;
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
    for (HighsInt i = lpos + 1; i < bestlCplusend; ++i) {
      u[i] -= upperl;
      m[i] -= mlu;
    }
  }

  HighsInt cplussize = a.size();

  assert(mu > 10 * feastol);

  double mudival = static_cast<double>(mu / al);
  double eta = ceil(mudival);
  HighsCDouble r = mu - floor(mudival) * static_cast<HighsCDouble>(al);
  // we multiply with r and it is important that it does not flip the sign
  // so we safe guard against tiny numerical errors here
  if (r < 0) r = 0;

  HighsCDouble ulminusetaplusone =
      static_cast<HighsCDouble>(upperl) - eta + 1.0;
  HighsCDouble cplusthreshold = ulminusetaplusone * al;

  HighsInt kmin = floor(eta - upperl - 0.5);

  auto phi_l = [&](double a) {
    assert(a < 0);

    int64_t k = std::min(static_cast<int64_t>(a / static_cast<double>(al)),
                         int64_t{-1});

    for (; k >= kmin; --k) {
      if (a >= k * al + r) {
        assert(a < (k + 1) * al);
        return static_cast<double>(a - (k + 1) * r);
      }

      if (a >= k * al) {
        assert(a < k * al + r);
        return static_cast<double>(k * (al - r));
      }
    }

    assert(a <= -lambda + epsilon);
    return static_cast<double>(kmin * (al - r));
  };

  int64_t kmax = floor(upperl - eta + 0.5);

  auto gamma_l = [&](double z) {
    assert(z > 0);
    for (HighsInt i = 0; i < cplussize; ++i) {
      HighsInt upperi = static_cast<HighsInt>(cut.vars[cover[i]].upper);

      for (HighsInt h = 0; h <= upperi; ++h) {
        HighsCDouble mih = m[i] + h * a[i];
        HighsCDouble uih = u[i] + h;
        HighsCDouble mihplusdeltai = mih + a[i] - cplusthreshold;
        if (z <= mihplusdeltai) {
          assert(mih <= z + epsilon);
          return static_cast<double>(uih * ulminusetaplusone * (al - r));
        }

        int64_t k = static_cast<int64_t>(
                        static_cast<double>((z - mihplusdeltai) / al)) -
                    1;
        for (; k <= kmax; ++k) {
          if (z <= mihplusdeltai + k * al + r) {
            assert(mihplusdeltai + k * al <= z);
            return static_cast<double>((uih * ulminusetaplusone + k) *
                                       (al - r));
          }

          if (z <= mihplusdeltai + (k + 1) * al) {
            assert(mihplusdeltai + k * al + r <= z);
            return static_cast<double>((uih * ulminusetaplusone) * (al - r) +
                                       z - mih - a[i] + cplusthreshold -
                                       (k + 1) * r);
          }
        }
      }
    }

    int64_t p =
        static_cast<int64_t>(static_cast<double>((z - m[cplussize]) / al)) - 1;
    for (;; ++p) {
      if (z <= m[cplussize] + p * al + r) {
        assert(m[cplussize] + p * al <= z);
        return static_cast<double>((u[cplussize] * ulminusetaplusone + p) *
                                   (al - r));
      }

      if (z <= m[cplussize] + (p + 1) * al) {
        assert(m[cplussize] + p * al + r <= z);
        return static_cast<double>((u[cplussize] * ulminusetaplusone) *
                                       (al - r) +
                                   z - m[cplussize] - (p + 1) * r);
      }
    }
  };

  cut.rhs = (static_cast<HighsCDouble>(upperl) - eta) * r - lambda;
  integralSupport = true;
  integralCoefficients = false;
  for (HighsInt i = 0; i != numVars; ++i) {
    TransformedVar& var = cut.vars[i];
    if (var.val == 0.0) continue;
    if (!var.integral) {
      if (var.val < 0.0)
        integralSupport = false;
      else
        var.val = 0.0;
      continue;
    }

    if (coverflag[i]) {
      var.val = -phi_l(-var.val);
      cut.rhs += static_cast<HighsCDouble>(var.upper) * var.val;
    } else {
      var.val = gamma_l(var.val);
    }
  }

  return true;
}

static double fast_floor(double x) {
  return static_cast<int64_t>(x) - (x < static_cast<int64_t>(x));
}

bool HighsCutGeneration::cmirCutGenerationHeuristic(double minEfficacy,
                                                    bool onlyInitialCMIRScale) {
  using std::abs;
  using std::floor;
  using std::max;
  using std::sqrt;

  // initialize indicators
  integralSupport = false;
  integralCoefficients = false;

  HighsCDouble continuouscontribution = 0.0;
  HighsCDouble continuoussqrnorm = 0.0;

  const HighsInt numVars = cut.vars.size();
  deltas.clear();
  deltas.reserve(numVars + 3);
  integerinds.clear();
  integerinds.reserve(numVars);

  double maxabsdelta = 0.0;
  constexpr double maxCMirScale = 1e6;
  constexpr double f0min = 0.005;
  constexpr double f0max = 0.995;

  for (HighsInt i = 0; i != numVars; ++i) {
    TransformedVar& var = cut.vars[i];
    if (var.integral) {
      integerinds.push_back(i);

      if (var.upper < 2 * var.solval) cut.flipComplementation(var);

      if (onlyInitialCMIRScale) continue;

      if (var.solval > feastol) {
        double delta = abs(var.val);
        if (delta <= 1e-4 || delta == maxabsdelta) continue;
        maxabsdelta = max(maxabsdelta, delta);
        deltas.push_back(delta);
      }
    } else {
      updateViolationAndNorm(var, var.val, continuouscontribution,
                             continuoussqrnorm);
    }
  }

  if (continuoussqrnorm == 0 && deltas.size() > 1) {
    double intScale = HighsIntegers::integralScale(deltas, feastol, kHighsTiny);

    if (intScale != 0.0 && intScale <= 1e4) {
      double scalrhs = static_cast<double>(cut.rhs) * intScale;
      double downrhs = fast_floor(scalrhs);

      double f0 = scalrhs - downrhs;
      if (f0 >= f0min && f0 <= f0max) deltas.push_back(1.0 / intScale);
    }
  }

  deltas.push_back(std::min(1.0, initialScale));
  if (!onlyInitialCMIRScale)
    deltas.push_back(maxabsdelta + std::min(1.0, initialScale));

  if (deltas.empty()) return false;

  pdqsort(deltas.begin(), deltas.end());
  double curdelta = deltas[0];
  for (size_t i = 1; i < deltas.size(); ++i) {
    if (deltas[i] - curdelta <= 10 * feastol)
      deltas[i] = 0.0;
    else
      curdelta = deltas[i];
  }

  deltas.erase(std::remove(deltas.begin(), deltas.end(), 0.0), deltas.end());

  // initialize best delta and efficacy
  double bestdelta = -1;
  double bestefficacy = minEfficacy;

  // lambda for computing the coefficient of an integer variable
  auto computeIntegralCoef = [&](double val, const HighsCDouble& f0,
                                 const HighsCDouble& oneoveroneminusf0,
                                 const HighsCDouble& scale) {
    HighsCDouble scalaj = val * scale;
    HighsCDouble downaj = floor(scalaj + kHighsTiny);
    HighsCDouble fj = scalaj - downaj;
    return downaj + max(HighsCDouble{0.0}, fj - f0) * oneoveroneminusf0;
  };

  // lambda for computing rhs and fractional part
  auto computeRhs = [&](const HighsCDouble& scale, HighsCDouble& downrhs,
                        HighsCDouble& f0, HighsCDouble& oneoveroneminusf0) {
    HighsCDouble scalrhs = cut.rhs * scale;
    downrhs = floor(scalrhs);
    f0 = scalrhs - downrhs;
    oneoveroneminusf0 = 1.0 / (1.0 - f0);
  };

  // lambda for trying out a particular delta and/or complementation
  auto isCutBetter = [&](double delta, double& efficacy) {
    efficacy = -kHighsInf;

    HighsCDouble scale = 1.0 / static_cast<HighsCDouble>(delta);
    HighsCDouble downrhs, f0, oneoveroneminusf0;
    computeRhs(scale, downrhs, f0, oneoveroneminusf0);
    if (f0 < f0min || f0 > f0max || oneoveroneminusf0 > maxCMirScale)
      return false;

    HighsCDouble contscale = scale * oneoveroneminusf0;
    HighsCDouble sqrnorm = contscale * contscale * continuoussqrnorm;
    HighsCDouble viol = contscale * continuouscontribution - downrhs;

    for (HighsInt j : integerinds) {
      updateViolationAndNorm(
          cut.vars[j],
          computeIntegralCoef(cut.vars[j].val, f0, oneoveroneminusf0, scale),
          viol, sqrnorm);
    }

    if (sqrnorm == 0.0) return false;

    efficacy = static_cast<double>(viol / sqrt(sqrnorm));
    return efficacy > bestefficacy;
  };

  for (double delta : deltas) {
    double efficacy;
    if (isCutBetter(delta, efficacy)) {
      bestdelta = delta;
      bestefficacy = efficacy;
    }
  }

  if (bestdelta == -1) return false;

  /* try if multiplying best delta by 2 4 or 8 gives a better efficacy */
  for (HighsInt k = 1; !onlyInitialCMIRScale && k <= 3; ++k) {
    double delta = bestdelta * (1 << k);
    double efficacy;

    if (isCutBetter(delta, efficacy)) {
      bestdelta = delta;
      bestefficacy = efficacy;
    }
  }

  assert(bestdelta != -1);

  // try to flip complementation of integers to increase efficacy
  for (HighsInt k : integerinds) {
    TransformedVar& var = cut.vars[k];
    if (var.upper == kHighsInf) continue;
    if (var.solval <= feastol) continue;

    cut.flipComplementation(var);

    double efficacy;

    if (isCutBetter(bestdelta, efficacy)) {
      bestefficacy = efficacy;
    } else {
      cut.flipComplementation(var);
    }
  }

  HighsCDouble scale = 1.0 / static_cast<HighsCDouble>(bestdelta);
  HighsCDouble downrhs, f0, oneoveroneminusf0;
  computeRhs(scale, downrhs, f0, oneoveroneminusf0);

  cut.rhs = downrhs * bestdelta;
  integralSupport = true;
  integralCoefficients = false;
  for (TransformedVar& var : cut.vars) {
    if (var.val == 0.0) continue;
    if (!var.integral) {
      if (var.val > 0.0)
        var.val = 0.0;
      else {
        var.val = static_cast<double>(var.val * oneoveroneminusf0);
        integralSupport = false;
      }
    } else {
      var.val = static_cast<double>(
          computeIntegralCoef(var.val, f0, oneoveroneminusf0, scale) *
          bestdelta);
    }
  }

#ifndef NDEBUG
  // Check if the computed cut has the correct efficacy
  {
    HighsCDouble checkviol = -cut.rhs;
    HighsCDouble checknorm = 0.0;
    for (const TransformedVar& var : cut.vars) {
      if (var.val == 0.0) continue;
      updateViolationAndNorm(var, var.val, checkviol, checknorm);
    }
    if (checknorm != 0.0) {
      double checkefficacy = static_cast<double>(checkviol / sqrt(checknorm));
      // the efficacy can become infinite if the cut 0 <= -1 is derived
      assert(abs(checkefficacy - bestefficacy) < 0.001 ||
             abs(checkefficacy) >= 1e30);
    }
  }
#endif

  return true;
}

// returns a power of two to scale the inequality with such that the given value
// is scaled to about one; the scaling is exact
static double computeScale(double val) {
  int expshift = 0;
  std::frexp(val, &expshift);
  expshift = -expshift;

  // don't scale the coefficients by more than +1024 (violations can increase)
  expshift = std::min(10, expshift);

  return std::ldexp(1.0, expshift);
}

bool HighsCutGeneration::postprocessCut(Row& row,
                                        const HighsDomain& globaldom) {
  std::vector<HighsInt>& inds = row.inds;
  std::vector<double>& vals = row.vals;
  HighsCDouble& rhs = row.rhs;

  // right hand sides slightly below zero are likely due to numerical errors and
  // can cause numerical troubles with scaling, so set them to zero
  if (rhs < 0 && rhs > -epsilon) rhs = 0;

  if (integralSupport && integralCoefficients) {
    // if the cut is known to be integral no postprocessing is needed and we
    // simply remove zero coefficients
    row.removeZeros();
    return true;
  }

  // determine maximal absolute coefficient
  double maxAbsValue = 0.0;
  for (double val : vals) maxAbsValue = std::max(std::abs(val), maxAbsValue);

  // determine minimal allowed coefficient
  double minCoefficientValue = 100 * feastol * std::max(maxAbsValue, 1e-3);

  // remove small coefficients and check whether the remaining support is
  // integral
  integralSupport = true;
  for (HighsInt i = static_cast<HighsInt>(vals.size()) - 1; i >= 0; --i) {
    if (vals[i] == 0) continue;
    if (std::abs(vals[i]) <= minCoefficientValue) {
      if (vals[i] < 0) {
        double ub = globaldom.col_upper_[inds[i]];
        if (ub == kHighsInf)
          return false;
        else
          rhs -= static_cast<HighsCDouble>(ub) * vals[i];
      } else {
        double lb = globaldom.col_lower_[inds[i]];
        if (lb == -kHighsInf)
          return false;
        else
          rhs -= static_cast<HighsCDouble>(lb) * vals[i];
      }

      vals[i] = 0.0;
      continue;
    }

    if (integralSupport && !lpRelaxation.isColIntegral(inds[i]))
      integralSupport = false;
  }

  // remove zeros in place
  row.removeZeros();

  if (vals.empty()) return false;

  if (integralSupport) {
    // integral support -> determine scale to make all coefficients integral
    double intscale = HighsIntegers::integralScale(vals, feastol, epsilon);

    bool scaleSmallestValToOne = true;

    if (intscale != 0.0 && intscale * std::max(1.0, maxAbsValue) <=
                               static_cast<double>(uint64_t{1} << 52)) {
      // A scale to make all value integral was found. The scale is only
      // rejected if it is in a range where not all integral values are
      // representable in double precision anymore or cannot be correctly
      // rounded by adding 0.5 and casting to an int. The latter starts at 2^52
      // + 1 since adding 0.5 will round upwards to the next even number for
      // that magnitude. So the largest acceptable value is 2^52 and when at
      // most that value we want to always use the scale to adjust the
      // coefficients and right hand side for numerical safety reasons. If the
      // resulting integral values are too large, however, we scale the cut down
      // by shifting the exponent.
      rhs.renormalize();
      rhs *= intscale;
      maxAbsValue = HighsIntegers::nearestInteger(maxAbsValue * intscale);
      for (HighsInt i = 0; i != static_cast<HighsInt>(vals.size()); ++i) {
        HighsCDouble scaleval = intscale * static_cast<HighsCDouble>(vals[i]);
        double intval =
            HighsIntegers::nearestInteger(static_cast<double>(scaleval));
        double delta = static_cast<double>(scaleval - intval);

        vals[i] = intval;

        // if the coefficient would be strengthened by rounding, we add the
        // upperbound constraint to make it exactly integral instead and
        // therefore weaken the right hand side
        if (delta < 0.0) {
          double ub = globaldom.col_upper_[inds[i]];
          if (ub == kHighsInf) return false;

          rhs -= static_cast<HighsCDouble>(delta) * ub;
        } else {
          double lb = globaldom.col_lower_[inds[i]];
          if (lb == -kHighsInf) return false;

          rhs -= static_cast<HighsCDouble>(delta) * lb;
        }
      }

      // finally we can round down the right hand side. Therefore in most cases
      // small errors for which the upper bound constraints where used and the
      // right hand side was weakened, do not weaken the final cut.
      rhs = floor(rhs + feastol);

      if (intscale * maxAbsValue * feastol < 0.5) {
        // integral scale leads to small enough values, accept scale
        scaleSmallestValToOne = false;
        integralCoefficients = true;
      }
    }

    if (scaleSmallestValToOne) {
      // integral scale lead to very large coefficient values. We now shift the
      // exactly integral values down such that the smallest coefficient is
      // around 1
      double minAbsValue = kHighsInf;
      for (double val : vals)
        minAbsValue = std::min(std::abs(val), minAbsValue);

      row.scale(computeScale(minAbsValue - epsilon));
    }
  } else {
    // the support is not integral, scale cut to have the largest coefficient
    // around 1.0
    row.scale(computeScale(maxAbsValue - epsilon));
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
  HighsInt numZeros = 0;

  const HighsInt numVars = cut.vars.size();
  HighsCDouble maxact = -feastol;
  double maxAbsVal = 0;
  for (const TransformedVar& var : cut.vars)
    maxAbsVal = std::max(std::abs(var.val), maxAbsVal);

  initialScale = computeScale(maxAbsVal);
  cut.scale(initialScale);

  for (TransformedVar& var : cut.vars) {
    // we do not want to have integral variables with small coefficients as this
    // may lead to numerical instabilities during cut generation
    // Therefore we relax integral variables with small coefficients to
    // continuous ones because they still might have a non-negligible
    // contribution e.g. when they come from integral rows with coefficients on
    // the larger side. When we relax them to continuous variables they will be
    // complemented so that their solution value is closest to zero and then
    // will be relaxed if their value is positive or their maximal contribution
    // is below feasibility tolerance.
    var.integral = var.integral && std::abs(var.val) > 10 * feastol;

    if (!var.integral) {
      // complement non-integer variable (cmir separation heuristic complements
      // integral variables in the same way)
      if (var.upper < 2 * var.solval) cut.flipComplementation(var);

      // relax positive continuous variables and those with small contributions
      if (var.val > 0 || std::abs(var.val) * var.upper <= 10 * feastol) {
        // printf("remove: val = %g  upper = %g\n", var.val, var.upper);
        if (var.val < 0) {
          if (var.upper == kHighsInf) return false;
          cut.rhs -= static_cast<HighsCDouble>(var.upper) * var.val;
        }

        ++numZeros;
        var.val = 0.0;
        continue;
      }

      hasContinuous = true;
    } else {
      if (var.upper == kHighsInf) {
        hasUnboundedInts = true;
        hasGeneralInts = true;
      } else if (var.upper != 1.0) {
        hasGeneralInts = true;
      }

      if (var.val > 0) maxact += static_cast<HighsCDouble>(var.upper) * var.val;
    }
  }

  HighsInt maxLen = 100 + 0.15 * (lpRelaxation.numCols());

  if (numVars - numZeros > maxLen) {
    HighsInt numCancel = numVars - numZeros - maxLen;
    std::vector<HighsInt> cancelNzs;

    for (HighsInt i = 0; i != numVars; ++i) {
      const TransformedVar& var = cut.vars[i];
      double cancelSlack = var.val > 0 ? var.solval : var.flippedSolval;
      if (cancelSlack <= feastol) cancelNzs.push_back(i);
    }

    if (static_cast<HighsInt>(cancelNzs.size()) < numCancel) return false;
    if (static_cast<HighsInt>(cancelNzs.size()) > numCancel)
      std::partial_sort(cancelNzs.begin(), cancelNzs.begin() + numCancel,
                        cancelNzs.end(), [&](HighsInt a, HighsInt b) {
                          return std::abs(cut.vars[a].val) <
                                 std::abs(cut.vars[b].val);
                        });

    for (HighsInt i = 0; i < numCancel; ++i) {
      TransformedVar& var = cut.vars[cancelNzs[i]];

      if (var.val < 0) {
        cut.rhs -= static_cast<HighsCDouble>(var.upper) * var.val;
      } else
        maxact -= static_cast<HighsCDouble>(var.upper) * var.val;

      var.val = 0.0;
    }

    numZeros += numCancel;
  }

  if (numZeros != 0) {
    // remove zeros in place
    for (HighsInt i = numVars - 1; i >= 0; --i) {
      if (cut.vars[i].val == 0.0) {
        cut.vars[i] = cut.vars.back();
        cut.vars.pop_back();
        if (--numZeros == 0) break;
      }
    }
  }

  return maxact > cut.rhs;
}

#if 0
static void checkNumerics(const double* vals, HighsInt len, double rhs) {
  double maxAbsCoef = 0.0;
  double minAbsCoef = kHighsInf;
  HighsCDouble sqrnorm = 0;
  for (HighsInt i = 0; i < len; ++i) {
    sqrnorm += static_cast<HighsCDouble>(vals[i]) * vals[i];
    maxAbsCoef = std::max(std::abs(vals[i]), maxAbsCoef);
    minAbsCoef = std::min(std::abs(vals[i]), minAbsCoef);
  }

  double norm = static_cast<double>(sqrt(sqrnorm));

  // printf("length: %" HIGHSINT_FORMAT
  //       ", minCoef: %g, maxCoef, %g, norm %g, rhs: %g, dynamism=%g\n",
  //       len, minAbsCoef, maxAbsCoef, norm, rhs, maxAbsCoef / minAbsCoef);
}
#endif

bool HighsCutGeneration::generateCut(HighsTransformedLp& transLp,
                                     std::vector<HighsInt>& inds_,
                                     std::vector<double>& vals_, double& rhs_,
                                     bool onlyInitialCMIRScale) {
#if 0
  if (vals_.size() > 1) {
    bool intsPositive = true;
    if (!transLp.transform(Row{inds_, vals_, rhs_}, cut, intsPositive))
      return false;

    bool hasUnboundedInts = false;
    bool hasGeneralInts = false;
    bool hasContinuous = false;
    // printf("before preprocessing of base inequality:\n");
    std::vector<HighsInt> indsCheck_;
    std::vector<double> valsCheck_;
    Row checkRow{indsCheck_, valsCheck_, 0.0};
    cut.store(checkRow);
    checkNumerics(valsCheck_.data(), valsCheck_.size(),
                  static_cast<double>(checkRow.rhs));
    if (!preprocessBaseInequality(hasUnboundedInts, hasGeneralInts,
                                  hasContinuous))
      return false;
    cut.store(checkRow);
    // printf("after preprocessing of base inequality:\n");
    checkNumerics(valsCheck_.data(), valsCheck_.size(),
                  static_cast<double>(checkRow.rhs));

    if (!transLp.untransform(checkRow)) return false;

    // printf("after untransform of base inequality:\n");
    checkNumerics(valsCheck_.data(), valsCheck_.size(),
                  static_cast<double>(checkRow.rhs));

    // finally check whether the cut is violated
    lpRelaxation.getMipSolver().mipdata_->debugSolution.checkCut(
        indsCheck_.data(), valsCheck_.data(), indsCheck_.size(),
        static_cast<double>(checkRow.rhs));
  }
#endif

  // the base inequality, which is overwritten by the cut
  Row row{inds_, vals_, rhs_};
  bool intsPositive = true;
  if (!transLp.transform(row, cut, intsPositive)) return false;

  bool hasUnboundedInts = false;
  bool hasGeneralInts = false;
  bool hasContinuous = false;
  if (!preprocessBaseInequality(hasUnboundedInts, hasGeneralInts,
                                hasContinuous))
    return false;

  // it can happen that there is an unbounded integer variable during the
  // transform call so that the integers are not transformed to positive values.
  // Now the call to preprocessBaseInequality may have removed the unbounded
  // integer, e.g. due to a small coefficient value, so that we can still use
  // the lifted inequalities instead of cmir. We need to make sure, however,
  // that the cut values are transformed to positive coefficients first, which
  // we do below.
  if (!hasUnboundedInts && !intsPositive) {
    for (TransformedVar& var : cut.vars) {
      if (var.val > 0 || !var.integral) continue;

      cut.flipComplementation(var);
    }
  }

  // try to generate a cut
  if (!tryGenerateCut(hasUnboundedInts, hasGeneralInts, hasContinuous,
                      10 * feastol, onlyInitialCMIRScale))
    return false;

  // remove the complementation if exists
  cut.removeComplementation();

  // store the cut in the given vectors and remove zeros
  cut.store(row);
  row.removeZeros();

  // transform the cut back into the original space, i.e. remove the bound
  // substitution and replace implicit slack variables
  if (!transLp.untransform(row)) return false;

  lpRelaxation.getMipSolver().mipdata_->debugSolution.checkCut(
      inds_.data(), vals_.data(), inds_.size(), static_cast<double>(row.rhs));
  // apply cut postprocessing including scaling and removal of small
  // coefficients
  if (!postprocessCut(row, transLp.getGlobaldom())) return false;
  rhs_ = static_cast<double>(row.rhs);

  lpRelaxation.getMipSolver().mipdata_->debugSolution.checkCut(
      inds_.data(), vals_.data(), inds_.size(), rhs_);

  // finally determine the violation of the cut in the original space
  if (row.violation(lpRelaxation.getSolution().col_value) <= 10 * feastol)
    return false;

  transLp.getGlobaldom().tightenCoefficients(inds_.data(), vals_.data(),
                                             inds_.size(), rhs_);

  // if the cut is violated by a small factor above the feasibility
  // tolerance, add it to the cutpool
  HighsInt cutindex = cutpool.addCut(lpRelaxation.getMipSolver(), inds_.data(),
                                     vals_.data(), inds_.size(), rhs_,
                                     integralSupport && integralCoefficients);

  // only return true if cut was accepted by the cutpool, i.e. not a duplicate
  // of a cut already in the pool
  return cutindex != -1;
}

bool HighsCutGeneration::generateConflict(const HighsDomain& localdomain,
                                          const HighsDomain& globaldom,
                                          std::vector<HighsInt>& proofinds,
                                          std::vector<double>& proofvals,
                                          double& proofrhs) {
  // the proof, which is overwritten by the cut
  Row row{proofinds, proofvals, proofrhs};

  lpRelaxation.getMipSolver().mipdata_->debugSolution.checkCut(
      proofinds.data(), proofvals.data(), proofinds.size(), proofrhs);

  cut.load(row);

  HighsCDouble activity = 0.0;
  for (TransformedVar& var : cut.vars) {
    HighsInt col = var.col;

    var.upper = globaldom.col_upper_[col] - globaldom.col_lower_[col];
    var.integral = lpRelaxation.isColIntegral(col);

    var.solval =
        var.val < 0
            ? std::min(globaldom.col_upper_[col], localdomain.col_upper_[col])
            : std::max(globaldom.col_lower_[col], localdomain.col_lower_[col]);
    if (var.val < 0 && globaldom.col_upper_[col] != kHighsInf) {
      cut.rhs -= static_cast<HighsCDouble>(globaldom.col_upper_[col]) * var.val;
      var.val = -var.val;
      var.complemented = true;

      var.solval = globaldom.col_upper_[col] - var.solval;
    } else {
      cut.rhs -= static_cast<HighsCDouble>(globaldom.col_lower_[col]) * var.val;
      var.complemented = false;
      var.solval = var.solval - globaldom.col_lower_[col];
    }

    activity += static_cast<HighsCDouble>(var.solval) * var.val;
    var.flippedSolval = var.upper - var.solval;
  }

  if (activity > cut.rhs) {
    double solScale = static_cast<double>(cut.rhs / activity);
    for (TransformedVar& var : cut.vars) {
      var.solval *= solScale;
      var.flippedSolval = var.upper - var.solval;
    }
  }

  bool hasUnboundedInts = false;
  bool hasGeneralInts = false;
  bool hasContinuous = false;

  if (!preprocessBaseInequality(hasUnboundedInts, hasGeneralInts,
                                hasContinuous))
    return false;

  // try to generate a cut
  if (!tryGenerateCut(hasUnboundedInts, hasGeneralInts, hasContinuous, feastol,
                      false, false, false))
    return false;

  // remove the complementation
  for (TransformedVar& var : cut.vars) {
    if (var.complemented) {
      cut.rhs -=
          static_cast<HighsCDouble>(globaldom.col_upper_[var.col]) * var.val;
      var.val = -var.val;
    } else
      cut.rhs +=
          static_cast<HighsCDouble>(globaldom.col_lower_[var.col]) * var.val;
  }

  // store the cut in the given vectors
  cut.store(row);

  // apply cut postprocessing including scaling and removal of small
  // coefficients
  if (!postprocessCut(row, globaldom)) return false;

  proofrhs = static_cast<double>(row.rhs);

  bool cutintegral = integralSupport && integralCoefficients;

  globaldom.tightenCoefficients(proofinds.data(), proofvals.data(),
                                proofinds.size(), proofrhs);

  HighsInt cutindex = cutpool.addCut(
      lpRelaxation.getMipSolver(), proofinds.data(), proofvals.data(),
      proofinds.size(), proofrhs, cutintegral, true, true, true);

  // only return true if cut was accepted by the cutpool, i.e. not a duplicate
  // of a cut already in the pool
  return cutindex != -1;
}

bool HighsCutGeneration::finalizeAndAddCut(const HighsDomain& globaldom,
                                           Row& row) {
  // remove zeros in place and check whether the support is integral
  row.removeZeros();
  integralSupport = true;
  integralCoefficients = false;
  for (HighsInt col : row.inds)
    integralSupport &= lpRelaxation.isColIntegral(col);

  lpRelaxation.getMipSolver().mipdata_->debugSolution.checkCut(
      row.inds.data(), row.vals.data(), row.inds.size(),
      static_cast<double>(row.rhs));
  // apply cut postprocessing including scaling and removal of small
  // coefficients
  if (!postprocessCut(row, globaldom)) return false;
  double rhs = static_cast<double>(row.rhs);

  lpRelaxation.getMipSolver().mipdata_->debugSolution.checkCut(
      row.inds.data(), row.vals.data(), row.inds.size(), rhs);

  // finally determine the violation of the cut in the original space
  if (row.violation(lpRelaxation.getSolution().col_value) <= 10 * feastol)
    return false;

  globaldom.tightenCoefficients(row.inds.data(), row.vals.data(),
                                row.inds.size(), rhs);

  // if the cut is violated by a small factor above the feasibility
  // tolerance, add it to the cutpool
  HighsInt cutindex = cutpool.addCut(
      lpRelaxation.getMipSolver(), row.inds.data(), row.vals.data(),
      row.inds.size(), rhs, integralSupport && integralCoefficients);

  // only return true if cut was accepted by the cutpool, i.e. not a duplicate
  // of a cut already in the pool
  return cutindex != -1;
}

template <typename T>
void HighsCutGeneration::updateViolationAndNorm(const TransformedVar& var, T aj,
                                                HighsCDouble& violation,
                                                HighsCDouble& norm) const {
  // update violation
  violation += static_cast<HighsCDouble>(aj) * var.solval;

  // update norm
  if (aj > 0 && var.solval <= feastol) return;
  if (aj < 0 && var.flippedSolval <= feastol) return;
  norm += static_cast<HighsCDouble>(aj) * static_cast<HighsCDouble>(aj);
}

bool HighsCutGeneration::tryGenerateCut(bool hasUnboundedInts,
                                        bool hasGeneralInts, bool hasContinuous,
                                        double minEfficacy,
                                        bool onlyInitialCMIRScale,
                                        bool allowRejectCut, bool lpSol) {
  // use cmir if there are unbounded integer variables
  if (hasUnboundedInts)
    return cmirCutGenerationHeuristic(minEfficacy, onlyInitialCMIRScale);

  // 0. Save data before determining cover and applying lifting functions
  tmpCut = cut;

  // 1. Determine a cover, cover does not need to be minimal as neither of
  //    the lifting functions have minimality of the cover as necessary facet
  //    condition
  bool success = false;
  bool saveIntegralSupport = false;
  bool saveIntegralCoefficients = false;
  do {
    if (!determineCover(lpSol)) break;

    // 2. use superadditive lifting function depending on structure of base
    //    inequality:
    //    We have 3 lifting functions available for pure binary knapsack sets,
    //    for mixed-binary knapsack sets and for mixed integer knapsack sets.
    if (!hasContinuous && !hasGeneralInts) {
      separateLiftedKnapsackCover();
      success = true;
    } else if (hasGeneralInts) {
      success = separateLiftedMixedIntegerCover();
    } else {
      assert(hasContinuous);
      assert(!hasGeneralInts);
      success = separateLiftedMixedBinaryCover();
    }
  } while (false);

  double minMirEfficacy = minEfficacy;
  if (success) {
    // save data that might otherwise be overwritten when calling the cmir
    // separator
    saveIntegralSupport = integralSupport;
    saveIntegralCoefficients = integralCoefficients;

    // compute violation and squared norm
    HighsCDouble violation = -cut.rhs;
    HighsCDouble sqrnorm = 0.0;
    for (const TransformedVar& var : cut.vars)
      updateViolationAndNorm(var, var.val, violation, sqrnorm);

    // compute efficacy (distance cut off)
    double efficacy = static_cast<double>(violation / sqrt(sqrnorm));
    if (allowRejectCut && efficacy <= minEfficacy) {
      // reject cut
      success = false;
    } else {
      // accept cut and increase minimum efficiency requirement for cmir cut
      minMirEfficacy += efficacy;
    }
  }

  // restore the base inequality for the cmir separator and keep the lifted cut
  // in tmpCut
  std::swap(cut, tmpCut);

  if (cmirCutGenerationHeuristic(minMirEfficacy, onlyInitialCMIRScale)) {
    // take the cmir cut as it is better
    return true;
  } else if (success) {
    // take the previous lifted cut as cmir could not improve; note that the
    // complementation and the solution values are restored as well since the
    // cmir separator may have flipped the complementation
    std::swap(cut, tmpCut);
    // restore indicators
    integralSupport = saveIntegralSupport;
    integralCoefficients = saveIntegralCoefficients;
    return true;
  } else
    // neither cmir nor lifted cut successful
    return false;
}
