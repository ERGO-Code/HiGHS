/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_UTIL_INTEGERS_H_
#define HIGHS_UTIL_INTEGERS_H_

#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "util/HighsCDouble.h"
#include "util/HighsInt.h"

class HighsIntegers {
 public:
  static int64_t gcd(int64_t a, int64_t b) {
    int64_t h;
    if (a < 0) a = -a;
    if (b < 0) b = -b;

    if (a == 0) return b;
    if (b == 0) return a;

    do {
      h = a % b;
      a = b;
      b = h;
    } while (b != 0);

    return a;
  }

  // computes a rational approximation with given maximal denominator
  static int64_t denominator(double x, double eps, int64_t maxdenom) {
    int64_t ai = (int64_t)x;
    int64_t m[] = {ai, 1, 1, 0};

    HighsCDouble xi = x;
    HighsCDouble fraction = xi - double(ai);

    while (fraction > eps) {
      xi = 1.0 / fraction;
      if (double(xi) > double(int64_t{1} << 53)) break;

      ai = (int64_t)(double)xi;
      int64_t t = m[2] * ai + m[3];
      if (t > maxdenom) break;

      m[3] = m[2];
      m[2] = t;

      t = m[0] * ai + m[1];
      m[1] = m[0];
      m[0] = t;

      fraction = xi - ai;
    }

    ai = (maxdenom - m[3]) / m[2];
    m[1] += m[0] * ai;
    m[3] += m[2] * ai;

    double x0 = m[0] / (double)m[2];
    double x1 = m[1] / (double)m[3];
    x = std::abs(x);
    double err0 = std::abs(x - x0);
    double err1 = std::abs(x - x1);

    if (err0 < err1) return m[2];
    return m[3];
  }

  static double integralScale(const std::vector<double>& vals, double deltadown,
                              double deltaup) {
    if (vals.empty()) return 0.0;

    double minval = *std::min_element(
        vals.begin(), vals.end(),
        [](double a, double b) { return std::abs(a) < std::abs(b); });
    HighsInt numVals = vals.size();

    int expshift;

    // to cover many small denominators at once use a denominator of 75 * 2^n
    // with n-3 being large enough so that the smallest value is not below 0.5.
    std::frexp(minval, &expshift);
    expshift = std::max(-expshift, 0) + 3;

    uint64_t denom = 75 << expshift;
    HighsCDouble startdenom = denom;
    // now check if the values are integral and if not compute a common
    // denominator for their remaining fraction
    HighsCDouble val = startdenom * vals[0];
    HighsCDouble downval = floor(val + deltaup);
    HighsCDouble fraction = val - downval;

    if (fraction > deltadown) {
      // use a continued fraction algorithm to compute small missing
      // denominators for the remaining fraction
      denom *= denominator(double(fraction), deltaup, 1000);
      val = denom * vals[0];
      downval = floor(val + deltaup);
      fraction = val - downval;

      // if this is not sufficient for reaching integrality, we stop here
      if (fraction > deltadown) return 0.0;
    }

    uint64_t currgcd = (uint64_t)std::abs(double(downval));

    for (HighsInt i = 1; i != numVals; ++i) {
      val = denom * HighsCDouble(vals[i]);
      downval = floor(val + deltaup);
      fraction = val - downval;

      if (fraction > deltadown) {
        val = startdenom * vals[i];
        fraction = val - floor(val);
        denom *= denominator(double(fraction), deltaup, 1000);
        val = denom * vals[i];
        downval = floor(val + deltaup);
        fraction = val - downval;

        if (fraction > deltadown) return 0.0;
      }

      if (currgcd != 1) {
        currgcd = gcd(currgcd, (int64_t) double(downval));

        // if the denominator is large, divide by the current gcd to prevent
        // unecessary overflows
        if (denom > std::numeric_limits<unsigned int>::max()) {
          denom /= currgcd;
          currgcd = 1;
        }
      }
    }

    return denom / (double)currgcd;
  }
};

#endif
