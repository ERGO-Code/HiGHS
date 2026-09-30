#ifndef PDLP_HIPDLP_BACKEND_CUDA_CERTIFICATE_SCREEN_HPP_
#define PDLP_HIPDLP_BACKEND_CUDA_CERTIFICATE_SCREEN_HPP_

#include <cfloat>
#include <cmath>
#include <limits>

#ifdef __CUDACC__
#define HIPDLP_CERT_INLINE __host__ __device__ inline
#else
#define HIPDLP_CERT_INLINE inline
#endif

namespace hipdlp {
namespace backend {
namespace cuda {
namespace certificate {

struct Interval {
  double lower;
  double upper;
};

struct NormBounds {
  double scale;
  double lower_square;
  double upper_square;
};

struct SumBounds {
  Interval value;
  double absolute_sum;
};

// Rows longer than this are split across warps in the row-activity pass.
const int kSegmentLength = 1024;

struct ScreenPartial {
  SumBounds objective;
  SumBounds row_support;
  double primal_violation;
  double farkas_violation;
  int unsafe;
};

struct ScreenFlags {
  bool primal_ray_possible;
  bool farkas_possible;
};

HIPDLP_CERT_INLINE double addDown(const double left, const double right) {
#ifdef __CUDA_ARCH__
  return __dadd_rd(left, right);
#else
  if (left == 0.0) return right;
  if (right == 0.0) return left;
  volatile double rounded = left + right;
  return std::nextafter(rounded, -HUGE_VAL);
#endif
}

HIPDLP_CERT_INLINE double addUp(const double left, const double right) {
#ifdef __CUDA_ARCH__
  return __dadd_ru(left, right);
#else
  if (left == 0.0) return right;
  if (right == 0.0) return left;
  volatile double rounded = left + right;
  return std::nextafter(rounded, HUGE_VAL);
#endif
}

HIPDLP_CERT_INLINE double multiplyDown(const double left, const double right) {
#ifdef __CUDA_ARCH__
  return __dmul_rd(left, right);
#else
  if (left == 0.0 || right == 0.0) return 0.0;
  volatile double rounded = left * right;
  return std::nextafter(rounded, -HUGE_VAL);
#endif
}

HIPDLP_CERT_INLINE double multiplyUp(const double left, const double right) {
#ifdef __CUDA_ARCH__
  return __dmul_ru(left, right);
#else
  if (left == 0.0 || right == 0.0) return 0.0;
  volatile double rounded = left * right;
  return std::nextafter(rounded, HUGE_VAL);
#endif
}

HIPDLP_CERT_INLINE double divideDown(const double left, const double right) {
#ifdef __CUDA_ARCH__
  return __ddiv_rd(left, right);
#else
  if (left == 0.0 && right > 0.0) return 0.0;
  volatile double rounded = left / right;
  return std::nextafter(rounded, -HUGE_VAL);
#endif
}

HIPDLP_CERT_INLINE double divideUp(const double left, const double right) {
#ifdef __CUDA_ARCH__
  return __ddiv_ru(left, right);
#else
  if (left == 0.0 && right > 0.0) return 0.0;
  volatile double rounded = left / right;
  return std::nextafter(rounded, HUGE_VAL);
#endif
}

HIPDLP_CERT_INLINE double sqrtDown(const double value) {
#ifdef __CUDA_ARCH__
  return __dsqrt_rd(value);
#else
  if (value == 0.0) return 0.0;
  return std::nextafter(std::sqrt(value), 0.0);
#endif
}

HIPDLP_CERT_INLINE double sqrtUp(const double value) {
#ifdef __CUDA_ARCH__
  return __dsqrt_ru(value);
#else
  if (value == 0.0) return 0.0;
  return std::nextafter(std::sqrt(value), HUGE_VAL);
#endif
}

HIPDLP_CERT_INLINE bool finiteValue(const double value) {
  return value >= -DBL_MAX && value <= DBL_MAX;
}

HIPDLP_CERT_INLINE bool finite(const Interval value) {
  return finiteValue(value.lower) && finiteValue(value.upper) &&
         value.lower <= value.upper;
}

HIPDLP_CERT_INLINE bool positiveFinite(const Interval value) {
  return finite(value) && value.lower > 0.0;
}

HIPDLP_CERT_INLINE Interval unknown() { return {-HUGE_VAL, HUGE_VAL}; }

HIPDLP_CERT_INLINE double maximum(const double left, const double right) {
  return left > right ? left : right;
}

HIPDLP_CERT_INLINE double minimum(const double left, const double right) {
  return left < right ? left : right;
}

HIPDLP_CERT_INLINE double absoluteLower(const Interval value) {
  if (value.lower > 0.0) return value.lower;
  if (value.upper < 0.0) return -value.upper;
  return 0.0;
}

HIPDLP_CERT_INLINE double absoluteUpper(const Interval value) {
  return maximum(fabs(value.lower), fabs(value.upper));
}

HIPDLP_CERT_INLINE Interval add(const Interval left, const Interval right) {
  if (!finite(left) || !finite(right)) return unknown();
  return {addDown(left.lower, right.lower), addUp(left.upper, right.upper)};
}

HIPDLP_CERT_INLINE Interval multiply(const Interval value,
                                     const double factor) {
  if (!finite(value) || !finiteValue(factor)) return unknown();
  if (factor >= 0.0)
    return {multiplyDown(value.lower, factor), multiplyUp(value.upper, factor)};
  return {multiplyDown(value.upper, factor), multiplyUp(value.lower, factor)};
}

HIPDLP_CERT_INLINE Interval divide(const Interval value,
                                   const Interval denominator) {
  if (!finite(value) || !positiveFinite(denominator)) return unknown();
  return {divideDown(value.lower, value.lower >= 0.0 ? denominator.upper
                                                     : denominator.lower),
          divideUp(value.upper,
                   value.upper >= 0.0 ? denominator.lower : denominator.upper)};
}

HIPDLP_CERT_INLINE Interval recoverProduct(const double value,
                                           const double scale) {
  if (!finiteValue(value) || !finiteValue(scale) || scale <= 0.0)
    return unknown();
  return multiply({value, value}, scale);
}

HIPDLP_CERT_INLINE Interval recoverQuotient(const double value,
                                            const double scale) {
  return divide({value, value}, {scale, scale});
}

HIPDLP_CERT_INLINE bool safeRecovery(const Interval recovered,
                                     const double source) {
  return finite(recovered) && (source == 0.0 || absoluteLower(recovered) > 0.0);
}

HIPDLP_CERT_INLINE void mergeNorm(NormBounds& left, const NormBounds& right) {
  if (!finiteValue(left.scale) || !finiteValue(right.scale) ||
      !finiteValue(left.upper_square) || !finiteValue(right.upper_square)) {
    left.scale = HUGE_VAL;
    return;
  }
  if (right.scale == 0.0) return;
  if (left.scale == 0.0) {
    left = right;
    return;
  }
  const double scale = maximum(left.scale, right.scale);
  const double left_lower = maximum(0.0, divideDown(left.scale, scale));
  const double left_upper = divideUp(left.scale, scale);
  const double right_lower = maximum(0.0, divideDown(right.scale, scale));
  const double right_upper = divideUp(right.scale, scale);
  left.lower_square = maximum(
      0.0, addDown(multiplyDown(multiplyDown(left.lower_square, left_lower),
                                left_lower),
                   multiplyDown(multiplyDown(right.lower_square, right_lower),
                                right_lower)));
  left.upper_square = addUp(
      multiplyUp(multiplyUp(left.upper_square, left_upper), left_upper),
      multiplyUp(multiplyUp(right.upper_square, right_upper), right_upper));
  left.scale = scale;
}

HIPDLP_CERT_INLINE void addNorm(NormBounds& norm, const Interval value) {
  if (!finite(value)) {
    norm.scale = HUGE_VAL;
    return;
  }
  NormBounds term = {};
  term.scale = absoluteUpper(value);
  if (term.scale == 0.0) return;
  const double lower =
      maximum(0.0, divideDown(absoluteLower(value), term.scale));
  term.lower_square = maximum(0.0, multiplyDown(lower, lower));
  term.upper_square = 1.0;
  mergeNorm(norm, term);
}

HIPDLP_CERT_INLINE double errorGamma(const double operations) {
  if (!finiteValue(operations) || operations < 0.0) return HUGE_VAL;
  const double product = multiplyUp(operations, DBL_EPSILON / 2.0);
  const double denominator = addDown(1.0, -product);
  return denominator > 0.0 ? divideUp(product, denominator) : HUGE_VAL;
}

HIPDLP_CERT_INLINE double subnormalError(const double operations,
                                         const double gamma) {
  const double smallest = DBL_MIN * DBL_EPSILON;
  return multiplyUp(multiplyUp(operations, smallest), addUp(1.0, gamma));
}

HIPDLP_CERT_INLINE Interval hostNormBounds(const NormBounds& norm,
                                           const double hypot_count) {
  if (!finiteValue(norm.scale) || !finiteValue(norm.upper_square))
    return unknown();
  if (norm.scale == 0.0) return {0.0, 0.0};
  const double operations = multiplyUp(8.0, addUp(hypot_count, 1.0));
  const double gamma = errorGamma(operations);
  if (!finiteValue(gamma) || gamma >= 1.0) return unknown();
  const double absolute = subnormalError(operations, gamma);
  const double lower = multiplyDown(norm.scale, sqrtDown(norm.lower_square));
  const double upper = multiplyUp(norm.scale, sqrtUp(norm.upper_square));
  return {maximum(0.0, addDown(multiplyDown(lower, addDown(1.0, -gamma)),
                               -absolute)),
          addUp(multiplyUp(upper, addUp(1.0, gamma)), absolute)};
}

HIPDLP_CERT_INLINE void addSum(SumBounds& sum, const Interval term) {
  sum.value = add(sum.value, term);
  sum.absolute_sum = addUp(sum.absolute_sum, absoluteUpper(term));
}

HIPDLP_CERT_INLINE void mergeSum(SumBounds& left, const SumBounds& right) {
  left.value = add(left.value, right.value);
  left.absolute_sum = addUp(left.absolute_sum, right.absolute_sum);
}

HIPDLP_CERT_INLINE Interval serialSumBounds(const SumBounds& sum,
                                            const double additions) {
  if (!finite(sum.value) || !finiteValue(sum.absolute_sum)) return unknown();
  if (sum.absolute_sum == 0.0) return {0.0, 0.0};
  const double gamma = errorGamma(additions);
  if (!finiteValue(gamma)) return unknown();
  const double error = addUp(multiplyUp(gamma, sum.absolute_sum),
                             subnormalError(additions, gamma));
  return {addDown(sum.value.lower, -error), addUp(sum.value.upper, error)};
}

HIPDLP_CERT_INLINE double violationLower(const Interval value,
                                         const double lower, const double upper,
                                         const bool multiplier) {
  if (!finite(value)) return 0.0;
  const bool has_lower = finiteValue(lower);
  const bool has_upper = finiteValue(upper);
  if (has_lower && has_upper) return multiplier ? 0.0 : absoluteLower(value);
  if (!has_lower && !has_upper) return multiplier ? absoluteLower(value) : 0.0;
  if (has_lower) return maximum(0.0, -value.upper);
  return maximum(0.0, value.lower);
}

HIPDLP_CERT_INLINE bool rejectViolation(const double lower,
                                        const Interval verifier_norm,
                                        const double tolerance) {
  return positiveFinite(verifier_norm) && finiteValue(tolerance) &&
         tolerance > 0.0 && divideDown(lower, verifier_norm.upper) > tolerance;
}

HIPDLP_CERT_INLINE bool rejectObjective(const Interval objective,
                                        const Interval verifier_norm,
                                        const double sense_sign,
                                        const double tolerance) {
  const Interval direction = divide(objective, verifier_norm);
  if (!finite(direction) || !finiteValue(tolerance) || tolerance <= 0.0)
    return false;
  const double margin = sense_sign > 0.0 ? -direction.lower : direction.upper;
  return margin <= tolerance;
}

HIPDLP_CERT_INLINE Interval rowSupport(const Interval multiplier,
                                       const double lower, const double upper) {
  if (!finite(multiplier)) return unknown();
  if (multiplier.lower == 0.0 && multiplier.upper == 0.0) return {0.0, 0.0};
  if (multiplier.lower >= 0.0) return multiply(multiplier, lower);
  if (multiplier.upper <= 0.0) return multiply(multiplier, upper);
  const Interval positive = multiply({0.0, multiplier.upper}, lower);
  const Interval negative = multiply({multiplier.lower, 0.0}, upper);
  if (!finite(positive) || !finite(negative)) return unknown();
  return {minimum(positive.lower, negative.lower),
          maximum(positive.upper, negative.upper)};
}

HIPDLP_CERT_INLINE void mergeScreen(ScreenPartial& left,
                                    const ScreenPartial& right) {
  mergeSum(left.objective, right.objective);
  mergeSum(left.row_support, right.row_support);
  left.primal_violation =
      maximum(left.primal_violation, right.primal_violation);
  left.farkas_violation =
      maximum(left.farkas_violation, right.farkas_violation);
  left.unsafe |= right.unsafe;
}

HIPDLP_CERT_INLINE ScreenFlags
finishScreen(const ScreenPartial& partial, const Interval primal_verifier_norm,
             const Interval dual_verifier_norm, const double num_cols,
             const double num_rows, const double sense_sign,
             const bool column_bounds_straddle_zero, const bool unsafe,
             const double tolerance) {
  ScreenFlags flags = {true, true};
  if (unsafe || partial.unsafe || !positiveFinite(primal_verifier_norm) ||
      !positiveFinite(dual_verifier_norm) || !finiteValue(tolerance) ||
      tolerance <= 0.0 || (sense_sign != 1.0 && sense_sign != -1.0))
    return flags;
  flags.primal_ray_possible =
      !rejectViolation(partial.primal_violation, primal_verifier_norm,
                       tolerance) &&
      !rejectObjective(serialSumBounds(partial.objective, num_cols),
                       primal_verifier_norm, sense_sign, tolerance);
  flags.farkas_possible =
      !rejectViolation(partial.farkas_violation, dual_verifier_norm, tolerance);
  if (column_bounds_straddle_zero) {
    const Interval support = serialSumBounds(partial.row_support, num_rows);
    if (finite(support) && support.upper <= 0.0) flags.farkas_possible = false;
  }
  return flags;
}

}  // namespace certificate
}  // namespace cuda
}  // namespace backend
}  // namespace hipdlp

#undef HIPDLP_CERT_INLINE

#endif