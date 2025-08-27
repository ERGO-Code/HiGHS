#include "VectorOperations.h"

#include <cassert>
#include <cmath>

namespace hipo {

void vectorAdd(std::vector<double>& v1, const std::vector<double>& v2,
               double beta, double alpha) {
  for (Int i = 0; i < v1.size(); ++i) {
    v1[i] = alpha * v1[i] + beta * v2[i];
  }
}

void vectorAdd(std::vector<double>& v1, const double alpha) {
  for (Int i = 0; i < v1.size(); ++i) {
    v1[i] += alpha;
  }
}

void vectorMultiply(std::vector<double>& v1, const std::vector<double>& v2,
                    double alpha, double beta) {
  for (Int i = 0; i < v1.size(); ++i) {
    v1[i] = alpha * v1[i] * v2[i] + beta;
  }
}

void vectorAddMult(std::vector<double>& v1, const std::vector<double>& v2,
                   const std::vector<double>& v3, double beta, double alpha) {
  for (Int i = 0; i < v1.size(); ++i) {
    v1[i] = alpha * v1[i] + beta * v2[i] * v3[i];
  }
}

void vectorDivide(std::vector<double>& v1, const std::vector<double>& v2) {
  for (Int i = 0; i < v1.size(); ++i) {
    v1[i] /= v2[i];
  }
}

void vectorScale(std::vector<double>& v1, double alpha) {
  for (Int i = 0; i < v1.size(); ++i) {
    v1[i] *= alpha;
  }
}

double dotProd(const std::vector<double>& v1, const std::vector<double>& v2) {
  double result{};
  for (Int i = 0; i < v1.size(); ++i) {
    result += v1[i] * v2[i];
  }
  return result;
}

double norm2(const std::vector<double>& x) {
  double norm{};
  for (Int i = 0; i < x.size(); ++i) {
    norm += (x[i] * x[i]);
  }
  return std::sqrt(norm);
}

double norm2(const std::vector<double>& x0, const std::vector<double>& x1) {
  double norm{};
  for (Int i = 0; i < x0.size(); ++i) {
    norm += (x0[i] * x0[i]);
  }
  for (Int i = 0; i < x1.size(); ++i) {
    norm += (x1[i] * x1[i]);
  }
  return std::sqrt(norm);
}

double infNorm(const std::vector<double>& x) {
  double norm{};
  for (Int i = 0; i < x.size(); ++i) {
    norm = std::max(norm, std::fabs(x[i]));
  }
  return norm;
}

double infNorm(const std::vector<double>& x0, const std::vector<double>& x1) {
  double norm{};
  for (Int i = 0; i < x0.size(); ++i) {
    norm = std::max(norm, std::fabs(x0[i]));
  }
  for (Int i = 0; i < x1.size(); ++i) {
    norm = std::max(norm, std::fabs(x1[i]));
  }
  return norm;
}

double infNormDiff(const std::vector<double>& x, const std::vector<double>& y) {
  assert(x.size() == y.size());
  double inf_norm_diff = 0;
  for (Int i = 0; i < Int(x.size()); i++) {
    double diff = std::fabs(x[i] - y[i]);
    inf_norm_diff = std::max(diff, inf_norm_diff);
  }
  return inf_norm_diff;
}

bool isNanVector(const std::vector<double>& x) {
  for (Int i = 0; i < x.size(); ++i) {
    if (std::isnan(x[i])) return true;
  }
  return false;
}

bool isInfVector(const std::vector<double>& x) {
  for (Int i = 0; i < x.size(); ++i) {
    if (std::isinf(x[i])) return true;
  }
  return false;
}

}  // namespace hipo
