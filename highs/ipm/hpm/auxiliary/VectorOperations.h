#ifndef HIGHSPM_VECTOR_OPERATIONS_H
#define HIGHSPM_VECTOR_OPERATIONS_H

#include <vector>

#include "ipm/hpm/auxiliary/IntConfig.h"

namespace highspm {

// =======================================================================
// COMPONENT-WISE VECTOR OPERATIONS
// =======================================================================

// alpha * v1[i] + beta * v2[i]
void vectorAdd(std::vector<double>& v1, const std::vector<double>& v2,
               double beta = 1.0, double alpha = 1.0);

// v1[i] + alpha
void vectorAdd(std::vector<double>& v1, const double alpha);

// alpha * v1[i] * v2[i] + beta
void vectorMultiply(std::vector<double>& v1, const std::vector<double>& v2,
                    double beta = 1.0, double alpha = 0.0);

// alpha * v1[i] + beta * v2[i] * v3[i]
void vectorAddMult(std::vector<double>& v1, const std::vector<double>& v2,
                   const std::vector<double>& v3, double beta = 1.0,
                   double alpha = 1.0);

// v1[i] / v2[i]
void vectorDivide(std::vector<double>& v1, const std::vector<double>& v2);

// v1[i] * alpha
void vectorScale(std::vector<double>& v1, double alpha);

// =======================================================================

// scalar product
double dotProd(const std::vector<double>& v1, const std::vector<double>& v2);

// Euclidean norm
double norm2(const std::vector<double>& x);
double norm2(const std::vector<double>& x0, const std::vector<double>& x1);

double infNorm(const std::vector<double>& x);
double infNorm(const std::vector<double>& x0, const std::vector<double>& x1);

// Infinity norm of the difference of two vectors
double infNormDiff(const std::vector<double>& x, const std::vector<double>& y);

// check for NaN
bool isNanVector(const std::vector<double>& x);

// check for Inf
bool isInfVector(const std::vector<double>& x);

}  // namespace highspm

#endif