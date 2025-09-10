/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/hipdlp/linalg.hpp
 * @brief
 */
#ifndef LINALG_HPP
#define LINALG_HPP

#include <vector>

#include "Highs.h"

namespace linalg {
double project_box(double x, double l, double u);
double project_non_negative(double y);

// Function to compute A*x for a given HighsLp and vector x
void Ax(const HighsLp& lp, const std::vector<double>& x,
        std::vector<double>& result);

// Function to compute A^T*y for a given HighsLp and vector y
void ATy(const HighsLp& lp, const std::vector<double>& y,
         std::vector<double>& result);

double nrm2(const std::vector<double>& vec);
void scale(std::vector<double>& vec, double factor);

void normalize(std::vector<double>& vec);

double dot(const std::vector<double>& a, const std::vector<double>& b);

double diffTwoNorm(const std::vector<double>& v1,
                   const std::vector<double>& v2);

// General norm functions
double vector_norm(const std::vector<double>& vec, double p = 2.0);
double vector_norm(const double* values, size_t size, double p = 2.0);

// LP-specific norm calculations
double compute_cost_norm(const HighsLp& lp, double p = 2.0);
double compute_rhs_norm(const HighsLp& lp, double p = 2.0);

// Matrix column/row norm calculations
std::vector<double> compute_column_norms(
    const HighsLp& lp, double p = std::numeric_limits<double>::infinity());
std::vector<double> compute_row_norms(
    const HighsLp& lp, double p = std::numeric_limits<double>::infinity());
}  // namespace linalg

#endif  // LINALG_HPP
