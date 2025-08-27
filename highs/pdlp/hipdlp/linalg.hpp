/*
 * @Author: Zhou Yanyu（周妍妤） 47125824+Yanyu000@users.noreply.github.com
 * @Date: 2025-07-14 12:06:25
 * @LastEditors: Zhou Yanyu（周妍妤） 47125824+Yanyu000@users.noreply.github.com
 * @LastEditTime: 2025-08-05 14:47:53
 * @FilePath: /cupdlp-CPP/include/linalg.hpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#ifndef LINALG_HPP
#define LINALG_HPP

#include <vector>
#include "Highs.h"

namespace linalg {
    double project_box(double x, double l, double u);  
    double project_non_negative(double y);

    // Function to compute A*x for a given HighsLp and vector x
    void Ax(const HighsLp &lp, const std::vector<double> &x, std::vector<double> &result);

    // Function to compute A^T*y for a given HighsLp and vector y
    void ATy(const HighsLp &lp, const std::vector<double> &y, std::vector<double> &result);

    double nrm2(const std::vector<double>& vec);
    void scale(std::vector<double>& vec, double factor);

    void normalize(std::vector<double>& vec);
    
    double dot(const std::vector<double>& a, const std::vector<double>& b);

    double diffTwoNorm(const std::vector<double>& v1, const std::vector<double>& v2);

    // General norm functions
    double vector_norm(const std::vector<double>& vec, double p = 2.0);
    double vector_norm(const double* values, size_t size, double p = 2.0);
    
    // LP-specific norm calculations
    double compute_cost_norm(const HighsLp& lp, double p = 2.0);
    double compute_rhs_norm(const HighsLp& lp, double p = 2.0);
    
    // Matrix column/row norm calculations
    std::vector<double> compute_column_norms(const HighsLp& lp, double p = std::numeric_limits<double>::infinity());
    std::vector<double> compute_row_norms(const HighsLp& lp, double p = std::numeric_limits<double>::infinity());
} // namespace linalg

#endif // LINALG_HPP