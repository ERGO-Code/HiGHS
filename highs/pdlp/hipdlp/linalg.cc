/*
 * @Author: Zhou Yanyu（周妍妤） 47125824+Yanyu000@users.noreply.github.com
 * @Date: 2025-07-09 14:54:26
 * @LastEditors: Zhou Yanyu（周妍妤） 47125824+Yanyu000@users.noreply.github.com
 * @LastEditTime: 2025-08-05 14:47:19
 * @FilePath: /cupdlp-CPP/src/linalg.cpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include "linalg.hpp"
#include "Highs.h"

namespace linalg {

    double project_box(double x, double l, double u){
        if (x < l) {
            return l;
        } else if (x > u) {
            return u;
        } else {
            return x;
        }
    }

    double project_non_negative(double x){
        return (x < 0.0) ? 0.0 : x;
    }

    void Ax(const HighsLp &lp, const std::vector<double> &x, std::vector<double> & result){
        for (HighsInt i = 0; i < lp.num_row_; ++i) {
            result[i] = 0.0;
        }
        for (HighsInt col = 0; col < lp.num_col_; ++col) {  // Loop over columns
            for (HighsInt k_ptr = lp.a_matrix_.start_[col]; k_ptr < lp.a_matrix_.start_[col + 1]; ++k_ptr) {
                HighsInt row = lp.a_matrix_.index_[k_ptr];  // row index
                double a_val = lp.a_matrix_.value_[k_ptr];
                result[row] += a_val * x[col];              // A[row][col] * x[col]
            }
        }
    }

    void ATy(const HighsLp &lp, const std::vector<double> &y, std::vector<double> & result){
        for (HighsInt i = 0; i < lp.num_col_; ++i) {
            result[i] = 0.0;
        }
        for (HighsInt col = 0; col < lp.num_col_; ++col) {  // Loop over columns
            for (HighsInt k_ptr = lp.a_matrix_.start_[col]; k_ptr < lp.a_matrix_.start_[col + 1]; ++k_ptr) {
                HighsInt row = lp.a_matrix_.index_[k_ptr];  // row index
                double a_val = lp.a_matrix_.value_[k_ptr];  // A[row][col]
                result[col] += a_val * y[row];              // A^T[col][row] * y[row]
            }
        }
    }

    double nrm2(const std::vector<double>& vec) {
        double sum_sq = 0.0;
        for (double val : vec) {
            sum_sq += val * val;
        }
        return std::sqrt(sum_sq);
    }

    // ADD THE IMPLEMENTATION FOR scale:
    void scale(std::vector<double>& vec, double factor) {
        for (size_t i = 0; i < vec.size(); ++i) {
            vec[i] *= factor;
        }
    }

    void normalize(std::vector<double>& vec){
        double norm = nrm2(vec);
        if (norm > 0.0) {
            scale(vec, 1.0 / norm);
        } else {
            // If the vector is zero, we can choose to leave it as is or set it to a default value
            // Here we choose to leave it unchanged
        }
    }

    double dot(const std::vector<double>& a, const std::vector<double>& b) {
        if (a.size() != b.size()) {
            throw std::invalid_argument("Vectors must be of the same size");
        }
        double result = 0.0;
        for (size_t i = 0; i < a.size(); ++i) {
            result += a[i] * b[i];
        }
        return result;
    }

    // Computes the L2 norm of the difference between two vectors (v1 - v2).
    double diffTwoNorm(const std::vector<double>& v1, const std::vector<double>& v2) {
        double norm_sq = 0.0;
        if (v1.size() != v2.size()) {
            // Handle error: vectors must have the same dimension.
            return -1.0; 
        }
        for (size_t i = 0; i < v1.size(); ++i) {
            double diff = v1[i] - v2[i];
            norm_sq += diff * diff;
        }
        return std::sqrt(norm_sq);
    }

    double vector_norm(const std::vector<double>& vec, double p) {
        if (std::isinf(p)) {
            // Infinity norm
            double max_val = 0.0;
            for (double val : vec) {
                max_val = std::max(max_val, std::abs(val));
            }
            return max_val;
        } else if (p == 2.0) {
            // L2 norm (use existing optimized function)
            return nrm2(vec);
        } else if (p == 1.0) {
            // L1 norm
            double sum = 0.0;
            for (double val : vec) {
                sum += std::abs(val);
            }
            return sum;
        } else {
            // General p-norm
            double sum = 0.0;
            for (double val : vec) {
                sum += std::pow(std::abs(val), p);
            }
            return std::pow(sum, 1.0 / p);
        }
    }
    
    // General vector norm (for raw array)
    double vector_norm(const double* values, size_t size, double p) {
        if (std::isinf(p)) {
            // Infinity norm
            double max_val = 0.0;
            for (size_t i = 0; i < size; ++i) {
                max_val = std::max(max_val, std::abs(values[i]));
            }
            return max_val;
        } else if (p == 2.0) {
            // L2 norm
            double sum_sq = 0.0;
            for (size_t i = 0; i < size; ++i) {
                sum_sq += values[i] * values[i];
            }
            return std::sqrt(sum_sq);
        } else if (p == 1.0) {
            // L1 norm
            double sum = 0.0;
            for (size_t i = 0; i < size; ++i) {
                sum += std::abs(values[i]);
            }
            return sum;
        } else {
            // General p-norm
            double sum = 0.0;
            for (size_t i = 0; i < size; ++i) {
                sum += std::pow(std::abs(values[i]), p);
            }
            return std::pow(sum, 1.0 / p);
        }
    }

    double compute_cost_norm(const HighsLp& lp, double p) {
        return vector_norm(lp.col_cost_, p);
    }
    
    // Compute norm of RHS vector (only for finite values)
    double compute_rhs_norm(const HighsLp& lp, double p) {
        std::vector<double> finite_rhs;
        finite_rhs.reserve(lp.num_row_);
        
        for (HighsInt i = 0; i < lp.num_row_; ++i) {
            if (lp.row_lower_[i] > -kHighsInf && lp.row_lower_[i] < kHighsInf) {
                finite_rhs.push_back(lp.row_lower_[i]);
            }
        }
        
        return finite_rhs.empty() ? 0.0 : vector_norm(finite_rhs, p);
    }
    
    // Compute column norms of the constraint matrix
    std::vector<double> compute_column_norms(const HighsLp& lp, double p) {
        std::vector<double> col_norms(lp.num_col_, 0.0);
        
        for (HighsInt col = 0; col < lp.num_col_; ++col) {
            HighsInt start = lp.a_matrix_.start_[col];
            HighsInt end = lp.a_matrix_.start_[col + 1];
            
            if (start < end) {
                col_norms[col] = vector_norm(&lp.a_matrix_.value_[start], end - start, p);
            }
        }
        
        return col_norms;
    }
    
    // Compute row norms of the constraint matrix
    std::vector<double> compute_row_norms(const HighsLp& lp, double p) {
        std::vector<double> row_norms(lp.num_row_, 0.0);
        
        if (std::isinf(p)) {
            // Infinity norm - find max absolute value in each row
            for (HighsInt col = 0; col < lp.num_col_; ++col) {
                for (HighsInt el = lp.a_matrix_.start_[col]; 
                     el < lp.a_matrix_.start_[col + 1]; ++el) {
                    HighsInt row = lp.a_matrix_.index_[el];
                    double abs_val = std::abs(lp.a_matrix_.value_[el]);
                    row_norms[row] = std::max(row_norms[row], abs_val);
                }
            }
        } else if (p == 2.0) {
            // L2 norm - sum of squares
            for (HighsInt col = 0; col < lp.num_col_; ++col) {
                for (HighsInt el = lp.a_matrix_.start_[col]; 
                     el < lp.a_matrix_.start_[col + 1]; ++el) {
                    HighsInt row = lp.a_matrix_.index_[el];
                    double val = lp.a_matrix_.value_[el];
                    row_norms[row] += val * val;
                }
            }
            for (HighsInt row = 0; row < lp.num_row_; ++row) {
                row_norms[row] = std::sqrt(row_norms[row]);
            }
        } else if (p == 1.0) {
            // L1 norm - sum of absolute values
            for (HighsInt col = 0; col < lp.num_col_; ++col) {
                for (HighsInt el = lp.a_matrix_.start_[col]; 
                     el < lp.a_matrix_.start_[col + 1]; ++el) {
                    HighsInt row = lp.a_matrix_.index_[el];
                    row_norms[row] += std::abs(lp.a_matrix_.value_[el]);
                }
            }
        } else {
            // General p-norm
            for (HighsInt col = 0; col < lp.num_col_; ++col) {
                for (HighsInt el = lp.a_matrix_.start_[col]; 
                     el < lp.a_matrix_.start_[col + 1]; ++el) {
                    HighsInt row = lp.a_matrix_.index_[el];
                    row_norms[row] += std::pow(std::abs(lp.a_matrix_.value_[el]), p);
                }
            }
            for (HighsInt row = 0; row < lp.num_row_; ++row) {
                if (row_norms[row] > 0.0) {
                    row_norms[row] = std::pow(row_norms[row], 1.0 / p);
                }
            }
        }
        
        return row_norms;
    }
}
