#ifndef HIPO_MODEL_H
#define HIPO_MODEL_H

#include <limits>
#include <string>
#include <vector>

#include "CurtisReidScaling.h"
#include "LogHighs.h"
#include "ipm/hipo/auxiliary/IntConfig.h"
#include "ipm/ipx/lp_solver.h"
#include "lp_data/HighsLp.h"
#include "util/HighsSparseMatrix.h"

namespace hipo {

// Optimization problem:
// min  c^T * x
// s.t. A * x = b
//      x - xl = lower
//      x + xu = upper
//      xl, xu >= 0
//
// A is of size num_con x num_var, stored in CSC format using ptr, rows, vals.
//
// See Schork, Gondzio "Implementation of an interior point method with basis
// preconditioning", Math. Prog. Comput. 12, 2020
//

class Model {
 private:
  // data of original problem
  Int n_orig_{};
  Int m_orig_{};
  const HighsLp* lp_orig_ = nullptr;
  double offset_;

  // data of reformulated problem
  Int n_{};
  Int m_{};
  std::vector<double> c_{};
  std::vector<double> b_{};
  std::vector<double> lower_{};
  std::vector<double> upper_{};
  HighsSparseMatrix A_{};
  std::vector<char> constraints_{};
  Int num_dense_cols_{};
  double max_col_density_{};

  double norm_unscaled_rhs_, norm_scaled_rhs_, norm_unscaled_obj_,
      norm_scaled_obj_;

  bool ready_ = false;

  // coefficients for scaling
  std::vector<double> colscale_{};
  std::vector<double> rowscale_{};
  Int CG_iter_scaling_{};

  // information about empty rows, for postprocessing
  std::vector<Int> rows_shift_{};
  Int empty_rows_{};

  // norms of rows and cols of A
  std::vector<double> one_norm_cols_, one_norm_rows_, inf_norm_cols_,
      inf_norm_rows_;

  void reformulate();
  void scale();
  void preprocess();
  void denseColumns();
  Int checkData(const Int num_var, const Int num_con,
                const std::vector<double>& obj, const std::vector<double>& rhs,
                const std::vector<double>& lower,
                const std::vector<double>& upper, const std::vector<Int>& A_ptr,
                const std::vector<Int>& A_rows,
                const std::vector<double>& A_vals,
                const std::vector<char>& constraints) const;
  void computeNorms();

 public:
  // Initialise the model
  Int init(const HighsLp& lp);

  // Print information of model
  void print(const LogHighs& log) const;

  void postprocess(std::vector<double>& slack, std::vector<double>& y) const;

  // Unscale a given solution
  void unscale(std::vector<double>& x, std::vector<double>& xl,
               std::vector<double>& xu, std::vector<double>& slack,
               std::vector<double>& y, std::vector<double>& zl,
               std::vector<double>& zu) const;
  void unscale(std::vector<double>& x, std::vector<double>& slack,
               std::vector<double>& y, std::vector<double>& z) const;

  double normScaledRhs() const { return norm_scaled_rhs_; }
  double normScaledObj() const { return norm_scaled_obj_; }
  double normUnscaledObj() const { return norm_unscaled_obj_; }
  double normUnscaledRhs() const { return norm_unscaled_rhs_; }

  // multiply by A or A^T without slacks
  void multWithoutSlack(double alpha, const std::vector<double>& x,
                        std::vector<double>& y, bool trans = false) const;

  // Check if variable has finite lower/upper bound
  bool hasLb(Int j) const { return std::isfinite(lower_[j]); }
  bool hasUb(Int j) const { return std::isfinite(upper_[j]); }

  Int m() const { return m_; }
  Int n() const { return n_; }
  Int n_orig() const { return n_orig_; }
  Int m_orig() const { return m_orig_; }
  const HighsSparseMatrix& A() const { return A_; }
  const std::vector<double>& b() const { return b_; }
  const std::vector<double>& c() const { return c_; }
  double lb(Int i) const { return lower_[i]; }
  double ub(Int i) const { return upper_[i]; }
  char constraint(Int i) const { return constraints_[i]; }
  double colScale(Int i) const { return colscale_[i]; }
  double rowScale(Int i) const { return rowscale_[i]; }
  bool ready() const { return ready_; }
  bool scaled() const { return colscale_.size() > 0; }
  double offset() const { return offset_; }
  double maxColDensity() const { return max_col_density_; }
  Int numDenseCols() const { return num_dense_cols_; }
  double oneNormRows(Int i) const { return one_norm_rows_[i]; }
  double oneNormCols(Int i) const { return one_norm_cols_[i]; }
  double infNormRows(Int i) const { return inf_norm_rows_[i]; }
  double infNormCols(Int i) const { return inf_norm_cols_[i]; }

  Int loadIntoIpx(ipx::LpSolver& lps) const;
};

}  // namespace hipo

#endif
