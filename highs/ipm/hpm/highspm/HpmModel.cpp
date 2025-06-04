#include "HpmModel.h"

#include "HpmConst.h"
#include "HpmStatus.h"
#include "ipm/hpm/auxiliary/HpmLog.h"

namespace highspm {

Int HpmModel::init(const Int num_var, const Int num_con, const double* obj,
                   const double* rhs, const double* lower, const double* upper,
                   const Int* A_ptr, const Int* A_rows, const double* A_vals,
                   const char* constraints, double offset) {
  // copy the input into the model

  if (checkData(num_var, num_con, obj, rhs, lower, upper, A_ptr, A_rows, A_vals,
                constraints))
    return kIpmStatusBadModel;

  n_orig_ = num_var;
  m_orig_ = num_con;
  c_orig_ = obj;
  b_orig_ = rhs;
  lower_orig_ = lower;
  upper_orig_ = upper;
  A_ptr_orig_ = A_ptr;
  A_rows_orig_ = A_rows;
  A_vals_orig_ = A_vals;
  constraints_orig_ = constraints;
  offset_ = offset;

  n_ = num_var;
  m_ = num_con;
  c_ = std::vector<double>(obj, obj + n_);
  b_ = std::vector<double>(rhs, rhs + m_);
  lower_ = std::vector<double>(lower, lower + n_);
  upper_ = std::vector<double>(upper, upper + n_);

  Int Annz = A_ptr[n_];
  A_.num_col_ = n_;
  A_.num_row_ = m_;
  A_.start_ = std::vector<Int>(A_ptr, A_ptr + n_ + 1);
  A_.index_ = std::vector<Int>(A_rows, A_rows + Annz);
  A_.value_ = std::vector<double>(A_vals, A_vals + Annz);

  constraints_ = std::vector<char>(constraints, constraints + m_);

  preprocess();
  scale();
  reformulate();
  denseColumns();

  ready_ = true;

  return 0;
}

Int HpmModel::checkData(const Int num_var, const Int num_con, const double* obj,
                        const double* rhs, const double* lower,
                        const double* upper, const Int* A_ptr,
                        const Int* A_rows, const double* A_vals,
                        const char* constraints) const {
  // Check if model provided by the user is ok.
  // Return kIpmStatusBadModel if something is wrong.

  // Pointers are valid
  if (!obj || !rhs || !lower || !upper || !A_ptr || !A_rows || !A_vals ||
      !constraints)
    return kIpmStatusBadModel;

  // Dimensions are valid
  if (num_var <= 0 || num_con < 0) return kIpmStatusBadModel;

  // Vectors are valid
  for (Int i = 0; i < num_var; ++i)
    if (!std::isfinite(obj[i])) return kIpmStatusBadModel;
  for (Int i = 0; i < num_con; ++i)
    if (!std::isfinite(rhs[i])) return kIpmStatusBadModel;
  for (Int i = 0; i < num_var; ++i) {
    if (!std::isfinite(lower[i]) && lower[i] != -INFINITY)
      return kIpmStatusBadModel;
    if (!std::isfinite(upper[i]) && upper[i] != INFINITY)
      return kIpmStatusBadModel;
    if (lower[i] > upper[i]) return kIpmStatusBadModel;
  }
  for (Int i = 0; i < num_con; ++i)
    if (constraints[i] != '<' && constraints[i] != '=' && constraints[i] != '>')
      return kIpmStatusBadModel;

  // Matrix is valid
  for (Int i = 0; i < A_ptr[num_var]; ++i)
    if (!std::isfinite(A_vals[i])) return kIpmStatusBadModel;

  return 0;
}

void HpmModel::preprocess() {
  // Perform some basic preprocessing, in case the problem is run without
  // presolve

  // ==========================================
  // Remove empty rows
  // ==========================================

  // find empty rows
  std::vector<Int> entries_per_row(m_, 0);
  for (Int col = 0; col < n_; ++col) {
    for (Int el = A_.start_[col]; el < A_.start_[col + 1]; ++el) {
      const Int row = A_.index_[el];
      ++entries_per_row[row];
    }
  }

  rows_shift_.assign(m_, 0);
  empty_rows_ = 0;
  for (Int i = 0; i < m_; ++i) {
    if (entries_per_row[i] == 0) {
      // count number of empty rows
      ++empty_rows_;

      // count how many empty rows there are before a given row
      for (Int j = i + 1; j < m_; ++j) ++rows_shift_[j];

      rows_shift_[i] = -1;
    }
  }

  if (empty_rows_ > 0) {
    // shift each row index by the number of empty rows before it
    for (Int col = 0; col < n_; ++col) {
      for (Int el = A_.start_[col]; el < A_.start_[col + 1]; ++el) {
        const Int row = A_.index_[el];
        A_.index_[el] -= rows_shift_[row];
      }
    }
    A_.num_row_ -= empty_rows_;

    // shift entries in b and constraints
    for (Int i = 0; i < m_; ++i) {
      // ignore entries to be removed
      if (rows_shift_[i] == -1) continue;

      Int shifted_pos = i - rows_shift_[i];
      b_[shifted_pos] = b_[i];
      constraints_[shifted_pos] = constraints_[i];
    }
    b_.resize(A_.num_row_);
    constraints_.resize(A_.num_row_);

    m_ = A_.num_row_;
  }
}

void HpmModel::postprocess(std::vector<double>& slack,
                           std::vector<double>& y) const {
  // Add Lagrange multiplier for empty rows that were removed
  // Add slack for constraints that were removed

  if (empty_rows_ == 0) return;

  std::vector<double> new_y(rows_shift_.size(), 0.0);
  std::vector<double> new_slack(rows_shift_.size(), 0.0);

  // position to read from y and slack
  Int pos{};

  for (Int i = 0; i < rows_shift_.size(); ++i) {
    // ignore shift of empty rows, they will receive a value of 0
    if (rows_shift_[i] == -1) continue;

    // re-align value of y and slack, considering empty rows
    new_y[pos + rows_shift_[i]] = y[pos];
    new_slack[pos + rows_shift_[i]] = slack[pos];
    ++pos;
  }

  y = std::move(new_y);
  slack = std::move(new_slack);
}

void HpmModel::reformulate() {
  // put the model into correct formulation

  Int Annz = A_.numNz();

  for (Int i = 0; i < m_; ++i) {
    if (constraints_[i] != '=') {
      // inequality constraint, add slack variable

      ++n_;

      // lower/upper bound for new slack
      if (constraints_[i] == '>') {
        lower_.push_back(-kHighsInf);
        upper_.push_back(0.0);
      } else {
        lower_.push_back(0.0);
        upper_.push_back(kHighsInf);
      }

      // cost for new slack
      c_.push_back(0.0);

      // add column of identity to A_
      std::vector<Int> temp_ind{i};
      std::vector<double> temp_val{1.0};
      A_.addVec(1, temp_ind.data(), temp_val.data());

      // set scaling to 1
      if (scaled()) colscale_.push_back(1.0);
    }
  }
}

void HpmModel::print() const {
  std::stringstream log_stream;

  log_stream << textline("Rows:") << sci(m_, 0, 1) << '\n';
  log_stream << textline("Cols:") << sci(n_, 0, 1) << '\n';
  log_stream << textline("Nnz:") << sci(A_.numNz(), 0, 1) << '\n';
  if (num_dense_cols_ > 0)
    log_stream << textline("Dense cols:") << integer(num_dense_cols_, 0)
               << '\n';
  if (empty_rows_ > 0)
    log_stream << "Removed " << empty_rows_ << " empty rows\n";

  // compute max and min entry of A in absolute value
  double Amin = kHighsInf;
  double Amax = 0.0;
  for (Int col = 0; col < A_.num_col_; ++col) {
    for (Int el = A_.start_[col]; el < A_.start_[col + 1]; ++el) {
      double val = std::abs(A_.value_[el]);
      if (val != 0.0) {
        Amin = std::min(Amin, val);
        Amax = std::max(Amax, val);
      }
    }
  }
  if (std::isinf(Amin)) Amin = 0.0;

  // compute max and min entry of c
  double cmin = kHighsInf;
  double cmax = 0.0;
  for (Int i = 0; i < n_; ++i) {
    if (c_[i] != 0.0) {
      cmin = std::min(cmin, std::abs(c_[i]));
      cmax = std::max(cmax, std::abs(c_[i]));
    }
  }
  if (std::isinf(cmin)) cmin = 0.0;

  // compute max and min entry of b
  double bmin = kHighsInf;
  double bmax = 0.0;
  for (Int i = 0; i < m_; ++i) {
    if (b_[i] != 0.0) {
      bmin = std::min(bmin, std::abs(b_[i]));
      bmax = std::max(bmax, std::abs(b_[i]));
    }
  }
  if (std::isinf(bmin)) bmin = 0.0;

  // compute max and min for bounds
  double boundmin = kHighsInf;
  double boundmax = 0.0;
  for (Int i = 0; i < n_; ++i) {
    if (lower_[i] != 0.0 && std::isfinite(lower_[i])) {
      boundmin = std::min(boundmin, std::abs(lower_[i]));
      boundmax = std::max(boundmax, std::abs(lower_[i]));
    }
    if (upper_[i] != 0.0 && std::isfinite(upper_[i])) {
      boundmin = std::min(boundmin, std::abs(upper_[i]));
      boundmax = std::max(boundmax, std::abs(upper_[i]));
    }
  }
  if (std::isinf(boundmin)) boundmin = 0.0;

  // compute max and min scaling
  double scalemin = kHighsInf;
  double scalemax = 0.0;
  if (scaled()) {
    for (Int i = 0; i < n_; ++i) {
      scalemin = std::min(scalemin, colscale_[i]);
      scalemax = std::max(scalemax, colscale_[i]);
    }
    for (Int i = 0; i < m_; ++i) {
      scalemin = std::min(scalemin, rowscale_[i]);
      scalemax = std::max(scalemax, rowscale_[i]);
    }
  }
  if (std::isinf(scalemin)) scalemin = 0.0;

  // print ranges
  log_stream << textline("Range of A:") << "[" << sci(Amin, 5, 1) << ", "
             << sci(Amax, 5, 1) << "], ratio ";
  if (Amin != 0.0)
    log_stream << sci(Amax / Amin, 0, 1) << '\n';
  else
    log_stream << "-\n";

  log_stream << textline("Range of b:") << "[" << sci(bmin, 5, 1) << ", "
             << sci(bmax, 5, 1) << "], ratio ";
  if (bmin != 0.0)
    log_stream << sci(bmax / bmin, 0, 1) << '\n';
  else
    log_stream << "-\n";

  log_stream << textline("Range of c:") << "[" << sci(cmin, 5, 1) << ", "
             << sci(cmax, 5, 1) << "], ratio ";
  if (cmin != 0.0)
    log_stream << sci(cmax / cmin, 0, 1) << '\n';
  else
    log_stream << "-\n";

  log_stream << textline("Range of bounds:") << "[" << sci(boundmin, 5, 1)
             << ", " << sci(boundmax, 5, 1) << "], ratio ";
  if (boundmin != 0.0)
    log_stream << sci(boundmax / boundmin, 0, 1) << '\n';
  else
    log_stream << "-\n";

  log_stream << textline("Scaling coefficients:") << "[" << sci(scalemin, 5, 1)
             << ", " << sci(scalemax, 5, 1) << "], ratio ";
  if (scalemin != 0.0)
    log_stream << sci(scalemax / scalemin, 0, 1) << '\n';
  else
    log_stream << "-\n";

  Log::print(log_stream);
}

void HpmModel::scale() {
  // Apply Curtis-Reid scaling and scale the problem accordingly

  // check if scaling is needed
  bool need_scaling = false;
  for (Int col = 0; col < n_; ++col) {
    for (Int el = A_.start_[col]; el < A_.start_[col + 1]; ++el) {
      if (std::abs(A_.value_[el]) != 1.0) {
        need_scaling = true;
        break;
      }
    }
  }

  if (!need_scaling) return;

  // *********************************************************************
  // Compute scaling
  // *********************************************************************
  // Transformation:
  // A -> R * A * C
  // b -> R * b
  // c -> C * c
  // x -> C^-1 * x
  // y -> R^-1 * y
  // z -> C * z
  // where R is row scaling, C is col scaling.

  // Compute exponents for CR scaling of matrix A
  std::vector<Int> colexp(n_);
  std::vector<Int> rowexp(m_);
  CurtisReidScaling(A_.start_, A_.index_, A_.value_, rowexp, colexp);

  // Compute scaling from exponents
  colscale_.resize(n_);
  rowscale_.resize(m_);
  for (Int i = 0; i < n_; ++i) colscale_[i] = std::ldexp(1.0, colexp[i]);
  for (Int i = 0; i < m_; ++i) rowscale_[i] = std::ldexp(1.0, rowexp[i]);

  // *********************************************************************
  // Apply scaling
  // *********************************************************************

  // Column has been scaled up by colscale_[col], so cost is scaled up and
  // bounds are scaled down
  for (Int col = 0; col < n_; ++col) {
    c_[col] *= colscale_[col];
    lower_[col] /= colscale_[col];
    upper_[col] /= colscale_[col];
  }

  // Row has been scaled up by rowscale_[row], so b is scaled up
  for (Int row = 0; row < m_; ++row) b_[row] *= rowscale_[row];

  // Each entry of the matrix is scaled by the corresponding row and col factor
  for (Int col = 0; col < n_; ++col) {
    for (Int el = A_.start_[col]; el < A_.start_[col + 1]; ++el) {
      Int row = A_.index_[el];
      A_.value_[el] *= rowscale_[row];
      A_.value_[el] *= colscale_[col];
    }
  }
}

void HpmModel::unscale(std::vector<double>& x, std::vector<double>& xl,
                       std::vector<double>& xu, std::vector<double>& slack,
                       std::vector<double>& y, std::vector<double>& zl,
                       std::vector<double>& zu) const {
  // Undo the scaling with internal format

  if (scaled()) {
    for (Int i = 0; i < n_orig_; ++i) {
      x[i] *= colscale_[i];
      xl[i] *= colscale_[i];
      xu[i] *= colscale_[i];
      zl[i] /= colscale_[i];
      zu[i] /= colscale_[i];
    }
    for (Int i = 0; i < m_; ++i) {
      y[i] *= rowscale_[i];
      slack[i] /= rowscale_[i];
    }
  }

  // set variables that were ignored
  for (Int i = 0; i < n_orig_; ++i) {
    if (!hasLb(i)) {
      xl[i] = kHighsInf;
      zl[i] = 0.0;
    }
    if (!hasUb(i)) {
      xu[i] = kHighsInf;
      zu[i] = 0.0;
    }
  }
}

void HpmModel::unscale(std::vector<double>& x, std::vector<double>& slack,
                       std::vector<double>& y, std::vector<double>& z) const {
  // Undo the scaling with format for crossover

  if (scaled()) {
    for (Int i = 0; i < n_orig_; ++i) {
      x[i] *= colscale_[i];
      z[i] /= colscale_[i];
    }
    for (Int i = 0; i < m_; ++i) {
      y[i] *= rowscale_[i];
      slack[i] /= rowscale_[i];
    }
  }
}

void HpmModel::denseColumns() {
  // Compute the maximum density of any column of A and count the number of
  // dense columns.

  max_col_density_ = 0.0;
  num_dense_cols_ = 0;
  for (Int col = 0; col < n_; ++col) {
    Int col_nz = A_.start_[col + 1] - A_.start_[col];
    double col_density = (double)col_nz / m_;
    max_col_density_ = std::max(max_col_density_, col_density);
    if (A_.num_row_ > kMinRowsForDensity && col_density > kDenseColThresh)
      ++num_dense_cols_;
  }
}

double HpmModel::normScaledRhs() const {
  double norm_rhs = infNorm(b_);
  for (double d : lower_)
    if (std::isfinite(d)) norm_rhs = std::max(norm_rhs, std::abs(d));
  for (double d : upper_)
    if (std::isfinite(d)) norm_rhs = std::max(norm_rhs, std::abs(d));
  return norm_rhs;
}

double HpmModel::normScaledObj() const { return infNorm(c_); }

double HpmModel::normUnscaledObj() const {
  double norm_obj = 0.0;
  for (Int i = 0; i < n_; ++i) {
    double val = std::abs(c_[i]);
    if (scaled()) val /= colscale_[i];
    norm_obj = std::max(norm_obj, val);
  }
  return norm_obj;
}

double HpmModel::normUnscaledRhs() const {
  double norm_rhs = 0.0;
  for (Int i = 0; i < m_; ++i) {
    double val = std::abs(b_[i]);
    if (scaled()) val /= rowscale_[i];
    norm_rhs = std::max(norm_rhs, val);
  }
  for (Int i = 0; i < n_; ++i) {
    if (std::isfinite(lower_[i])) {
      double val = std::abs(lower_[i]);
      if (scaled()) val *= colscale_[i];
      norm_rhs = std::max(norm_rhs, val);
    }
    if (std::isfinite(upper_[i])) {
      double val = std::abs(upper_[i]);
      if (scaled()) val *= colscale_[i];
      norm_rhs = std::max(norm_rhs, val);
    }
  }
  return norm_rhs;
}

Int HpmModel::loadIntoIpx(ipx::LpSolver& lps) const {
  Int load_status = lps.LoadModel(
      n_orig_, offset_, c_orig_, lower_orig_, upper_orig_, m_orig_, A_ptr_orig_,
      A_rows_orig_, A_vals_orig_, b_orig_, constraints_orig_);

  return load_status;
}

void HpmModel::multWithoutSlack(double alpha, const std::vector<double>& x,
                                std::vector<double>& y, bool trans) const {
  assert(x.size() == trans ? m_ : n_orig_);
  assert(y.size() == trans ? n_orig_ : m_);

  if (trans) {
    for (Int col = 0; col < n_orig_; ++col) {
      for (Int el = A_.start_[col]; el < A_.start_[col + 1]; ++el) {
        y[col] += alpha * A_.value_[el] * x[A_.index_[el]];
      }
    }
  } else {
    for (Int col = 0; col < n_orig_; ++col) {
      for (Int el = A_.start_[col]; el < A_.start_[col + 1]; ++el) {
        y[A_.index_[el]] += alpha * A_.value_[el] * x[col];
      }
    }
  }
}

}  // namespace highspm