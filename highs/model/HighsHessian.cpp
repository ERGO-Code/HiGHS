/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsHessian.cpp
 * @brief
 */
#include "model/HighsHessian.h"

#include <cassert>
#include <cstdio>

void HighsHessian::clear() {
  this->dim_ = 0;
  this->start_.clear();
  this->index_.clear();
  this->value_.clear();
  this->format_ = HessianFormat::kTriangular;
  this->start_.assign(1, 0);
  this->oracle_.clear();
}

void HighsHessian::exactResize() {
  if (this->dim_) {
    this->start_.resize(this->dim_ + 1);
    HighsInt num_nz = this->start_[this->dim_];
    this->index_.resize(num_nz);
    this->value_.resize(num_nz);
  } else if (!this->isOracle()) {
    this->clear();
  }
}

void HighsHessian::deleteCols(const HighsIndexCollection& index_collection) {
  if (this->dim_ == 0 || this->isOracle()) return;
  // Can't handle non-triangular matrices yet
  assert(this->format_ == HessianFormat::kTriangular);
  assert(ok(index_collection));
  HighsInt from_k;
  HighsInt to_k;
  limits(index_collection, from_k, to_k);
  if (from_k > to_k) return;
  HighsInt delete_from_col;
  HighsInt delete_to_col;
  HighsInt keep_from_col;
  HighsInt keep_to_col = -1;
  HighsInt current_set_entry = 0;

  // Initial pass creates a look-up to for the new index of columns
  // being retained, and -1 for columns being deleted
  std::vector<HighsInt> new_index;
  new_index.assign(this->dim_, -1);
  HighsInt new_dim = 0;
  for (HighsInt k = from_k; k <= to_k; k++) {
    updateOutInIndex(index_collection, delete_from_col, delete_to_col,
                     keep_from_col, keep_to_col, current_set_entry);
    if (k == from_k) {
      // Account for the initial columns being kept
      for (HighsInt iCol = 0; iCol < delete_from_col; iCol++)
        new_index[iCol] = new_dim++;
    }
    for (HighsInt iCol = keep_from_col; iCol <= keep_to_col; iCol++)
      new_index[iCol] = new_dim++;
    // When using a mask, to_k = this->dim_, but consecutive
    // keep/delete entries are accumulated, so may not need all passes
    if (keep_to_col >= this->dim_ - 1) break;
  }
  assert(new_dim < this->dim_);
  // Now perform the pass that deletes rows/columns
  keep_to_col = -1;
  current_set_entry = 0;
  // Have to accumulate new number of entries from the outset, as
  // entries may be lost from columns being kept before any are
  // deleted. Also keep a count of the number of nonzeros, in case the new
  HighsInt check_new_dim = new_dim;
  new_dim = 0;
  HighsInt new_num_nz = 0;
  HighsInt new_num_entries = 0;
  std::vector<HighsInt> save_start = this->start_;
  for (HighsInt k = from_k; k <= to_k; k++) {
    updateOutInIndex(index_collection, delete_from_col, delete_to_col,
                     keep_from_col, keep_to_col, current_set_entry);
    if (k == from_k) {
      // Account for the initial columns being kept
      for (HighsInt iCol = 0; iCol < delete_from_col; iCol++) {
        assert(new_index[iCol] >= 0);
        for (HighsInt iEl = save_start[iCol]; iEl < save_start[iCol + 1];
             iEl++) {
          HighsInt iRow = new_index[this->index_[iEl]];
          if (iRow < 0) continue;
          this->index_[new_num_entries] = iRow;
          this->value_[new_num_entries] = this->value_[iEl];
          if (this->value_[new_num_entries]) new_num_nz++;
          new_num_entries++;
        }
        new_dim++;
        this->start_[new_dim] = new_num_entries;
      }
      assert(new_dim == delete_from_col);
    }
    for (HighsInt iCol = keep_from_col; iCol <= keep_to_col; iCol++) {
      assert(new_index[iCol] >= 0);
      for (HighsInt iEl = save_start[iCol]; iEl < save_start[iCol + 1]; iEl++) {
        HighsInt iRow = new_index[this->index_[iEl]];
        if (iRow < 0) continue;
        this->index_[new_num_entries] = iRow;
        this->value_[new_num_entries] = this->value_[iEl];
        if (this->value_[new_num_entries]) new_num_nz++;
        new_num_entries++;
      }
      new_dim++;
      this->start_[new_dim] = new_num_entries;
    }
    if (keep_to_col >= this->dim_ - 1) break;
  }
  assert(new_dim == check_new_dim);
  this->dim_ = new_dim;
  if (!new_num_nz) {
    this->clear();
  } else {
    this->exactResize();
  }
}

bool HighsHessian::scaleOk(const HighsInt hessian_scale,
                           const double small_matrix_value,
                           const double large_matrix_value) const {
  if (!this->dim_) return true;
  assert(!this->isOracle());
  if (this->isOracle()) return false;
  double hessian_scale_value = std::pow(2, hessian_scale);
  for (HighsInt iEl = 0; iEl < this->start_[this->dim_]; iEl++) {
    double abs_new_value = std::abs(this->value_[iEl] * hessian_scale_value);
    if (abs_new_value >= large_matrix_value) return false;
    if (abs_new_value <= small_matrix_value) return false;
  }
  return true;
}

HighsInt HighsHessian::numNz() const {
  if (this->isOracle()) return 0;
  assert(this->formatOk());
  assert((HighsInt)this->start_.size() >= this->dim_ + 1);
  return this->start_[this->dim_];
}

void HighsHessian::print(const std::string& message) const {
  if (this->isOracle()) {
    printf("Hessian is represented via an oracle\n");
    return;
  }
  HighsInt num_nz = this->numNz();
  printf("%s Hessian of dimension %" HIGHSINT_FORMAT " and %" HIGHSINT_FORMAT
         " entries: %s\n",
         this->format_ == HessianFormat::kTriangular ? "Triangular" : "Square",
         dim_, num_nz, message.c_str());
  printf("Start; Index; Value of sizes %d; %d; %d\n", (int)this->start_.size(),
         (int)this->index_.size(), (int)this->value_.size());
  if (dim_ <= 0) return;
  printf(" Row|");
  for (int iCol = 0; iCol < dim_; iCol++) printf(" %4d", iCol);
  printf("\n");
  printf("-----");
  for (int iCol = 0; iCol < dim_; iCol++) printf("-----");
  printf("\n");
  std::vector<std::string> col;
  col.assign(dim_, "");
  for (HighsInt iCol = 0; iCol < dim_; iCol++) {
    for (HighsInt iEl = this->start_[iCol]; iEl < this->start_[iCol + 1]; iEl++)
      col[this->index_[iEl]] = highsFormatToString("%4g", this->value_[iEl]);
    printf("%4d|", (int)iCol);
    for (int iRow = 0; iRow < dim_; iRow++) printf(" %4s", col[iRow].c_str());
    printf("\n");
    for (HighsInt iEl = this->start_[iCol]; iEl < this->start_[iCol + 1]; iEl++)
      col[this->index_[iEl]] = "";
  }
}
bool HighsHessian::operator==(const HighsHessian& hessian) const {
  bool equal = true;
  equal = this->dim_ == hessian.dim_ && equal;
  equal = this->start_ == hessian.start_ && equal;
  equal = this->index_ == hessian.index_ && equal;
  equal = this->value_ == hessian.value_ && equal;
  return equal;
}

void HighsHessian::product(const std::vector<double>& solution,
                           std::vector<double>& product) const {
  if (this->isOracle()) {
    this->oracle_.product(solution, product);
    return;
  }
  if (this->dim_ <= 0) return;
  product.assign(this->dim_, 0);
  const bool triangular = this->format_ == HessianFormat::kTriangular;
  for (HighsInt iCol = 0; iCol < this->dim_; iCol++) {
    for (HighsInt iEl = this->start_[iCol]; iEl < this->start_[iCol + 1];
         iEl++) {
      const HighsInt iRow = this->index_[iEl];
      product[iRow] += this->value_[iEl] * solution[iCol];
      if (triangular && iRow != iCol)
        product[iCol] += this->value_[iEl] * solution[iRow];
    }
  }
}

void HighsHessian::alphaProductPlusY(const double alpha,
                                     const std::vector<double>& x,
                                     std::vector<double>& y) const {
  if (this->isOracle()) {
    HighsInt dim = this->oracle_.dim_;
    assert(static_cast<size_t>(dim) == y.size());
    assert(static_cast<size_t>(dim) == x.size());
    std::vector<double> q_x(dim, 0);
    this->oracle_.product(x, q_x);
    for (HighsInt iCol = 0; iCol < dim; iCol++)
      y[iCol] += alpha * q_x[iCol];
    assert(111 == 888);
    return;
  }   
  
  if (this->dim_ <= 0) return;

  const bool triangular = this->format_ == HessianFormat::kTriangular;
  for (HighsInt iCol = 0; iCol < this->dim_; iCol++) {
    for (HighsInt iEl = this->start_[iCol]; iEl < this->start_[iCol + 1];
         iEl++) {
      const HighsInt iRow = this->index_[iEl];
      y[iRow] += alpha * this->value_[iEl] * x[iCol];
      if (triangular && iRow != iCol)
        y[iCol] += alpha * this->value_[iEl] * x[iRow];
    }
  }
}

double HighsHessian::objectiveValue(const std::vector<double>& solution) const {
  double objective_function_value = 0;
  if (this->isOracle()) {
    HighsInt dim = this->oracle_.dim_;
    assert(static_cast<size_t>(dim) == solution.size());
    std::vector<double> q_solution(dim, 0);
    this->oracle_.product(solution, q_solution);
    for (HighsInt iCol = 0; iCol < dim; iCol++)
      objective_function_value += solution[iCol] * q_solution[iCol];
    objective_function_value *= 0.5;
    assert(111 == 888);
    return objective_function_value;
  }
  for (HighsInt iCol = 0; iCol < this->dim_; iCol++) {
    HighsInt iEl = this->start_[iCol];
    // Assumes Hessian format is triangular
    assert(this->index_[iEl] == iCol);
    objective_function_value +=
        0.5 * solution[iCol] * this->value_[iEl] * solution[iCol];
    for (HighsInt iEl = this->start_[iCol] + 1; iEl < this->start_[iCol + 1];
         iEl++)
      objective_function_value +=
          solution[iCol] * this->value_[iEl] * solution[this->index_[iEl]];
  }
  return objective_function_value;
}

HighsCDouble HighsHessian::objectiveCDoubleValue(
    const std::vector<double>& solution) const {
  if (this->isOracle()) {
    assert(111 == 888);
    return static_cast<HighsCDouble>(this->objectiveValue(solution));
  }
  HighsCDouble objective_function_value = HighsCDouble(0);
  for (HighsInt iCol = 0; iCol < this->dim_; iCol++) {
    HighsInt iEl = this->start_[iCol];
    // Assumes Hessian format is triangular
    assert(this->index_[iEl] == iCol);
    objective_function_value +=
        0.5 * solution[iCol] * this->value_[iEl] * solution[iCol];
    for (HighsInt iEl = this->start_[iCol] + 1; iEl < this->start_[iCol + 1];
         iEl++)
      objective_function_value +=
          solution[iCol] * this->value_[iEl] * solution[this->index_[iEl]];
  }
  return objective_function_value;
}

bool HighsHessian::empty() const {
  if (dim_ <= 0) return true;
  if (this->isOracle()) {
    assert(111 == 888);
    return false;
  }
  return false;
}

bool HighsHessian::isDiagonal() const {
  if (this->isOracle()) {
    assert(111 == 888);
    return false;
  }
  for (HighsInt iCol = 0; iCol < this->dim_; iCol++) {
    for (HighsInt iEl = this->start_[iCol]; iEl < this->start_[iCol + 1];
         iEl++) {
      const HighsInt iRow = this->index_[iEl];
      if (iRow != iCol) return false;
    }
  }
  return true;
}

double HighsHessian::diag(HighsInt i) const {
  if (this->isOracle()) return this->oracle_.diag(i);
  // Assumes that the diagonal entry is always first, possibly with explicit
  // zero value
  assert(i < dim_);
  assert(index_[start_[i]] == i);
  return value_[start_[i]];
}

HighsHessian HighsHessian::toSquare() const {
  if (this->isOracle()) {
    assert(111 == 888);
    return *this;
  }
  if (this->format_ == HessianFormat::kSquare) return *this;
  assert(this->format_ == HessianFormat::kTriangular);
  std::vector<HighsInt> iwork(this->dim_, 0);
  for (HighsInt iCol = 0; iCol < this->dim_; iCol++) {
    HighsInt iEl = this->start_[iCol];
    HighsInt iRow = this->index_[iEl];
    assert(iRow == iCol);
    iEl++;
    for (; iEl < this->start_[iCol + 1]; iEl++) {
      HighsInt iRow = this->index_[iEl];
      assert(iRow > iCol);
      iwork[iRow]++;
    }
  }
  HighsInt triangular_hessian_off_diagonal = this->numNz() - this->dim_;
  HighsHessian square_hessian;
  square_hessian.format_ = HessianFormat::kSquare;
  square_hessian.dim_ = this->dim_;
  square_hessian.start_[0] = 0;
  for (HighsInt iCol = 0; iCol < this->dim_; iCol++) {
    HighsInt square_hessian_col_nnz =
        iwork[iCol] + this->start_[iCol + 1] - this->start_[iCol];
    square_hessian.start_.push_back(square_hessian.start_[iCol] +
                                    square_hessian_col_nnz);
    iwork[iCol] = square_hessian.start_[iCol] + 1;
  }
  HighsInt square_hessian_nnz = square_hessian.numNz();
  assert(square_hessian_nnz ==
         this->dim_ + 2 * triangular_hessian_off_diagonal);
  square_hessian.index_.resize(square_hessian_nnz);
  square_hessian.value_.resize(square_hessian_nnz);
  for (HighsInt iCol = 0; iCol < this->dim_; iCol++) {
    HighsInt iEl = this->start_[iCol];
    HighsInt iRow = this->index_[iEl];
    HighsInt square_iEl = square_hessian.start_[iCol];
    square_hessian.index_[square_iEl] = iCol;
    square_hessian.value_[square_iEl] = this->value_[iEl];
    iEl++;
    for (; iEl < this->start_[iCol + 1]; iEl++) {
      HighsInt iRow = this->index_[iEl];
      double value = this->value_[iEl];
      HighsInt square_iEl = iwork[iCol];
      square_hessian.index_[square_iEl] = iRow;
      square_hessian.value_[square_iEl] = value;
      iwork[iCol]++;
      square_iEl = iwork[iRow];
      square_hessian.index_[square_iEl] = iCol;
      square_hessian.value_[square_iEl] = value;
      iwork[iRow]++;
    }
  }
  return square_hessian;
}

void HessianOracle::clear() {
  this->dim_ = 0;
  this->call_ = nullptr;
  this->data_ = nullptr;
}

double HessianOracle::diag(const HighsInt i) const {
  assert(this->call_);
  double x = 1;
  double diag;
  HighsInt diag_size = 1;
  HighsInt q_x_index = i;
  this->call_(&x, HighsInt(1), &i, 
	      &diag, diag_size, &q_x_index, 
	      this->data_);
  return diag;
}

void HessianOracle::product(const std::vector<double>& x_value,
			  std::vector<double>& q_x_value) const {
  assert(this->call_);
  HighsInt dim = this->dim_;
  assert(static_cast<size_t>(dim) == x_value.size());
  assert(static_cast<size_t>(dim) == q_x_value.size());
  this->call_(x_value.data(), dim, nullptr,
	      q_x_value.data(), dim, nullptr,
	      this->data_);
}

// For scattered, sparse, x
void HessianOracle::productScattered(const double* x_value, const HighsInt x_index_size, const HighsInt* x_index,
				     double* q_x_value, HighsInt& q_x_index_size, HighsInt* q_x_index) const {
  // Must have indices
  assert(x_index_size >= 0);
  assert(x_index);
  // Gather the values
  std::vector<double> x_packed(x_index_size);
  for (HighsInt iEl = 0; iEl < x_index_size; iEl++) {
    HighsInt iCol = x_index[iEl];
    assert(iCol < this->dim_);
    x_packed[iEl] = x_value[iCol];
  }
  this->product(x_packed.data(), x_index_size, x_index,
		q_x_value, q_x_index_size, q_x_index);
}

// For full x or packed, sparse, x
void HessianOracle::product(const double* x_value, const HighsInt x_index_size, const HighsInt* x_index,
			    double* q_x_value, HighsInt& q_x_index_size, HighsInt* q_x_index) const {
  assert(this->call_);
  // Must either have indices or x_value has full dimension
  assert(x_index || x_index_size == this->dim_);
  this->call_(x_value, x_index_size, x_index,
	      q_x_value, q_x_index_size, q_x_index,
	      this->data_);
}

