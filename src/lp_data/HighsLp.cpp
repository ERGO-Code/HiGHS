/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2024 by Julian Hall, Ivet Galabova,    */
/*    Leona Gottwald and Michael Feldmeier                               */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsLp.cpp
 * @brief
 */
#include "lp_data/HighsLp.h"

#include <cassert>

#include "lp_data/HighsLpUtils.h"
#include "util/HighsMatrixUtils.h"

bool HighsLp::isMip() const {
  HighsInt integrality_size = this->integrality_.size();
  if (integrality_size) {
    assert(integrality_size == this->num_col_);
    for (HighsInt iCol = 0; iCol < this->num_col_; iCol++)
      if (this->integrality_[iCol] != HighsVarType::kContinuous) return true;
  }
  return false;
}

bool HighsLp::hasInfiniteCost(const double infinite_cost) const {
  for (HighsInt iCol = 0; iCol < this->num_col_; iCol++) {
    if (this->col_cost_[iCol] >= infinite_cost) return true;
    if (this->col_cost_[iCol] <= -infinite_cost) return true;
  }
  return false;
}

bool HighsLp::hasSemiVariables() const {
  HighsInt integrality_size = this->integrality_.size();
  if (integrality_size) {
    assert(integrality_size == this->num_col_);
    for (HighsInt iCol = 0; iCol < this->num_col_; iCol++)
      if (this->integrality_[iCol] == HighsVarType::kSemiContinuous ||
          this->integrality_[iCol] == HighsVarType::kSemiInteger)
        return true;
  }
  return false;
}

bool HighsLp::hasMods() const {
  return this->mods_.save_non_semi_variable_index.size() > 0 ||
         this->mods_.save_inconsistent_semi_variable_index.size() > 0 ||
         this->mods_.save_relaxed_semi_variable_lower_bound_index.size() > 0 ||
         this->mods_.save_tightened_semi_variable_upper_bound_index.size() >
             0 ||
         this->mods_.save_inf_cost_variable_index.size() > 0;
}

bool HighsLp::needsMods(const double infinite_cost) const {
  assert(this->has_infinite_cost_ == this->hasInfiniteCost(infinite_cost));
  return this->has_infinite_cost_ || this->hasSemiVariables();
}

bool HighsLp::operator==(const HighsLp& lp) const {
  bool equal = equalButForNames(lp);
  equal = equalNames(lp) && equal;
  return equal;
}

bool HighsLp::equalButForNames(const HighsLp& lp) const {
  bool equal = equalButForScalingAndNames(lp);
  equal = equalScaling(lp) && equal;
  return equal;
}

bool HighsLp::equalButForScalingAndNames(const HighsLp& lp) const {
  bool equal_vectors = true;
  equal_vectors = this->num_col_ == lp.num_col_ && equal_vectors;
  equal_vectors = this->num_row_ == lp.num_row_ && equal_vectors;
  equal_vectors = this->sense_ == lp.sense_ && equal_vectors;
  equal_vectors = this->offset_ == lp.offset_ && equal_vectors;
  equal_vectors = this->model_name_ == lp.model_name_ && equal_vectors;
  equal_vectors = this->col_cost_ == lp.col_cost_ && equal_vectors;
  equal_vectors = this->col_upper_ == lp.col_upper_ && equal_vectors;
  equal_vectors = this->col_lower_ == lp.col_lower_ && equal_vectors;
  equal_vectors = this->row_upper_ == lp.row_upper_ && equal_vectors;
  equal_vectors = this->row_lower_ == lp.row_lower_ && equal_vectors;
#ifndef NDEBUG
  if (!equal_vectors) printf("HighsLp::equalButForNames: Unequal vectors\n");
#endif
  const bool equal_matrix = this->a_matrix_ == lp.a_matrix_;
#ifndef NDEBUG
  if (!equal_matrix) printf("HighsLp::equalButForNames: Unequal matrix\n");
#endif
  return equal_vectors && equal_matrix;
}

bool HighsLp::equalNames(const HighsLp& lp) const {
  bool equal = true;
  equal = this->objective_name_ == lp.objective_name_ && equal;
  equal = this->row_names_ == lp.row_names_ && equal;
  equal = this->col_names_ == lp.col_names_ && equal;
  return equal;
}

bool HighsLp::equalScaling(const HighsLp& lp) const {
  bool equal = true;
  equal = this->scale_.strategy == lp.scale_.strategy && equal;
  equal = this->scale_.has_scaling == lp.scale_.has_scaling && equal;
  equal = this->scale_.num_col == lp.scale_.num_col && equal;
  equal = this->scale_.num_row == lp.scale_.num_row && equal;
  equal = this->scale_.cost == lp.scale_.cost && equal;
  equal = this->scale_.col == lp.scale_.col && equal;
  equal = this->scale_.row == lp.scale_.row && equal;
#ifndef NDEBUG
  if (!equal) printf("HighsLp::equalScaling: Unequal scaling\n");
#endif
  return equal;
}

double HighsLp::objectiveValue(const std::vector<double>& solution) const {
  assert((int)solution.size() >= this->num_col_);
  double objective_function_value = this->offset_;
  for (HighsInt iCol = 0; iCol < this->num_col_; iCol++)
    objective_function_value += this->col_cost_[iCol] * solution[iCol];
  return objective_function_value;
}

HighsCDouble HighsLp::objectiveCDoubleValue(
    const std::vector<double>& solution) const {
  assert((int)solution.size() >= this->num_col_);
  HighsCDouble objective_function_value = this->offset_;
  for (HighsInt iCol = 0; iCol < this->num_col_; iCol++)
    objective_function_value += this->col_cost_[iCol] * solution[iCol];
  return objective_function_value;
}

void HighsLp::setMatrixDimensions() {
  this->a_matrix_.num_col_ = this->num_col_;
  this->a_matrix_.num_row_ = this->num_row_;
}

void HighsLp::resetScale() {
  // Should allow user-supplied scale to be retained
  //  const bool dimensions_ok =
  //    this->scale_.num_col_ == this->num_col_ &&
  //    this->scale_.num_row_ == this->num_row_;
  this->clearScale();
}

void HighsLp::setFormat(const MatrixFormat format) {
  this->a_matrix_.setFormat(format);
}

void HighsLp::exactResize() {
  this->col_cost_.resize(this->num_col_);
  this->col_lower_.resize(this->num_col_);
  this->col_upper_.resize(this->num_col_);
  this->row_lower_.resize(this->num_row_);
  this->row_upper_.resize(this->num_row_);
  this->a_matrix_.exactResize();

  if ((int)this->col_names_.size()) this->col_names_.resize(this->num_col_);
  if ((int)this->row_names_.size()) this->row_names_.resize(this->num_row_);
  if ((int)this->integrality_.size()) this->integrality_.resize(this->num_col_);
}

void HighsLp::clear() {
  this->num_col_ = 0;
  this->num_row_ = 0;

  this->col_cost_.clear();
  this->col_lower_.clear();
  this->col_upper_.clear();
  this->row_lower_.clear();
  this->row_upper_.clear();

  this->a_matrix_.clear();

  this->sense_ = ObjSense::kMinimize;
  this->offset_ = 0;

  this->model_name_ = "";
  this->objective_name_ = "";

  this->new_col_name_ix_ = 0;
  this->new_row_name_ix_ = 0;
  this->col_names_.clear();
  this->row_names_.clear();

  this->integrality_.clear();

  this->col_hash_.clear();
  this->row_hash_.clear();

  this->user_cost_scale_ = 0;
  this->user_bound_scale_ = 0;

  this->clearScale();
  this->is_scaled_ = false;
  this->is_moved_ = false;
  this->cost_row_location_ = -1;
  this->has_infinite_cost_ = false;
  this->stats_.clear();
  this->mods_.clear();
}

void HighsLp::clearScale() {
  this->scale_.strategy = kSimplexScaleStrategyOff;
  this->scale_.has_scaling = false;
  this->scale_.num_col = 0;
  this->scale_.num_row = 0;
  this->scale_.cost = 0;
  this->scale_.col.clear();
  this->scale_.row.clear();
}

void HighsLp::clearScaling() {
  this->unapplyScale();
  this->clearScale();
}

void HighsLp::applyScale() {
  // Ensure that any scaling is applied
  const HighsScale& scale = this->scale_;
  if (this->is_scaled_) {
    // Already scaled - so check that there is scaling and return
    assert(scale.has_scaling);
    return;
  }
  // No scaling currently applied
  this->is_scaled_ = false;
  if (scale.has_scaling) {
    // Apply the scaling to the bounds, costs and matrix, and record
    // that it has been applied
    for (HighsInt iCol = 0; iCol < this->num_col_; iCol++) {
      this->col_lower_[iCol] /= scale.col[iCol];
      this->col_upper_[iCol] /= scale.col[iCol];
      this->col_cost_[iCol] *= scale.col[iCol];
    }
    for (HighsInt iRow = 0; iRow < this->num_row_; iRow++) {
      this->row_lower_[iRow] *= scale.row[iRow];
      this->row_upper_[iRow] *= scale.row[iRow];
    }
    this->a_matrix_.applyScale(scale);
    this->is_scaled_ = true;
  }
}

void HighsLp::unapplyScale() {
  // Ensure that any scaling is not applied
  const HighsScale& scale = this->scale_;
  if (!this->is_scaled_) {
    // Not scaled so return
    return;
  }
  // Already scaled - so check that there is scaling
  assert(scale.has_scaling);
  // Unapply the scaling to the bounds, costs and matrix, and record
  // that it has been unapplied
  for (HighsInt iCol = 0; iCol < this->num_col_; iCol++) {
    this->col_lower_[iCol] *= scale.col[iCol];
    this->col_upper_[iCol] *= scale.col[iCol];
    this->col_cost_[iCol] /= scale.col[iCol];
  }
  for (HighsInt iRow = 0; iRow < this->num_row_; iRow++) {
    this->row_lower_[iRow] /= scale.row[iRow];
    this->row_upper_[iRow] /= scale.row[iRow];
  }
  this->a_matrix_.unapplyScale(scale);
  this->is_scaled_ = false;
}

void HighsLp::moveBackLpAndUnapplyScaling(HighsLp& lp) {
  assert(this->is_moved_ == true);
  *this = std::move(lp);
  this->unapplyScale();
  assert(this->is_moved_ == false);
}

bool HighsLp::userBoundScaleOk(const HighsInt user_bound_scale,
                               const double infinite_bound) const {
  const HighsInt dl_user_bound_scale =
      user_bound_scale - this->user_bound_scale_;
  if (!dl_user_bound_scale) return true;
  if (!boundScaleOk(this->col_lower_, this->col_upper_, dl_user_bound_scale,
                    infinite_bound))
    return false;
  return boundScaleOk(this->row_lower_, this->row_upper_, dl_user_bound_scale,
                      infinite_bound);
}

void HighsLp::userBoundScale(const HighsInt user_bound_scale) {
  const HighsInt dl_user_bound_scale =
      user_bound_scale - this->user_bound_scale_;
  if (!dl_user_bound_scale) return;
  double dl_user_bound_scale_value = std::pow(2, dl_user_bound_scale);
  for (HighsInt iCol = 0; iCol < this->num_col_; iCol++) {
    this->col_lower_[iCol] *= dl_user_bound_scale_value;
    this->col_upper_[iCol] *= dl_user_bound_scale_value;
  }
  for (HighsInt iRow = 0; iRow < this->num_row_; iRow++) {
    this->row_lower_[iRow] *= dl_user_bound_scale_value;
    this->row_upper_[iRow] *= dl_user_bound_scale_value;
  }
  // Record the current user bound scaling applied to the LP
  this->user_bound_scale_ = user_bound_scale;
}

bool HighsLp::userCostScaleOk(const HighsInt user_cost_scale,
                              const double infinite_cost) const {
  const HighsInt dl_user_cost_scale = user_cost_scale - this->user_cost_scale_;
  if (!dl_user_cost_scale) return true;
  return costScaleOk(this->col_cost_, dl_user_cost_scale, infinite_cost);
}

void HighsLp::userCostScale(const HighsInt user_cost_scale) {
  const HighsInt dl_user_cost_scale = user_cost_scale - this->user_cost_scale_;
  if (!dl_user_cost_scale) return;
  double dl_user_cost_scale_value = std::pow(2, dl_user_cost_scale);
  for (HighsInt iCol = 0; iCol < this->num_col_; iCol++)
    this->col_cost_[iCol] *= dl_user_cost_scale_value;
  this->user_cost_scale_ = user_cost_scale;
}

void HighsLp::addColNames(const std::string name, const HighsInt num_new_col) {
  // Don't add names if there are no columns, or if the names are
  // already incomplete
  if (this->num_col_ == 0) return;
  HighsInt col_names_size = this->col_names_.size();
  if (col_names_size < this->num_col_) return;
  if (!this->col_hash_.name2index.size())
    this->col_hash_.form(this->col_names_);
  // Handle the addition of user-defined names later
  assert(name == "");
  for (HighsInt iCol = this->num_col_; iCol < this->num_col_ + num_new_col;
       iCol++) {
    const std::string col_name =
        "col_ekk_" + std::to_string(this->new_col_name_ix_++);
    bool added = false;
    auto search = this->col_hash_.name2index.find(col_name);
    if (search == this->col_hash_.name2index.end()) {
      // Name not found in hash
      if (col_names_size == this->num_col_) {
        // No space (or name) for this col name
        this->col_names_.push_back(col_name);
        added = true;
      } else if (col_names_size > iCol) {
        // Space for this col name. Only add if name is blank
        if (this->col_names_[iCol] == "") {
          this->col_names_[iCol] = col_name;
          added = true;
        }
      }
    }
    if (added) {
      const bool duplicate =
          !this->col_hash_.name2index.emplace(col_name, iCol).second;
      assert(!duplicate);
      assert(this->col_names_[iCol] == col_name);
      assert(this->col_hash_.name2index.find(col_name)->second == iCol);
    } else {
      // Duplicate name or other failure
      this->col_hash_.name2index.clear();
      return;
    }
  }
}

void HighsLp::addRowNames(const std::string name, const HighsInt num_new_row) {
  // Don't add names if there are no rows, or if the names are already
  // incomplete
  if (this->num_row_ == 0) return;
  HighsInt row_names_size = this->row_names_.size();
  if (row_names_size < this->num_row_) return;
  if (!this->row_hash_.name2index.size())
    this->row_hash_.form(this->row_names_);
  // Handle the addition of user-defined names later
  assert(name == "");
  for (HighsInt iRow = this->num_row_; iRow < this->num_row_ + num_new_row;
       iRow++) {
    const std::string row_name =
        "row_ekk_" + std::to_string(this->new_row_name_ix_++);
    bool added = false;
    auto search = this->row_hash_.name2index.find(row_name);
    if (search == this->row_hash_.name2index.end()) {
      // Name not found in hash
      if (row_names_size == this->num_row_) {
        // No space (or name) for this row name
        this->row_names_.push_back(row_name);
        added = true;
      } else if (row_names_size > iRow) {
        // Space for this row name. Only add if name is blank
        if (this->row_names_[iRow] == "") {
          this->row_names_[iRow] = row_name;
          added = true;
        }
      }
    }
    if (added) {
      const bool duplicate =
          !this->row_hash_.name2index.emplace(row_name, iRow).second;
      assert(!duplicate);
      assert(this->row_names_[iRow] == row_name);
      assert(this->row_hash_.name2index.find(row_name)->second == iRow);
    } else {
      // Duplicate name or other failure
      this->row_hash_.name2index.clear();
      return;
    }
  }
}

void HighsLp::deleteColsFromVectors(
    HighsInt& new_num_col, const HighsIndexCollection& index_collection) {
  assert(ok(index_collection));
  HighsInt from_k;
  HighsInt to_k;
  limits(index_collection, from_k, to_k);
  // Initialise new_num_col in case none is removed due to from_k > to_k
  new_num_col = this->num_col_;
  if (from_k > to_k) return;

  HighsInt delete_from_col;
  HighsInt delete_to_col;
  HighsInt keep_from_col;
  HighsInt keep_to_col = -1;
  HighsInt current_set_entry = 0;

  HighsInt col_dim = this->num_col_;
  new_num_col = 0;
  bool have_names = (this->col_names_.size() != 0);
  bool have_integrality = (this->integrality_.size() != 0);
  for (HighsInt k = from_k; k <= to_k; k++) {
    updateOutInIndex(index_collection, delete_from_col, delete_to_col,
                     keep_from_col, keep_to_col, current_set_entry);
    // Account for the initial columns being kept
    if (k == from_k) new_num_col = delete_from_col;
    if (delete_to_col >= col_dim - 1) break;
    assert(delete_to_col < col_dim);
    for (HighsInt col = keep_from_col; col <= keep_to_col; col++) {
      this->col_cost_[new_num_col] = this->col_cost_[col];
      this->col_lower_[new_num_col] = this->col_lower_[col];
      this->col_upper_[new_num_col] = this->col_upper_[col];
      if (have_names) this->col_names_[new_num_col] = this->col_names_[col];
      if (have_integrality)
        this->integrality_[new_num_col] = this->integrality_[col];
      new_num_col++;
    }
    if (keep_to_col >= col_dim - 1) break;
  }
  this->col_cost_.resize(new_num_col);
  this->col_lower_.resize(new_num_col);
  this->col_upper_.resize(new_num_col);
  if (have_integrality) this->integrality_.resize(new_num_col);
  if (have_names) this->col_names_.resize(new_num_col);
}

void HighsLp::deleteRowsFromVectors(
    HighsInt& new_num_row, const HighsIndexCollection& index_collection) {
  assert(ok(index_collection));
  HighsInt from_k;
  HighsInt to_k;
  limits(index_collection, from_k, to_k);
  // Initialise new_num_row in case none is removed due to from_k > to_k
  new_num_row = this->num_row_;
  if (from_k > to_k) return;

  HighsInt delete_from_row;
  HighsInt delete_to_row;
  HighsInt keep_from_row;
  HighsInt keep_to_row = -1;
  HighsInt current_set_entry = 0;

  HighsInt row_dim = this->num_row_;
  new_num_row = 0;
  bool have_names = (HighsInt)this->row_names_.size() > 0;
  for (HighsInt k = from_k; k <= to_k; k++) {
    updateOutInIndex(index_collection, delete_from_row, delete_to_row,
                     keep_from_row, keep_to_row, current_set_entry);
    if (k == from_k) {
      // Account for the initial rows being kept
      new_num_row = delete_from_row;
    }
    if (delete_to_row >= row_dim - 1) break;
    assert(delete_to_row < row_dim);
    for (HighsInt row = keep_from_row; row <= keep_to_row; row++) {
      this->row_lower_[new_num_row] = this->row_lower_[row];
      this->row_upper_[new_num_row] = this->row_upper_[row];
      if (have_names) this->row_names_[new_num_row] = this->row_names_[row];
      new_num_row++;
    }
    if (keep_to_row >= row_dim - 1) break;
  }
  this->row_lower_.resize(new_num_row);
  this->row_upper_.resize(new_num_row);
  if (have_names) this->row_names_.resize(new_num_row);
}

void HighsLp::deleteCols(const HighsIndexCollection& index_collection) {
  HighsInt new_num_col;
  this->deleteColsFromVectors(new_num_col, index_collection);
  this->a_matrix_.deleteCols(index_collection);
  this->num_col_ = new_num_col;
}

void HighsLp::deleteRows(const HighsIndexCollection& index_collection) {
  HighsInt new_num_row;
  this->deleteRowsFromVectors(new_num_row, index_collection);
  this->a_matrix_.deleteRows(index_collection);
  this->num_row_ = new_num_row;
}

void HighsLpStats::clear() {
  valid = false;
  model = "Not set";                            // "";
  num_col = -kHighsIInf;                        // 0;
  num_row = -kHighsIInf;                        // 0;
  num_nz = -kHighsIInf;                         // 0;
  relative_max_cost_entry = -kHighsInf;         // 0;
  relative_num_equal_cost = -kHighsInf;         // 0;
  relative_num_inf_upper = -kHighsInf;          // 0;
  relative_num_equations = -kHighsInf;          // 0;
  relative_max_rhs_entry = -kHighsInf;          // 0;
  relative_num_equal_rhs = -kHighsInf;          // 0;
  a_matrix_density = -kHighsInf;                // 0;
  a_matrix_nz_per_col = -kHighsInf;             // 0;
  a_matrix_nz_per_row = -kHighsInf;             // 0;
  relative_max_matrix_entry = -kHighsInf;       // 0;
  relative_num_equal_a_matrix_nz = -kHighsInf;  // 0;
  relative_num_dense_row = -kHighsInf;          // 0;
}

void HighsLpStats::report(FILE* file, std::string message,
                          const HighsInt style) {
  if (style == HighsLpStatsReportPretty) {
    fprintf(file, "\nLP stats\n");
    if (message == "") {
      fprintf(file, "\n");
    } else {
      fprintf(file, " for %s\n", message.c_str());
    }
    fprintf(file, "   Valid             = %d\n", valid);
    fprintf(file, "   Model             = %s\n", model.c_str());
    fprintf(file, "   Number of columns = %d\n", num_col);
    fprintf(file, "   Number of rows =    %d\n", num_row);
    fprintf(file, "   Number of entries = %d\n", num_nz);
    fprintf(file, "   Relative maximum cost_entry =                     %g\n",
            relative_max_cost_entry);
    fprintf(file, "   Relative number of identical costs =              %g\n",
            relative_num_equal_cost);
    fprintf(file, "   Relative number of infinite column upper bounds = %g\n",
            relative_num_inf_upper);
    fprintf(file, "   Relative number of equations =                    %g\n",
            relative_num_equations);
    fprintf(file, "   Relative maximum rhs entry =                      %g\n",
            relative_max_rhs_entry);
    fprintf(file, "   Relative number of identical rhs entries =        %g\n",
            relative_num_equal_rhs);
    fprintf(file, "   Constraint matrix stats\n");
    fprintf(file, "      Density =                                     %g\n",
            a_matrix_density);
    fprintf(file, "      Nonzeros per column =                         %g\n",
            a_matrix_nz_per_col);
    fprintf(file, "      Nonzeros per row =                            %g\n",
            a_matrix_nz_per_row);
    fprintf(file, "      Relative maximum entry =                      %g\n",
            relative_max_matrix_entry);
    fprintf(file, "      Relative number of almost identical entries = %g\n",
            relative_num_equal_a_matrix_nz);
    fprintf(file, "      Relative number of dense rows =               %g\n",
            relative_num_dense_row);
  } else if (style == HighsLpStatsReportCsvHeader) {
    fprintf(file,
            "valid,model,col,row,nz,relative_max_cost_entry,relative_num_equal_"
            "cost,relative_num_inf_upper,relative_num_equations,relative_max_"
            "rhs_entry,relative_num_equal_rhs,a_matrix_density,a_matrix_nz_per_"
            "col,a_matrix_nz_per_row,relative_max_matrix_entry,relative_num_"
            "equal_a_matrix_nz,relative_num_dense_row");
  } else if (style == HighsLpStatsReportCsvData) {
    fprintf(file, "%d,%s,%d,%d,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g",
            int(this->valid), this->model.c_str(), int(this->num_col),
            int(this->num_row), int(this->num_nz), relative_max_cost_entry,
            relative_num_equal_cost, relative_num_inf_upper,
            relative_num_equations, relative_max_rhs_entry,
            relative_num_equal_rhs, a_matrix_density, a_matrix_nz_per_col,
            a_matrix_nz_per_row, relative_max_matrix_entry,
            relative_num_equal_a_matrix_nz, relative_num_dense_row);
  } else {
    fprintf(file, "Unknown LP stats report style of %d\n", int(style));
    assert(123 == 456);
  }
}

/*void reportValueCount(
    const std::vector<std::pair<double, HighsInt>> value_count,
    const double tolerance) {
  printf("Index              Value Count");
  if (tolerance > 0)
    printf(": %s %g\n", ": tolerance = ", tolerance);
  else
    printf("\n");

  for (HighsInt iX = 0; iX < HighsInt(value_count.size()); iX++)
    printf("   %2d %18.12g    %2d\n", int(iX), value_count[iX].first,
           int(value_count[iX].second));
}
*/

void HighsLp::stats() {
  const HighsInt num_nz = this->a_matrix_.numNz();
  stats_.model = this->model_name_;
  stats_.num_col = this->num_col_;
  stats_.num_row = this->num_row_;
  stats_.num_nz = num_nz;
  std::vector<std::pair<double, HighsInt>> value_count;
  double max_cost = 0;
  double min_cost = kHighsInf;
  value_count = valueCountSorted(this->col_cost_);
  reportValueCount(value_count, "Column cost");
  for (HighsInt iX = 0; iX < HighsInt(value_count.size()); iX++) {
    double cost = value_count[iX].first;
    if (cost == 0) continue;
    double abs_cost = std::fabs(cost);
    max_cost = std::max(abs_cost, max_cost);
    min_cost = std::min(abs_cost, min_cost);
  }
  // If there is a nonzero cost then min_cost and max_cost will both
  // be positive and finite
  assert(max_cost == 0 || (0 < min_cost && min_cost < kHighsInf));
  this->stats_.relative_max_cost_entry = max_cost > 0 ? max_cost / min_cost : 0;
  // Lambda for relative number of identical values
  auto numRelativeIdentical = [&](const HighsInt num_value) {
    if (num_value == 0) return 0.0;
    if (value_count.size() == 1) {
      // If all values are the same, the number of different values is
      // 1, so the relative number of different values is 0, so the
      // relative number of identical values is 1
      //
      // Yes, there's an inconsistency, in that the code below would
      // return 1 - 1/num_value < 1
      return 1.0;
    } else {
      // The relative number of identical values is (1 - relative
      // number of different values), and the relative number of
      // different values is
      //
      // (number of different values)/(number of values)
      //
      // Hence, if all values are different, the number of different
      // values is the number of values, so the relative number of
      // different values is 1, so the relative number of identical
      // values is 0
      return 1.0 - (1.0 * value_count.size()) / num_value;
    }
  };

  this->stats_.relative_num_equal_cost = numRelativeIdentical(this->num_col_);

  HighsInt num_inf_upper = 0;
  for (HighsInt iCol = 0; iCol < this->num_col_; iCol++)
    if (this->col_upper_[iCol] >= kHighsInf) num_inf_upper++;
  this->stats_.relative_num_inf_upper =
      this->num_col_ > 0 ? (1.0 * num_inf_upper) / this->num_col_ : 0;

  std::vector<double> rhs;
  double max_rhs = 0;
  double min_rhs = kHighsInf;
  HighsInt num_equations = 0;
  for (HighsInt iRow = 0; iRow < this->num_row_; iRow++) {
    double lower = this->row_lower_[iRow];
    double upper = this->row_upper_[iRow];
    if (lower == upper) {
      num_equations++;
      rhs.push_back(lower);
    } else {
      if (lower > -kHighsInf) {
        rhs.push_back(lower);
        if (lower) {
          double abs_rhs = std::fabs(lower);
          max_rhs = std::max(abs_rhs, max_rhs);
          min_rhs = std::min(abs_rhs, min_rhs);
        }
      }
      if (upper < kHighsInf) {
        rhs.push_back(upper);
        if (upper) {
          double abs_rhs = std::fabs(upper);
          max_rhs = std::max(abs_rhs, max_rhs);
          min_rhs = std::min(abs_rhs, min_rhs);
        }
      }
    }
  }
  value_count = valueCountSorted(rhs);
  HighsInt num_rhs = rhs.size();
  reportValueCount(value_count, "RHS");
  // If there is a nonzero rhs then min_rhs and max_rhs will both
  // be positive and finite
  assert(max_rhs == 0 || (0 < min_rhs && min_rhs < kHighsInf));
  this->stats_.relative_max_rhs_entry = max_rhs > 0 ? max_rhs / min_rhs : 0;
  this->stats_.relative_num_equations =
      this->num_row_ > 0 ? (1.0 * num_equations) / this->num_row_ : 0;

  this->stats_.relative_num_equal_rhs = numRelativeIdentical(num_rhs);

  const bool nontrivial_matrix = this->num_col_ > 0 && this->num_row_ > 0;
  this->stats_.a_matrix_density =
      nontrivial_matrix ? ((1.0 * num_nz) / this->num_col_) / this->num_row_
                        : 0;
  this->stats_.a_matrix_nz_per_col =
      nontrivial_matrix ? ((1.0 * num_nz) / this->num_col_) : 0;
  this->stats_.a_matrix_nz_per_row =
      nontrivial_matrix ? ((1.0 * num_nz) / this->num_row_) : 0;
  if (num_nz > 0) {
    const double value_cluster_size = 1e-4;
    value_count =
        valueCountSorted(this->a_matrix_.value_, true, value_cluster_size);
    reportValueCount(value_count, "Matrix", value_cluster_size);
    this->stats_.relative_num_equal_a_matrix_nz = numRelativeIdentical(num_nz);
  }
  double max_value = 0;
  double min_value = kHighsInf;
  for (HighsInt iX = 0; iX < HighsInt(value_count.size()); iX++) {
    double value = value_count[iX].first;
    if (value == 0) continue;
    double abs_value = std::fabs(value);
    max_value = std::max(abs_value, max_value);
    min_value = std::min(abs_value, min_value);
  }
  // If there is a nonzero value then min_value and max_value will both
  // be positive and finite
  assert(max_value == 0 || (0 < min_value && min_value < kHighsInf));
  this->stats_.relative_max_matrix_entry =
      max_value > 0 ? max_value / min_value : 0;

  const double kDenseRowTolerance = 0.5;
  const HighsInt dense_row_count = this->num_col_ * kDenseRowTolerance;
  std::vector<HighsInt> row_count;
  this->a_matrix_.getRowCounts(row_count);
  HighsInt num_dense_row = 0;
  for (HighsInt iRow = 0; iRow < this->num_row_; iRow++)
    if (row_count[iRow] >= dense_row_count) num_dense_row++;

  this->stats_.relative_num_dense_row =
      this->num_row_ > 0 ? (1.0 * num_dense_row) / this->num_row_ : 0;
  this->stats_.valid = true;
}

void HighsLp::unapplyMods() {
  // Restore any non-semi types
  const HighsInt num_non_semi = this->mods_.save_non_semi_variable_index.size();
  for (HighsInt k = 0; k < num_non_semi; k++) {
    HighsInt iCol = this->mods_.save_non_semi_variable_index[k];
    assert(this->integrality_[iCol] == HighsVarType::kContinuous ||
           this->integrality_[iCol] == HighsVarType::kInteger);
    if (this->integrality_[iCol] == HighsVarType::kContinuous) {
      this->integrality_[iCol] = HighsVarType::kSemiContinuous;
    } else {
      this->integrality_[iCol] = HighsVarType::kSemiInteger;
    }
  }
  // Restore any inconsistent semi variables
  const HighsInt num_inconsistent_semi =
      this->mods_.save_inconsistent_semi_variable_index.size();
  if (!num_inconsistent_semi) {
    assert(
        !this->mods_.save_inconsistent_semi_variable_lower_bound_value.size());
    assert(
        !this->mods_.save_inconsistent_semi_variable_upper_bound_value.size());
    assert(!this->mods_.save_inconsistent_semi_variable_type.size());
  }
  for (HighsInt k = 0; k < num_inconsistent_semi; k++) {
    HighsInt iCol = this->mods_.save_inconsistent_semi_variable_index[k];
    this->col_lower_[iCol] =
        this->mods_.save_inconsistent_semi_variable_lower_bound_value[k];
    this->col_upper_[iCol] =
        this->mods_.save_inconsistent_semi_variable_upper_bound_value[k];
    this->integrality_[iCol] =
        this->mods_.save_inconsistent_semi_variable_type[k];
  }
  // Restore any relaxed lower bounds
  std::vector<HighsInt>& relaxed_semi_variable_lower_index =
      this->mods_.save_relaxed_semi_variable_lower_bound_index;
  std::vector<double>& relaxed_semi_variable_lower_value =
      this->mods_.save_relaxed_semi_variable_lower_bound_value;
  const HighsInt num_lower_bound = relaxed_semi_variable_lower_index.size();
  if (!num_lower_bound) {
    assert(!relaxed_semi_variable_lower_value.size());
  }
  for (HighsInt k = 0; k < num_lower_bound; k++) {
    HighsInt iCol = relaxed_semi_variable_lower_index[k];
    assert(this->integrality_[iCol] == HighsVarType::kSemiContinuous ||
           this->integrality_[iCol] == HighsVarType::kSemiInteger);
    this->col_lower_[iCol] = relaxed_semi_variable_lower_value[k];
  }
  // Restore any tightened upper bounds
  std::vector<HighsInt>& tightened_semi_variable_upper_bound_index =
      this->mods_.save_tightened_semi_variable_upper_bound_index;
  std::vector<double>& tightened_semi_variable_upper_bound_value =
      this->mods_.save_tightened_semi_variable_upper_bound_value;
  const HighsInt num_upper_bound =
      tightened_semi_variable_upper_bound_index.size();
  if (!num_upper_bound) {
    assert(!tightened_semi_variable_upper_bound_value.size());
  }
  for (HighsInt k = 0; k < num_upper_bound; k++) {
    HighsInt iCol = tightened_semi_variable_upper_bound_index[k];
    assert(this->integrality_[iCol] == HighsVarType::kSemiContinuous ||
           this->integrality_[iCol] == HighsVarType::kSemiInteger);
    this->col_upper_[iCol] = tightened_semi_variable_upper_bound_value[k];
  }

  this->mods_.clear();
}

void HighsLpMods::clear() {
  this->save_non_semi_variable_index.clear();
  this->save_inconsistent_semi_variable_index.clear();
  this->save_inconsistent_semi_variable_lower_bound_value.clear();
  this->save_inconsistent_semi_variable_upper_bound_value.clear();
  this->save_inconsistent_semi_variable_type.clear();
  this->save_relaxed_semi_variable_lower_bound_index.clear();
  this->save_relaxed_semi_variable_lower_bound_value.clear();
  this->save_tightened_semi_variable_upper_bound_index.clear();
  this->save_tightened_semi_variable_upper_bound_value.clear();
  this->save_inf_cost_variable_index.clear();
  this->save_inf_cost_variable_cost.clear();
  this->save_inf_cost_variable_lower.clear();
  this->save_inf_cost_variable_upper.clear();
}

bool HighsLpMods::isClear() {
  if (this->save_non_semi_variable_index.size()) return false;
  if (this->save_inconsistent_semi_variable_index.size()) return false;
  if (this->save_inconsistent_semi_variable_lower_bound_value.size())
    return false;
  if (this->save_inconsistent_semi_variable_upper_bound_value.size())
    return false;
  if (this->save_inconsistent_semi_variable_type.size()) return false;
  if (this->save_relaxed_semi_variable_lower_bound_value.size()) return false;
  if (this->save_relaxed_semi_variable_lower_bound_value.size()) return false;
  if (this->save_tightened_semi_variable_upper_bound_index.size()) return false;
  if (this->save_tightened_semi_variable_upper_bound_value.size()) return false;
  return true;
}

void HighsNameHash::form(const std::vector<std::string>& name) {
  size_t num_name = name.size();
  this->clear();
  for (size_t index = 0; index < num_name; index++) {
    auto emplace_result = this->name2index.emplace(name[index], index);
    const bool duplicate = !emplace_result.second;
    if (duplicate) {
      // Find the original and mark it as duplicate
      auto& search = emplace_result.first;
      assert(int(search->second) < int(this->name2index.size()));
      search->second = kHashIsDuplicate;
    }
  }
}

bool HighsNameHash::hasDuplicate(const std::vector<std::string>& name) {
  size_t num_name = name.size();
  this->clear();
  bool has_duplicate = false;
  for (size_t index = 0; index < num_name; index++) {
    has_duplicate = !this->name2index.emplace(name[index], index).second;
    if (has_duplicate) break;
  }
  this->clear();
  return has_duplicate;
}

void HighsNameHash::update(int index, const std::string& old_name,
                           const std::string& new_name) {
  this->name2index.erase(old_name);
  auto emplace_result = this->name2index.emplace(new_name, index);
  const bool duplicate = !emplace_result.second;
  if (duplicate) {
    // Find the original and mark it as duplicate
    auto& search = emplace_result.first;
    assert(int(search->second) < int(this->name2index.size()));
    search->second = kHashIsDuplicate;
  }
}

void HighsNameHash::clear() { this->name2index.clear(); }
