/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/FilereaderLp.cpp
 * @brief
 */

#include "io/FilereaderLp.h"

#include <cstdarg>
#include <cstdio>
#include <exception>
#include <map>

#include "../extern/filereaderlp/reader.hpp"
#include "lp_data/HighsLpUtils.h"

FilereaderRetcode FilereaderLp::readModelFromFile(const HighsOptions& options,
                                                  const std::string filename,
                                                  HighsModel& model) {
  bool warning_issued = false;
  HighsLp& lp = model.lp_;
  HighsHessian& hessian = model.hessian_;
  try {
    Model m = readinstance(filename);

    if (!m.soss.empty()) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "SOS not supported by HiGHS\n");
      return FilereaderRetcode::kParserError;
    }

    // build variable index and gather variable information
    std::map<std::string, unsigned int> varindex;

    lp.num_col_ = m.variables.size();
    lp.num_row_ = m.constraints.size();
    lp.row_names_.resize(m.constraints.size());
    lp.integrality_.assign(lp.num_col_, HighsVarType::kContinuous);
    HighsInt num_continuous = 0;
    for (size_t i = 0; i < m.variables.size(); i++) {
      varindex[m.variables[i]->name] = i;
      lp.col_lower_.push_back(m.variables[i]->lowerbound);
      lp.col_upper_.push_back(m.variables[i]->upperbound);
      lp.col_names_.push_back(m.variables[i]->name);
      if (m.variables[i]->type == VariableType::BINARY ||
          m.variables[i]->type == VariableType::GENERAL) {
        lp.integrality_[i] = HighsVarType::kInteger;
      } else if (m.variables[i]->type == VariableType::SEMICONTINUOUS) {
        lp.integrality_[i] = HighsVarType::kSemiContinuous;
      } else if (m.variables[i]->type == VariableType::SEMIINTEGER) {
        lp.integrality_[i] = HighsVarType::kSemiInteger;
      } else {
        lp.integrality_[i] = HighsVarType::kContinuous;
        num_continuous++;
      }
    }
    // Clear lp.integrality_ if problem is pure LP
    if (static_cast<size_t>(num_continuous) == m.variables.size())
      lp.integrality_.clear();
    // get objective
    lp.objective_name_ = m.objective->name;
    // ToDo: Fix m.objective->offset and then use it here
    //
    lp.offset_ = m.objective->offset;
    lp.col_cost_.resize(lp.num_col_, 0.0);
    for (const auto& lt : m.objective->linterms) {
      lp.col_cost_[varindex[lt->var->name]] = lt->coef;
    }

    std::map<std::shared_ptr<Variable>, std::vector<std::shared_ptr<Variable>>>
        mat;
    std::map<std::shared_ptr<Variable>, std::vector<double>> mat2;
    for (std::shared_ptr<QuadTerm> qt : m.objective->quadterms) {
      if (qt->var1 != qt->var2) {
        mat[qt->var1].push_back(qt->var2);
        mat2[qt->var1].push_back(qt->coef / 2);
        mat[qt->var2].push_back(qt->var1);
        mat2[qt->var2].push_back(qt->coef / 2);
      } else {
        mat[qt->var1].push_back(qt->var2);
        mat2[qt->var1].push_back(qt->coef);
      }
    }

    // Determine whether there is a Hessian to set up by counting its
    // nonzero entries
    unsigned int qnnz = 0;
    for (std::shared_ptr<Variable> var : m.variables)
      for (size_t i = 0; i < mat[var].size(); i++)
        if (mat2[var][i]) qnnz++;
    if (qnnz) {
      hessian.dim_ = m.variables.size();
      qnnz = 0;
      // model_.hessian_ is initialised with start_[0] for fictitious
      // column 0, so have to clear this before pushing back start
      hessian.start_.clear();
      assert((int)hessian.start_.size() == 0);
      for (std::shared_ptr<Variable> var : m.variables) {
        hessian.start_.push_back(qnnz);
        for (size_t i = 0; i < mat[var].size(); i++) {
          double value = mat2[var][i];
          if (value) {
            hessian.index_.push_back(varindex[mat[var][i]->name]);
            hessian.value_.push_back(value);
            qnnz++;
          }
        }
      }
      hessian.start_.push_back(qnnz);
      hessian.format_ = HessianFormat::kSquare;
    } else {
      assert(hessian.dim_ == 0 && hessian.start_[0] == 0);
    }

    // handle constraints
    std::map<std::shared_ptr<Variable>, std::vector<unsigned int>>
        consofvarmap_index;
    std::map<std::shared_ptr<Variable>, std::vector<double>> consofvarmap_value;
    for (size_t i = 0; i < m.constraints.size(); i++) {
      std::shared_ptr<Constraint> con = m.constraints[i];
      lp.row_names_[i] = con->expr->name;
      for (size_t j = 0; j < con->expr->linterms.size(); j++) {
        std::shared_ptr<LinTerm> lt = con->expr->linterms[j];
        if (consofvarmap_index.count(lt->var) == 0) {
          consofvarmap_index[lt->var] = std::vector<unsigned int>();
          consofvarmap_value[lt->var] = std::vector<double>();
        }
        consofvarmap_index[lt->var].push_back(i);
        consofvarmap_value[lt->var].push_back(lt->coef);
      }

      lp.row_lower_.push_back(con->lowerbound);
      lp.row_upper_.push_back(con->upperbound);

      if (!con->expr->quadterms.empty()) {
        highsLogUser(options.log_options, HighsLogType::kError,
                     "Quadratic constraints not supported by HiGHS\n");
        return FilereaderRetcode::kParserError;
      }
    }

    // Check for empty row names, giving them a special name if possible
    bool highs_prefix_ok = true;
    bool used_highs_prefix = false;
    std::string highs_prefix = "HiGHS_R";
    for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
      // Look to see whether the name begins HiGHS_R
      if (strncmp(lp.row_names_[iRow].c_str(), highs_prefix.c_str(), 7) == 0) {
        printf("Name %s begins with \"HiGHS_R\"\n",
               lp.row_names_[iRow].c_str());
        highs_prefix_ok = false;
      } else if (lp.row_names_[iRow] == "") {
        // Make up a name beginning HiGHS_R
        lp.row_names_[iRow] = highs_prefix + std::to_string(iRow);
        used_highs_prefix = true;
      }
    }
    if (used_highs_prefix && !highs_prefix_ok) {
      // Have made up a name beginning HiGHS_R, but this occurs with
      // other "natural" rows, so abandon the row names
      lp.row_names_.clear();
      highsLogUser(options.log_options, HighsLogType::kWarning,
                   "Cannot create row name beginning \"HiGHS_R\" due to others "
                   "with same prefix: row names cleared\n");
    }

    HighsInt num_nz = 0;
    // lp.a_matrix_ is initialised with start_[0] for fictitious
    // column 0, so have to clear this before pushing back start
    lp.a_matrix_.start_.clear();
    assert((int)lp.a_matrix_.start_.size() == 0);
    for (const auto& var : m.variables) {
      lp.a_matrix_.start_.push_back(num_nz);
      for (size_t j = 0; j < consofvarmap_index[var].size(); j++) {
        double value = consofvarmap_value[var][j];
        if (value) {
          lp.a_matrix_.index_.push_back(consofvarmap_index[var][j]);
          lp.a_matrix_.value_.push_back(value);
          num_nz++;
        }
      }
    }
    lp.a_matrix_.start_.push_back(num_nz);
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    lp.sense_ = m.sense == ObjectiveSense::MIN ? ObjSense::kMinimize
                                               : ObjSense::kMaximize;
    // In a .lp file, more than one term involving the same variable
    // can appear in a constraint, resulting in repeated row indices
    // in a column
    HighsSparseMatrix& matrix = lp.a_matrix_;
    assert(matrix.isColwise());
    num_nz = 0;
    std::vector<double> column(lp.num_row_, 0);
    std::vector<HighsInt> nz_count(lp.num_row_, 0);
    std::vector<HighsInt> zero_count(lp.num_row_, 0);
    HighsInt sum_num_duplicate = 0;
    HighsInt sum_num_zero = 0;
    HighsInt sum_cancellation = 0;
    HighsInt num_report = 0;
    HighsInt max_num_report = 10;
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
      // Save the start, since this will be reduced if there are
      // repeated row indices in a column
      const HighsInt from_el = matrix.start_[iCol];
      for (HighsInt iEl = from_el; iEl < matrix.start_[iCol + 1]; iEl++) {
        // Add in the value to zero or any previous nonzero in this
        // row
        HighsInt iRow = matrix.index_[iEl];
        double value = matrix.value_[iEl];
        if (value) {
          column[iRow] += value;
          nz_count[iRow]++;
        } else {
          zero_count[iRow]++;
        }
      }
      // Pass through the column again, storing and then zeroing the
      // entries in column - both to eliminate the duplicate and
      // ensure that column is zeroed for the next matrix column.
      matrix.start_[iCol] = num_nz;
      for (HighsInt iEl = from_el; iEl < matrix.start_[iCol + 1]; iEl++) {
        HighsInt iRow = matrix.index_[iEl];
        if (column[iRow]) {
          assert(num_nz <= iEl);
          matrix.index_[num_nz] = iRow;
          matrix.value_[num_nz] = column[iRow];
          num_nz++;
        }
        // Report explicit zeros and/or sum to a nonzero value
        HighsInt num_ocurrence = zero_count[iRow] + nz_count[iRow];
        if (num_ocurrence > 1) {
          if (nz_count[iRow] > 1) {
            if (num_report < max_num_report)
              highsLogUser(options.log_options, HighsLogType::kWarning,
                           "Column %d (name \"%s\") occurs %d times in row %d "
                           "(name \"%s\"): values summed to %g\n",
                           int(iCol), lp.col_names_[iCol].c_str(),
                           int(num_ocurrence), int(iRow),
                           lp.row_names_[iRow].c_str(), column[iRow]);
            num_report++;
          }
          if (zero_count[iRow] > 0) {
            if (num_report < max_num_report)
              highsLogUser(options.log_options, HighsLogType::kWarning,
                           "Column %d (name \"%s\") contains %d explicit zero "
                           "coefficient%s in row %d (name \"%s\")\n",
                           int(iCol), lp.col_names_[iCol].c_str(),
                           int(zero_count[iRow]),
                           zero_count[iRow] > 1 ? "s" : "", int(iRow),
                           lp.row_names_[iRow].c_str());
            num_report++;
          }
          sum_num_duplicate += (num_ocurrence - 1);
          sum_num_zero += zero_count[iRow];
          if (column[iRow] == 0 && nz_count[iRow] > 0) sum_cancellation++;
        }
        zero_count[iRow] = 0;
        nz_count[iRow] = 0;
        column[iRow] = 0;
      }
      for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
        assert(zero_count[iRow] == 0);
        assert(nz_count[iRow] == 0);
        assert(column[iRow] == 0);
      }
    }
    matrix.start_[lp.num_col_] = num_nz;
    warning_issued = sum_num_duplicate > 0 || sum_num_zero > 0;
    HighsInt num_report_skipped = num_report - max_num_report;
    if (num_report_skipped > 0)
      highsLogUser(options.log_options, HighsLogType::kInfo,
                   "Skipped %d further warning%s of this kind\n",
                   int(num_report_skipped), num_report_skipped > 1 ? "s" : "");

    if (sum_num_duplicate > 0)
      highsLogUser(options.log_options, HighsLogType::kWarning,
                   "lp file contains %d repeated variable%s in constraints: "
                   "summing them yielded %d cancellation%s\n",
                   int(sum_num_duplicate), sum_num_duplicate > 1 ? "s" : "",
                   int(sum_cancellation),
                   (sum_cancellation == 0 || sum_cancellation > 1) ? "s" : "");
    if (sum_num_zero > 0)
      highsLogUser(options.log_options, HighsLogType::kWarning,
                   "lp file contains %d explicit zero%s\n", int(sum_num_zero),
                   sum_num_zero > 1 ? "s" : "");

  } catch (std::invalid_argument& ex) {
    // lpassert in extern/filereaderlp/def.hpp throws
    // std::invalid_argument whatever the error. Hence, unless
    // something is done specially - here or elsewhere -
    // FilereaderRetcode::kParserError will be returned.
    //
    // This is misleading when the file isn't found, as it's not a
    // parser error
    FILE* file = fopen(filename.c_str(), "r");
    if (file == nullptr) return FilereaderRetcode::kFileNotFound;
    fclose(file);
    return FilereaderRetcode::kParserError;
  }
  lp.ensureColwise();
  return warning_issued ? FilereaderRetcode::kWarning : FilereaderRetcode::kOk;
}

void FilereaderLp::writeToFile(FILE* file, const char* format, ...) {
  va_list argptr;
  va_start(argptr, format);
  std::array<char, LP_MAX_LINE_LENGTH + 1> stringbuffer = {};
  HighsInt tokenlength =
      vsnprintf(stringbuffer.data(), stringbuffer.size(), format, argptr);
  va_end(argptr);
  if (this->linelength + tokenlength >= LP_MAX_LINE_LENGTH) {
    fprintf(file, "\n");
    fprintf(file, "%s", stringbuffer.data());
    this->linelength = tokenlength;
  } else {
    fprintf(file, "%s", stringbuffer.data());
    this->linelength += tokenlength;
  }
}

void FilereaderLp::writeToFileLineEnd(FILE* file) {
  fprintf(file, "\n");
  this->linelength = 0;
}

void FilereaderLp::writeToFileValue(FILE* file, const double value,
                                    const bool force_plus) {
  // As for writeModelAsMps
  if (force_plus) {
    this->writeToFile(file, " %+.15g", value);
  } else {
    this->writeToFile(file, " %.15g", value);
  }
}

void FilereaderLp::writeToFileVar(FILE* file, const std::string var_name) {
  this->writeToFile(file, " %s", var_name.c_str());
}

void FilereaderLp::writeToFileMatrixRow(FILE* file, const HighsInt iRow,
                                        const HighsSparseMatrix ar_matrix,
                                        const std::vector<string> col_names) {
  assert(ar_matrix.isRowwise());

  for (HighsInt iEl = ar_matrix.start_[iRow]; iEl < ar_matrix.start_[iRow + 1];
       iEl++) {
    HighsInt iCol = ar_matrix.index_[iEl];
    double coef = ar_matrix.value_[iEl];
    this->writeToFileValue(file, coef);
    this->writeToFileVar(file, col_names[iCol]);
  }
}

HighsStatus FilereaderLp::writeModelToFile(const HighsOptions& options,
                                           const std::string filename,
                                           const HighsModel& model) {
  const HighsLp& lp = model.lp_;

  const bool ok_names = lp.okNames();
  assert(ok_names);
  if (!ok_names) return HighsStatus::kError;

  // Create a row-wise copy of the matrix
  HighsSparseMatrix ar_matrix = lp.a_matrix_;
  ar_matrix.ensureRowwise();

  FILE* file = fopen(filename.c_str(), "w");

  // write comment at the start of the file
  this->writeToFile(file, "\\ %s", LP_COMMENT_FILESTART);
  this->writeToFileLineEnd(file);

  // write objective
  this->writeToFile(file, "%s",
                    lp.sense_ == ObjSense::kMinimize ? "min" : "max");
  this->writeToFileLineEnd(file);
  this->writeToFile(file, " obj:");
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    double coef = lp.col_cost_[iCol];
    if (coef != 0.0) {
      this->writeToFileValue(file, coef);
      this->writeToFileVar(file, lp.col_names_[iCol]);
    }
  }
  this->writeToFile(file,
                    " ");  // ToDo Unnecessary, but only to give empty diff
  if (model.isQp()) {
    this->writeToFile(file, "+ [");
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
      for (HighsInt iEl = model.hessian_.start_[iCol];
           iEl < model.hessian_.start_[iCol + 1]; iEl++) {
        HighsInt iRow = model.hessian_.index_[iEl];
        if (iCol <= iRow) {
          double coef = model.hessian_.value_[iEl];
          if (iCol != iRow) coef *= 2;
          if (coef != 0.0) {
            this->writeToFileValue(file, coef);
            this->writeToFileVar(file, lp.col_names_[iCol]);
            this->writeToFile(file, " *");
            this->writeToFileVar(file, lp.col_names_[iRow]);
          }
        }
      }
    }
    this->writeToFile(file,
                      "  ]/2 ");  // ToDo Surely needs only to be one space
  }
  double coef = lp.offset_;
  if (coef != 0) this->writeToFileValue(file, coef);
  this->writeToFileLineEnd(file);

  // write constraint section, lower & upper bounds are one constraint
  // each
  this->writeToFile(file, "st");
  this->writeToFileLineEnd(file);
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    if (lp.row_lower_[iRow] == lp.row_upper_[iRow]) {
      // Equality constraint
      this->writeToFileVar(file, lp.row_names_[iRow]);
      this->writeToFile(file, ":");
      this->writeToFileMatrixRow(file, iRow, ar_matrix, lp.col_names_);
      this->writeToFile(file, " =");
      this->writeToFileValue(file, lp.row_lower_[iRow], true);
      this->writeToFileLineEnd(file);
    } else {
      // Need to distinguish the names when writing out boxed
      // constraint row as two single-sided constraints
      const bool boxed =
          lp.row_lower_[iRow] > -kHighsInf && lp.row_upper_[iRow] < kHighsInf;
      if (lp.row_lower_[iRow] > -kHighsInf) {
        // Has a lower bound
        this->writeToFileVar(file, lp.row_names_[iRow]);
        if (boxed) {
          this->writeToFile(file, "lo:");
        } else {
          this->writeToFile(file, ":");
        }
        this->writeToFileMatrixRow(file, iRow, ar_matrix, lp.col_names_);
        this->writeToFile(file, " >=");
        this->writeToFileValue(file, lp.row_lower_[iRow], true);
        this->writeToFileLineEnd(file);
      }
      if (lp.row_upper_[iRow] < kHighsInf) {
        // Has an upper bound
        this->writeToFileVar(file, lp.row_names_[iRow]);
        if (boxed) {
          this->writeToFile(file, "up:");
        } else {
          this->writeToFile(file, ":");
        }
        this->writeToFileMatrixRow(file, iRow, ar_matrix, lp.col_names_);
        this->writeToFile(file, " <=");
        this->writeToFileValue(file, lp.row_upper_[iRow], true);
        this->writeToFileLineEnd(file);
      }
    }
  }

  // write bounds section
  this->writeToFile(file, "bounds");
  this->writeToFileLineEnd(file);
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    const bool default_bounds =
        lp.col_lower_[iCol] == 0 && lp.col_upper_[iCol] == kHighsInf;
    if (default_bounds) continue;
    if (lp.col_lower_[iCol] <= -kHighsInf && lp.col_upper_[iCol] >= kHighsInf) {
      // Free variable
      this->writeToFileVar(file, lp.col_names_[iCol]);
      this->writeToFile(file, " free");
    } else if (lp.col_lower_[iCol] == lp.col_upper_[iCol]) {
      // Fixed variable
      this->writeToFileVar(file, lp.col_names_[iCol]);
      this->writeToFile(file, " =");
      this->writeToFileValue(file, lp.col_upper_[iCol], false);
    } else {
      assert(!default_bounds);
      // Non-default bound
      if (lp.col_lower_[iCol] != 0) {
        // Nonzero lower bound
        this->writeToFileValue(file, lp.col_lower_[iCol], false);
        this->writeToFile(file, " <=");
      }
      this->writeToFileVar(file, lp.col_names_[iCol]);
      if (lp.col_upper_[iCol] < kHighsInf) {
        // Finite upper bound
        this->writeToFile(file, " <=");
        this->writeToFileValue(file, lp.col_upper_[iCol], false);
      }
    }
    this->writeToFileLineEnd(file);
  }
  if (lp.integrality_.size() > 0) {
    // write binary section
    this->writeToFile(file, "bin");
    this->writeToFileLineEnd(file);
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
      if (lp.integrality_[iCol] == HighsVarType::kInteger) {
        if (lp.col_lower_[iCol] == 0.0 && lp.col_upper_[iCol] == 1.0) {
          this->writeToFileVar(file, lp.col_names_[iCol]);
          this->writeToFileLineEnd(file);
        }
      }
    }

    // write general section
    this->writeToFile(file, "gen");
    this->writeToFileLineEnd(file);
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
      if (lp.integrality_[iCol] == HighsVarType::kInteger) {
        if (lp.col_lower_[iCol] != 0.0 || lp.col_upper_[iCol] != 1.0) {
          this->writeToFileVar(file, lp.col_names_[iCol]);
          this->writeToFileLineEnd(file);
        }
      }
    }

    // write semi section
    this->writeToFile(file, "semi");
    this->writeToFileLineEnd(file);
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
      if (lp.integrality_[iCol] == HighsVarType::kSemiContinuous ||
          lp.integrality_[iCol] == HighsVarType::kSemiInteger) {
        this->writeToFileVar(file, lp.col_names_[iCol]);
        this->writeToFileLineEnd(file);
      }
    }
  }
  // write end
  this->writeToFile(file, "end");
  this->writeToFileLineEnd(file);

  fclose(file);
  return HighsStatus::kOk;
}
