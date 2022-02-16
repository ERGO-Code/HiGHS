/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2022 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Leona Gottwald and Michael    */
/*    Feldmeier                                                          */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/FilereaderLp.cpp
 * @brief
 */

#include "io/FilereaderLp.h"

#include <cstdarg>
#include <exception>
#include <map>

#include "../extern/filereaderlp/reader.hpp"
#include "lp_data/HighsLpUtils.h"

FilereaderRetcode FilereaderLp::readModelFromFile(const HighsOptions& options,
                                                  const std::string filename,
                                                  HighsModel& model) {
  HighsLp& lp = model.lp_;
  HighsHessian& hessian = model.hessian_;
  try {
    Model m = readinstance(filename);

    // build variable index and gather variable information
    std::map<std::string, unsigned int> varindex;

    lp.num_col_ = m.variables.size();
    lp.num_row_ = m.constraints.size();
    lp.integrality_.assign(lp.num_col_, HighsVarType::kContinuous);
    HighsInt num_continuous = 0;
    for (HighsUInt i = 0; i < m.variables.size(); i++) {
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
    if (num_continuous == m.variables.size()) lp.integrality_.clear();
    // get objective
    if (m.objective->offset) {
      highsLogUser(options.log_options, HighsLogType::kWarning,
                   "Ignoring m.objective->offset = %g\n", m.objective->offset);
      lp.offset_ = 0;  // m.objective->offset;
    }
    lp.col_cost_.resize(lp.num_col_, 0.0);
    for (HighsUInt i = 0; i < m.objective->linterms.size(); i++) {
      std::shared_ptr<LinTerm> lt = m.objective->linterms[i];
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
        hessian.dim_++;
      }
    }

    unsigned int qnnz = 0;
    // model_.hessian_ is initialised with start_[0] for fictitious
    // column 0, so have to clear this before pushing back start
    hessian.start_.clear();
    assert((int)hessian.start_.size() == 0);
    for (std::shared_ptr<Variable> var : m.variables) {
      hessian.start_.push_back(qnnz);

      for (unsigned int i = 0; i < mat[var].size(); i++) {
        hessian.index_.push_back(varindex[mat[var][i]->name]);
        hessian.value_.push_back(mat2[var][i]);
        qnnz++;
      }
    }
    hessian.start_.push_back(qnnz);
    hessian.format_ = HessianFormat::kSquare;

    // handle constraints
    std::map<std::shared_ptr<Variable>, std::vector<unsigned int>>
        consofvarmap_index;
    std::map<std::shared_ptr<Variable>, std::vector<double>> consofvarmap_value;
    for (HighsUInt i = 0; i < m.constraints.size(); i++) {
      std::shared_ptr<Constraint> con = m.constraints[i];
      for (HighsUInt j = 0; j < con->expr->linterms.size(); j++) {
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
    }

    HighsInt nz = 0;
    // lp.a_matrix_ is initialised with start_[0] for fictitious
    // column 0, so have to clear this before pushing back start
    lp.a_matrix_.start_.clear();
    assert((int)lp.a_matrix_.start_.size() == 0);
    for (HighsInt i = 0; i < lp.num_col_; i++) {
      std::shared_ptr<Variable> var = m.variables[i];
      lp.a_matrix_.start_.push_back(nz);
      for (HighsUInt j = 0; j < consofvarmap_index[var].size(); j++) {
        lp.a_matrix_.index_.push_back(consofvarmap_index[var][j]);
        lp.a_matrix_.value_.push_back(consofvarmap_value[var][j]);
        nz++;
      }
    }
    lp.a_matrix_.start_.push_back(nz);
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    lp.sense_ = m.sense == ObjectiveSense::MIN ? ObjSense::kMinimize
                                               : ObjSense::kMaximize;
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
  return FilereaderRetcode::kOk;
}

void FilereaderLp::writeToFile(FILE* file, const char* format, ...) {
  va_list argptr;
  va_start(argptr, format);
  char stringbuffer[LP_MAX_LINE_LENGTH + 1];
  HighsInt tokenlength = vsprintf(stringbuffer, format, argptr);
  if (this->linelength + tokenlength >= LP_MAX_LINE_LENGTH) {
    fprintf(file, "\n");
    fprintf(file, "%s", stringbuffer);
    this->linelength = tokenlength;
  } else {
    fprintf(file, "%s", stringbuffer);
    this->linelength += tokenlength;
  }
}

void FilereaderLp::writeToFileLineend(FILE* file) {
  fprintf(file, "\n");
  this->linelength = 0;
}

HighsStatus FilereaderLp::writeModelToFile(const HighsOptions& options,
                                           const std::string filename,
                                           const HighsModel& model) {
  const HighsLp& lp = model.lp_;
  assert(lp.a_matrix_.isColwise());
  FILE* file = fopen(filename.c_str(), "w");

  // write comment at the start of the file
  this->writeToFile(file, "\\ %s", LP_COMMENT_FILESTART);
  this->writeToFileLineend(file);

  // write objective
  this->writeToFile(file, "%s",
                    lp.sense_ == ObjSense::kMinimize ? "min" : "max");
  this->writeToFileLineend(file);
  this->writeToFile(file, " obj: ");
  for (HighsInt i = 0; i < lp.num_col_; i++) {
    this->writeToFile(file, "%+g x%" HIGHSINT_FORMAT " ", lp.col_cost_[i],
                      (i + 1));
  }
  if (model.isQp()) {
    this->writeToFile(file, "+ [ ");
    for (HighsInt col = 0; col < lp.num_col_; col++) {
      for (HighsInt i = model.hessian_.start_[col];
           i < model.hessian_.start_[col + 1]; i++) {
        if (col <= model.hessian_.index_[i]) {
          this->writeToFile(
              file, "%+g x%" HIGHSINT_FORMAT " * x%" HIGHSINT_FORMAT " ",
              model.hessian_.value_[i], col, model.hessian_.index_[i]);
        }
      }
    }
    this->writeToFile(file, " ]/2 ");
  }
  this->writeToFileLineend(file);

  // write constraint section, lower & upper bounds are one constraint
  // each
  this->writeToFile(file, "st");
  this->writeToFileLineend(file);
  for (HighsInt row = 0; row < lp.num_row_; row++) {
    if (lp.row_lower_[row] == lp.row_upper_[row]) {
      // equality constraint
      this->writeToFile(file, " con%" HIGHSINT_FORMAT ": ", row + 1);
      for (HighsInt var = 0; var < lp.num_col_; var++) {
        for (HighsInt idx = lp.a_matrix_.start_[var];
             idx < lp.a_matrix_.start_[var + 1]; idx++) {
          if (lp.a_matrix_.index_[idx] == row) {
            this->writeToFile(file, "%+g x%" HIGHSINT_FORMAT " ",
                              lp.a_matrix_.value_[idx], var + 1);
          }
        }
      }
      this->writeToFile(file, "= %+g", lp.row_lower_[row]);
      this->writeToFileLineend(file);
    } else {
      if (lp.row_lower_[row] > -kHighsInf) {
        // has a lower bounds
        this->writeToFile(file, " con%" HIGHSINT_FORMAT "lo: ", row + 1);
        for (HighsInt var = 0; var < lp.num_col_; var++) {
          for (HighsInt idx = lp.a_matrix_.start_[var];
               idx < lp.a_matrix_.start_[var + 1]; idx++) {
            if (lp.a_matrix_.index_[idx] == row) {
              this->writeToFile(file, "%+g x%" HIGHSINT_FORMAT " ",
                                lp.a_matrix_.value_[idx], var + 1);
            }
          }
        }
        this->writeToFile(file, ">= %+g", lp.row_lower_[row]);
        this->writeToFileLineend(file);
      } else if (lp.row_upper_[row] < kHighsInf) {
        // has an upper bounds
        this->writeToFile(file, " con%" HIGHSINT_FORMAT "up: ", row + 1);
        for (HighsInt var = 0; var < lp.num_col_; var++) {
          for (HighsInt idx = lp.a_matrix_.start_[var];
               idx < lp.a_matrix_.start_[var + 1]; idx++) {
            if (lp.a_matrix_.index_[idx] == row) {
              this->writeToFile(file, "%+g x%" HIGHSINT_FORMAT " ",
                                lp.a_matrix_.value_[idx], var + 1);
            }
          }
        }
        this->writeToFile(file, "<= %+g", lp.row_upper_[row]);
        this->writeToFileLineend(file);
      } else {
        // constraint has infinite lower & upper bounds so not a proper
        // constraint, does not get written
      }
    }
  }

  // write bounds section
  this->writeToFile(file, "bounds");
  this->writeToFileLineend(file);
  for (HighsInt i = 0; i < lp.num_col_; i++) {
    // if both lower/upper bound are +/-infinite: [name] free
    if (lp.col_lower_[i] > -kHighsInf && lp.col_upper_[i] < kHighsInf) {
      this->writeToFile(file, " %+g <= x%" HIGHSINT_FORMAT " <= %+g",
                        lp.col_lower_[i], i + 1, lp.col_upper_[i]);
      this->writeToFileLineend(file);
    } else if (lp.col_lower_[i] <= -kHighsInf && lp.col_upper_[i] < kHighsInf) {
      this->writeToFile(file, " -inf <= x%" HIGHSINT_FORMAT " <= %+g", i + 1,
                        lp.col_upper_[i]);
      this->writeToFileLineend(file);

    } else if (lp.col_lower_[i] > -kHighsInf && lp.col_upper_[i] >= kHighsInf) {
      this->writeToFile(file, " %+g <= x%" HIGHSINT_FORMAT " <= +inf",
                        lp.col_lower_[i], i + 1);
      this->writeToFileLineend(file);
    } else {
      this->writeToFile(file, " x%" HIGHSINT_FORMAT " free", i + 1);
      this->writeToFileLineend(file);
    }
  }

  if (lp.integrality_.size() > 0) {
    // write binary section
    this->writeToFile(file, "bin");
    this->writeToFileLineend(file);
    for (HighsInt i = 0; i < lp.num_col_; i++) {
      if (lp.integrality_[i] == HighsVarType::kInteger ||
          lp.integrality_[i] == HighsVarType::kSemiInteger) {
        if (lp.col_lower_[i] == 0.0 && lp.col_upper_[i] == 1.0) {
          this->writeToFile(file, " x%" HIGHSINT_FORMAT, i + 1);
          this->writeToFileLineend(file);
        }
      }
    }

    // write general section
    this->writeToFile(file, "gen");
    this->writeToFileLineend(file);
    for (HighsInt i = 0; i < lp.num_col_; i++) {
      if (lp.integrality_[i] == HighsVarType::kInteger ||
          lp.integrality_[i] == HighsVarType::kSemiInteger) {
        if (lp.col_lower_[i] != 0.0 || lp.col_upper_[i] != 1.0) {
          this->writeToFile(file, " x%" HIGHSINT_FORMAT, i + 1);
          this->writeToFileLineend(file);
        }
      }
    }

    // write semi section
    this->writeToFile(file, "semi");
    this->writeToFileLineend(file);
    for (HighsInt i = 0; i < lp.num_col_; i++) {
      if (lp.integrality_[i] == HighsVarType::kSemiContinuous ||
          lp.integrality_[i] == HighsVarType::kSemiInteger) {
        if (lp.col_lower_[i] != 0.0 || lp.col_upper_[i] != 1.0) {
          this->writeToFile(file, " x%" HIGHSINT_FORMAT, i + 1);
          this->writeToFileLineend(file);
        }
      }
    }
  }
  // write end
  this->writeToFile(file, "end");
  this->writeToFileLineend(file);

  fclose(file);
  return HighsStatus::kOk;
}
