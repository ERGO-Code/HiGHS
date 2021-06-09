/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/FilereaderLp.cpp
 * @brief
 */

#include "io/FilereaderLp.h"

#include <cstdarg>
#include <exception>
#include <map>

#include "../external/filereaderlp/reader.hpp"
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

    lp.numCol_ = m.variables.size();
    lp.numRow_ = m.constraints.size();
    for (HighsUInt i = 0; i < m.variables.size(); i++) {
      varindex[m.variables[i]->name] = i;
      lp.colLower_.push_back(m.variables[i]->lowerbound);
      lp.colUpper_.push_back(m.variables[i]->upperbound);
      lp.col_names_.push_back(m.variables[i]->name);
    }

    // get objective
    lp.offset_ = m.objective->offset;
    lp.colCost_.resize(lp.numCol_, 0.0);
    for (HighsUInt i = 0; i < m.objective->linterms.size(); i++) {
      std::shared_ptr<LinTerm> lt = m.objective->linterms[i];
      lp.colCost_[varindex[lt->var->name]] = lt->coef;
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
    for (std::shared_ptr<Variable> var : m.variables) {
      hessian.q_start_.push_back(qnnz);

      for (unsigned int i = 0; i < mat[var].size(); i++) {
        hessian.q_index_.push_back(varindex[mat[var][i]->name]);
        hessian.q_value_.push_back(mat2[var][i]);
        qnnz++;
      }
    }
    hessian.q_start_.push_back(qnnz);

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

      lp.rowLower_.push_back(con->lowerbound);
      lp.rowUpper_.push_back(con->upperbound);
    }

    HighsInt nz = 0;
    for (HighsInt i = 0; i < lp.numCol_; i++) {
      std::shared_ptr<Variable> var = m.variables[i];
      lp.Astart_.push_back(nz);
      for (HighsUInt j = 0; j < consofvarmap_index[var].size(); j++) {
        lp.Aindex_.push_back(consofvarmap_index[var][j]);
        lp.Avalue_.push_back(consofvarmap_value[var][j]);
        nz++;
      }
    }
    lp.Astart_.push_back(nz);
    lp.orientation_ = MatrixOrientation::kColwise;
    lp.sense_ = m.sense == ObjectiveSense::MIN ? ObjSense::kMinimize
                                               : ObjSense::kMaximize;
  } catch (std::invalid_argument& ex) {
    return FilereaderRetcode::kParserError;
  }
  if (setOrientation(lp) != HighsStatus::kOk)
    return FilereaderRetcode::kParserError;
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
  assert(lp.orientation_ != MatrixOrientation::kRowwise);
  FILE* file = fopen(filename.c_str(), "w");

  // write comment at the start of the file
  this->writeToFile(file, "\\ %s", LP_COMMENT_FILESTART);
  this->writeToFileLineend(file);

  // write objective
  this->writeToFile(file, "%s",
                    lp.sense_ == ObjSense::kMinimize ? "min" : "max");
  this->writeToFileLineend(file);
  this->writeToFile(file, " obj: ");
  for (HighsInt i = 0; i < lp.numCol_; i++) {
    this->writeToFile(file, "%+g x%" HIGHSINT_FORMAT " ", lp.colCost_[i],
                      (i + 1));
  }
  this->writeToFileLineend(file);

  // write constraint section, lower & upper bounds are one constraint
  // each
  this->writeToFile(file, "st");
  this->writeToFileLineend(file);
  for (HighsInt row = 0; row < lp.numRow_; row++) {
    if (lp.rowLower_[row] == lp.rowUpper_[row]) {
      // equality constraint
      this->writeToFile(file, " con%" HIGHSINT_FORMAT ": ", row + 1);
      for (HighsInt var = 0; var < lp.numCol_; var++) {
        for (HighsInt idx = lp.Astart_[var]; idx < lp.Astart_[var + 1]; idx++) {
          if (lp.Aindex_[idx] == row) {
            this->writeToFile(file, "%+g x%" HIGHSINT_FORMAT " ",
                              lp.Avalue_[idx], var + 1);
          }
        }
      }
      this->writeToFile(file, "= %+g", lp.rowLower_[row]);
      this->writeToFileLineend(file);
    } else {
      if (lp.rowLower_[row] > -kHighsInf) {
        // has a lower bounds
        this->writeToFile(file, " con%" HIGHSINT_FORMAT "lo: ", row + 1);
        for (HighsInt var = 0; var < lp.numCol_; var++) {
          for (HighsInt idx = lp.Astart_[var]; idx < lp.Astart_[var + 1];
               idx++) {
            if (lp.Aindex_[idx] == row) {
              this->writeToFile(file, "%+g x%" HIGHSINT_FORMAT " ",
                                lp.Avalue_[idx], var + 1);
            }
          }
        }
        this->writeToFile(file, ">= %+g", lp.rowLower_[row]);
        this->writeToFileLineend(file);
      } else if (lp.rowUpper_[row] < kHighsInf) {
        // has an upper bounds
        this->writeToFile(file, " con%" HIGHSINT_FORMAT "up: ", row + 1);
        for (HighsInt var = 0; var < lp.numCol_; var++) {
          for (HighsInt idx = lp.Astart_[var]; idx < lp.Astart_[var + 1];
               idx++) {
            if (lp.Aindex_[idx] == row) {
              this->writeToFile(file, "%+g x%" HIGHSINT_FORMAT " ",
                                lp.Avalue_[idx], var + 1);
            }
          }
        }
        this->writeToFile(file, "<= %+g", lp.rowUpper_[row]);
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
  for (HighsInt i = 0; i < lp.numCol_; i++) {
    // if both lower/upper bound are +/-infinite: [name] free
    if (lp.colLower_[i] > -kHighsInf && lp.colUpper_[i] < kHighsInf) {
      this->writeToFile(file, " %+g <= x%" HIGHSINT_FORMAT " <= %+g",
                        lp.colLower_[i], i + 1, lp.colUpper_[i]);
      this->writeToFileLineend(file);
    } else if (lp.colLower_[i] <= -kHighsInf && lp.colUpper_[i] < kHighsInf) {
      this->writeToFile(file, " -inf <= x%" HIGHSINT_FORMAT " <= %+g", i + 1,
                        lp.colUpper_[i]);
      this->writeToFileLineend(file);

    } else if (lp.colLower_[i] > -kHighsInf && lp.colUpper_[i] >= kHighsInf) {
      this->writeToFile(file, " %+g <= x%" HIGHSINT_FORMAT " <= +inf",
                        lp.colLower_[i], i + 1);
      this->writeToFileLineend(file);
    } else {
      this->writeToFile(file, " x%" HIGHSINT_FORMAT " free", i + 1);
      this->writeToFileLineend(file);
    }
  }

  // write binary section
  this->writeToFile(file, "bin");
  this->writeToFileLineend(file);

  // write general section
  this->writeToFile(file, "gen");
  this->writeToFileLineend(file);

  // write semi section
  this->writeToFile(file, "semi");
  this->writeToFileLineend(file);

  // write end
  this->writeToFile(file, "end");
  this->writeToFileLineend(file);

  fclose(file);
  return HighsStatus::kOk;
}
