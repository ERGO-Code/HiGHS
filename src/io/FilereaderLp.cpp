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

#include "../extern/filereaderlp/reader.hpp"
#include "lp_data/HighsLpUtils.h"

FilereaderRetcode FilereaderLp::readModelFromFile(const HighsOptions& options,
                                                  const std::string filename,
                                                  HighsLp& model) {
  try {
    Model m = readinstance(filename);

    // build variable index and gather variable information
    std::map<std::string, unsigned int> varindex;

    model.numCol_ = m.variables.size();
    model.numRow_ = m.constraints.size();
    for (HighsUInt i = 0; i < m.variables.size(); i++) {
      varindex[m.variables[i]->name] = i;
      model.colLower_.push_back(m.variables[i]->lowerbound);
      model.colUpper_.push_back(m.variables[i]->upperbound);
      model.col_names_.push_back(m.variables[i]->name);
    }

    // get objective
    model.offset_ = m.objective->offset;
    model.colCost_.resize(model.numCol_, 0.0);
    for (HighsUInt i = 0; i < m.objective->linterms.size(); i++) {
      std::shared_ptr<LinTerm> lt = m.objective->linterms[i];
      model.colCost_[varindex[lt->var->name]] = lt->coef;
    }

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

      model.rowLower_.push_back(con->lowerbound);
      model.rowUpper_.push_back(con->upperbound);
    }

    HighsInt nz = 0;
    for (HighsInt i = 0; i < model.numCol_; i++) {
      std::shared_ptr<Variable> var = m.variables[i];
      model.Astart_.push_back(nz);
      for (HighsUInt j = 0; j < consofvarmap_index[var].size(); j++) {
        model.Aindex_.push_back(consofvarmap_index[var][j]);
        model.Avalue_.push_back(consofvarmap_value[var][j]);
        nz++;
      }
    }
    model.Astart_.push_back(nz);
    model.sense_ = m.sense == ObjectiveSense::MIN ? ObjSense::kMinimize
                                                  : ObjSense::kMaximize;
  } catch (std::invalid_argument& ex) {
    return FilereaderRetcode::kParserError;
  }
  setOrientation(model);
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
                                           const HighsLp& model) {
  assert(model.orientation_ != MatrixOrientation::kRowwise);
  FILE* file = fopen(filename.c_str(), "w");

  // write comment at the start of the file
  this->writeToFile(file, "\\ %s", LP_COMMENT_FILESTART);
  this->writeToFileLineend(file);

  // write objective
  this->writeToFile(file, "%s",
                    model.sense_ == ObjSense::kMinimize ? "min" : "max");
  this->writeToFileLineend(file);
  this->writeToFile(file, " obj: ");
  for (HighsInt i = 0; i < model.numCol_; i++) {
    this->writeToFile(file, "%+g x%" HIGHSINT_FORMAT " ", model.colCost_[i],
                      (i + 1));
  }
  this->writeToFileLineend(file);

  // write constraint section, lower & upper bounds are one constraint
  // each
  this->writeToFile(file, "st");
  this->writeToFileLineend(file);
  for (HighsInt row = 0; row < model.numRow_; row++) {
    if (model.rowLower_[row] == model.rowUpper_[row]) {
      // equality constraint
      this->writeToFile(file, " con%" HIGHSINT_FORMAT ": ", row + 1);
      for (HighsInt var = 0; var < model.numCol_; var++) {
        for (HighsInt idx = model.Astart_[var]; idx < model.Astart_[var + 1];
             idx++) {
          if (model.Aindex_[idx] == row) {
            this->writeToFile(file, "%+g x%" HIGHSINT_FORMAT " ",
                              model.Avalue_[idx], var + 1);
          }
        }
      }
      this->writeToFile(file, "= %+g", model.rowLower_[row]);
      this->writeToFileLineend(file);
    } else {
      if (model.rowLower_[row] > -kHighsInf) {
        // has a lower bounds
        this->writeToFile(file, " con%" HIGHSINT_FORMAT "lo: ", row + 1);
        for (HighsInt var = 0; var < model.numCol_; var++) {
          for (HighsInt idx = model.Astart_[var]; idx < model.Astart_[var + 1];
               idx++) {
            if (model.Aindex_[idx] == row) {
              this->writeToFile(file, "%+g x%" HIGHSINT_FORMAT " ",
                                model.Avalue_[idx], var + 1);
            }
          }
        }
        this->writeToFile(file, ">= %+g", model.rowLower_[row]);
        this->writeToFileLineend(file);
      } else if (model.rowUpper_[row] < kHighsInf) {
        // has an upper bounds
        this->writeToFile(file, " con%" HIGHSINT_FORMAT "up: ", row + 1);
        for (HighsInt var = 0; var < model.numCol_; var++) {
          for (HighsInt idx = model.Astart_[var]; idx < model.Astart_[var + 1];
               idx++) {
            if (model.Aindex_[idx] == row) {
              this->writeToFile(file, "%+g x%" HIGHSINT_FORMAT " ",
                                model.Avalue_[idx], var + 1);
            }
          }
        }
        this->writeToFile(file, "<= %+g", model.rowUpper_[row]);
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
  for (HighsInt i = 0; i < model.numCol_; i++) {
    // if both lower/upper bound are +/-infinite: [name] free
    if (model.colLower_[i] > -kHighsInf && model.colUpper_[i] < kHighsInf) {
      this->writeToFile(file, " %+g <= x%" HIGHSINT_FORMAT " <= %+g",
                        model.colLower_[i], i + 1, model.colUpper_[i]);
      this->writeToFileLineend(file);
    } else if (model.colLower_[i] <= -kHighsInf &&
               model.colUpper_[i] < kHighsInf) {
      this->writeToFile(file, " -inf <= x%" HIGHSINT_FORMAT " <= %+g", i + 1,
                        model.colUpper_[i]);
      this->writeToFileLineend(file);

    } else if (model.colLower_[i] > -kHighsInf &&
               model.colUpper_[i] >= kHighsInf) {
      this->writeToFile(file, " %+g <= x%" HIGHSINT_FORMAT " <= +inf",
                        model.colLower_[i], i + 1);
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
