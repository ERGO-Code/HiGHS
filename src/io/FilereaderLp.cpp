/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/FilereaderLp.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#include "io/FilereaderLp.h"

#include <cstdarg>
#include <exception>
#include <map>

// #include "lp_data/HConst.h"
// #include "util/stringutil.h"

#include "../external/filereaderlp/reader.hpp"

FilereaderLp::FilereaderLp() {
}

void emptyTokenQueue(std::list<LpToken*>& list) {
  while (list.size() > 0) {
    LpToken* token = list.front();
    list.pop_front();
    delete token;
  }
}

FilereaderLp::~FilereaderLp() {
}

//   lp->numCol_ = this->variables.size();
//   lp->numRow_ = this->linearConstraints.size();

//   lp->sense_ = this->objSense;

//   // determine order of variables
//   HighsVar** variables = new HighsVar*[lp->numCol_];
//   for (int i = 0; i < lp->numCol_; i++) {
//     HighsVar* front = this->variables.front();
//     this->variables.pop_front();
//     this->variables.push_back(front);
//     variables[i] = front;
//     lp->colCost_.push_back(front->obj);
//     lp->colLower_.push_back(front->lowerBound);
//     lp->colUpper_.push_back(front->upperBound);
//   }

//   // determine order of constraints
//   HighsLinearCons** constraints = new HighsLinearCons*[lp->numRow_];
//   for (int i = 0; i < lp->numRow_; i++) {
//     HighsLinearCons* front = this->linearConstraints.front();
//     this->linearConstraints.pop_front();
//     this->linearConstraints.push_back(front);
//     constraints[i] = front;
//     lp->rowLower_.push_back(front->lowerBound);
//     lp->rowUpper_.push_back(front->upperBound);
//   }

//   // handle constraints
//   lp->Astart_.clear();
//   lp->Astart_.push_back(0);
//   for (int var = 0; var < lp->numCol_; var++) {
//     VarConsCoefsMap::iterator iter =
//         this->variableConstraintCoefficientMap.find(variables[var]);
//     if (iter != this->variableConstraintCoefficientMap.end()) {
//       std::list<HighsLinearConsCoef*>* coefs = iter->second;
//       int numberOfCoefficients = coefs->size();

//       lp->Astart_.push_back(lp->Astart_[var] + numberOfCoefficients);

//       for (int coef = 0; coef < numberOfCoefficients; coef++) {
//         HighsLinearConsCoef* front = coefs->front();
//         coefs->pop_front();
//         coefs->push_back(front);
//         lp->Avalue_.push_back(front->coef);
//         CoefConsMap::iterator it = this->coefficientConstraintMap.find(front);
//         if (it != this->coefficientConstraintMap.end()) {
//           // find index of constraint
//           HighsCons* currentCons = it->second;
//           for (int cons = 0; cons < lp->numRow_; cons++) {
//             if (constraints[cons] == currentCons) {
//               lp->Aindex_.push_back(cons);
//               break;
//             }
//           }
//         } else {
//           // ERROR
//         }
//       }
//     }
//   }

//   delete[] variables;
//   delete[] constraints;



FilereaderRetcode FilereaderLp::readModelFromFile(const HighsOptions& options,
                                                  HighsLp& model) {
  try {
    Model m = readinstance(options.model_file);
    // transform Model to HighsLp

    // build variable index and gather variable information
    std::map<std::shared_ptr<Variable>, unsigned int> varindex;


    model.numCol_ = m.variables.size();
    model.numRow_ = m.constraints.size();
    for (int i=0; i<m.variables.size(); i++) {
      varindex[m.variables[i]] = i;
      model.colLower_.push_back(m.variables[i]->lowerbound);
      model.colUpper_.push_back(m.variables[i]->upperbound);
      model.col_names_.push_back(m.variables[i]->name);
    }

    // get objective
    model.offset_ = m.objective->offset;
    for (int i=0; i<m.objective->linterms.size(); i++) {
      std::shared_ptr<LinTerm> lt = m.objective->linterms[i];
      model.colCost_[varindex[lt->var]] = lt->coef;
    }

    //handle constraints
    std::map<std::shared_ptr<Variable>, std::vector<unsigned int>> consofvarmap_index;
    std::map<std::shared_ptr<Variable>, std::vector<double>> consofvarmap_value;
    for(int i=0; i<m.constraints.size(); i++) {
      std::shared_ptr<Constraint> con = m.constraints[i];
      for(int j=0; j<con->expr->linterms.size(); j++) {
        std::shared_ptr<LinTerm> lt = con->expr->linterms[j];
        if (consofvarmap_index.count(lt->var) == 0) {
           consofvarmap_index[lt->var] = std::vector<unsigned int>();
           consofvarmap_value[lt->var] = std::vector<double>();
        }
        consofvarmap_index[lt->var].push_back(i);
        consofvarmap_value[lt->var].push_back(lt->coef);
      }

      model.rowLower_[i] = con->lowerbound;
      model.rowUpper_[i] = con->upperbound;
    }

    int nz = 0;
    for(int i=0; i<model.numCol_; i++) {
      std::shared_ptr<Variable> var = m.variables[i];
      model.Astart_.push_back(nz);
      for(int j=0; j<consofvarmap_index[var].size(); j++) {
        model.Aindex_[nz] = consofvarmap_index[var][j];
        model.Avalue_[nz] = consofvarmap_value[var][j];
        nz++;
      }
    }
    model.Astart_.push_back(nz);
    
  } catch(std::invalid_argument ex) {
     return FilereaderRetcode::PARSERERROR;
  }
  return FilereaderRetcode::OK;
}

void FilereaderLp::writeToFile(FILE* file, const char* format, ...) {
  va_list argptr;
  va_start(argptr, format);
  char stringbuffer[LP_MAX_LINE_LENGTH+1];
  int tokenlength = vsprintf(stringbuffer, format, argptr);
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
                                           const char* filename,
                                           HighsLp& model) {
  FILE* file = fopen(filename, "w");

  // write comment at the start of the file
  this->writeToFile(file, "\\ %s", LP_COMMENT_FILESTART);
  this->writeToFileLineend(file);

  // write objective
  this->writeToFile(file, "%s", model.sense_ == ObjSense::MINIMIZE
                              ? LP_KEYWORD_MIN[0]
                              : LP_KEYWORD_MAX[0]);
  this->writeToFileLineend(file);
  this->writeToFile(file, " obj: ");
  for (int i = 0; i < model.numCol_; i++) {
    this->writeToFile(file, "%+g x%d ", model.colCost_[i], (i + 1));
  }
  this->writeToFileLineend(file);

  // write constraint section, lower & upper bounds are one constraint each
  this->writeToFile(file, "%s", LP_KEYWORD_ST[2]);
  this->writeToFileLineend(file);
  for (int row = 0; row < model.numRow_; row++) {
    if (model.rowLower_[row] == model.rowUpper_[row]) {
      // equality constraint
      this->writeToFile(file, " con%d: ", row + 1);
      for (int var = 0; var < model.numCol_; var++) {
        for (int idx = model.Astart_[var]; idx < model.Astart_[var + 1];
             idx++) {
          if (model.Aindex_[idx] == row) {
            this->writeToFile(file, "%+g x%d ", model.Avalue_[idx], var + 1);
          }
        }
      }
      this->writeToFile(file, "= %+g", model.rowLower_[row]);
      this->writeToFileLineend(file);
    } else {
      if (model.rowLower_[row] >= -10E10) {
        // has a lower bounds
        this->writeToFile(file, " con%dlo: ", row + 1);
        for (int var = 0; var < model.numCol_; var++) {
          for (int idx = model.Astart_[var]; idx < model.Astart_[var + 1];
               idx++) {
            if (model.Aindex_[idx] == row) {
              this->writeToFile(file, "%+g x%d ", model.Avalue_[idx], var + 1);
            }
          }
        }
        this->writeToFile(file, ">= %+g", model.rowLower_[row]);
        this->writeToFileLineend(file);
      } else if (model.rowUpper_[row] <= 10E10) {
        // has an upper bounds
        this->writeToFile(file, " con%dup: ", row + 1);
        for (int var = 0; var < model.numCol_; var++) {
          for (int idx = model.Astart_[var]; idx < model.Astart_[var + 1];
               idx++) {
            if (model.Aindex_[idx] == row) {
              this->writeToFile(file, "%+g x%d ", model.Avalue_[idx], var + 1);
            }
          }
        }
        this->writeToFile(file, "<= %+g", model.rowLower_[row]);
        this->writeToFileLineend(file);
      } else {
        // constraint has infinite lower & upper bounds so not a proper
        // constraint, does not get written
      }
    }
  }

  // write bounds section
  this->writeToFile(file, "%s", LP_KEYWORD_BOUNDS[0]);
  this->writeToFileLineend(file);
  for (int i = 0; i < model.numCol_; i++) {
    // if both lower/upper bound are +/-infinite: [name] free
    if (model.colLower_[i] > -HIGHS_CONST_INF &&
        model.colUpper_[i] < HIGHS_CONST_INF) {
      this->writeToFile(file, " %+g <= x%d <= %+g", model.colLower_[i], i + 1,
                        model.colUpper_[i]);
      this->writeToFileLineend(file);
    } else if (model.colLower_[i] <= -HIGHS_CONST_INF &&
               model.colUpper_[i] < HIGHS_CONST_INF) {
      this->writeToFile(file, " -inf <= x%d <= %+g", i + 1, model.colUpper_[i]);
      this->writeToFileLineend(file);

    } else if (model.colLower_[i] > -HIGHS_CONST_INF &&
               model.colUpper_[i] >= HIGHS_CONST_INF) {
      this->writeToFile(file, " %+g <= x%d <= +inf", model.colLower_[i], i + 1);
      this->writeToFileLineend(file);
    } else {
      this->writeToFile(file, " x%d %s", i + 1, LP_KEYWORD_FREE[0]);
      this->writeToFileLineend(file);
    }
  }

  // write binary section
  this->writeToFile(file, "%s", LP_KEYWORD_BIN[0]);
  this->writeToFileLineend(file);

  // write general section
  this->writeToFile(file, "%s", LP_KEYWORD_GEN[0]);
  this->writeToFileLineend(file);

  // write semi section
  this->writeToFile(file, "%s", LP_KEYWORD_SEMI[1]);
  this->writeToFileLineend(file);

  // write end
  this->writeToFile(file, "%s", LP_KEYWORD_END[0]);
  this->writeToFileLineend(file);

  fclose(file);
  return HighsStatus::OK;
}
