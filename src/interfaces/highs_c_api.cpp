/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "highs_c_api.h"

#include "Highs.h"

HighsInt Highs_call(HighsInt numcol, HighsInt numrow, HighsInt numnz,
                    double* colcost, double* collower, double* colupper,
                    double* rowlower, double* rowupper, HighsInt* astart,
                    HighsInt* aindex, double* avalue, double* colvalue,
                    double* coldual, double* rowvalue, double* rowdual,
                    HighsInt* colbasisstatus, HighsInt* rowbasisstatus,
                    int* modelstatus) {
  Highs highs;

  HighsInt status =
      Highs_passLp(&highs, numcol, numrow, numnz, colcost, collower, colupper,
                   rowlower, rowupper, astart, aindex, avalue);
  if (status != 0) {
    return status;
  }

  status = (HighsInt)highs.run();

  if (status == 0) {
    HighsSolution solution;
    HighsBasis basis;
    solution = highs.getSolution();
    basis = highs.getBasis();
    *modelstatus = (HighsInt)highs.getModelStatus();

    for (HighsInt i = 0; i < numcol; i++) {
      colvalue[i] = solution.col_value[i];
      coldual[i] = solution.col_dual[i];

      colbasisstatus[i] = (HighsInt)basis.col_status[i];
    }

    for (HighsInt i = 0; i < numrow; i++) {
      rowvalue[i] = solution.row_value[i];
      rowdual[i] = solution.row_dual[i];

      rowbasisstatus[i] = (HighsInt)basis.row_status[i];
    }
  }

  return status;
}

void* Highs_create() { return new Highs(); }

void Highs_destroy(void* highs) { delete (Highs*)highs; }

HighsInt Highs_run(void* highs) { return (HighsInt)((Highs*)highs)->run(); }

HighsInt Highs_readModel(void* highs, const char* filename) {
  return (HighsInt)((Highs*)highs)->readModel(std::string(filename));
}

HighsInt Highs_writeModel(void* highs, const char* filename) {
  return (HighsInt)((Highs*)highs)->writeModel(std::string(filename));
}

HighsInt Highs_writeSolution(void* highs, const char* filename) {
  return (HighsInt)((Highs*)highs)->writeSolution(std::string(filename));
}

HighsInt Highs_passLp(void* highs, const HighsInt numcol, const HighsInt numrow,
                      const HighsInt numnz, const double* colcost,
                      const double* collower, const double* colupper,
                      const double* rowlower, const double* rowupper,
                      const HighsInt* astart, const HighsInt* aindex,
                      const double* avalue) {
  return (HighsInt)((Highs*)highs)
      ->passModel(numcol, numrow, numnz, colcost, collower, colupper, rowlower,
                  rowupper, astart, aindex, avalue);
}

HighsInt Highs_clearModel(void* highs) {
  return (HighsInt)((Highs*)highs)->clearModel();
}

HighsInt Highs_runQuiet(void* highs) {
  return (HighsInt)((Highs*)highs)->setHighsOptionValue("output_flag", false);
}

HighsInt Highs_setHighsLogfile(void* highs, void* logfile) {
  return (HighsInt)((Highs*)highs)->setHighsOptionValue("output_flag", false);
}

HighsInt Highs_setHighsOutput(void* highs, void* outputfile) {
  return (HighsInt)((Highs*)highs)->setHighsOptionValue("output_flag", false);
}

HighsInt Highs_setHighsBoolOptionValue(void* highs, const char* option,
                                       const HighsInt value) {
  return (HighsInt)((Highs*)highs)
      ->setHighsOptionValue(std::string(option), (bool)value);
}

HighsInt Highs_setHighsIntOptionValue(void* highs, const char* option,
                                      const HighsInt value) {
  return (HighsInt)((Highs*)highs)
      ->setHighsOptionValue(std::string(option), value);
}

HighsInt Highs_setHighsDoubleOptionValue(void* highs, const char* option,
                                         const double value) {
  return (HighsInt)((Highs*)highs)
      ->setHighsOptionValue(std::string(option), value);
}

HighsInt Highs_setHighsStringOptionValue(void* highs, const char* option,
                                         const char* value) {
  return (HighsInt)((Highs*)highs)
      ->setHighsOptionValue(std::string(option), std::string(value));
}

HighsInt Highs_setHighsOptionValue(void* highs, const char* option,
                                   const char* value) {
  return (HighsInt)((Highs*)highs)
      ->setHighsOptionValue(std::string(option), std::string(value));
}

HighsInt Highs_getHighsBoolOptionValue(void* highs, const char* option,
                                       HighsInt* value) {
  bool v;
  HighsInt retcode =
      (HighsInt)((Highs*)highs)->getHighsOptionValue(std::string(option), v);
  *value = (HighsInt)v;
  return retcode;
}

HighsInt Highs_getHighsIntOptionValue(void* highs, const char* option,
                                      HighsInt* value) {
  return (HighsInt)((Highs*)highs)
      ->getHighsOptionValue(std::string(option), *value);
}

HighsInt Highs_getHighsDoubleOptionValue(void* highs, const char* option,
                                         double* value) {
  return (HighsInt)((Highs*)highs)
      ->getHighsOptionValue(std::string(option), *value);
}

HighsInt Highs_getHighsStringOptionValue(void* highs, const char* option,
                                         char* value) {
  std::string v;
  HighsInt retcode =
      (HighsInt)((Highs*)highs)->getHighsOptionValue(std::string(option), v);
  strcpy(value, v.c_str());
  return retcode;
}

HighsInt Highs_getHighsOptionType(void* highs, const char* option,
                                  HighsInt* type) {
  HighsOptionType t;
  HighsInt retcode =
      (HighsInt)((Highs*)highs)->getHighsOptionType(std::string(option), t);
  *type = (HighsInt)t;
  return retcode;
}

HighsInt Highs_resetHighsOptions(void* highs) {
  return (HighsInt)((Highs*)highs)->resetHighsOptions();
}

HighsInt Highs_getHighsIntInfoValue(void* highs, const char* info,
                                    HighsInt* value) {
  return (HighsInt)((Highs*)highs)->getHighsInfoValue(info, *value);
}

HighsInt Highs_getHighsDoubleInfoValue(void* highs, const char* info,
                                       double* value) {
  return (HighsInt)((Highs*)highs)->getHighsInfoValue(info, *value);
}

void Highs_getSolution(void* highs, double* colvalue, double* coldual,
                       double* rowvalue, double* rowdual) {
  HighsSolution solution = ((Highs*)highs)->getSolution();

  for (HighsInt i = 0; i < (HighsInt)solution.col_value.size(); i++) {
    colvalue[i] = solution.col_value[i];
  }

  for (HighsInt i = 0; i < (HighsInt)solution.col_dual.size(); i++) {
    coldual[i] = solution.col_dual[i];
  }

  for (HighsInt i = 0; i < (HighsInt)solution.row_value.size(); i++) {
    rowvalue[i] = solution.row_value[i];
  }

  for (HighsInt i = 0; i < (HighsInt)solution.row_dual.size(); i++) {
    rowdual[i] = solution.row_dual[i];
  }
}

void Highs_getBasis(void* highs, HighsInt* colstatus, HighsInt* rowstatus) {
  HighsBasis basis = ((Highs*)highs)->getBasis();
  for (HighsInt i = 0; i < (HighsInt)basis.col_status.size(); i++) {
    colstatus[i] = (HighsInt)basis.col_status[i];
  }

  for (HighsInt i = 0; i < (HighsInt)basis.row_status.size(); i++) {
    rowstatus[i] = (HighsInt)basis.row_status[i];
  }
}

HighsInt Highs_getModelStatus(void* highs, const HighsInt scaled_model) {
  return (HighsInt)((Highs*)highs)->getModelStatus(scaled_model);
}

HighsInt Highs_getDualRay(void* highs, HighsInt* has_dual_ray,
                          double* dual_ray_value) {
  bool v;
  HighsInt retcode = (HighsInt)((Highs*)highs)->getDualRay(v, dual_ray_value);
  *has_dual_ray = (HighsInt)v;
  return retcode;
}

HighsInt Highs_getPrimalRay(void* highs, HighsInt* has_primal_ray,
                            double* primal_ray_value) {
  bool v;
  HighsInt retcode =
      (HighsInt)((Highs*)highs)->getPrimalRay(v, primal_ray_value);
  *has_primal_ray = (HighsInt)v;
  return retcode;
}

double Highs_getObjectiveValue(void* highs) {
  return ((Highs*)highs)->getObjectiveValue();
}

HighsInt Highs_getIterationCount(void* highs) {
  return Highs_getSimplexIterationCount(highs);
}

HighsInt Highs_getSimplexIterationCount(void* highs) {
  return (HighsInt)((Highs*)highs)->getSimplexIterationCount();
}

HighsInt Highs_getBasicVariables(void* highs, HighsInt* basic_variables) {
  return (HighsInt)((Highs*)highs)->getBasicVariables(basic_variables);
}

HighsInt Highs_getBasisInverseRow(void* highs, const HighsInt row,
                                  double* row_vector, HighsInt* row_num_nz,
                                  HighsInt* row_indices) {
  return (HighsInt)((Highs*)highs)
      ->getBasisInverseRow(row, row_vector, row_num_nz, row_indices);
}

HighsInt Highs_getBasisInverseCol(void* highs, const HighsInt col,
                                  double* col_vector, HighsInt* col_num_nz,
                                  HighsInt* col_indices) {
  return (HighsInt)((Highs*)highs)
      ->getBasisInverseCol(col, col_vector, col_num_nz, col_indices);
}

HighsInt Highs_getBasisSolve(void* highs, const double* rhs,
                             double* solution_vector, HighsInt* solution_num_nz,
                             HighsInt* solution_indices) {
  return (HighsInt)((Highs*)highs)
      ->getBasisSolve(rhs, solution_vector, solution_num_nz, solution_indices);
}

HighsInt Highs_getBasisTransposeSolve(void* highs, const double* rhs,
                                      double* solution_vector,
                                      HighsInt* solution_nz,
                                      HighsInt* solution_indices) {
  return (HighsInt)((Highs*)highs)
      ->getBasisTransposeSolve(rhs, solution_vector, solution_nz,
                               solution_indices);
}

HighsInt Highs_getReducedRow(void* highs, const HighsInt row,
                             double* row_vector, HighsInt* row_num_nz,
                             HighsInt* row_indices) {
  return (HighsInt)((Highs*)highs)
      ->getReducedRow(row, row_vector, row_num_nz, row_indices);
}

HighsInt Highs_getReducedColumn(void* highs, const HighsInt col,
                                double* col_vector, HighsInt* col_num_nz,
                                HighsInt* col_indices) {
  return (HighsInt)((Highs*)highs)
      ->getReducedColumn(col, col_vector, col_num_nz, col_indices);
}

HighsInt Highs_setBasis(void* highs, const HighsInt* colstatus,
                        const HighsInt* rowstatus) {
  HighsBasis basis;
  const HighsInt num_col = Highs_getNumCols(highs);
  basis.col_status.resize(num_col);
  for (HighsInt i = 0; i < num_col; i++) {
    if (colstatus[i] == (HighsInt)HighsBasisStatus::LOWER) {
      basis.col_status[i] = HighsBasisStatus::LOWER;
    } else if (colstatus[i] == (HighsInt)HighsBasisStatus::BASIC) {
      basis.col_status[i] = HighsBasisStatus::BASIC;
    } else if (colstatus[i] == (HighsInt)HighsBasisStatus::UPPER) {
      basis.col_status[i] = HighsBasisStatus::UPPER;
    } else if (colstatus[i] == (HighsInt)HighsBasisStatus::ZERO) {
      basis.col_status[i] = HighsBasisStatus::ZERO;
    } else if (colstatus[i] == (HighsInt)HighsBasisStatus::NONBASIC) {
      basis.col_status[i] = HighsBasisStatus::NONBASIC;
    } else {
      return (HighsInt)HighsStatus::Error;
    }
  }
  const HighsInt num_row = Highs_getNumRows(highs);
  basis.row_status.resize(num_row);
  for (HighsInt i = 0; i < num_row; i++) {
    if (rowstatus[i] == (HighsInt)HighsBasisStatus::LOWER) {
      basis.row_status[i] = HighsBasisStatus::LOWER;
    } else if (rowstatus[i] == (HighsInt)HighsBasisStatus::BASIC) {
      basis.row_status[i] = HighsBasisStatus::BASIC;
    } else if (rowstatus[i] == (HighsInt)HighsBasisStatus::UPPER) {
      basis.row_status[i] = HighsBasisStatus::UPPER;
    } else if (rowstatus[i] == (HighsInt)HighsBasisStatus::ZERO) {
      basis.row_status[i] = HighsBasisStatus::ZERO;
    } else if (rowstatus[i] == (HighsInt)HighsBasisStatus::NONBASIC) {
      basis.row_status[i] = HighsBasisStatus::NONBASIC;
    } else {
      return (HighsInt)HighsStatus::Error;
    }
  }
  return (HighsInt)((Highs*)highs)->setBasis(basis);
}

HighsInt Highs_setLogicalBasis(void* highs) {
  return (HighsInt)((Highs*)highs)->setBasis();
}

double Highs_getHighsRunTime(void* highs) {
  return (double)((Highs*)highs)->getHighsRunTime();
}

HighsInt Highs_addRow(void* highs, const double lower, const double upper,
                      const HighsInt num_new_nz, const HighsInt* indices,
                      const double* values) {
  return ((Highs*)highs)->addRow(lower, upper, num_new_nz, indices, values);
}

HighsInt Highs_addRows(void* highs, const HighsInt num_new_row,
                       const double* lower, const double* upper,
                       const HighsInt num_new_nz, const HighsInt* starts,
                       const HighsInt* indices, const double* values) {
  return ((Highs*)highs)
      ->addRows(num_new_row, lower, upper, num_new_nz, starts, indices, values);
}

HighsInt Highs_addCol(void* highs, const double cost, const double lower,
                      const double upper, const HighsInt num_new_nz,
                      const HighsInt* indices, const double* values) {
  return ((Highs*)highs)
      ->addCol(cost, lower, upper, num_new_nz, indices, values);
}

HighsInt Highs_addCols(void* highs, const HighsInt num_new_col,
                       const double* costs, const double* lower,
                       const double* upper, const HighsInt num_new_nz,
                       const HighsInt* starts, const HighsInt* indices,
                       const double* values) {
  return ((Highs*)highs)
      ->addCols(num_new_col, costs, lower, upper, num_new_nz, starts, indices,
                values);
}

HighsInt Highs_changeObjectiveSense(void* highs, const HighsInt sense) {
  ObjSense pass_sense = ObjSense::MINIMIZE;
  if (sense == (HighsInt)ObjSense::MAXIMIZE) pass_sense = ObjSense::MAXIMIZE;
  return ((Highs*)highs)->changeObjectiveSense(pass_sense);
}

HighsInt Highs_changeColsCost(void* highs, const HighsInt col,
                              const double cost) {
  return ((Highs*)highs)->changeColCost(col, cost);
}

HighsInt Highs_changeColsCostByRange(void* highs, const HighsInt from_col,
                                     const HighsInt to_col,
                                     const double* cost) {
  return ((Highs*)highs)->changeColsCost(from_col, to_col, cost);
}

HighsInt Highs_changeColsCostBySet(void* highs, const HighsInt num_set_entries,
                                   const HighsInt* set, const double* cost) {
  return ((Highs*)highs)->changeColsCost(num_set_entries, set, cost);
}

HighsInt Highs_changeColsCostByMask(void* highs, const HighsInt* mask,
                                    const double* cost) {
  return ((Highs*)highs)->changeColsCost(mask, cost);
}

HighsInt Highs_changeColBounds(void* highs, const HighsInt col,
                               const double lower, const double upper) {
  return ((Highs*)highs)->changeColBounds(col, lower, upper);
}

HighsInt Highs_changeColsBoundsByRange(void* highs, const HighsInt from_col,
                                       const HighsInt to_col,
                                       const double* lower,
                                       const double* upper) {
  return ((Highs*)highs)->changeColsBounds(from_col, to_col, lower, upper);
}

HighsInt Highs_changeColsBoundsBySet(void* highs,
                                     const HighsInt num_set_entries,
                                     const HighsInt* set, const double* lower,
                                     const double* upper) {
  return ((Highs*)highs)->changeColsBounds(num_set_entries, set, lower, upper);
}

HighsInt Highs_changeColsBoundsByMask(void* highs, const HighsInt* mask,
                                      const double* lower,
                                      const double* upper) {
  return ((Highs*)highs)->changeColsBounds(mask, lower, upper);
}

HighsInt Highs_changeRowBounds(void* highs, const HighsInt row,
                               const double lower, const double upper) {
  return ((Highs*)highs)->changeRowBounds(row, lower, upper);
}

HighsInt Highs_changeRowsBoundsBySet(void* highs,
                                     const HighsInt num_set_entries,
                                     const HighsInt* set, const double* lower,
                                     const double* upper) {
  return ((Highs*)highs)->changeRowsBounds(num_set_entries, set, lower, upper);
}

HighsInt Highs_changeRowsBoundsByMask(void* highs, const HighsInt* mask,
                                      const double* lower,
                                      const double* upper) {
  return ((Highs*)highs)->changeRowsBounds(mask, lower, upper);
}

HighsInt Highs_changeCoeff(void* highs, const HighsInt row, const HighsInt col,
                           const double value) {
  return ((Highs*)highs)->changeCoeff(row, col, value);
}

HighsInt Highs_getObjectiveSense(void* highs, HighsInt* sense) {
  ObjSense get_sense;
  HighsInt status = ((Highs*)highs)->getObjectiveSense(get_sense);
  *sense = (HighsInt)get_sense;
  return status;
}

HighsInt Highs_getColsByRange(void* highs, const HighsInt from_col,
                              const HighsInt to_col, HighsInt* num_col,
                              double* costs, double* lower, double* upper,
                              HighsInt* num_nz, HighsInt* matrix_start,
                              HighsInt* matrix_index, double* matrix_value) {
  HighsInt numcol, numnz;
  HighsInt status =
      ((Highs*)highs)
          ->getCols(from_col, to_col, numcol, costs, lower, upper, numnz,
                    matrix_start, matrix_index, matrix_value);
  *num_col = numcol;
  *num_nz = numnz;
  return status;
}

HighsInt Highs_getColsBySet(void* highs, const HighsInt num_set_entries,
                            const HighsInt* set, HighsInt* num_col,
                            double* costs, double* lower, double* upper,
                            HighsInt* num_nz, HighsInt* matrix_start,
                            HighsInt* matrix_index, double* matrix_value) {
  HighsInt numcol, numnz;
  HighsInt status =
      ((Highs*)highs)
          ->getCols(num_set_entries, set, numcol, costs, lower, upper, numnz,
                    matrix_start, matrix_index, matrix_value);
  *num_col = numcol;
  *num_nz = numnz;
  return status;
}

HighsInt Highs_getColsByMask(void* highs, const HighsInt* mask,
                             HighsInt* num_col, double* costs, double* lower,
                             double* upper, HighsInt* num_nz,
                             HighsInt* matrix_start, HighsInt* matrix_index,
                             double* matrix_value) {
  HighsInt numcol, numnz;
  HighsInt status = ((Highs*)highs)
                        ->getCols(mask, numcol, costs, lower, upper, numnz,
                                  matrix_start, matrix_index, matrix_value);
  *num_col = numcol;
  *num_nz = numnz;
  return status;
}

HighsInt Highs_getRowsByRange(void* highs, const HighsInt from_row,
                              const HighsInt to_row, HighsInt* num_row,
                              double* lower, double* upper, HighsInt* num_nz,
                              HighsInt* matrix_start, HighsInt* matrix_index,
                              double* matrix_value) {
  HighsInt numrow, numnz;
  HighsInt status = ((Highs*)highs)
                        ->getRows(from_row, to_row, numrow, lower, upper, numnz,
                                  matrix_start, matrix_index, matrix_value);
  *num_row = numrow;
  *num_nz = numnz;
  return status;
}

HighsInt Highs_getRowsBySet(void* highs, const HighsInt num_set_entries,
                            const HighsInt* set, HighsInt* num_row,
                            double* lower, double* upper, HighsInt* num_nz,
                            HighsInt* matrix_start, HighsInt* matrix_index,
                            double* matrix_value) {
  HighsInt numrow, numnz;
  HighsInt status =
      ((Highs*)highs)
          ->getRows(num_set_entries, set, numrow, lower, upper, numnz,
                    matrix_start, matrix_index, matrix_value);
  *num_row = numrow;
  *num_nz = numnz;
  return status;
}

HighsInt Highs_getRowsByMask(void* highs, const HighsInt* mask,
                             HighsInt* num_row, double* lower, double* upper,
                             HighsInt* num_nz, HighsInt* matrix_start,
                             HighsInt* matrix_index, double* matrix_value) {
  HighsInt numrow, numnz;
  HighsInt status = ((Highs*)highs)
                        ->getRows(mask, numrow, lower, upper, numnz,
                                  matrix_start, matrix_index, matrix_value);
  *num_row = numrow;
  *num_nz = numnz;
  return status;
}

HighsInt Highs_deleteColsByRange(void* highs, const HighsInt from_col,
                                 const HighsInt to_col) {
  return ((Highs*)highs)->deleteCols(from_col, to_col);
}

HighsInt Highs_deleteColsBySet(void* highs, const HighsInt num_set_entries,
                               const HighsInt* set) {
  return ((Highs*)highs)->deleteCols(num_set_entries, set);
}

HighsInt Highs_deleteColsByMask(void* highs, HighsInt* mask) {
  return ((Highs*)highs)->deleteCols(mask);
}

HighsInt Highs_deleteRowsByRange(void* highs, const HighsInt from_row,
                                 const HighsInt to_row) {
  return ((Highs*)highs)->deleteRows(from_row, to_row);
}

HighsInt Highs_deleteRowsBySet(void* highs, const HighsInt num_set_entries,
                               const HighsInt* set) {
  return ((Highs*)highs)->deleteRows(num_set_entries, set);
}

HighsInt Highs_deleteRowsByMask(void* highs, HighsInt* mask) {
  return ((Highs*)highs)->deleteRows(mask);
}

double Highs_getHighsInfinity(void* highs) {
  return ((Highs*)highs)->getHighsInfinity();
}

HighsInt Highs_getNumCols(void* highs) {
  return ((Highs*)highs)->getLp().numCol_;
}

HighsInt Highs_getNumRows(void* highs) {
  return ((Highs*)highs)->getLp().numRow_;
}

HighsInt Highs_getNumNz(void* highs) {
  HighsInt numCol = Highs_getNumCols(highs);
  if (numCol <= 0) return 0;
  return ((Highs*)highs)->getLp().Astart_[numCol];
}

const char* Highs_highsModelStatusToChar(void* highs,
                                         HighsInt int_highs_model_status) {
  const char* illegal_highs_model_status = "Model status out of range";
  if (int_highs_model_status <
          (HighsInt)HighsModelStatus::HIGHS_MODEL_STATUS_MIN ||
      int_highs_model_status >
          (HighsInt)HighsModelStatus::HIGHS_MODEL_STATUS_MAX)
    return illegal_highs_model_status;
  return ((Highs*)highs)
      ->highsModelStatusToString(
          static_cast<HighsModelStatus>(int_highs_model_status))
      .c_str();
}

const char* Highs_primalDualStatusToChar(void* highs,
                                         HighsInt int_primal_dual_status) {
  const char* illegal_primal_dual_status = "Primal/Dual status out of range";
  if (int_primal_dual_status < PrimalDualStatus::STATUS_MIN ||
      int_primal_dual_status > PrimalDualStatus::STATUS_MAX)
    return illegal_primal_dual_status;
  return ((Highs*)highs)
      ->primalDualStatusToString(int_primal_dual_status)
      .c_str();
}
