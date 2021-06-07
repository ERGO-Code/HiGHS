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
#include "highs_c_api.h"

#include "Highs.h"

HighsInt Highs_lpCall(const HighsInt numcol, const HighsInt numrow,
                      const HighsInt numnz, const HighsInt rowwise,
                      const HighsInt sense, const double offset,
                      const double* colcost, const double* collower,
                      const double* colupper, const double* rowlower,
                      const double* rowupper, const HighsInt* astart,
                      const HighsInt* aindex, const double* avalue,
                      double* colvalue, double* coldual, double* rowvalue,
                      double* rowdual, HighsInt* colbasisstatus,
                      HighsInt* rowbasisstatus, HighsInt* modelstatus) {
  Highs highs;
  highs.setOptionValue("output_flag", false);
  HighsStatus status = highs.passModel(
      numcol, numrow, numnz, (bool)rowwise, sense, offset, colcost, collower,
      colupper, rowlower, rowupper, astart, aindex, avalue);
  if (status != HighsStatus::kOk) return (HighsInt)status;

  status = highs.run();

  if (status == HighsStatus::kOk) {
    HighsSolution solution;
    HighsBasis basis;
    solution = highs.getSolution();
    basis = highs.getBasis();
    *modelstatus = (HighsInt)highs.getModelStatus();
    const HighsInfo& info = highs.getInfo();
    const bool has_value =
        info.primal_solution_status != SolutionStatus::kSolutionStatusNone;
    const bool has_dual =
        info.dual_solution_status != SolutionStatus::kSolutionStatusNone;
    const bool has_basis = basis.valid;

    for (HighsInt i = 0; i < numcol; i++) {
      if (has_value) colvalue[i] = solution.col_value[i];
      if (has_dual) coldual[i] = solution.col_dual[i];
      if (has_basis) colbasisstatus[i] = (HighsInt)basis.col_status[i];
    }

    for (HighsInt i = 0; i < numrow; i++) {
      if (has_value) rowvalue[i] = solution.row_value[i];
      if (has_dual) rowdual[i] = solution.row_dual[i];
      if (has_basis) rowbasisstatus[i] = (HighsInt)basis.row_status[i];
    }
  }

  return (HighsInt)status;
}

HighsInt Highs_mipCall(const HighsInt numcol, const HighsInt numrow,
                       const HighsInt numnz, const HighsInt rowwise,
                       const HighsInt sense, const double offset,
                       const double* colcost, const double* collower,
                       const double* colupper, const double* rowlower,
                       const double* rowupper, const HighsInt* astart,
                       const HighsInt* aindex, const double* avalue,
                       const HighsInt* integrality, double* colvalue,
                       double* rowvalue, HighsInt* modelstatus) {
  Highs highs;
  highs.setOptionValue("output_flag", false);
  HighsStatus status = highs.passModel(
      numcol, numrow, numnz, (bool)rowwise, sense, offset, colcost, collower,
      colupper, rowlower, rowupper, astart, aindex, avalue, integrality);
  if (status != HighsStatus::kOk) return (HighsInt)status;

  status = highs.run();

  if (status == HighsStatus::kOk) {
    HighsSolution solution;
    solution = highs.getSolution();
    *modelstatus = (HighsInt)highs.getModelStatus();
    const bool has_value = highs.getInfo().primal_solution_status !=
                           SolutionStatus::kSolutionStatusNone;

    if (has_value) {
      for (HighsInt i = 0; i < numcol; i++) colvalue[i] = solution.col_value[i];
      for (HighsInt i = 0; i < numrow; i++) rowvalue[i] = solution.row_value[i];
    }
  }

  return (HighsInt)status;
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

HighsInt Highs_writeSolutionPretty(void* highs, const char* filename) {
  return (HighsInt)((Highs*)highs)->writeSolution(std::string(filename), true);
}

HighsInt Highs_passLp(void* highs, const HighsInt numcol, const HighsInt numrow,
                      const HighsInt numnz, const HighsInt rowwise,
                      const HighsInt sense, const double offset,
                      const double* colcost, const double* collower,
                      const double* colupper, const double* rowlower,
                      const double* rowupper, const HighsInt* astart,
                      const HighsInt* aindex, const double* avalue) {
  const bool bool_rowwise = rowwise;
  return (HighsInt)((Highs*)highs)
      ->passModel(numcol, numrow, numnz, bool_rowwise, sense, offset, colcost,
                  collower, colupper, rowlower, rowupper, astart, aindex,
                  avalue);
}

HighsInt Highs_passMip(void* highs, const HighsInt numcol,
                       const HighsInt numrow, const HighsInt numnz,
                       const HighsInt rowwise, const HighsInt sense,
                       const double offset, const double* colcost,
                       const double* collower, const double* colupper,
                       const double* rowlower, const double* rowupper,
                       const HighsInt* astart, const HighsInt* aindex,
                       const double* avalue, const HighsInt* integrality) {
  const bool bool_rowwise = rowwise;
  return (HighsInt)((Highs*)highs)
      ->passModel(numcol, numrow, numnz, bool_rowwise, sense, offset, colcost,
                  collower, colupper, rowlower, rowupper, astart, aindex,
                  avalue, integrality);
}

HighsInt Highs_passModel(void* highs, const HighsInt numcol,
                         const HighsInt numrow, const HighsInt numnz,
                         const HighsInt hessian_num_nz, const HighsInt rowwise,
                         const HighsInt sense, const double offset,
                         const double* colcost, const double* collower,
                         const double* colupper, const double* rowlower,
                         const double* rowupper, const HighsInt* astart,
                         const HighsInt* aindex, const double* avalue,
                         const HighsInt* qstart, const HighsInt* qindex,
                         const double* qvalue, const HighsInt* integrality) {
  const bool bool_rowwise = rowwise;
  return (HighsInt)((Highs*)highs)
      ->passModel(numcol, numrow, numnz, hessian_num_nz, bool_rowwise, sense,
                  offset, colcost, collower, colupper, rowlower, rowupper,
                  astart, aindex, avalue, qstart, qindex, qvalue, integrality);
}

HighsInt Highs_clearModel(void* highs) {
  return (HighsInt)((Highs*)highs)->clearModel();
}

HighsInt Highs_setBoolOptionValue(void* highs, const char* option,
                                  const HighsInt value) {
  return (HighsInt)((Highs*)highs)
      ->setOptionValue(std::string(option), (bool)value);
}

HighsInt Highs_setIntOptionValue(void* highs, const char* option,
                                 const HighsInt value) {
  return (HighsInt)((Highs*)highs)->setOptionValue(std::string(option), value);
}

HighsInt Highs_setDoubleOptionValue(void* highs, const char* option,
                                    const double value) {
  return (HighsInt)((Highs*)highs)->setOptionValue(std::string(option), value);
}

HighsInt Highs_setStringOptionValue(void* highs, const char* option,
                                    const char* value) {
  return (HighsInt)((Highs*)highs)
      ->setOptionValue(std::string(option), std::string(value));
}

HighsInt Highs_setOptionValue(void* highs, const char* option,
                              const char* value) {
  return (HighsInt)((Highs*)highs)
      ->setOptionValue(std::string(option), std::string(value));
}

HighsInt Highs_getBoolOptionValue(void* highs, const char* option,
                                  HighsInt* value) {
  bool v;
  HighsInt retcode =
      (HighsInt)((Highs*)highs)->getOptionValue(std::string(option), v);
  *value = (HighsInt)v;
  return retcode;
}

HighsInt Highs_getIntOptionValue(void* highs, const char* option,
                                 HighsInt* value) {
  return (HighsInt)((Highs*)highs)->getOptionValue(std::string(option), *value);
}

HighsInt Highs_getDoubleOptionValue(void* highs, const char* option,
                                    double* value) {
  return (HighsInt)((Highs*)highs)->getOptionValue(std::string(option), *value);
}

HighsInt Highs_getStringOptionValue(void* highs, const char* option,
                                    char* value) {
  std::string v;
  memset(value, 0, 7);
  HighsInt retcode =
      (HighsInt)((Highs*)highs)->getOptionValue(std::string(option), v);
  strcpy(value, v.c_str());
  return retcode;
}

HighsInt Highs_getOptionType(void* highs, const char* option, HighsInt* type) {
  HighsOptionType t;
  HighsInt retcode =
      (HighsInt)((Highs*)highs)->getOptionType(std::string(option), t);
  *type = (HighsInt)t;
  return retcode;
}

HighsInt Highs_resetOptions(void* highs) {
  return (HighsInt)((Highs*)highs)->resetOptions();
}

HighsInt Highs_getIntInfoValue(void* highs, const char* info, HighsInt* value) {
  return (HighsInt)((Highs*)highs)->getInfoValue(info, *value);
}

HighsInt Highs_getDoubleInfoValue(void* highs, const char* info,
                                  double* value) {
  return (HighsInt)((Highs*)highs)->getInfoValue(info, *value);
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

HighsInt Highs_getModelStatus(void* highs) {
  return (HighsInt)((Highs*)highs)->getModelStatus();
}

HighsInt Highs_getScaledModelStatus(void* highs) {
  return (HighsInt)((Highs*)highs)->getModelStatus(true);
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
  if (num_col > 0) {
    basis.col_status.resize(num_col);
    for (HighsInt i = 0; i < num_col; i++) {
      if (colstatus[i] == (HighsInt)HighsBasisStatus::kLower) {
        basis.col_status[i] = HighsBasisStatus::kLower;
      } else if (colstatus[i] == (HighsInt)HighsBasisStatus::kBasic) {
        basis.col_status[i] = HighsBasisStatus::kBasic;
      } else if (colstatus[i] == (HighsInt)HighsBasisStatus::kUpper) {
        basis.col_status[i] = HighsBasisStatus::kUpper;
      } else if (colstatus[i] == (HighsInt)HighsBasisStatus::kZero) {
        basis.col_status[i] = HighsBasisStatus::kZero;
      } else if (colstatus[i] == (HighsInt)HighsBasisStatus::kNonbasic) {
        basis.col_status[i] = HighsBasisStatus::kNonbasic;
      } else {
        return (HighsInt)HighsStatus::kError;
      }
    }
  }
  const HighsInt num_row = Highs_getNumRows(highs);
  if (num_row > 0) {
    basis.row_status.resize(num_row);
    for (HighsInt i = 0; i < num_row; i++) {
      if (rowstatus[i] == (HighsInt)HighsBasisStatus::kLower) {
        basis.row_status[i] = HighsBasisStatus::kLower;
      } else if (rowstatus[i] == (HighsInt)HighsBasisStatus::kBasic) {
        basis.row_status[i] = HighsBasisStatus::kBasic;
      } else if (rowstatus[i] == (HighsInt)HighsBasisStatus::kUpper) {
        basis.row_status[i] = HighsBasisStatus::kUpper;
      } else if (rowstatus[i] == (HighsInt)HighsBasisStatus::kZero) {
        basis.row_status[i] = HighsBasisStatus::kZero;
      } else if (rowstatus[i] == (HighsInt)HighsBasisStatus::kNonbasic) {
        basis.row_status[i] = HighsBasisStatus::kNonbasic;
      } else {
        return (HighsInt)HighsStatus::kError;
      }
    }
  }
  return (HighsInt)((Highs*)highs)->setBasis(basis);
}

HighsInt Highs_setLogicalBasis(void* highs) {
  return (HighsInt)((Highs*)highs)->setBasis();
}

double Highs_getRunTime(void* highs) {
  return (double)((Highs*)highs)->getRunTime();
}

HighsInt Highs_addRow(void* highs, const double lower, const double upper,
                      const HighsInt num_new_nz, const HighsInt* indices,
                      const double* values) {
  return (HighsInt)((Highs*)highs)
      ->addRow(lower, upper, num_new_nz, indices, values);
}

HighsInt Highs_addRows(void* highs, const HighsInt num_new_row,
                       const double* lower, const double* upper,
                       const HighsInt num_new_nz, const HighsInt* starts,
                       const HighsInt* indices, const double* values) {
  return (HighsInt)((Highs*)highs)
      ->addRows(num_new_row, lower, upper, num_new_nz, starts, indices, values);
}

HighsInt Highs_addCol(void* highs, const double cost, const double lower,
                      const double upper, const HighsInt num_new_nz,
                      const HighsInt* indices, const double* values) {
  return (HighsInt)((Highs*)highs)
      ->addCol(cost, lower, upper, num_new_nz, indices, values);
}

HighsInt Highs_addCols(void* highs, const HighsInt num_new_col,
                       const double* costs, const double* lower,
                       const double* upper, const HighsInt num_new_nz,
                       const HighsInt* starts, const HighsInt* indices,
                       const double* values) {
  return (HighsInt)((Highs*)highs)
      ->addCols(num_new_col, costs, lower, upper, num_new_nz, starts, indices,
                values);
}

HighsInt Highs_changeObjectiveSense(void* highs, const HighsInt sense) {
  ObjSense pass_sense = ObjSense::kMinimize;
  if (sense == (HighsInt)ObjSense::kMaximize) pass_sense = ObjSense::kMaximize;
  return (HighsInt)((Highs*)highs)->changeObjectiveSense(pass_sense);
}

HighsInt Highs_changeColIntegrality(void* highs, const HighsInt col,
                                    const HighsInt integrality) {
  return (HighsInt)((Highs*)highs)
      ->changeColIntegrality(col, (HighsVarType)integrality);
}

HighsInt Highs_changeColsIntegralityByRange(void* highs,
                                            const HighsInt from_col,
                                            const HighsInt to_col,
                                            const HighsInt* integrality) {
  vector<HighsVarType> pass_integrality;
  HighsInt num_ix = to_col - from_col + 1;
  if (num_ix > 0) {
    pass_integrality.resize(num_ix);
    for (HighsInt ix = 0; ix < num_ix; ix++) {
      pass_integrality[ix] = (HighsVarType)integrality[ix];
    }
  }
  return (HighsInt)((Highs*)highs)
      ->changeColsIntegrality(from_col, to_col, &pass_integrality[0]);
}

HighsInt Highs_changeColsIntegralityBySet(void* highs,
                                          const HighsInt num_set_entries,
                                          const HighsInt* set,
                                          const HighsInt* integrality) {
  vector<HighsVarType> pass_integrality;
  if (num_set_entries > 0) {
    pass_integrality.resize(num_set_entries);
    for (HighsInt ix = 0; ix < num_set_entries; ix++) {
      pass_integrality[ix] = (HighsVarType)integrality[ix];
    }
  }
  return (HighsInt)((Highs*)highs)
      ->changeColsIntegrality(num_set_entries, set, &pass_integrality[0]);
}

HighsInt Highs_changeColsIntegralityByMask(void* highs, const HighsInt* mask,
                                           const HighsInt* integrality) {
  const HighsInt num_col = Highs_getNumCols(highs);
  vector<HighsVarType> pass_integrality;
  if (num_col > 0) {
    pass_integrality.resize(num_col);
    for (HighsInt iCol = 0; iCol < num_col; iCol++) {
      pass_integrality[iCol] = (HighsVarType)integrality[iCol];
    }
  }
  return (HighsInt)((Highs*)highs)
      ->changeColsIntegrality(mask, &pass_integrality[0]);
}

HighsInt Highs_changeColCost(void* highs, const HighsInt col,
                             const double cost) {
  return (HighsInt)((Highs*)highs)->changeColCost(col, cost);
}

HighsInt Highs_changeColsCostByRange(void* highs, const HighsInt from_col,
                                     const HighsInt to_col,
                                     const double* cost) {
  return (HighsInt)((Highs*)highs)->changeColsCost(from_col, to_col, cost);
}

HighsInt Highs_changeColsCostBySet(void* highs, const HighsInt num_set_entries,
                                   const HighsInt* set, const double* cost) {
  return (HighsInt)((Highs*)highs)->changeColsCost(num_set_entries, set, cost);
}

HighsInt Highs_changeColsCostByMask(void* highs, const HighsInt* mask,
                                    const double* cost) {
  return (HighsInt)((Highs*)highs)->changeColsCost(mask, cost);
}

HighsInt Highs_changeColBounds(void* highs, const HighsInt col,
                               const double lower, const double upper) {
  return (HighsInt)((Highs*)highs)->changeColBounds(col, lower, upper);
}

HighsInt Highs_changeColsBoundsByRange(void* highs, const HighsInt from_col,
                                       const HighsInt to_col,
                                       const double* lower,
                                       const double* upper) {
  return (HighsInt)((Highs*)highs)
      ->changeColsBounds(from_col, to_col, lower, upper);
}

HighsInt Highs_changeColsBoundsBySet(void* highs,
                                     const HighsInt num_set_entries,
                                     const HighsInt* set, const double* lower,
                                     const double* upper) {
  return (HighsInt)((Highs*)highs)
      ->changeColsBounds(num_set_entries, set, lower, upper);
}

HighsInt Highs_changeColsBoundsByMask(void* highs, const HighsInt* mask,
                                      const double* lower,
                                      const double* upper) {
  return (HighsInt)((Highs*)highs)->changeColsBounds(mask, lower, upper);
}

HighsInt Highs_changeRowBounds(void* highs, const HighsInt row,
                               const double lower, const double upper) {
  return (HighsInt)((Highs*)highs)->changeRowBounds(row, lower, upper);
}

HighsInt Highs_changeRowsBoundsBySet(void* highs,
                                     const HighsInt num_set_entries,
                                     const HighsInt* set, const double* lower,
                                     const double* upper) {
  return (HighsInt)((Highs*)highs)
      ->changeRowsBounds(num_set_entries, set, lower, upper);
}

HighsInt Highs_changeRowsBoundsByMask(void* highs, const HighsInt* mask,
                                      const double* lower,
                                      const double* upper) {
  return (HighsInt)((Highs*)highs)->changeRowsBounds(mask, lower, upper);
}

HighsInt Highs_changeCoeff(void* highs, const HighsInt row, const HighsInt col,
                           const double value) {
  return (HighsInt)((Highs*)highs)->changeCoeff(row, col, value);
}

HighsInt Highs_getObjectiveSense(void* highs, HighsInt* sense) {
  ObjSense get_sense;
  HighsStatus status = ((Highs*)highs)->getObjectiveSense(get_sense);
  *sense = (HighsInt)get_sense;
  return (HighsInt)status;
}

HighsInt Highs_getColsByRange(void* highs, const HighsInt from_col,
                              const HighsInt to_col, HighsInt* num_col,
                              double* costs, double* lower, double* upper,
                              HighsInt* num_nz, HighsInt* matrix_start,
                              HighsInt* matrix_index, double* matrix_value) {
  HighsInt numcol, numnz;
  HighsStatus status =
      ((Highs*)highs)
          ->getCols(from_col, to_col, numcol, costs, lower, upper, numnz,
                    matrix_start, matrix_index, matrix_value);
  *num_col = numcol;
  *num_nz = numnz;
  return (HighsInt)status;
}

HighsInt Highs_getColsBySet(void* highs, const HighsInt num_set_entries,
                            const HighsInt* set, HighsInt* num_col,
                            double* costs, double* lower, double* upper,
                            HighsInt* num_nz, HighsInt* matrix_start,
                            HighsInt* matrix_index, double* matrix_value) {
  HighsInt numcol, numnz;
  HighsStatus status =
      ((Highs*)highs)
          ->getCols(num_set_entries, set, numcol, costs, lower, upper, numnz,
                    matrix_start, matrix_index, matrix_value);
  *num_col = numcol;
  *num_nz = numnz;
  return (HighsInt)status;
}

HighsInt Highs_getColsByMask(void* highs, const HighsInt* mask,
                             HighsInt* num_col, double* costs, double* lower,
                             double* upper, HighsInt* num_nz,
                             HighsInt* matrix_start, HighsInt* matrix_index,
                             double* matrix_value) {
  HighsInt numcol, numnz;
  HighsStatus status = ((Highs*)highs)
                           ->getCols(mask, numcol, costs, lower, upper, numnz,
                                     matrix_start, matrix_index, matrix_value);
  *num_col = numcol;
  *num_nz = numnz;
  return (HighsInt)status;
}

HighsInt Highs_getRowsByRange(void* highs, const HighsInt from_row,
                              const HighsInt to_row, HighsInt* num_row,
                              double* lower, double* upper, HighsInt* num_nz,
                              HighsInt* matrix_start, HighsInt* matrix_index,
                              double* matrix_value) {
  HighsInt numrow, numnz;
  HighsStatus status =
      ((Highs*)highs)
          ->getRows(from_row, to_row, numrow, lower, upper, numnz, matrix_start,
                    matrix_index, matrix_value);
  *num_row = numrow;
  *num_nz = numnz;
  return (HighsInt)status;
}

HighsInt Highs_getRowsBySet(void* highs, const HighsInt num_set_entries,
                            const HighsInt* set, HighsInt* num_row,
                            double* lower, double* upper, HighsInt* num_nz,
                            HighsInt* matrix_start, HighsInt* matrix_index,
                            double* matrix_value) {
  HighsInt numrow, numnz;
  HighsStatus status =
      ((Highs*)highs)
          ->getRows(num_set_entries, set, numrow, lower, upper, numnz,
                    matrix_start, matrix_index, matrix_value);
  *num_row = numrow;
  *num_nz = numnz;
  return (HighsInt)status;
}

HighsInt Highs_getRowsByMask(void* highs, const HighsInt* mask,
                             HighsInt* num_row, double* lower, double* upper,
                             HighsInt* num_nz, HighsInt* matrix_start,
                             HighsInt* matrix_index, double* matrix_value) {
  HighsInt numrow, numnz;
  HighsStatus status = ((Highs*)highs)
                           ->getRows(mask, numrow, lower, upper, numnz,
                                     matrix_start, matrix_index, matrix_value);
  *num_row = numrow;
  *num_nz = numnz;
  return (HighsInt)status;
}

HighsInt Highs_deleteColsByRange(void* highs, const HighsInt from_col,
                                 const HighsInt to_col) {
  return (HighsInt)((Highs*)highs)->deleteCols(from_col, to_col);
}

HighsInt Highs_deleteColsBySet(void* highs, const HighsInt num_set_entries,
                               const HighsInt* set) {
  return (HighsInt)((Highs*)highs)->deleteCols(num_set_entries, set);
}

HighsInt Highs_deleteColsByMask(void* highs, HighsInt* mask) {
  return (HighsInt)((Highs*)highs)->deleteCols(mask);
}

HighsInt Highs_deleteRowsByRange(void* highs, const HighsInt from_row,
                                 const HighsInt to_row) {
  return (HighsInt)((Highs*)highs)->deleteRows(from_row, to_row);
}

HighsInt Highs_deleteRowsBySet(void* highs, const HighsInt num_set_entries,
                               const HighsInt* set) {
  return (HighsInt)((Highs*)highs)->deleteRows(num_set_entries, set);
}

HighsInt Highs_deleteRowsByMask(void* highs, HighsInt* mask) {
  return (HighsInt)((Highs*)highs)->deleteRows(mask);
}

HighsInt Highs_scaleCol(void* highs, const HighsInt col,
                        const double scaleval) {
  return (HighsInt)((Highs*)highs)->scaleCol(col, scaleval);
}

HighsInt Highs_scaleRow(void* highs, const HighsInt row,
                        const double scaleval) {
  return (HighsInt)((Highs*)highs)->scaleRow(row, scaleval);
}

double Highs_getInfinity(void* highs) { return ((Highs*)highs)->getInfinity(); }

HighsInt Highs_getNumCols(void* highs) { return ((Highs*)highs)->getNumCols(); }

HighsInt Highs_getNumRows(void* highs) { return ((Highs*)highs)->getNumRows(); }

HighsInt Highs_getNumNz(void* highs) { return ((Highs*)highs)->getNumNz(); }

HighsInt Highs_getHessianNumNz(void* highs) {
  return ((Highs*)highs)->getHessianNumNz();
}

void Highs_getModel(void* highs, const HighsInt orientation, HighsInt* numcol,
                    HighsInt* numrow, HighsInt* numnz, HighsInt* hessian_num_nz,
                    HighsInt* sense, double* offset, double* colcost,
                    double* collower, double* colupper, double* rowlower,
                    double* rowupper, HighsInt* astart, HighsInt* aindex,
                    double* avalue, HighsInt* qstart, HighsInt* qindex,
                    double* qvalue, HighsInt* integrality) {
  const HighsModel& model = ((Highs*)highs)->getModel();
  const HighsLp& lp = model.lp_;
  const HighsHessian& hessian = model.hessian_;
  ObjSense obj_sense = ObjSense::kMinimize;
  *sense = (HighsInt)obj_sense;
  *offset = lp.offset_;
  *numcol = lp.numCol_;
  *numrow = lp.numRow_;
  if (*numcol > 0) {
    memcpy(colcost, &lp.colCost_[0], *numcol * sizeof(double));
    memcpy(collower, &lp.colLower_[0], *numcol * sizeof(double));
    memcpy(colupper, &lp.colUpper_[0], *numcol * sizeof(double));
  }
  if (*numrow > 0) {
    memcpy(rowlower, &lp.rowLower_[0], *numrow * sizeof(double));
    memcpy(rowupper, &lp.rowUpper_[0], *numrow * sizeof(double));
  }

  // Save the original orientation so that it is recovered
  MatrixOrientation original_orientation = lp.orientation_;
  // Determine the desired orientation and number of start entries to
  // be copied
  MatrixOrientation desired_orientation = MatrixOrientation::kColwise;
  HighsInt num_start_entries = *numcol;
  if (orientation == (HighsInt)MatrixOrientation::kRowwise) {
    desired_orientation = MatrixOrientation::kRowwise;
    num_start_entries = *numrow;
  }
  // Ensure the desired orientation
  ((Highs*)highs)->setMatrixOrientation(desired_orientation);

  if (*numcol > 0 && *numrow > 0) {
    memcpy(astart, &lp.Astart_[0], num_start_entries * sizeof(HighsInt));
    *numnz = lp.Astart_[*numcol];
    memcpy(aindex, &lp.Aindex_[0], *numnz * sizeof(HighsInt));
    memcpy(avalue, &lp.Avalue_[0], *numnz * sizeof(double));
  }
  if (hessian.dim_ > 0) {
    memcpy(qstart, &hessian.q_start_[0], *numcol * sizeof(HighsInt));
    *hessian_num_nz = hessian.q_start_[*numcol];
    memcpy(qindex, &hessian.q_index_[0], *hessian_num_nz * sizeof(HighsInt));
    memcpy(qvalue, &hessian.q_value_[0], *hessian_num_nz * sizeof(double));
  }
  if ((HighsInt)lp.integrality_.size()) {
    for (int iCol = 0; iCol < *numcol; iCol++)
      integrality[iCol] = (HighsInt)lp.integrality_[iCol];
  }
  // Restore the original orientation
  ((Highs*)highs)->setMatrixOrientation(original_orientation);
}

// Fails on Windows and MacOS since string_model_status is destroyed
// after the method returns, so what's returned is a pointer to
// something that no longer exists.
//
// const char* Highs_modelStatusToChar(void* highs, HighsInt int_model_status) {
//   const char* illegal_model_status = "Model status out of range";
//   if (int_model_status < (HighsInt)HighsModelStatus::kMin ||
//       int_model_status > (HighsInt)HighsModelStatus::kMax)
//     return illegal_model_status;
//   const std::string string_model_status =
//       ((Highs*)highs)
//           ->modelStatusToString(
//               static_cast<HighsModelStatus>(int_model_status));
//   return string_model_status.c_str();
// }

// Fails on Windows and MacOS since string_solution_status is
// destroyed after the method returns, so what's returned is a pointer
// to something that no longer exists.
//
// const char* Highs_solutionStatusToChar(void* highs,
//                                        HighsInt int_solution_status) {
//   const char* illegal_solution_status = "Solution status out of range";
//   if (int_solution_status < kSolutionStatusMin ||
//       int_solution_status > kSolutionStatusMax)
//     return illegal_solution_status;
//   const std::string string_solution_status =
//       ((Highs*)highs)->solutionStatusToString(int_solution_status);
//   return string_solution_status.c_str();
// }

// *********************
// * Deprecated methods*
// *********************

HighsInt Highs_call(const HighsInt numcol, const HighsInt numrow,
                    const HighsInt numnz, const double* colcost,
                    const double* collower, const double* colupper,
                    const double* rowlower, const double* rowupper,
                    const HighsInt* astart, const HighsInt* aindex,
                    const double* avalue, double* colvalue, double* coldual,
                    double* rowvalue, double* rowdual, HighsInt* colbasisstatus,
                    HighsInt* rowbasisstatus, HighsInt* modelstatus) {
  printf(
      "Method Highs_call is deprecated: alternative method is Highs_lpCall\n");
  const HighsInt rowwise = 0;
  const HighsInt sense = 1;
  const double offset = 0;
  return Highs_lpCall(numcol, numrow, numnz, rowwise, sense, offset, colcost,
                      collower, colupper, rowlower, rowupper, astart, aindex,
                      avalue, colvalue, coldual, rowvalue, rowdual,
                      colbasisstatus, rowbasisstatus, modelstatus);
}

HighsInt Highs_runQuiet(void* highs) {
  ((Highs*)highs)->deprecationMessage("Highs_runQuiet", "None");
  return (HighsInt)((Highs*)highs)->setOptionValue("output_flag", false);
}

HighsInt Highs_setHighsLogfile(void* highs, void* logfile) {
  ((Highs*)highs)->deprecationMessage("Highs_setHighsLogfile", "None");
  return (HighsInt)((Highs*)highs)->setOptionValue("output_flag", false);
}

HighsInt Highs_setHighsOutput(void* highs, void* outputfile) {
  ((Highs*)highs)->deprecationMessage("Highs_setHighsOutput", "None");
  return (HighsInt)((Highs*)highs)->setOptionValue("output_flag", false);
}

HighsInt Highs_getIterationCount(void* highs) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_getIterationCount", "Highs_getIntInfoValue");
  return (HighsInt)((Highs*)highs)->getInfo().simplex_iteration_count;
}

HighsInt Highs_getSimplexIterationCount(void* highs) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_getSimplexIterationCount",
                           "Highs_getIntInfoValue");
  return (HighsInt)((Highs*)highs)->getInfo().simplex_iteration_count;
}

HighsInt Highs_setHighsBoolOptionValue(void* highs, const char* option,
                                       const HighsInt value) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_setHighsBoolOptionValue",
                           "Highs_setBoolOptionValue");
  return Highs_setBoolOptionValue(highs, option, value);
}

HighsInt Highs_setHighsIntOptionValue(void* highs, const char* option,
                                      const HighsInt value) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_setHighsIntOptionValue",
                           "Highs_setIntOptionValue");
  return Highs_setIntOptionValue(highs, option, value);
}

HighsInt Highs_setHighsDoubleOptionValue(void* highs, const char* option,
                                         const double value) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_setHighsDoubleOptionValue",
                           "Highs_setDoubleOptionValue");
  return Highs_setDoubleOptionValue(highs, option, value);
}

HighsInt Highs_setHighsStringOptionValue(void* highs, const char* option,
                                         const char* value) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_setHighsStringOptionValue",
                           "Highs_setStringOptionValue");
  return Highs_setStringOptionValue(highs, option, value);
}

HighsInt Highs_setHighsOptionValue(void* highs, const char* option,
                                   const char* value) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_setHighsOptionValue", "Highs_setOptionValue");
  return Highs_setOptionValue(highs, option, value);
}

HighsInt Highs_getHighsBoolOptionValue(void* highs, const char* option,
                                       HighsInt* value) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_getHighsBoolOptionValue",
                           "Highs_getBoolOptionValue");
  return Highs_getBoolOptionValue(highs, option, value);
}

HighsInt Highs_getHighsIntOptionValue(void* highs, const char* option,
                                      HighsInt* value) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_getHighsIntOptionValue",
                           "Highs_getIntOptionValue");
  return Highs_getIntOptionValue(highs, option, value);
}

HighsInt Highs_getHighsDoubleOptionValue(void* highs, const char* option,
                                         double* value) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_getHighsDoubleOptionValue",
                           "Highs_getDoubleOptionValue");
  return Highs_getDoubleOptionValue(highs, option, value);
}

HighsInt Highs_getHighsStringOptionValue(void* highs, const char* option,
                                         char* value) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_getHighsStringOptionValue",
                           "Highs_getStringOptionValue");
  return Highs_getStringOptionValue(highs, option, value);
}

HighsInt Highs_getHighsOptionType(void* highs, const char* option,
                                  HighsInt* type) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_getHighsOptionType", "Highs_getOptionType");
  return Highs_getOptionType(highs, option, type);
}

HighsInt Highs_resetHighsOptions(void* highs) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_resetHighsOptions", "Highs_resetOptions");
  return Highs_resetOptions(highs);
}

HighsInt Highs_getHighsIntInfoValue(void* highs, const char* info,
                                    HighsInt* value) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_getHighsIntInfoValue",
                           "Highs_getIntInfoValue");
  return Highs_getIntInfoValue(highs, info, value);
}

HighsInt Highs_getHighsDoubleInfoValue(void* highs, const char* info,
                                       double* value) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_getHighsDoubleInfoValue",
                           "Highs_getDoubleInfoValue");
  return Highs_getDoubleInfoValue(highs, info, value);
}

double Highs_getHighsRunTime(void* highs) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_getHighsRunTime", "Highs_getRunTime");
  return Highs_getRunTime(highs);
}

double Highs_getHighsInfinity(void* highs) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_getHighsInfinity", "Highs_getInfinity");
  return Highs_getInfinity(highs);
}

// const char* Highs_highsModelStatusToChar(void* highs,
//                                          HighsInt int_model_status) {
//   return Highs_modelStatusToChar(highs, int_model_status);
// }
