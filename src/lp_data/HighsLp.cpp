/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsLp.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "lp_data/HighsLp.h"
#include "lp_data/HConst.h"
#include "io/HighsIO.h" // For HighsModelStatusToString and HighsModelStatusReport

bool isSolutionConsistent(const HighsLp& lp, const HighsSolution& solution) {
  if (solution.col_dual.size() == (size_t)lp.numCol_ ||
      solution.col_value.size() == (size_t)lp.numCol_ ||
      solution.row_dual.size() == (size_t)lp.numRow_ ||
      solution.row_value.size() == (size_t)lp.numRow_)
    return true;
  return false;
}

// Return a string representation of HighsModelStatus.
std::string HighsModelStatusToString(HighsModelStatus model_status) {

  switch (model_status) {
  case HighsModelStatus::NOTSET:
      return "Not Set";
      break;
  case HighsModelStatus::LOAD_ERROR:
      return "Load error";
      break;
  case HighsModelStatus::MODEL_ERROR:
      return "Model error";
      break;
  case HighsModelStatus::MODEL_EMPTY:
      return "Model empty";
      break;
  case HighsModelStatus::PRESOLVE_ERROR:
      return "Presolve error";
      break;
  case HighsModelStatus::SOLVE_ERROR:
      return "Solve error";
      break;
  case HighsModelStatus::POSTSOLVE_ERROR:
      return "Postsolve error";
      break;
  case HighsModelStatus::PRIMAL_FEASIBLE:
      return "Primal feasible";
      break;
  case HighsModelStatus::DUAL_FEASIBLE:
      return "Dual feasible";
      break;
  case HighsModelStatus::PRIMAL_INFEASIBLE:
      return "Primal infeasible";
      break;
  case HighsModelStatus::PRIMAL_UNBOUNDED:
      return "Primal unbounded";
      break;
  case HighsModelStatus::OPTIMAL:
      return "Optimal";
      break;
  case HighsModelStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND:
      return "Reached dual objective upper bound";
      break;
  case HighsModelStatus::REACHED_TIME_LIMIT:
      return "Reached time limit";
      break;
  case HighsModelStatus::REACHED_ITERATION_LIMIT:
      return "Reached iteration limit";
      break;
    default:
      return "Status toString() not implemented.";
      break;
  }
  return "";
}

// Report a HighsModelStatus.
void HighsModelStatusReport(const char* message, HighsModelStatus model_status) {
  HighsLogMessage(HighsMessageType::INFO, "%s: HighsModelStatus = %d - %s\n",
                  message, (int)model_status, HighsModelStatusToString(model_status).c_str());
}

