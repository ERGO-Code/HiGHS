/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsUtils.cpp
 * @brief Class-independent utilities for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#include "lp_data/HighsModelUtils.h"

#include <algorithm>
#include <vector>

#include "HConfig.h"
#include "io/HighsIO.h"
#include "lp_data/HConst.h"
#include "util/HighsUtils.h"

void analyseModelBounds(const HighsLogOptions& log_options, const char* message,
                        HighsInt numBd, const std::vector<double>& lower,
                        const std::vector<double>& upper) {
  if (numBd == 0) return;
  HighsInt numFr = 0;
  HighsInt numLb = 0;
  HighsInt numUb = 0;
  HighsInt numBx = 0;
  HighsInt numFx = 0;
  for (HighsInt ix = 0; ix < numBd; ix++) {
    if (highs_isInfinity(-lower[ix])) {
      // Infinite lower bound
      if (highs_isInfinity(upper[ix])) {
        // Infinite lower bound and infinite upper bound: Fr
        numFr++;
      } else {
        // Infinite lower bound and   finite upper bound: Ub
        numUb++;
      }
    } else {
      // Finite lower bound
      if (highs_isInfinity(upper[ix])) {
        // Finite lower bound and infinite upper bound: Lb
        numLb++;
      } else {
        // Finite lower bound and   finite upper bound:
        if (lower[ix] < upper[ix]) {
          // Distinct finite bounds: Bx
          numBx++;
        } else {
          // Equal finite bounds: Fx
          numFx++;
        }
      }
    }
  }
  highsLogDev(log_options, HighsLogType::INFO,
              "Analysing %" HIGHSINT_FORMAT " %s bounds\n", numBd, message);
  if (numFr > 0)
    highsLogDev(log_options, HighsLogType::INFO,
                "   Free:  %7" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT "%%)\n",
                numFr, (100 * numFr) / numBd);
  if (numLb > 0)
    highsLogDev(log_options, HighsLogType::INFO,
                "   LB:    %7" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT "%%)\n",
                numLb, (100 * numLb) / numBd);
  if (numUb > 0)
    highsLogDev(log_options, HighsLogType::INFO,
                "   UB:    %7" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT "%%)\n",
                numUb, (100 * numUb) / numBd);
  if (numBx > 0)
    highsLogDev(log_options, HighsLogType::INFO,
                "   Boxed: %7" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT "%%)\n",
                numBx, (100 * numBx) / numBd);
  if (numFx > 0)
    highsLogDev(log_options, HighsLogType::INFO,
                "   Fixed: %7" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT "%%)\n",
                numFx, (100 * numFx) / numBd);
  highsLogDev(log_options, HighsLogType::INFO,
              "grep_CharMl,%s,Free,LB,UB,Boxed,Fixed\n", message);
  highsLogDev(log_options, HighsLogType::INFO,
              "grep_CharMl,%" HIGHSINT_FORMAT ",%" HIGHSINT_FORMAT
              ",%" HIGHSINT_FORMAT ",%" HIGHSINT_FORMAT ",%" HIGHSINT_FORMAT
              ",%" HIGHSINT_FORMAT "\n",
              numBd, numFr, numLb, numUb, numBx, numFx);
}

std::string statusToString(const HighsBasisStatus status, const double lower,
                           const double upper) {
  switch (status) {
    case HighsBasisStatus::LOWER:
      if (lower == upper) {
        return "FX";
      } else {
        return "LB";
      }
      break;
    case HighsBasisStatus::BASIC:
      return "BS";
      break;
    case HighsBasisStatus::UPPER:
      return "UB";
      break;
    case HighsBasisStatus::ZERO:
      return "FR";
      break;
    case HighsBasisStatus::NONBASIC:
      return "NB";
      break;
  }
  return "";
}

void writeModelBoundSol(FILE* file, const bool columns, const HighsInt dim,
                        const std::vector<double>& lower,
                        const std::vector<double>& upper,
                        const std::vector<std::string>& names,
                        const std::vector<double>& primal,
                        const std::vector<double>& dual,
                        const std::vector<HighsBasisStatus>& status) {
  const bool have_names = names.size() > 0;
  const bool have_basis = status.size() > 0;
  const bool have_primal = primal.size() > 0;
  const bool have_dual = dual.size() > 0;
  std::string var_status_string;
  if (columns) {
    fprintf(file, "Columns\n");
  } else {
    fprintf(file, "Rows\n");
  }
  fprintf(
      file,
      "    Index Status        Lower        Upper       Primal         Dual");
  if (have_names) {
    fprintf(file, "  Name\n");
  } else {
    fprintf(file, "\n");
  }
  for (HighsInt ix = 0; ix < dim; ix++) {
    if (have_basis) {
      var_status_string = statusToString(status[ix], lower[ix], upper[ix]);
    } else {
      var_status_string = "";
    }
    fprintf(file, "%9" HIGHSINT_FORMAT "   %4s %12g %12g", ix,
            var_status_string.c_str(), lower[ix], upper[ix]);
    if (have_primal) {
      fprintf(file, " %12g", primal[ix]);
    } else {
      fprintf(file, "             ");
    }
    if (have_dual) {
      fprintf(file, " %12g", dual[ix]);
    } else {
      fprintf(file, "             ");
    }
    if (have_names) {
      fprintf(file, "  %-s\n", names[ix].c_str());
    } else {
      fprintf(file, "\n");
    }
  }
}

bool namesWithSpaces(const HighsInt num_name,
                     const std::vector<std::string>& names, const bool report) {
  bool names_with_spaces = false;
  for (HighsInt ix = 0; ix < num_name; ix++) {
    HighsInt space_pos = names[ix].find(" ");
    if (space_pos >= 0) {
      if (report)
        printf(
            "Name |%s| contains a space character in position %" HIGHSINT_FORMAT
            "\n",
            names[ix].c_str(), space_pos);
      names_with_spaces = true;
    }
  }
  return names_with_spaces;
}

HighsInt maxNameLength(const HighsInt num_name,
                       const std::vector<std::string>& names) {
  HighsInt max_name_length = 0;
  for (HighsInt ix = 0; ix < num_name; ix++)
    max_name_length = std::max((HighsInt)names[ix].length(), max_name_length);
  return max_name_length;
}

HighsStatus normaliseNames(const HighsLogOptions& log_options,
                           const std::string name_type, const HighsInt num_name,
                           std::vector<std::string>& names,
                           HighsInt& max_name_length) {
  // Record the desired maximum name length
  HighsInt desired_max_name_length = max_name_length;
  // First look for empty names
  HighsInt num_empty_name = 0;
  std::string name_prefix = name_type.substr(0, 1);
  bool names_with_spaces = false;
  for (HighsInt ix = 0; ix < num_name; ix++) {
    if ((HighsInt)names[ix].length() == 0) num_empty_name++;
  }
  // If there are no empty names - in which case they will all be
  // replaced - find the maximum name length
  if (!num_empty_name) max_name_length = maxNameLength(num_name, names);
  bool construct_names =
      num_empty_name || max_name_length > desired_max_name_length;
  if (construct_names) {
    // Construct names, either because they are empty names, or
    // because the existing names are too long

    highsLogUser(log_options, HighsLogType::WARNING,
                 "There are empty or excessively-long %s names: using "
                 "constructed names with prefix %s\n",
                 name_type.c_str(), name_prefix.c_str());
    for (HighsInt ix = 0; ix < num_name; ix++)
      names[ix] = name_prefix + std::to_string(ix);
  } else {
    // Using original names, so look to see whether there are names with spaces
    names_with_spaces = namesWithSpaces(num_name, names);
  }
  // Find the final maximum name length
  max_name_length = maxNameLength(num_name, names);
  // Can't have names with spaces and more than 8 characters
  if (max_name_length > 8 && names_with_spaces) return HighsStatus::Error;
  if (construct_names) return HighsStatus::Warning;
  return HighsStatus::OK;
}

HighsBasisStatus checkedVarHighsNonbasicStatus(
    const HighsBasisStatus ideal_status, const double lower,
    const double upper) {
  HighsBasisStatus checked_status;
  if (ideal_status == HighsBasisStatus::LOWER ||
      ideal_status == HighsBasisStatus::ZERO) {
    // Looking to give status LOWER or ZERO
    if (highs_isInfinity(-lower)) {
      // Lower bound is infinite
      if (highs_isInfinity(upper)) {
        // Upper bound is infinite
        checked_status = HighsBasisStatus::ZERO;
      } else {
        // Upper bound is finite
        checked_status = HighsBasisStatus::UPPER;
      }
    } else {
      checked_status = HighsBasisStatus::LOWER;
    }
  } else {
    // Looking to give status UPPER
    if (highs_isInfinity(upper)) {
      // Upper bound is infinite
      if (highs_isInfinity(-lower)) {
        // Lower bound is infinite
        checked_status = HighsBasisStatus::ZERO;
      } else {
        // Upper bound is finite
        checked_status = HighsBasisStatus::LOWER;
      }
    } else {
      checked_status = HighsBasisStatus::UPPER;
    }
  }
  return checked_status;
}

// Return a string representation of PrimalDualStatus
std::string utilPrimalDualStatusToString(const HighsInt primal_dual_status) {
  switch (primal_dual_status) {
    case kHighsPrimalDualStatusNotset:
      return "Not set";
      break;
    case kHighsPrimalDualStatusNoSolution:
      return "No solution";
      break;
    case kHighsPrimalDualStatusUnknown:
      return "Point of unknown feasibility";
      break;
    case kHighsPrimalDualStatusInfeasiblePoint:
      return "Infeasible point";
      break;
    case kHighsPrimalDualStatusFeasiblePoint:
      return "Feasible point";
      break;
    default:
#ifdef HiGHSDEV
      printf("Primal/dual status %" HIGHSINT_FORMAT " not recognised\n",
             primal_dual_status);
#endif
      return "Unrecognised primal/dual status";
      break;
  }
  return "";
}

// Return a string representation of HighsModelStatus.
std::string utilModelStatusToString(const HighsModelStatus model_status) {
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
    case HighsModelStatus::PRESOLVE_ERROR:
      return "Presolve error";
      break;
    case HighsModelStatus::SOLVE_ERROR:
      return "Solve error";
      break;
    case HighsModelStatus::POSTSOLVE_ERROR:
      return "Postsolve error";
      break;
    case HighsModelStatus::MODEL_EMPTY:
      return "Model empty";
      break;
    case HighsModelStatus::PRIMAL_INFEASIBLE:
      return "Infeasible";  //"Primal infeasible";
      break;
    case HighsModelStatus::PRIMAL_UNBOUNDED:
      return "Unbounded";  //"Primal unbounded";
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
    case HighsModelStatus::PRIMAL_DUAL_INFEASIBLE:
      return "Primal and dual infeasible";
      break;
    case HighsModelStatus::PRIMAL_INFEASIBLE_OR_UNBOUNDED:
      return "Primal infeasible or unbounded";
      break;
    case HighsModelStatus::DUAL_INFEASIBLE:
      return "Dual infeasible";
      break;
    default:
#ifdef HiGHSDEV
      printf("HiGHS model status %" HIGHSINT_FORMAT " not recognised\n",
             (HighsInt)model_status);
#endif
      return "Unrecognised HiGHS model status";
      break;
  }
  return "";
}

void zeroHighsIterationCounts(HighsIterationCounts& iteration_counts) {
  iteration_counts.simplex = 0;
  iteration_counts.ipm = 0;
  iteration_counts.crossover = 0;
}

void zeroHighsIterationCounts(HighsInfo& info) {
  info.simplex_iteration_count = 0;
  info.ipm_iteration_count = 0;
  info.crossover_iteration_count = 0;
  info.mip_node_count = -1;
}

void copyHighsIterationCounts(const HighsIterationCounts& iteration_counts,
                              HighsInfo& info) {
  info.simplex_iteration_count = iteration_counts.simplex;
  info.ipm_iteration_count = iteration_counts.ipm;
  info.crossover_iteration_count = iteration_counts.crossover;
  info.mip_node_count = -1;
}

void copyHighsIterationCounts(const HighsInfo& info,
                              HighsIterationCounts& iteration_counts) {
  iteration_counts.simplex = info.simplex_iteration_count;
  iteration_counts.ipm = info.ipm_iteration_count;
  iteration_counts.crossover = info.crossover_iteration_count;
}

// Deduce the HighsStatus value corresponding to a HighsModelStatus value.
HighsStatus highsStatusFromHighsModelStatus(HighsModelStatus model_status) {
  switch (model_status) {
    case HighsModelStatus::NOTSET:
      return HighsStatus::Error;
    case HighsModelStatus::LOAD_ERROR:
      return HighsStatus::Error;
    case HighsModelStatus::MODEL_ERROR:
      return HighsStatus::Error;
    case HighsModelStatus::PRESOLVE_ERROR:
      return HighsStatus::Error;
    case HighsModelStatus::SOLVE_ERROR:
      return HighsStatus::Error;
    case HighsModelStatus::POSTSOLVE_ERROR:
      return HighsStatus::Error;
    case HighsModelStatus::MODEL_EMPTY:
      return HighsStatus::OK;
    case HighsModelStatus::OPTIMAL:
      return HighsStatus::OK;
    case HighsModelStatus::PRIMAL_INFEASIBLE:
      return HighsStatus::OK;
    case HighsModelStatus::PRIMAL_INFEASIBLE_OR_UNBOUNDED:
      return HighsStatus::OK;
    case HighsModelStatus::PRIMAL_UNBOUNDED:
      return HighsStatus::OK;
    case HighsModelStatus::PRIMAL_DUAL_INFEASIBLE:
      return HighsStatus::OK;
    case HighsModelStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND:
      return HighsStatus::OK;
    case HighsModelStatus::REACHED_TIME_LIMIT:
      return HighsStatus::Warning;
    case HighsModelStatus::REACHED_ITERATION_LIMIT:
      return HighsStatus::Warning;
    case HighsModelStatus::DUAL_INFEASIBLE:
      return HighsStatus::Warning;
    default:
      return HighsStatus::Error;
  }
}
