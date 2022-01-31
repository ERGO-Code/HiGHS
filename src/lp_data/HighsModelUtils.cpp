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
/**@file lp_data/HighsUtils.cpp
 * @brief Class-independent utilities for HiGHS
 */

#include "lp_data/HighsModelUtils.h"

#include <algorithm>
#include <sstream>
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
  highsLogDev(log_options, HighsLogType::kInfo,
              "Analysing %" HIGHSINT_FORMAT " %s bounds\n", numBd, message);
  if (numFr > 0)
    highsLogDev(log_options, HighsLogType::kInfo,
                "   Free:  %7" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT "%%)\n",
                numFr, (100 * numFr) / numBd);
  if (numLb > 0)
    highsLogDev(log_options, HighsLogType::kInfo,
                "   LB:    %7" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT "%%)\n",
                numLb, (100 * numLb) / numBd);
  if (numUb > 0)
    highsLogDev(log_options, HighsLogType::kInfo,
                "   UB:    %7" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT "%%)\n",
                numUb, (100 * numUb) / numBd);
  if (numBx > 0)
    highsLogDev(log_options, HighsLogType::kInfo,
                "   Boxed: %7" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT "%%)\n",
                numBx, (100 * numBx) / numBd);
  if (numFx > 0)
    highsLogDev(log_options, HighsLogType::kInfo,
                "   Fixed: %7" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT "%%)\n",
                numFx, (100 * numFx) / numBd);
  highsLogDev(log_options, HighsLogType::kInfo,
              "grep_CharMl,%s,Free,LB,UB,Boxed,Fixed\n", message);
  highsLogDev(log_options, HighsLogType::kInfo,
              "grep_CharMl,%" HIGHSINT_FORMAT ",%" HIGHSINT_FORMAT
              ",%" HIGHSINT_FORMAT ",%" HIGHSINT_FORMAT ",%" HIGHSINT_FORMAT
              ",%" HIGHSINT_FORMAT "\n",
              numBd, numFr, numLb, numUb, numBx, numFx);
}

std::string statusToString(const HighsBasisStatus status, const double lower,
                           const double upper) {
  switch (status) {
    case HighsBasisStatus::kLower:
      if (lower == upper) {
        return "FX";
      } else {
        return "LB";
      }
      break;
    case HighsBasisStatus::kBasic:
      return "BS";
      break;
    case HighsBasisStatus::kUpper:
      return "UB";
      break;
    case HighsBasisStatus::kZero:
      return "FR";
      break;
    case HighsBasisStatus::kNonbasic:
      return "NB";
      break;
  }
  return "";
}

std::string typeToString(const HighsVarType type) {
  switch (type) {
    case HighsVarType::kContinuous:
      return "Continuous";
    case HighsVarType::kInteger:
      return "Integer   ";
    case HighsVarType::kSemiContinuous:
      return "Semi-conts";
    case HighsVarType::kSemiInteger:
      return "Semi-int  ";
    case HighsVarType::kImplicitInteger:
      return "ImpliedInt";
  }
  return "";
}

void writeModelBoundSolution(
    FILE* file, const bool columns, const HighsInt dim,
    const std::vector<double>& lower, const std::vector<double>& upper,
    const std::vector<std::string>& names, const bool have_primal,
    const std::vector<double>& primal, const bool have_dual,
    const std::vector<double>& dual, const bool have_basis,
    const std::vector<HighsBasisStatus>& status,
    const HighsVarType* integrality) {
  const bool have_names = names.size() > 0;
  if (have_names) assert((int)names.size() >= dim);
  if (have_primal) assert((int)primal.size() >= dim);
  if (have_dual) assert((int)dual.size() >= dim);
  if (have_basis) assert((int)status.size() >= dim);
  const bool have_integrality = integrality != NULL;
  std::string var_status_string;
  if (columns) {
    fprintf(file, "Columns\n");
  } else {
    fprintf(file, "Rows\n");
  }
  fprintf(
      file,
      "    Index Status        Lower        Upper       Primal         Dual");
  if (have_integrality) fprintf(file, "  Type      ");
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
    if (have_integrality)
      fprintf(file, "  %s", typeToString(integrality[ix]).c_str());
    if (have_names) {
      fprintf(file, "  %-s\n", names[ix].c_str());
    } else {
      fprintf(file, "\n");
    }
  }
}

void writeModelSolution(FILE* file, const HighsLp& lp,
                        const HighsSolution& solution, const HighsInfo& info) {
  const bool have_col_names = lp.col_names_.size() > 0;
  const bool have_row_names = lp.row_names_.size() > 0;
  const bool have_primal = solution.value_valid;
  const bool have_dual = solution.dual_valid;
  std::stringstream ss;
  if (have_col_names) assert((int)lp.col_names_.size() >= lp.num_col_);
  if (have_row_names) assert((int)lp.row_names_.size() >= lp.num_row_);
  if (have_primal) {
    assert((int)solution.col_value.size() >= lp.num_col_);
    assert((int)solution.row_value.size() >= lp.num_row_);
    assert(info.primal_solution_status != kSolutionStatusNone);
  }
  if (have_dual) {
    assert((int)solution.col_dual.size() >= lp.num_col_);
    assert((int)solution.row_dual.size() >= lp.num_row_);
    assert(info.dual_solution_status != kSolutionStatusNone);
  }
  const double double_tolerance = 1e-13;
  fprintf(file, "\n# Primal solution values\n");
  if (!have_primal || info.primal_solution_status == kSolutionStatusNone) {
    fprintf(file, "None\n");
  } else {
    if (info.primal_solution_status == kSolutionStatusFeasible) {
      fprintf(file, "Feasible\n");
    } else {
      assert(info.primal_solution_status == kSolutionStatusInfeasible);
      fprintf(file, "Infeasible\n");
    }
    HighsCDouble objective_function_value = lp.offset_;
    for (HighsInt i = 0; i < lp.num_col_; ++i)
      objective_function_value += lp.col_cost_[i] * solution.col_value[i];
    std::array<char, 32> objStr =
        highsDoubleToString((double)objective_function_value, double_tolerance);
    fprintf(file, "Objective %s\n", objStr.data());
    fprintf(file, "# Columns %" HIGHSINT_FORMAT "\n", lp.num_col_);
    for (HighsInt ix = 0; ix < lp.num_col_; ix++) {
      std::array<char, 32> valStr =
          highsDoubleToString(solution.col_value[ix], double_tolerance);
      // Create a column name
      ss.str(std::string());
      ss << "C" << ix;
      const std::string name = have_col_names ? lp.col_names_[ix] : ss.str();
      fprintf(file, "%-s %s\n", name.c_str(), valStr.data());
    }
    fprintf(file, "# Rows %" HIGHSINT_FORMAT "\n", lp.num_row_);
    for (HighsInt ix = 0; ix < lp.num_row_; ix++) {
      std::array<char, 32> valStr =
          highsDoubleToString(solution.row_value[ix], double_tolerance);
      // Create a row name
      ss.str(std::string());
      ss << "R" << ix;
      const std::string name = have_row_names ? lp.row_names_[ix] : ss.str();
      fprintf(file, "%-s %s\n", name.c_str(), valStr.data());
    }
  }
  fprintf(file, "\n# Dual solution values\n");
  if (!have_dual || info.dual_solution_status == kSolutionStatusNone) {
    fprintf(file, "None\n");
  } else {
    if (info.dual_solution_status == kSolutionStatusFeasible) {
      fprintf(file, "Feasible\n");
    } else {
      assert(info.dual_solution_status == kSolutionStatusInfeasible);
      fprintf(file, "Infeasible\n");
    }
    fprintf(file, "# Columns %" HIGHSINT_FORMAT "\n", lp.num_col_);
    for (HighsInt ix = 0; ix < lp.num_col_; ix++) {
      std::array<char, 32> valStr =
          highsDoubleToString(solution.col_dual[ix], double_tolerance);
      ss.str(std::string());
      ss << "C" << ix;
      const std::string name = have_col_names ? lp.col_names_[ix] : ss.str();
      fprintf(file, "%-s %s\n", name.c_str(), valStr.data());
    }
    fprintf(file, "# Rows %" HIGHSINT_FORMAT "\n", lp.num_row_);
    for (HighsInt ix = 0; ix < lp.num_row_; ix++) {
      std::array<char, 32> valStr =
          highsDoubleToString(solution.row_dual[ix], double_tolerance);
      ss.str(std::string());
      ss << "R" << ix;
      const std::string name = have_row_names ? lp.row_names_[ix] : ss.str();
      fprintf(file, "%-s %s\n", name.c_str(), valStr.data());
    }
  }
}

bool hasNamesWithSpaces(const HighsLogOptions& log_options,
                        const HighsInt num_name,
                        const std::vector<std::string>& names) {
  HighsInt num_names_with_spaces = 0;
  for (HighsInt ix = 0; ix < num_name; ix++) {
    HighsInt space_pos = names[ix].find(" ");
    if (space_pos >= 0) {
      if (num_names_with_spaces == 0) {
        highsLogDev(
            log_options, HighsLogType::kInfo,
            "Name |%s| contains a space character in position %" HIGHSINT_FORMAT
            "\n",
            names[ix].c_str(), space_pos);
        num_names_with_spaces++;
      }
    }
  }
  if (num_names_with_spaces)
    highsLogDev(log_options, HighsLogType::kInfo,
                "There are %" HIGHSINT_FORMAT " names with spaces\n",
                num_names_with_spaces);
  return num_names_with_spaces > 0;
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

    highsLogUser(log_options, HighsLogType::kWarning,
                 "There are empty or excessively-long %s names: using "
                 "constructed names with prefix %s\n",
                 name_type.c_str(), name_prefix.c_str());
    for (HighsInt ix = 0; ix < num_name; ix++)
      names[ix] = name_prefix + std::to_string(ix);
  } else {
    // Using original names, so look to see whether there are names with spaces
    names_with_spaces = hasNamesWithSpaces(log_options, num_name, names);
  }
  // Find the final maximum name length
  max_name_length = maxNameLength(num_name, names);
  // Can't have names with spaces and more than 8 characters
  if (max_name_length > 8 && names_with_spaces) return HighsStatus::kError;
  if (construct_names) return HighsStatus::kWarning;
  return HighsStatus::kOk;
}

HighsBasisStatus checkedVarHighsNonbasicStatus(
    const HighsBasisStatus ideal_status, const double lower,
    const double upper) {
  HighsBasisStatus checked_status;
  if (ideal_status == HighsBasisStatus::kLower ||
      ideal_status == HighsBasisStatus::kZero) {
    // Looking to give status LOWER or ZERO
    if (highs_isInfinity(-lower)) {
      // Lower bound is infinite
      if (highs_isInfinity(upper)) {
        // Upper bound is infinite
        checked_status = HighsBasisStatus::kZero;
      } else {
        // Upper bound is finite
        checked_status = HighsBasisStatus::kUpper;
      }
    } else {
      checked_status = HighsBasisStatus::kLower;
    }
  } else {
    // Looking to give status UPPER
    if (highs_isInfinity(upper)) {
      // Upper bound is infinite
      if (highs_isInfinity(-lower)) {
        // Lower bound is infinite
        checked_status = HighsBasisStatus::kZero;
      } else {
        // Upper bound is finite
        checked_status = HighsBasisStatus::kLower;
      }
    } else {
      checked_status = HighsBasisStatus::kUpper;
    }
  }
  return checked_status;
}

// Return a string representation of SolutionStatus
std::string utilSolutionStatusToString(const HighsInt solution_status) {
  switch (solution_status) {
    case kSolutionStatusNone:
      return "None";
      break;
    case kSolutionStatusInfeasible:
      return "Infeasible";
      break;
    case kSolutionStatusFeasible:
      return "Feasible";
      break;
    default:
      assert(1 == 0);
      return "Unrecognised solution status";
  }
}

// Return a string representation of HighsBasisStatus
std::string utilBasisStatusToString(const HighsBasisStatus basis_status) {
  switch (basis_status) {
    case HighsBasisStatus::kLower:
      return "At lower/fixed bound";
      break;
    case HighsBasisStatus::kBasic:
      return "Basic";
      break;
    case HighsBasisStatus::kUpper:
      return "At upper bound";
      break;
    case HighsBasisStatus::kZero:
      return "Free at zero";
      break;
    case HighsBasisStatus::kNonbasic:
      return "Nonbasic";
      break;
    default:
      assert(1 == 0);
      return "Unrecognised solution status";
  }
}

// Return a string representation of basis validity
std::string utilBasisValidityToString(const HighsInt basis_validity) {
  if (basis_validity) {
    return "Valid";
  } else {
    return "Not valid";
  }
}

// Return a string representation of HighsModelStatus.
std::string utilModelStatusToString(const HighsModelStatus model_status) {
  switch (model_status) {
    case HighsModelStatus::kNotset:
      return "Not Set";
      break;
    case HighsModelStatus::kLoadError:
      return "Load error";
      break;
    case HighsModelStatus::kModelError:
      return "Model error";
      break;
    case HighsModelStatus::kPresolveError:
      return "Presolve error";
      break;
    case HighsModelStatus::kSolveError:
      return "Solve error";
      break;
    case HighsModelStatus::kPostsolveError:
      return "Postsolve error";
      break;
    case HighsModelStatus::kModelEmpty:
      return "Empty";
      break;
    case HighsModelStatus::kOptimal:
      return "Optimal";
      break;
    case HighsModelStatus::kInfeasible:
      return "Infeasible";
      break;
    case HighsModelStatus::kUnboundedOrInfeasible:
      return "Primal infeasible or unbounded";
      break;
    case HighsModelStatus::kUnbounded:
      return "Unbounded";
      break;
    case HighsModelStatus::kObjectiveBound:
      return "Bound on objective reached";
      break;
    case HighsModelStatus::kObjectiveTarget:
      return "Target for objective reached";
      break;
    case HighsModelStatus::kTimeLimit:
      return "Time limit reached";
      break;
    case HighsModelStatus::kIterationLimit:
      return "Iteration limit reached";
      break;
    case HighsModelStatus::kUnknown:
      return "Unknown";
      break;
    default:
      assert(1 == 0);
      return "Unrecognised HiGHS model status";
  }
}

// Deduce the HighsStatus value corresponding to a HighsModelStatus value.
HighsStatus highsStatusFromHighsModelStatus(HighsModelStatus model_status) {
  switch (model_status) {
    case HighsModelStatus::kNotset:
      return HighsStatus::kError;
    case HighsModelStatus::kLoadError:
      return HighsStatus::kError;
    case HighsModelStatus::kModelError:
      return HighsStatus::kError;
    case HighsModelStatus::kPresolveError:
      return HighsStatus::kError;
    case HighsModelStatus::kSolveError:
      return HighsStatus::kError;
    case HighsModelStatus::kPostsolveError:
      return HighsStatus::kError;
    case HighsModelStatus::kModelEmpty:
      return HighsStatus::kOk;
    case HighsModelStatus::kOptimal:
      return HighsStatus::kOk;
    case HighsModelStatus::kInfeasible:
      return HighsStatus::kOk;
    case HighsModelStatus::kUnboundedOrInfeasible:
      return HighsStatus::kOk;
    case HighsModelStatus::kUnbounded:
      return HighsStatus::kOk;
    case HighsModelStatus::kObjectiveBound:
      return HighsStatus::kOk;
    case HighsModelStatus::kObjectiveTarget:
      return HighsStatus::kOk;
    case HighsModelStatus::kTimeLimit:
      return HighsStatus::kWarning;
    case HighsModelStatus::kIterationLimit:
      return HighsStatus::kWarning;
    case HighsModelStatus::kUnknown:
      return HighsStatus::kWarning;
    default:
      return HighsStatus::kError;
  }
}
