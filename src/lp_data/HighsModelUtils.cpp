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
/**@file lp_data/HighsUtils.cpp
 * @brief Class-independent utilities for HiGHS
 */

#include "lp_data/HighsModelUtils.h"

#include <algorithm>
#include <vector>

#include "HConfig.h"
#include "io/HighsIO.h"
#include "lp_data/HConst.h"
#include "util/HighsUtils.h"

HighsStatus assessMatrixDimensions(const HighsLogOptions& log_options,
                                   const std::string matrix_name,
                                   const HighsInt num_vec,
                                   const vector<HighsInt>& matrix_start,
                                   const vector<HighsInt>& matrix_index,
                                   const vector<double>& matrix_value) {
  HighsStatus return_status = HighsStatus::kOk;
  // Use error_found to track whether an error has been found in multiple tests
  bool error_found = false;
  // Assess main dimensions
  bool legal_num_vec = num_vec >= 0;
  if (!legal_num_vec) {
    highsLogUser(log_options, HighsLogType::kError,
                 "%s matrix has illegal number of vectors = %" HIGHSINT_FORMAT
                 "\n",
                 matrix_name.c_str(), num_vec);
    error_found = true;
  }
  HighsInt matrix_start_size = matrix_start.size();
  bool legal_matrix_start_size = false;
  // Don't expect the matrix_start_size to be legal if there are no vectors
  if (num_vec > 0) {
    legal_matrix_start_size = matrix_start_size >= num_vec + 1;
    if (!legal_matrix_start_size) {
      highsLogUser(log_options, HighsLogType::kError,
                   "%s matrix has illegal start vector size = %" HIGHSINT_FORMAT
                   " < %" HIGHSINT_FORMAT "\n",
                   matrix_name.c_str(), matrix_start_size, num_vec + 1);
      error_found = true;
    }
  }
  if (matrix_start_size > 0) {
    // Check whether the first start is zero
    if (matrix_start[0]) {
      highsLogUser(log_options, HighsLogType::kWarning,
                   "%s matrix start vector does not begin with 0\n",
                   matrix_name.c_str());
      error_found = true;
    }
  }
  // Possibly check the sizes of the index and value vectors. Can only
  // do this with the number of nonzeros, and this is only known if
  // the start vector has a legal size. Setting num_nz = 0 otherwise
  // means that all tests pass, as they just check that the sizes of
  // the index and value vectors are non-negative.
  HighsInt num_nz = 0;
  if (legal_matrix_start_size) num_nz = matrix_start[num_vec];
  bool legal_num_nz = num_nz >= 0;
  if (!legal_num_nz) {
    highsLogUser(log_options, HighsLogType::kError,
                 "%s matrix has illegal number of nonzeros = %" HIGHSINT_FORMAT
                 "\n",
                 matrix_name.c_str(), num_nz);
    error_found = true;
  } else {
    HighsInt matrix_index_size = matrix_index.size();
    HighsInt matrix_value_size = matrix_value.size();
    bool legal_matrix_index_size = matrix_index_size >= num_nz;
    bool legal_matrix_value_size = matrix_value_size >= num_nz;
    if (!legal_matrix_index_size) {
      highsLogUser(log_options, HighsLogType::kError,
                   "%s matrix has illegal index vector size = %" HIGHSINT_FORMAT
                   " < %" HIGHSINT_FORMAT "\n",
                   matrix_name.c_str(), matrix_index_size, num_nz);
      error_found = true;
    }
    if (!legal_matrix_value_size) {
      highsLogUser(log_options, HighsLogType::kError,
                   "%s matrix has illegal value vector size = %" HIGHSINT_FORMAT
                   " < %" HIGHSINT_FORMAT "\n",
                   matrix_name.c_str(), matrix_value_size, num_nz);
      error_found = true;
    }
  }
  assert(!error_found);
  if (error_found)
    return_status = HighsStatus::kError;
  else
    return_status = HighsStatus::kOk;
  return return_status;
}

HighsStatus assessMatrix(const HighsLogOptions& log_options,
                         const std::string matrix_name, const HighsInt vec_dim,
                         const HighsInt num_vec, vector<HighsInt>& matrix_start,
                         vector<HighsInt>& matrix_index,
                         vector<double>& matrix_value,
                         const double small_matrix_value,
                         const double large_matrix_value) {
  if (assessMatrixDimensions(log_options, matrix_name, num_vec, matrix_start,
                             matrix_index, matrix_value) == HighsStatus::kError)
    return HighsStatus::kError;
  const HighsInt num_nz = matrix_start[num_vec];
  if (num_vec <= 0) return HighsStatus::kOk;
  if (num_nz <= 0) return HighsStatus::kOk;

  HighsStatus return_status = HighsStatus::kOk;
  bool error_found = false;
  bool warning_found = false;

  // Assess the starts
  // Set up previous_start for a fictitious previous empty packed vector
  HighsInt previous_start = matrix_start[0];
  for (HighsInt ix = 0; ix < num_vec; ix++) {
    HighsInt this_start = matrix_start[ix];
    bool this_start_too_small = this_start < previous_start;
    if (this_start_too_small) {
      highsLogUser(log_options, HighsLogType::kError,
                   "%s matrix packed vector %" HIGHSINT_FORMAT
                   " has illegal start of %" HIGHSINT_FORMAT
                   " < %" HIGHSINT_FORMAT
                   " = "
                   "previous start\n",
                   matrix_name.c_str(), ix, this_start, previous_start);
      return HighsStatus::kError;
    }
    bool this_start_too_big = this_start > num_nz;
    if (this_start_too_big) {
      highsLogUser(log_options, HighsLogType::kError,
                   "%s matrix packed vector %" HIGHSINT_FORMAT
                   " has illegal start of %" HIGHSINT_FORMAT
                   " > %" HIGHSINT_FORMAT
                   " = "
                   "number of nonzeros\n",
                   matrix_name.c_str(), ix, this_start, num_nz);
      return HighsStatus::kError;
    }
  }

  // Assess the indices and values
  // Count the number of acceptable indices/values
  HighsInt num_new_nz = 0;
  HighsInt num_small_values = 0;
  double max_small_value = 0;
  double min_small_value = kHighsInf;
  // Set up a zeroed vector to detect duplicate indices
  vector<HighsInt> check_vector;
  if (vec_dim > 0) check_vector.assign(vec_dim, 0);
  for (HighsInt ix = 0; ix < num_vec; ix++) {
    HighsInt from_el = matrix_start[ix];
    HighsInt to_el = matrix_start[ix + 1];
    // Account for any index-value pairs removed so far
    matrix_start[ix] = num_new_nz;
    for (HighsInt el = from_el; el < to_el; el++) {
      // Check the index
      HighsInt component = matrix_index[el];
      // Check that the index is non-negative
      bool legal_component = component >= 0;
      if (!legal_component) {
        highsLogUser(log_options, HighsLogType::kError,
                     "%s matrix packed vector %" HIGHSINT_FORMAT
                     ", entry %" HIGHSINT_FORMAT
                     ", is illegal index %" HIGHSINT_FORMAT "\n",
                     matrix_name.c_str(), ix, el, component);
        return HighsStatus::kError;
      }
      // Check that the index does not exceed the vector dimension
      legal_component = component < vec_dim;
      if (!legal_component) {
        highsLogUser(log_options, HighsLogType::kError,
                     "%s matrix packed vector %" HIGHSINT_FORMAT
                     ", entry %" HIGHSINT_FORMAT
                     ", is illegal index "
                     "%12" HIGHSINT_FORMAT " >= %" HIGHSINT_FORMAT
                     " = vector dimension\n",
                     matrix_name.c_str(), ix, el, component, vec_dim);
        return HighsStatus::kError;
      }
      // Check that the index has not already ocurred
      legal_component = check_vector[component] == 0;
      if (!legal_component) {
        highsLogUser(log_options, HighsLogType::kError,
                     "%s matrix packed vector %" HIGHSINT_FORMAT
                     ", entry %" HIGHSINT_FORMAT
                     ", is duplicate index %" HIGHSINT_FORMAT "\n",
                     matrix_name.c_str(), ix, el, component);
        return HighsStatus::kError;
      }
      // Indicate that the index has occurred
      check_vector[component] = 1;
      // Check the value
      double abs_value = fabs(matrix_value[el]);
      /*
      // Check that the value is not zero
      bool zero_value = abs_value == 0;
      if (zero_value) {
        highsLogUser(log_options, HighsLogType::kError,
                        "%s matrix packed vector %" HIGHSINT_FORMAT ", entry %"
      HIGHSINT_FORMAT ", is zero\n", matrix_name.c_str(),  ix, el); return
      HighsStatus::kError;
      }
      */
      // Check that the value is not too large
      bool large_value = abs_value > large_matrix_value;
      if (large_value) {
        highsLogUser(
            log_options, HighsLogType::kError,
            "%s matrix packed vector %" HIGHSINT_FORMAT
            ", entry %" HIGHSINT_FORMAT ", is large value |%g| >= %g\n",
            matrix_name.c_str(), ix, el, abs_value, large_matrix_value);
        return HighsStatus::kError;
      }
      bool ok_value = abs_value > small_matrix_value;
      if (!ok_value) {
        if (max_small_value < abs_value) max_small_value = abs_value;
        if (min_small_value > abs_value) min_small_value = abs_value;
        num_small_values++;
      }
      if (ok_value) {
        // Shift the index and value of the OK entry to the new
        // position in the index and value vectors, and increment
        // the new number of nonzeros
        matrix_index[num_new_nz] = matrix_index[el];
        matrix_value[num_new_nz] = matrix_value[el];
        num_new_nz++;
      } else {
        // Zero the check_vector entry since the small value
        // _hasn't_ occurred
        check_vector[component] = 0;
      }
    }
    // Zero check_vector
    for (HighsInt el = matrix_start[ix]; el < num_new_nz; el++)
      check_vector[matrix_index[el]] = 0;
#ifdef HiGHSDEV
    // NB This is very expensive so shouldn't be true
    const bool check_check_vector = false;
    if (check_check_vector) {
      // Check zeroing of check vector
      for (HighsInt component = 0; component < vec_dim; component++) {
        if (check_vector[component]) error_found = true;
      }
      if (error_found)
        highsLogUser(log_options, HighsLogType::kError,
                     "assessMatrix: check_vector not zeroed\n");
    }
#endif
  }
  if (num_small_values) {
    highsLogUser(log_options, HighsLogType::kWarning,
                 "%s matrix packed vector contains %" HIGHSINT_FORMAT
                 " |values| in [%g, %g] "
                 "less than %g: ignored\n",
                 matrix_name.c_str(), num_small_values, min_small_value,
                 max_small_value, small_matrix_value);
    warning_found = true;
  }
  matrix_start[num_vec] = num_new_nz;
  if (error_found)
    return_status = HighsStatus::kError;
  else if (warning_found)
    return_status = HighsStatus::kWarning;
  else
    return_status = HighsStatus::kOk;

  return return_status;
}

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

void writeModelBoundSol(FILE* file, const bool columns, const HighsInt dim,
                        const std::vector<double>& lower,
                        const std::vector<double>& upper,
                        const std::vector<std::string>& names,
                        const std::vector<double>& primal,
                        const std::vector<double>& dual,
                        const std::vector<HighsBasisStatus>& status) {
  const bool have_names = names.size() > 0;
  const bool have_primal = primal.size() > 0;
  const bool have_dual = dual.size() > 0;
  const bool have_basis = status.size() > 0;
  if (have_names) assert((int)names.size() >= dim);
  if (have_primal) assert((int)primal.size() >= dim);
  if (have_dual) assert((int)dual.size() >= dim);
  if (have_basis) assert((int)status.size() >= dim);
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

    highsLogUser(log_options, HighsLogType::kWarning,
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
      return "Model empty";
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
      return "Reached objective bound";
      break;
    case HighsModelStatus::kObjectiveTarget:
      return "Reached objective target";
      break;
    case HighsModelStatus::kTimeLimit:
      return "Reached time limit";
      break;
    case HighsModelStatus::kIterationLimit:
      return "Reached iteration limit";
      break;
    case HighsModelStatus::kUnknown:
      return "Unknown";
      break;
    default:
      assert(1 == 0);
      return "Unrecognised HiGHS model status";
  }
}

void zeroHighsIterationCounts(HighsIterationCounts& iteration_counts) {
  iteration_counts.simplex = 0;
  iteration_counts.ipm = 0;
  iteration_counts.crossover = 0;
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
