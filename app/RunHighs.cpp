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
/**@file ../app/RunHighs.cpp
 * @brief HiGHS main
 */
#include "Highs.h"
#include "lp_data/HighsRuntimeOptions.h"

void printHighsVersionCopyright(const HighsLogOptions& log_options);
void reportModelStatsOrError(const HighsLogOptions& log_options,
                             const HighsStatus read_status,
                             const HighsModel& model);

int main(int argc, char** argv) {
  // Create the Highs instance
  Highs highs;
  
  // Load user options.
  HighsOptions options;
  //  highsOpenLogFile(options, "Highs.log");

  printHighsVersionCopyright(options.log_options);

  std::string model_file;
  bool options_ok = loadOptions(argc, argv, options, model_file);
  if (!options_ok) return 0;
  //
  // Pass the option seetings to HiGHS
  highs.passOptions(options);
  //
  // Load the model from model_file
  HighsStatus read_status = highs.readModel(model_file);
  reportModelStatsOrError(options.log_options, read_status, highs.getModel());
  if (read_status == HighsStatus::kError)
    return 1;  // todo: change to read error
  //
  // Solve the model
  HighsStatus run_status = highs.run();

  // Possibly compute the ranging information
  if (options.ranging == kHighsOnString) highs.getRanging();

  // Possibly write the solution to a file
  if (options.write_solution_to_file)
    highs.writeSolution(options.solution_file, options.write_solution_style);

  return (int)run_status;
}

void printHighsVersionCopyright(const HighsLogOptions& log_options) {
  highsLogUser(log_options, HighsLogType::kInfo,
               "Running HiGHS %" HIGHSINT_FORMAT ".%" HIGHSINT_FORMAT
               ".%" HIGHSINT_FORMAT " [date: %s, git hash: %s]\n",
               HIGHS_VERSION_MAJOR, HIGHS_VERSION_MINOR, HIGHS_VERSION_PATCH,
               HIGHS_COMPILATION_DATE, HIGHS_GITHASH);
  highsLogUser(log_options, HighsLogType::kInfo,
               "Copyright (c) 2021 ERGO-Code under MIT licence terms\n");
}

void reportModelStatsOrError(const HighsLogOptions& log_options,
                             const HighsStatus read_status,
                             const HighsModel& model) {
  const HighsLp& lp = model.lp_;
  const HighsHessian& hessian = model.hessian_;
  if (read_status == HighsStatus::kError) {
    highsLogUser(log_options, HighsLogType::kInfo, "Error loading file\n");
  } else {
    HighsInt num_integer = 0;
    HighsInt num_semi_continuous = 0;
    HighsInt num_semi_integer = 0;
    for (HighsUInt i = 0; i < lp.integrality_.size(); i++) {
      switch (lp.integrality_[i]) {
        case HighsVarType::kInteger:
          num_integer++;
          break;
        case HighsVarType::kSemiContinuous:
          num_semi_continuous++;
          break;
        case HighsVarType::kSemiInteger:
          num_semi_integer++;
          break;
        default:
          break;
      }
    }
    std::string problem_type;
    const bool non_continuous =
        num_integer + num_semi_continuous + num_semi_integer;
    if (hessian.dim_) {
      if (non_continuous) {
        problem_type = "MIQP";
      } else {
        problem_type = "QP  ";
      }
    } else {
      if (non_continuous) {
        problem_type = "MIP ";
      } else {
        problem_type = "LP  ";
      }
    }
    const HighsInt a_num_nz = lp.a_matrix_.numNz();
    HighsInt q_num_nz = hessian.numNz();
    if (*log_options.log_dev_level) {
      highsLogDev(log_options, HighsLogType::kInfo, "%4s      : %s\n",
                  problem_type.c_str(), lp.model_name_.c_str());
      highsLogDev(log_options, HighsLogType::kInfo,
                  "Rows      : %" HIGHSINT_FORMAT "\n", lp.num_row_);
      highsLogDev(log_options, HighsLogType::kInfo,
                  "Cols      : %" HIGHSINT_FORMAT "\n", lp.num_col_);
      if (q_num_nz) {
        highsLogDev(log_options, HighsLogType::kInfo,
                    "Matrix Nz : %" HIGHSINT_FORMAT "\n", a_num_nz);
        highsLogDev(log_options, HighsLogType::kInfo,
                    "Hessian Nz: %" HIGHSINT_FORMAT "\n", q_num_nz);
      } else {
        highsLogDev(log_options, HighsLogType::kInfo,
                    "Nonzeros  : %" HIGHSINT_FORMAT "\n", a_num_nz);
      }
      if (num_integer)
        highsLogDev(log_options, HighsLogType::kInfo,
                    "Integer  : %" HIGHSINT_FORMAT "\n", num_integer);
      if (num_semi_continuous)
        highsLogDev(log_options, HighsLogType::kInfo,
                    "SemiConts: %" HIGHSINT_FORMAT "\n", num_semi_continuous);
      if (num_semi_integer)
        highsLogDev(log_options, HighsLogType::kInfo,
                    "SemiInt  : %" HIGHSINT_FORMAT "\n", num_semi_integer);
    } else {
      highsLogUser(log_options, HighsLogType::kInfo, "%s",
                   problem_type.c_str());
      if (lp.model_name_.length())
        highsLogUser(log_options, HighsLogType::kInfo, " %s",
                     lp.model_name_.c_str());
      highsLogUser(log_options, HighsLogType::kInfo,
                   " has %" HIGHSINT_FORMAT " rows; %" HIGHSINT_FORMAT " cols",
                   lp.num_row_, lp.num_col_);
      if (q_num_nz) {
        highsLogUser(log_options, HighsLogType::kInfo,
                     "; %" HIGHSINT_FORMAT " matrix nonzeros", a_num_nz);
        highsLogUser(log_options, HighsLogType::kInfo,
                     "; %" HIGHSINT_FORMAT " Hessian nonzeros", q_num_nz);
      } else {
        highsLogUser(log_options, HighsLogType::kInfo,
                     "; %" HIGHSINT_FORMAT " nonzeros", a_num_nz);
      }
      if (num_integer)
        highsLogUser(log_options, HighsLogType::kInfo,
                     "; %" HIGHSINT_FORMAT " integer variables", num_integer);
      if (num_semi_continuous)
        highsLogUser(log_options, HighsLogType::kInfo,
                     "; %" HIGHSINT_FORMAT " semi-continuous variables",
                     num_semi_continuous);
      if (num_semi_integer)
        highsLogUser(log_options, HighsLogType::kInfo,
                     "; %" HIGHSINT_FORMAT " semi-integer variables",
                     num_semi_integer);
      highsLogUser(log_options, HighsLogType::kInfo, "\n");
    }
  }
}
