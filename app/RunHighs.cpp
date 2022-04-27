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
/**@file ../app/RunHighs.cpp
 * @brief HiGHS main
 */
#include "Highs.h"
//#include "io/HighsIO.h"
#include "lp_data/HighsRuntimeOptions.h"

void reportModelStatsOrError(const HighsLogOptions& log_options,
                             const HighsStatus read_status,
                             const HighsModel& model);

int main(int argc, char** argv) {
  // Create the Highs instance
  Highs highs;
  const HighsOptions& options = highs.getOptions();
  const HighsLogOptions& log_options = options.log_options;

  // Load user options
  std::string model_file;
  HighsOptions loaded_options;
  // Set "HiGHS.log" as the default log_file for the app so that
  // log_file has this value if it isn't set in the file
  loaded_options.log_file = "HiGHS.log";
  // When loading the options file, any messages are reported using
  // the default HighsLogOptions
  if (!loadOptions(log_options, argc, argv, loaded_options, model_file))
    return (int)HighsStatus::kError;
  // Open the app log file - unless output_flag is false, to avoid
  // creating an empty file. It does nothing if its name is "".
  if (loaded_options.output_flag) highs.openLogFile(loaded_options.log_file);

  // Pass the option settings to HiGHS. Only error-checking produces
  // output, but values are checked in loadOptions, so it's safe to
  // call this first so that printHighsVersionCopyright uses reporting
  // settings defined in any options file.
  highs.passOptions(loaded_options);

  // Load the model from model_file
  HighsStatus read_status = highs.readModel(model_file);
  reportModelStatsOrError(log_options, read_status, highs.getModel());
  if (read_status == HighsStatus::kError) return (int)read_status;

  // Solve the model
  HighsStatus run_status = highs.run();
  if (run_status == HighsStatus::kError) return (int)run_status;

  // Possibly compute the ranging information
  if (options.ranging == kHighsOnString) highs.getRanging();

  // Possibly write the solution to a file
  if (options.write_solution_to_file)
    highs.writeSolution(options.solution_file, options.write_solution_style);

  // Possibly write the model to a file
  if (options.write_model_to_file) {
    HighsStatus write_model_status = highs.writeModel(options.write_model_file);
    if (write_model_status == HighsStatus::kError)
      return (int)write_model_status;  // todo: change to write model error
  }
  return (int)run_status;
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
