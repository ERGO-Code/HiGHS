/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/FilereaderMps.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "io/FilereaderMps.h"

#include "io/HMPSIO.h"
#include "io/HMpsFF.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelUtils.h"

using free_format_parser::HMpsFF;

FilereaderRetcode FilereaderMps::readModelFromFile(const HighsOptions& options,
                                                   HighsLp& model) {
  const std::string filename = options.model_file;

  // if free format parser
  // Parse file and return status.
  if (options.mps_parser_type_free) {
    HMpsFF parser{};
    if (options.time_limit < HIGHS_CONST_INF && options.time_limit > 0)
      parser.time_limit = options.time_limit;

    FreeFormatParserReturnCode result =
        parser.loadProblem(options.log_options, filename, model);
    switch (result) {
      case FreeFormatParserReturnCode::SUCCESS:
	setOrientation(model);
        return FilereaderRetcode::OK;
      case FreeFormatParserReturnCode::PARSERERROR:
        return FilereaderRetcode::PARSERERROR;
      case FreeFormatParserReturnCode::FILENOTFOUND:
        return FilereaderRetcode::FILENOTFOUND;
      case FreeFormatParserReturnCode::FIXED_FORMAT:
        highsLogUser(options.log_options, HighsLogType::WARNING,
                     "Free format reader has detected row/col names with "
                     "spaces: switching to fixed format parser\n");
        break;
      case FreeFormatParserReturnCode::TIMEOUT:
        highsLogUser(options.log_options, HighsLogType::WARNING,
                     "Free format reader reached time_limit while parsing "
                     "the input file\n");
        return FilereaderRetcode::TIMEOUT;
    }
  }

  // else use fixed format parser
  FilereaderRetcode return_code = readMPS(
      options.log_options, filename, -1, -1, model.numRow_, model.numCol_,
      model.sense_, model.offset_, model.Astart_, model.Aindex_, model.Avalue_,
      model.colCost_, model.colLower_, model.colUpper_, model.rowLower_,
      model.rowUpper_, model.integrality_, model.col_names_, model.row_names_,
      options.keep_n_rows);
  if (return_code == FilereaderRetcode::OK) setOrientation(model);
  if (namesWithSpaces(model.numCol_, model.col_names_)) {
    highsLogUser(options.log_options, HighsLogType::WARNING,
                 "Model has column names with spaces\n");
#ifdef HiGHSDEV
    namesWithSpaces(model.numCol_, model.col_names_, true);
#endif
  }
  if (namesWithSpaces(model.numRow_, model.row_names_)) {
    highsLogUser(options.log_options, HighsLogType::WARNING,
                 "Model has row names with spaces\n");
#ifdef HiGHSDEV
    namesWithSpaces(model.numRow_, model.row_names_, true);
#endif
  }
  return return_code;
}

HighsStatus FilereaderMps::writeModelToFile(const HighsOptions& options,
                                            const std::string filename,
                                            const HighsLp& model) {
  assert(model.orientation_ != MatrixOrientation::ROWWISE);
  return writeLpAsMPS(options, filename, model);
}
