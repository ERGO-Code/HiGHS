/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
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

FilereaderRetcode FilereaderMps::readModelFromFile(const HighsOptions& options,
                                                   HighsLp& model) {
  int status = 1;
  const char* filename = options.filename.c_str();

  // if free format parser
  // Parse file and return status.
  if (options.mps_parser_type == HighsMpsParserType::free) {
    HMpsFF parser{};
    FreeFormatParserReturnCode result = parser.loadProblem(filename, model);
    switch (result) {
      case FreeFormatParserReturnCode::SUCCESS:
        return FilereaderRetcode::OKAY;
      case FreeFormatParserReturnCode::PARSERERROR:
        return FilereaderRetcode::PARSERERROR;
      case FreeFormatParserReturnCode::FILENOTFOUND:
        return FilereaderRetcode::FILENOTFOUND;
      case FreeFormatParserReturnCode::FIXED_FORMAT:
        HighsPrintMessage(ML_DETAILED | ML_VERBOSE, "%s %s\n",
                          "Whitespaces encountered in row / col name.",
                          "Switching to fixed format parser.");
        break;
    }
  }

  // else use fixed format parser
  status = readMPS(
      filename, -1, -1, model.numRow_, model.numCol_, model.numInt_,
      model.sense_, model.offset_, model.Astart_, model.Aindex_, model.Avalue_,
      model.colCost_, model.colLower_, model.colUpper_, model.rowLower_,
      model.rowUpper_, model.integrality_, model.col_names_, model.row_names_);

  if (status) return FilereaderRetcode::PARSERERROR;
  return FilereaderRetcode::OKAY;
}

FilereaderRetcode FilereaderMps::writeModelToFile(const char* filename,
                                                  HighsLp& model) {
  int objsense = 1;
  double objoffset = 0;
  writeMPS(filename, model.numRow_, model.numCol_, model.numInt_, objsense,
           objoffset, model.Astart_, model.Aindex_, model.Avalue_,
           model.colCost_, model.colLower_, model.colUpper_, model.rowLower_,
           model.rowUpper_, model.integrality_, model.col_names_,
           model.row_names_);
  return FilereaderRetcode::OKAY;
}

FilereaderRetcode FilereaderMps::readModelFromFile(const char* filename,
                                                   HighsModelBuilder& model) {
  if (filename) {
  }  // surpress warning.
  if (model.getNumberOfVariables() > 0) {
  }  // surpress warning.
  return FilereaderRetcode::PARSERERROR;
}
