/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/LoadProblem.h
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef IO_LOAD_PROBLEM_H_
#define IO_LOAD_PROBLEM_H_

#include "FilereaderLp.h"
#include "FilereaderMps.h"
#include "HModel.h"

// Parses the file in options.filename using the parser specified in
// options.parser
HighsInputStatus loadLpFromFile(const HighsOptions& options, HighsLp& lp) {
  // Check if file exists
  if (options.fileName && access(options.fileName, F_OK) == -1) {
    return HighsInputStatus::FileNotFound;
  } else if (!options.fileName) {
    return HighsInputStatus::FileNotFound;
  }

  // todo: check file name extenion // until parsers work with HighsLp
  HModel model_in;
  int RtCd = model_in.load_fromMPS(options.fileName);

  lp = HModelToHighsLp(model_in);
  lp.nnz_ = lp.Avalue_.size();

  // make sure old tests pass before you start work on the
  // parsers. Then remove traces of read_fromMPS from below and replace the code
  // above with
  // if (mps) use FilereaderMps

  // else if (lp) use FilereaderLp

  return checkLp(lp); 
}

#endif
