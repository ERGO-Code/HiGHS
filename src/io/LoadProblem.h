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

#include "FilereaderMps.h"
#include "HModel.h"

// Parses the file in options.filename using the parser specified in
// options.parser
HighsInputStatus loadLpFromFile(const HighsOptions& options, HighsLp& lp) {
  // Check if file exists
  if (options.filenames.size() == 0 || access(options.filenames.c_str(), F_OK) == -1)
    return HighsInputStatus::FileNotFound;

  // if (mps) use FilereaderMps
  FilereaderMps reader;
  reader.readModelFromFile(options.filenames.c_str(), lp);
  lp.nnz_ = lp.Avalue_.size();

  // else if (lp) use FilereaderLp

  return checkLp(lp); 
}

#endif
