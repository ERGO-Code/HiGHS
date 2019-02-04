/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
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

#include "Filereader.h"
#include "HighsIO.h"

// Parses the file in options.filename using the parser specified in
// options.parser
HighsInputStatus loadLpFromFile(const HighsOptions& options, HighsLp& lp) {
  // Check if file exists
  if (options.filename.size() == 0 || access(options.filename.c_str(), F_OK) == -1)
    return HighsInputStatus::FileNotFound;


  Filereader* reader = Filereader::getFilereader(options.filename.c_str());
  FilereaderRetcode success =  reader->readModelFromFile(options.filename.c_str(), lp);
  delete reader;
  lp.nnz_ = lp.Avalue_.size();

  // Extract model name.
  std::string name = options.filename;
  std::size_t found = name.find_last_of("/\\");
  if (found < name.size())
    name = name.substr(found + 1);
  found = name.find_last_of(".");
  if (found < name.size())
    name.erase(found, name.size() - found);
  lp.model_name_ = name;

  if(success != FilereaderRetcode::OKAY) {
    HighsLogMessage(HighsMessageType::INFO, "Error when parsing file\n");
  }

  return checkLp(lp);
}

#endif
