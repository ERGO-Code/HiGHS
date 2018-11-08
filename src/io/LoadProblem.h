#ifndef IO_LOAD_PROBLEM_H_
#define IO_LOAD_PROBLEM_H_

#include "Filereader.h"

// Parses the file in options.filename using the parser specified in
// options.parser
HighsInputStatus loadLpFromFile(const HighsOptions& options, HighsLp& lp) {
  // Check if file exists
  if (options.fileName && access(options.fileName, F_OK) == -1) {
    return HighsInputStatus::FileNotFound;
  } else if (!options.fileName) {
    return HighsInputStatus::FileNotFound;
  }

  Filereader* filereader = Filereader::getFilereader(options.fileName);
  filereader->readModelFromFile(options.fileName, lp);
  delete filereader;
  exit(0);

  return checkLp(lp);
}

#endif