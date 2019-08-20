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

#include <sys/stat.h>
#include <sys/types.h>
#include <cstdio>
#include <fstream>

#include "io/Filereader.h"
#include "io/HighsIO.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsOptions.h"

// Parses the file in options.filename using the parser specified in
// options.parser
HighsStatus loadLpFromFile(const HighsOptions& options, HighsLp& lp) {
  if (options.filename.size() == 0) return HighsStatus::Error;

  // Make sure it is not a folder.

  struct stat info;
  const char* pathname = options.filename.c_str();
  if (stat(pathname, &info) != 0) {
    HighsPrintMessage(ML_ALWAYS, "Cannot access %s\n", pathname);
    return HighsStatus::Error;
  } else if (info.st_mode & S_IFDIR) {
    HighsPrintMessage(ML_ALWAYS, "%s is a directory. Please specify a file.\n",
                      pathname);
    return HighsStatus::Error;
  }

  Filereader* reader = Filereader::getFilereader(options.filename.c_str());
  FilereaderRetcode success = reader->readModelFromFile(options, lp);
  delete reader;

  switch (success) {
    case FilereaderRetcode::FILENOTFOUND:
      HighsPrintMessage(ML_ALWAYS, "File not found.\n");
      return HighsStatus::Error;
    case FilereaderRetcode::PARSERERROR:
      HighsPrintMessage(ML_ALWAYS, "Error when parsing file.\n");
      return HighsStatus::Error;
    default:
      break;
  }

  lp.nnz_ = lp.Avalue_.size();

  // Extract model name.
  std::string name = options.filename;
  std::size_t found = name.find_last_of("/\\");
  if (found < name.size()) name = name.substr(found + 1);
  found = name.find_last_of(".");
  if (found < name.size()) name.erase(found, name.size() - found);
  lp.model_name_ = name;

  //  HighsSetMessagelevel(HighsPrintMessageLevel::ML_ALWAYS); reportLp(lp, 1);
  //  return checkLp(lp);
  bool normalise = true;
  return assessLp(lp, options, normalise);
}

// For extended options to be parsed from a file. Assuming options file is
// specified.
bool loadOptionsFromFile(HighsOptions& options) {
  if (options.options_file.size() == 0) return false;

  string line, option, value;
  int line_count = 0;
  std::ifstream file(options.options_file);
  if (file.is_open()) {
    while (file.good()) {
      getline(file, line);
      line_count++;
      if (line.size() == 0 || line[0] == '#') continue;

      int equals = line.find_first_of("=");
      if (equals < 0 || equals >= (int)line.size() - 1) {
        HighsLogMessage(HighsMessageType::ERROR,
                        "Error on line %d of options file.", line_count);
        return false;
      }
      option = line.substr(0, equals);
      value = line.substr(equals + 1, line.size() - equals);
      OptionStatus status;
      HighsOptionType type;
      status = getOptionType(option, type);
      switch( type )
	{
	case HighsOptionType::BOOL:
	  printf("ERROR: No method to set options of type bool\n");
	  return false;
	case HighsOptionType::INT:
	  status = setOptionValue(options, option, atoi(value.c_str()));
	  break;
	case HighsOptionType::DOUBLE:
	  status = setOptionValue(options, option, atof(value.c_str()));
	  break;
	case HighsOptionType::STRING:
	  status = setOptionValue(options, option, value);
	  break;
	default:
	  printf("ERROR: No method to set option %s of unknown type %d\n", option.c_str(), (int)type);
	  return false;
	}
      if (status != OptionStatus::OK) return false;
    }
  } else {
    HighsLogMessage(HighsMessageType::ERROR, "Options file not found.");
    return false;
  }

  return true;
}
#endif
