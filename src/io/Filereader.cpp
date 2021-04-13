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

#include "io/Filereader.h"

#include "io/FilereaderEms.h"
#include "io/FilereaderLp.h"
#include "io/FilereaderMps.h"
#include "io/HighsIO.h"

static const std::string getFilenameExt(const std::string filename) {
  // Extract file name extension
  std::string name = filename;
  std::size_t found = name.find_last_of(".");
  if (found < name.size()) {
    name = name.substr(found + 1);
  } else {
    name = "";
  }
  return name;
}

Filereader* Filereader::getFilereader(const std::string filename) {
  Filereader* reader;
  const std::string extension = getFilenameExt(filename);
  if (extension.compare("mps") == 0) {
    reader = new FilereaderMps();
  } else if (extension.compare("lp") == 0) {
    reader = new FilereaderLp();
  } else if (extension.compare("ems") == 0) {
    reader = new FilereaderEms();
  } else {
    reader = NULL;
  }
  return reader;
}

void interpretFilereaderRetcode(const HighsLogOptions& log_options,
                                const std::string filename,
                                const FilereaderRetcode code) {
  switch (code) {
    case FilereaderRetcode::kOk:
      break;
    case FilereaderRetcode::kFileNotFound:
      highsLogUser(log_options, HighsLogType::kError, "File %s not found\n",
                   filename.c_str());
      break;
    case FilereaderRetcode::kParserError:
      highsLogUser(log_options, HighsLogType::kError,
                   "Parser error reading %s\n", filename.c_str());
      break;
    case FilereaderRetcode::kNotImplemented:
      highsLogUser(log_options, HighsLogType::kError,
                   "Parser not implemented for %s", filename.c_str());
      break;
    case FilereaderRetcode::kTimeout:
      highsLogUser(log_options, HighsLogType::kError,
                   "Parser reached timeout\n", filename.c_str());
      break;
  }
}

std::string extractModelName(const std::string filename) {
  // Extract model name
  std::string name = filename;
  std::size_t found = name.find_last_of("/\\");
  if (found < name.size()) name = name.substr(found + 1);
  found = name.find_last_of(".");
  if (found < name.size()) name.erase(found, name.size() - found);
  return name;
}
