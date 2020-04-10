/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "io/Filereader.h"

#include <cstring>  // For strrchr

#include "io/FilereaderEms.h"
#include "io/FilereaderLp.h"
#include "io/FilereaderMps.h"
#include "io/HighsIO.h"

static const char* getFilenameExt(const char* filename) {
  const char* dot = strrchr(filename, '.');
  if (!dot || dot == filename) return "";
  return dot + 1;
}

Filereader* Filereader::getFilereader(const char* filename) {
  Filereader* reader;
  const char* extension = getFilenameExt(filename);
  if (strcmp(extension, "mps") == 0) {
    reader = new FilereaderMps();
  } else if (strcmp(extension, "lp") == 0) {
    reader = new FilereaderLp();
  } else if (strcmp(extension, "ems") == 0) {
    reader = new FilereaderEms();
  } else {
    reader = NULL;
  }
  return reader;
}

void interpretFilereaderRetcode(FILE* logfile, const std::string filename,
                                const FilereaderRetcode code) {
  switch (code) {
    case FilereaderRetcode::OK:
      break;
    case FilereaderRetcode::FILENOTFOUND:
      HighsLogMessage(logfile, HighsMessageType::ERROR, "File %s not found",
                      filename.c_str());
      break;
    case FilereaderRetcode::PARSERERROR:
      HighsLogMessage(logfile, HighsMessageType::ERROR,
                      "Parser error reading %s", filename.c_str());
      break;
    case FilereaderRetcode::NOT_IMPLEMENTED:
      HighsLogMessage(logfile, HighsMessageType::ERROR,
                      "Parser not implemented for %s", filename.c_str());
      break;
  }
}

std::string extractModelName(const std::string filename) {
  // Extract model name.
  std::string name = filename;
  std::size_t found = name.find_last_of("/\\");
  if (found < name.size()) name = name.substr(found + 1);
  found = name.find_last_of(".");
  if (found < name.size()) name.erase(found, name.size() - found);
  return name;
}
