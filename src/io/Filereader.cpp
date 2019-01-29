/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "Filereader.h"
#include "FilereaderLp.h"
#include "FilereaderMps.h"
#include "FilereaderEms.h"

#include <string.h>
#include <stdexcept>

static __inline__ const char* getFilenameExt(const char* filename) {
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
    // use .mps filereader by default
    reader = new FilereaderMps();
  }
  return reader;
}