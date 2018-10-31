#ifndef IO_FILEREADER_H_
#define IO_FILEREADER_H_

#include "HighsLp.h"

enum class FilereaderRetcode { OKAY = 0, FILENOTFOUND = 1, PARSERERROR = 2 };

class Filereader {
 public:
  virtual FilereaderRetcode readModelFromFile(const char filename,
                                              HighsLp& model) = 0;
};

#endif