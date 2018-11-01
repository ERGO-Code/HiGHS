#ifndef FILEREADER_H
#define FILEREADER_H

#include <stdio.h>

#include "../lp_data/HighsLp.h"

enum class FilereaderRetcode { OKAY, FILENOTFOUND, PARSERERROR };

class Filereader {
 public:
  virtual FilereaderRetcode readModelFromFile(const char* filename,
                                              HighsLp& model);
  static void readLineFromFile(FILE* file, char* buffer, int buffersize);
};

#endif