#ifndef IO_FILEREADER_MPS_H_
#define IO_FILEREADER_MPS_H_

#include "Filereader.h"

class FilereaderMps : public Filereader {
 public:
  FilereaderRetcode readModelFromFile(const char filename, HighsLp& model);
};

#endif