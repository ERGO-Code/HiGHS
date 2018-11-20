#ifndef IO_FILEREADER_MPS_H_
#define IO_FILEREADER_MPS_H_

#include "Filereader.h"

#include "HMPSIO.h"
#ifdef Boost_FOUND
  #include "HMpsFF.h"
#endif

class FilereaderMps : public Filereader {
 public:
  FilereaderRetcode readModelFromFile(const char* filename, HighsLp& model);
  FilereaderRetcode writeModelToFile(const char* filename, HighsLp model);
};

#endif