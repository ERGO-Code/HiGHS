#ifndef FILEREADER_LP
#define FILEREADER_LP

#include "Filereader.h"

#define BUFFERSIZE 256;

class FilereaderLp : public Filereader {
 public:
  FilereaderRetcode readModelFromFile(const char* filename, HighsLp& model);

 private:
  void readObjectiveSense(FILE* file);
};

#endif