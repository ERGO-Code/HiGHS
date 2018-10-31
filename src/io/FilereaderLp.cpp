#include "FilereaderLp.h"

FilereaderRetcode FilereaderLp::readModelFromFile(const char* filename,
                                                  HighsLp& model) {
  char buffer[BUFFERSIZE];
  FILE* file = fopen(filename, "r");
  if (file == 0) {
    return FilereaderRetcode::FILENOTFOUND;
  }

  return FilereaderRetcode::OKAY;
}

void FilereaderLp::readObjectiveSense(FILE* file, char* buffer) {
  Filereader::readLineFromFile(file, buffer, BUFFERSIZE);
}