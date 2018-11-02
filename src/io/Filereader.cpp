#include "Filereader.h"
#include "FilereaderLp.h"
#include "FilereaderMps.h"

#include <string.h>
#include <stdexcept>

void Filereader::readLineFromFile(FILE* file, char* buffer, int buffersize) {
  fgets(buffer, buffersize, file);
}

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
    throw std::invalid_argument("No suitable filereader found.");
    reader = new FilereaderLp();
  } else if (strcmp(extension, "gz") == 0) {
    // unzip the file, remove ".gz" from filename, call this function again
  } else {
    throw std::invalid_argument("No suitable filereader found.");
  }
  return reader;
}