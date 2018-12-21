#include "Filereader.h"

#include <stdexcept>
#include <string.h>

#include "FilereaderLp.h"
#include "FilereaderMps.h"

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
    reader = new FilereaderLp();
  } else {
    // use .mps filereader by default
    reader = new FilereaderMps();
  }
  return reader;
}