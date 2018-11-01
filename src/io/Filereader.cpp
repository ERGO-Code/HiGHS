#include "Filereader.h"

void __inline__ Filereader::readLineFromFile(FILE* file, char* buffer, int buffersize) {
  fgets(buffer, buffersize, file);
}