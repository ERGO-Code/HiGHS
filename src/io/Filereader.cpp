#include "Filereader.h"

void Filereader::readLineFromFile(FILE* file, char* buffer, int buffersize) {
  fgets(buffer, buffersize, file);
}