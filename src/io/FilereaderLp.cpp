#include "FilereaderLp.h"

#include "../util/stringutil.h"

#include "HighsIO.h"

FilereaderRetcode FilereaderLp::readModelFromFile(const char* filename,
                                                  HighsLp& model) {
  char buffer[BUFFERSIZE];
  FILE* file = fopen(filename, "r");
  if (file == 0) {
    return FilereaderRetcode::FILENOTFOUND;
  }

  HighsPrintMessage(HighsMessageType::DEBUG, "Reading objective sense.. \n");
  readObjectiveSense(file, buffer);
  return FilereaderRetcode::OKAY;
}

void FilereaderLp::readNextInputLine(FILE* file, char* buffer) {
  do {
    this->readLineFromFile(file, buffer, BUFFERSIZE);
  } while (1 == 1);
}

static int isKeyword(const char* str, const char* const* keywords,
                     const int nkeywords) {
  int i;
  for (i = 0; i < nkeywords; i++) {
    if (strcmp(str, keywords[i]) == 0) {
      HighsPrintMessage(HighsMessageType::DEBUG, "Found keyword %s\n", str);
      return 1;
    }
  }
  return 0;
}

LpSectionKeyword FilereaderLp::parseStrKeyword(const char* str) {
  strTrim((char*)str);
  strToLower((char*)str);

  // min?
  if (isKeyword(str, LP_KEYWORD_MIN, LP_KEYWORD_MIN_N)) {
    return LpSectionKeyword::MIN;
  }

  // max?
  if (isKeyword(str, LP_KEYWORD_MAX, LP_KEYWORD_MAX_N)) {
    return LpSectionKeyword::MAX;
  }

  // st?
  if (isKeyword(str, LP_KEYWORD_ST, LP_KEYWORD_ST_N)) {
    return LpSectionKeyword::ST;
  }

  return LpSectionKeyword::NONE;
}

void FilereaderLp::readObjectiveSense(FILE* file, char* buffer) {
  Filereader::readLineFromFile(file, buffer, BUFFERSIZE);
  LpSectionKeyword keyword = parseStrKeyword(buffer);
}