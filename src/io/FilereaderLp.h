#ifndef IO_FILEREADER_LP_H_
#define IO_FILEREADER_LP_H_

#include "Filereader.h"

#define BUFFERSIZE 256

const char* const LP_KEYWORD_MIN[] = {"minimize", "min", "minimum"};
const char* const LP_KEYWORD_MAX[] = {"maximize", "max", "maximum"};
const char* const LP_KEYWORD_ST[] = {"subject to", "such that", "st", "s.t."};
const char* const LP_KEYWORD_BOUNDS[] = {"bounds", "bound"};
const char* const LP_KEYWORD_INF[] = {"infinity", "inf"};
const char* const LP_KEYWORD_GEN[] = {"general", "generals", "gen"};
const char* const LP_KEYWORD_BIN[] = {"binary", "binaries", "bin"};
const char* const LP_KEYWORD_SEMI[] = {"semi-continuous", "semi", "semis"};
const char* const LP_KEYWORD_SOS[] = {"sos"};
const char* const LP_KEYWORD_END[] = {"end"};

const int LP_KEYWORD_MIN_N = 3;
const int LP_KEYWORD_MAX_N = 3;
const int LP_KEYWORD_ST_N = 4;
const int LP_KEYWORD_BOUNDS_N = 2;
const int LP_KEYWORD_INF_N = 2;
const int LP_KEYWORD_GEN_N = 3;
const int LP_KEYWORD_BIN_N = 3;
const int LP_KEYWORD_SEMI_N = 3;
const int LP_KEYWORD_SOS_N = 1;
const int LP_KEYWORD_END_N = 1;

enum class LpSectionKeyword {
  NONE,
  MIN,
  MAX,
  ST,
  BOUNDS,
  INF,
  GEN,
  BIN,
  SEMI,
  SOS,
  END
};

const char* const LP_INDICATOR_COMMENT = "\\";

class HighsVar {
 public:
  double objectiveCoefficient;
  char* name;
  double lowerBound;
  double upperBound;
};

class HighsCons {
 public:
 // list of coefficients (HighsVar + coef)
  char* name;
  double lowerBound;
  double upperBound;
};

class FilereaderLp : public Filereader {
 public:
  FilereaderRetcode readModelFromFile(const char* filename, HighsLp& model);

 private:
  void readNextInputLine(FILE* file, char* buffer);

  LpSectionKeyword parseStrKeyword(const char* str);
  int isStrComment(const char* str);

  void readObjectiveSense(FILE* file, char* buffer);
  void readLinearConstraint(FILE* file, char* buffer);
};

#endif