#ifndef IO_FILEREADER_LP_H_
#define IO_FILEREADER_LP_H_

#include "Filereader.h"
#include "HighsModelBuilder.h"
#include "HighsIO.h"

#define BUFFERSIZE 561

const char* const LP_KEYWORD_MIN[] = {"minimize", "min", "minimum"};
const char* const LP_KEYWORD_MAX[] = {"maximize", "max", "maximum"};
const char* const LP_KEYWORD_ST[] = {"subject to", "such that", "st", "s.t."};
const char* const LP_KEYWORD_BOUNDS[] = {"bounds", "bound"};
const char* const LP_KEYWORD_INF[] = {"infinity", "inf"};
const char* const LP_KEYWORD_FREE[] = {"free"};
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
const int LP_KEYWORD_FREE_N = 1;
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
  FREE,
  GEN,
  BIN,
  SEMI,
  SOS,
  END
};

enum class LpComparisonIndicator { LEQ, L, EQ, G, GEQ };

enum LpTokenType {
  NONE,
  IDENTIFIER,
  KEYWORD,
  CONSTANT,
  SIGN,
  COLON,
  SENSE,
  LINEEND,
  FILEEND
};

const char * const LpTokenTypeString[] = {"NONE", "IDENTIFIER", "KEYWORD", "CONSTANT", "SIGN", "COLON", "SENSE", "LINEEND", "FILEEND"};

class LpToken {
 public:
  LpTokenType type;
  virtual void print() { HighsPrintMessage(HighsMessageType::INFO, "%s ", LpTokenTypeString[type]);}
};

class LpIdentifierToken : public LpToken {
 public:
  char identifier[BUFFERSIZE];
  void print() { }
};

class LpKeywordToken : public LpToken {
 public:
  LpSectionKeyword keyword;
};

class LpConstantToken : public LpToken {
 public:
  double constant;
};

class LpSignToken : public LpToken {
 public:
  int sign;
};

class LpComparisonToken : public LpToken {
 public:
  LpComparisonIndicator comparison;
};

const char LP_INDICATOR_COMMENT = '\\';
const char* const LP_INDICATOR_SPLIT = ":+-<>=[]*";

// reads .lp files according to https://www.ibm.com/support/knowledgecenter/SSSA5P_12.5.0/ilog.odms.cplex.help/CPLEX/FileFormats/topics/LP.html#File_formats_reference.uss_reffileformatscplex.162381__File_formats_reference.uss_reffileformatscplex.177305
// as it is not immedialy clear from the format, here are some of its limitations:
// 1) keywords ("min", "max", etc) can also be constraint/variable identifiers (judged from context)
// 2) quadratic constraints or objective funstions are not yet supported
// 3) integer, binary, sos, semi continuous variables are not yet supported
class FilereaderLp : public Filereader {
 public:
  FilereaderRetcode readModelFromFile(const char* filename, HighsLp& model);
  FilereaderLp();
  ~FilereaderLp();

 private:
  std::list<LpToken*> tokenQueue;
  FILE* file;
  char fileBuffer[BUFFERSIZE];
  char stringBuffer[BUFFERSIZE];
  char* readingPosition;

  bool isFileBufferFullyRead;

  bool tryReadNextToken();

  bool isKeyword(const char* str, const char* const* keywords,
                 const int nkeywords);
  bool tryParseSectionKeyword(const char* str, LpSectionKeyword* keyword);

};

#endif