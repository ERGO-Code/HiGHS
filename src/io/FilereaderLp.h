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

enum LpSectionKeyword {
  NOKEYWORD,
  MIN,
  MAX,
  ST,
  BOUNDS,
  GEN,
  BIN,
  SEMI,
  SOS,
  END
};

enum class LpSpecialKeyword {
  NONE,
  INF,
  FREE
};

enum class LpComparisonIndicator { LEQ, L, EQ, G, GEQ };

enum LpTokenType {
  NOTOKEN,
  IDENTIFIER,
  KEYWORD,
  SPECIAL,
  CONSTANT,
  SIGN,
  COLON,
  SENSE,
  LINEEND,
  FILEEND
};

const char * const LpTokenTypeString[] = {"NONE", "IDENTIFIER", "KEYWORD", "SPECIAL", "CONSTANT", "SIGN", "COLON", "SENSE", "LINEEND", "FILEEND"};
const char* const LpSectionKeywordString[] = {"NONE", "MIN", "MAX", "ST", "BOUNDS", "GEN", "BIN", "SEMI", "SOS", "END"};

class LpToken {
 public:
  LpTokenType type;
  virtual void print() { HighsPrintMessage(HighsMessageType::INFO, "%s ", LpTokenTypeString[type]);}
};

class LpIdentifierToken : public LpToken {
 public:
  char identifier[BUFFERSIZE];
};

class LpKeywordToken : public LpToken {
 public:
  LpSectionKeyword keyword;
};

class LpSpecialKeywordToken : public LpToken {
 public:
  LpSpecialKeyword keyword;
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
// 1.5) special keywords (inf, free) can not be constraint/variable identifiers
// 2) quadratic constraints or objective funstions are not yet supported //not necessary now
// 3) integer, binary, sos, semi continuous variables are not yet supported //not necessary now
// 4) Line continueing after keyword not yet supported //TODO
class FilereaderLp : public Filereader {
 public:
  FilereaderRetcode readModelFromFile(const char* filename, HighsLp& model);
  FilereaderLp();
  ~FilereaderLp();

 private:
  // list of all tokens after initial scan of file
  std::list<LpToken*> tokenQueue;

  // tokens split according to their section
  std::list<LpToken*> objectiveSection;
  std::list<LpToken*> constraintSection;
  std::list<LpToken*> boundsSection;
  std::list<LpToken*> binSection;
  std::list<LpToken*> generalSection;
  std::list<LpToken*> sosSection;
  std::list<LpToken*> semiSection;

  FILE* file;
  char fileBuffer[BUFFERSIZE];
  char stringBuffer[BUFFERSIZE];
  char* readingPosition;

  bool isFileBufferFullyRead;

  LpSectionKeyword getSectionTokens(LpSectionKeyword expectedSection);

  bool tryReadNextToken();

  bool isKeyword(const char* str, const char* const* keywords,
                 const int nkeywords);
  bool tryParseSectionKeyword(const char* str, LpSectionKeyword* keyword);
  bool tryParseSpecialKeyword(const char* str, LpSpecialKeyword* keyword);

};

#endif