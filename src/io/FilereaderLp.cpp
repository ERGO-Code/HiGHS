#include "FilereaderLp.h"

#include "../util/stringutil.h"

#include <assert.h>
#include "HighsIO.h"

FilereaderLp::FilereaderLp() {
  HighsPrintMessage(HighsMessageType::INFO, "Instantiated .lp file reader\n");
  this->isFileBufferFullyRead = true;  // nothing to read
}

FilereaderLp::~FilereaderLp() {
  HighsPrintMessage(HighsMessageType::INFO, "Destroyed .lp file reader\n");
}

FilereaderRetcode FilereaderLp::readModelFromFile(const char* filename,
                                                  HighsLp& model) {
  this->file = fopen(filename, "r");
  if (file == NULL) {
    return FilereaderRetcode::FILENOTFOUND;
  }
  HighsModelBuilder* modelBuilder = new HighsModelBuilder();

  // add extra new line
  LpToken* newToken = new LpToken();
  newToken->type = LpTokenType::LINEEND;
  this->tokenQueue.push_back(*newToken);

  // this->readFile(file, *modelBuilder);
  while (this->tryReadNextToken())
    ;
  
  while(!this->tokenQueue.empty()) {
    LpToken token = this->tokenQueue.front();
    this->tokenQueue.pop_front();
    token.print();
  }

  return FilereaderRetcode::OKAY;
}

bool FilereaderLp::tryReadNextToken() {
  LpToken previousToken = this->tokenQueue.back();
  if (this->isFileBufferFullyRead) {
    char* eof = fgets(this->fileBuffer, BUFFERSIZE, this->file);
    if (eof == NULL) {
      LpToken* newToken = new LpToken();
      newToken->type = LpTokenType::FILEEND;
      this->tokenQueue.push_back(*newToken);
      return false;
    }
    this->isFileBufferFullyRead = false;
    this->readingPosition = this->fileBuffer;
    HighsPrintMessage(HighsMessageType::INFO, "Now reading '%s'\n", this->fileBuffer);
  }
  bool isKeyword; LpSectionKeyword keyword; int charactersConsumedKeyword;
  bool isComment; 
  bool isIdentifier; char identifierBuffer[BUFFERSIZE]; int charactersConsumedIdentifier;
  bool isConstant; double constantBuffer; int charactersConsumedConstant;
  bool isSign;  int sign; int charactersConsumedSign;
  bool isComparison; LpComparisonIndicator comparision; int charactersConsumedComparision;
  bool isWhitespace; int charactersConsumedWhitespace;
  bool isLineEnd;

  int charactersConsumed;
  int nread;

  //zeroth, check if a comment starts at the current position
  if(*this->readingPosition == LP_INDICATOR_COMMENT) {
    HighsPrintMessage(HighsMessageType::INFO, "Ignored comment '%s'\n", this->readingPosition);
    *this->readingPosition = '\n';
    *(this->readingPosition+1) = '\0';
  }

  // first, check if the current position matches a keyword containing one space
  nread = sscanf(this->readingPosition, "%[^\t\n:+<>=-]%n", this->stringBuffer, &charactersConsumed); // probably won't work if the line continues after the keyword
  isKeyword = tryParseSectionKeyword(this->stringBuffer, &keyword);
  if(nread == 1 && isKeyword && previousToken.type == LpTokenType::LINEEND) {
    //check if previous token was a newline
    HighsPrintMessage(HighsMessageType::INFO, "Token was keyword: '%s'\n", this->stringBuffer);
    this->readingPosition += charactersConsumed;

    LpKeywordToken* newToken = new LpKeywordToken();
    newToken->keyword = keyword;
    newToken->type = LpTokenType::KEYWORD;
    this->tokenQueue.push_back(*newToken);
    return true;
  } 

  // check all other cases (which don't allow for tokens to contain a whitespace)
  nread = sscanf(this->readingPosition, "%lf%n", &constantBuffer, &charactersConsumedConstant);
  if(nread == 1) {
    HighsPrintMessage(HighsMessageType::INFO, "\tToken was constant '%lf'\n", constantBuffer);
      this->readingPosition += charactersConsumedConstant;

      LpConstantToken* newToken = new LpConstantToken();
      newToken->constant = constantBuffer;
      newToken->type = LpTokenType::CONSTANT;
      this->tokenQueue.push_back(*newToken);

      return true;

  }
  nread = sscanf(this->readingPosition, "%[^\t\n:+<>=\ -]%n", this->stringBuffer, &charactersConsumed);
  if (nread == 1) {
    HighsPrintMessage(HighsMessageType::INFO, "Next string to analyze for token: '%s'\n", this->stringBuffer);
    // can be a keyword, a constant, a name, a constant and a name
    isKeyword = tryParseSectionKeyword(this->stringBuffer, &keyword);
    if(isKeyword && previousToken.type == LpTokenType::LINEEND) {
      HighsPrintMessage(HighsMessageType::INFO, "\tToken was keyword '%s'\n", this->stringBuffer);
      this->readingPosition += charactersConsumed;

      LpKeywordToken* newToken = new LpKeywordToken();
      newToken->keyword = keyword;
      newToken->type = LpTokenType::KEYWORD;
      this->tokenQueue.push_back(*newToken);

      return true;
    } else {
      // check if identifier matches naming conventions
      HighsPrintMessage(HighsMessageType::INFO, "\tToken was identifier '%s'\n", this->stringBuffer);
      this->readingPosition += charactersConsumed;

      LpIdentifierToken* newToken = new LpIdentifierToken();
      strcpy(newToken->identifier, this->stringBuffer);
      newToken->type = LpTokenType::IDENTIFIER;
      this->tokenQueue.push_back(*newToken);

      return true;
    }
  } else {
    //no match, so not a more-than-one-character token

    //end of current buffer (i.e. end of line)
    if(*this->readingPosition == '\0') {
      HighsPrintMessage(HighsMessageType::INFO, "Symbol was string end\n");
      this->isFileBufferFullyRead = true;
      return true;
    }

    char symbol;
    nread = sscanf(this->readingPosition, "%c", &symbol);

    LpToken* newToken;
    switch(symbol){
      case '\n':
        HighsPrintMessage(HighsMessageType::INFO, "Symbol was newline\n");
        newToken = new LpToken();
        newToken->type = LpTokenType::LINEEND;
        this->tokenQueue.push_back(*newToken);
        break;
      case '\t':
      case ' ':
        HighsPrintMessage(HighsMessageType::INFO, "Symbol was whitespace\n");
        break;
      case '=':
        HighsPrintMessage(HighsMessageType::INFO, "Symbol was comparison indicator '%c'\n", symbol);
        newToken = new LpComparisonToken();
        newToken->type = LpTokenType::SENSE;
        ((LpComparisonToken*)newToken)->comparison = LpComparisonIndicator::EQ;
        this->tokenQueue.push_back(*newToken);
        break;
      case '>':
        HighsPrintMessage(HighsMessageType::INFO, "Symbol was comparison indicator '%c'\n", symbol);
        newToken = new LpComparisonToken();
        newToken->type = LpTokenType::SENSE;
        ((LpComparisonToken*)newToken)->comparison = LpComparisonIndicator::G;
        this->tokenQueue.push_back(*newToken);
        break;
      case '<':
        HighsPrintMessage(HighsMessageType::INFO, "Symbol was comparison indicator '%c'\n", symbol);
        newToken = new LpComparisonToken();
        newToken->type = LpTokenType::SENSE;
        ((LpComparisonToken*)newToken)->comparison = LpComparisonIndicator::L;
        this->tokenQueue.push_back(*newToken);
        break;
      case ':':
        newToken = new LpToken();
        newToken->type = LpTokenType::COLON;
        this->tokenQueue.push_back(*newToken);
        break;
      case '+':
        HighsPrintMessage(HighsMessageType::INFO, "Symbol was sign '%c'\n", symbol);
        newToken = new LpSignToken();
        newToken->type = LpTokenType::SIGN;
        ((LpSignToken*)newToken)->sign = +1;
        this->tokenQueue.push_back(*newToken);
        break;
      case '-':
        HighsPrintMessage(HighsMessageType::INFO, "Symbol was sign '%c'\n", symbol);
        newToken = new LpSignToken();
        newToken->type = LpTokenType::SIGN;
        ((LpSignToken*)newToken)->sign = -1;
        this->tokenQueue.push_back(*newToken);
        break;
      default:
      HighsPrintMessage(HighsMessageType::ERROR, "Symbol was unknown sign '%c'\n", symbol);
    }
    this->readingPosition++;
    return true;
  }
  
  //should never happen
  return false;
}



bool FilereaderLp::isKeyword(const char* str, const char* const* keywords,
                             const int nkeywords) {
  int i;
  for (i = 0; i < nkeywords; i++) {
    if (strcmp(str, keywords[i]) == 0) {
      return true;
    }
  }
  return false;
}

bool FilereaderLp::tryParseSectionKeyword(const char* str,
                                          LpSectionKeyword* keyword) {
  strTrim((char*)str);
  strToLower((char*)str);

  *keyword = LpSectionKeyword::NONE;

  // min?
  if (isKeyword(str, LP_KEYWORD_MIN, LP_KEYWORD_MIN_N)) {
    *keyword = LpSectionKeyword::MIN;
  }

  // max?
  if (isKeyword(str, LP_KEYWORD_MAX, LP_KEYWORD_MAX_N)) {
    *keyword = LpSectionKeyword::MAX;
  }

  // st?
  if (isKeyword(str, LP_KEYWORD_ST, LP_KEYWORD_ST_N)) {
    *keyword = LpSectionKeyword::ST;
  }

  // bounds?
  if (isKeyword(str, LP_KEYWORD_BOUNDS, LP_KEYWORD_BOUNDS_N)) {
    *keyword = LpSectionKeyword::BOUNDS;
  }

  // inf?
  if (isKeyword(str, LP_KEYWORD_INF, LP_KEYWORD_INF_N)) {
    *keyword = LpSectionKeyword::INF;
  }

  // free?
  if (isKeyword(str, LP_KEYWORD_FREE, LP_KEYWORD_FREE_N)) {
    *keyword = LpSectionKeyword::FREE;
  }

  // gen?
  if (isKeyword(str, LP_KEYWORD_GEN, LP_KEYWORD_GEN_N)) {
    *keyword = LpSectionKeyword::GEN;
  }

  // bin?
  if (isKeyword(str, LP_KEYWORD_BIN, LP_KEYWORD_BIN_N)) {
    *keyword = LpSectionKeyword::BIN;
  }

  // semi?
  if (isKeyword(str, LP_KEYWORD_SEMI, LP_KEYWORD_SEMI_N)) {
    *keyword = LpSectionKeyword::SEMI;
  }

  // sos?
  if (isKeyword(str, LP_KEYWORD_SOS, LP_KEYWORD_SOS_N)) {
    *keyword = LpSectionKeyword::SOS;
  }

  // end?
  if (isKeyword(str, LP_KEYWORD_END, LP_KEYWORD_END_N)) {
    *keyword = LpSectionKeyword::END;
  }

  return *keyword != LpSectionKeyword::NONE;
}
