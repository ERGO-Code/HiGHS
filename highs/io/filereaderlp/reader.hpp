#ifndef __READERLP_READER_HPP__
#define __READERLP_READER_HPP__

#include <string>

#include "io/filereaderlp/model.hpp"

const std::string kLpKeywordInf[] = {"infinity", "inf"};
const std::string kLpKeywordFree[] = {"free"};

const unsigned int kLpKeywordInfN = 2;
const unsigned int kLpKeywordFreeN = 1;

Model readinstance(std::string filename);

#endif
