#ifndef __READERLP_DEF_HPP__
#define __READERLP_DEF_HPP__

#include <stdexcept>
#include <string>

void inline lpassert(const bool condition, const std::string message = "") {
   if (!condition) {
     const std::string feedback = "LP file or format error" +
       message == "" ? "" : ": for " + message;
     printf("lpassert fails%s\n", (message == "" ? "" : (" for " + message).c_str()));
      throw std::invalid_argument(feedback);
   }
}

const std::string LP_KEYWORD_INF[] = {"infinity", "inf"};
const std::string LP_KEYWORD_FREE[] = {"free"};

const unsigned int LP_KEYWORD_INF_N = 2;
const unsigned int LP_KEYWORD_FREE_N = 1;

#endif
