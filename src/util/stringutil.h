#ifndef STRINGUTIL_H
#define STRINGUTIL_H

#include <ctype.h>
#include <string.h>
#include <string>

void strRemoveWhitespace(char *str);
char *strClone(const char *str);
int strIsWhitespace(const char *str);
void strToLower(char *str);
void strTrim(char *str);

std::string &ltrim(std::string &str, const std::string &chars = "\t\n\v\f\r ");
std::string &rtrim(std::string &str, const std::string &chars = "\t\n\v\f\r ");
std::string &trim(std::string &str, const std::string &chars = "\t\n\v\f\r ");

bool is_empty(std::string &str, const std::string &chars = "\t\n\v\f\r ");

std::string first_word(std::string &str, int start);
int first_word_end(std::string &str, int start);

#endif