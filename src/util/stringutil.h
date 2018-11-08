#ifndef STRINGUTIL_H
#define STRINGUTIL_H

#include <ctype.h>
#include <string.h>

void strRemoveWhitespace(char* str) {
  char* dest = str;
  do
    while (isspace(*str)) str++;
  while (*dest++ = *str++);
}

char* strClone(const char* str) {
  size_t n;
  n = strlen(str);

  char* cpy = new char[n];
  strcpy(cpy, str);
}

int strIsWhitespace(const char* str) {
  while (*str != '\0') {
    if (!isspace((unsigned char)*str)) {
      return 0;
    }
    str++;
  }
  return 1;
}

void strToLower(char* str) {
  int i;
  for (i = 0; str[i] != '\0'; i++) {
    str[i] = (char)tolower(str[i]);
  }
}

void strTrim(char* str) {
  int i;
  int begin = 0;
  int end = strlen(str) - 1;

  while (isspace((unsigned char)str[begin])) begin++;

  while ((end >= begin) && isspace((unsigned char)str[end])) end--;

  // Shift all characters back to the start of the string array.
  for (i = begin; i <= end; i++) str[i - begin] = str[i];

  str[i - begin] = '\0';  // Null terminate string.
}

#endif