#include <stringutil.h>

void strRemoveWhitespace(char *str) {
  char *dest = str;
  do
    while (isspace(*str))
      str++;
  while (*dest++ = *str++);
}

char *strClone(const char *str) {
  size_t n;
  n = strlen(str);

  char *cpy = new char[n + 1];
  strcpy(cpy, str);
  return cpy;
}

int strIsWhitespace(const char *str) {
  while (*str != '\0') {
    if (!isspace((unsigned char)*str)) {
      return 0;
    }
    str++;
  }
  return 1;
}

void strToLower(char *str) {
  int i;
  for (i = 0; str[i] != '\0'; i++) {
    str[i] = (char)tolower(str[i]);
  }
}

void strTrim(char *str) {
  int i;
  int begin = 0;
  int end = strlen(str) - 1;

  while (isspace((unsigned char)str[begin]))
    begin++;

  while ((end >= begin) && isspace((unsigned char)str[end]))
    end--;

  // Shift all characters back to the start of the string array.
  for (i = begin; i <= end; i++)
    str[i - begin] = str[i];

  str[i - begin] = '\0'; // Null terminate string.
}

std::string &ltrim(std::string &str, const std::string &chars) {
  str.erase(0, str.find_first_not_of(chars));
  return str;
}

std::string &rtrim(std::string &str, const std::string &chars) {
  str.erase(str.find_last_not_of(chars) + 1);
  return str;
}

std::string &trim(std::string &str, const std::string &chars) {
  return ltrim(rtrim(str, chars), chars);
}
