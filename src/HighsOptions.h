#ifndef HIGHSOPTIONS_H
#define HIGHSOPTIONS_H

#include <stdio.h>
#include <string.h>
#include <map>

struct char_cmp {
  bool operator()(const char* a, const char* b) const {
    return strcmp(a, b) < 0;
  }
};

std::map<char*, double, char_cmp> HighsDoubleOptions;
std::map<char*, int, char_cmp> HighsIntOptions;
std::map<char*, bool, char_cmp> HighsBoolOptions;
std::map<char*, FILE*, char_cmp> HighsFileOptions;

class HighsStringOptions {
 public:
  template <class T>
  static void getValue(char* key, T* ret) {
    if (key == NULL) {
      *ret = NULL;
    }
    if (typeid(T) == typeid(double)) {
      std::map<char*, double, char_cmp>::iterator iter =
          HighsDoubleOptions.find(key);
      if (iter != HighsDoubleOptions.end()) {
        *ret = *(T*)((void*)&iter->second);
      } else {
        *ret = NULL;
      }
    } else if (typeid(T) == typeid(int)) {
      std::map<char*, int, char_cmp>::iterator iter = HighsIntOptions.find(key);
      if (iter != HighsIntOptions.end()) {
        *ret = *(T*)((void*)&iter->second);
      } else {
        *ret = NULL;
      }
    } else if (typeid(T) == typeid(bool)) {
      std::map<char*, bool, char_cmp>::iterator iter =
          HighsBoolOptions.find(key);
      if (iter != HighsBoolOptions.end()) {
        *ret = *(T*)((void*)&iter->second);
      } else {
        *ret = NULL;
      }
    } else if (typeid(T) == typeid(FILE*)) {
      std::map<char*, FILE*, char_cmp>::iterator iter =
          HighsFileOptions.find(key);
      if (iter != HighsFileOptions.end()) {
        *ret = *(T*)((void*)&iter->second);
      } else {
        *ret = NULL;
      }
    }
  }

  template <class T>
  static void setValue(char* key, T value) {
    if (typeid(value) == typeid(double)) {
      std::map<char*, double, char_cmp>::iterator it =
          HighsDoubleOptions.find(key);
      if (it == HighsDoubleOptions.end()) {
        HighsDoubleOptions.insert(
            std::map<char*, double, char_cmp>::value_type(key, value));
      } else {
        it->second = value;
      }
    } else if (typeid(value) == typeid(int)) {
      std::map<char*, int, char_cmp>::iterator it = HighsIntOptions.find(key);
      if (it == HighsIntOptions.end()) {
        HighsIntOptions.insert(
            std::map<char*, int, char_cmp>::value_type(key, value));
      } else {
        it->second = value;
      }
    } else if (typeid(value) == typeid(bool)) {
      std::map<char*, bool, char_cmp>::iterator it = HighsBoolOptions.find(key);
      if (it == HighsBoolOptions.end()) {
        HighsBoolOptions.insert(
            std::map<char*, bool, char_cmp>::value_type(key, value));
      } else {
        it->second = value;
      }
    }
  }

  template <class T>
  static void setPtrValue(char* key, T value) {
    if (typeid(value) == typeid(FILE*)) {
      std::map<char*, FILE*, char_cmp>::iterator it =
          HighsFileOptions.find(key);
      if (it == HighsFileOptions.end()) {
        HighsFileOptions.insert(
            std::map<char*, FILE*, char_cmp>::value_type(key, value));
      } else {
        it->second = value;
      }
    }
  }
};
#endif