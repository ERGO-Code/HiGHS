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

class HighsStringOptions {
 private:
  static std::map<char*, double, char_cmp> doubleOptions;
  static std::map<char*, int, char_cmp> intOptions;
  static std::map<char*, bool, char_cmp> boolOptions;
  static std::map<char*, FILE*, char_cmp> fileOptions;

 public:
  template <class T>
  static void getValue(char* key, T* ret) {
    if (key == NULL) {
      *ret = NULL;
    }
    if (typeid(T) == typeid(double)) {
      std::map<char*, double, char_cmp>::iterator iter =
          HighsStringOptions::doubleOptions.find(key);
      if (iter != HighsStringOptions::doubleOptions.end()) {
        *ret = *(T*)((void*)&iter->second);
      } else {
        *ret = NULL;
      }
    } else if (typeid(T) == typeid(int)) {
      std::map<char*, int, char_cmp>::iterator iter =
          HighsStringOptions::intOptions.find(key);
      if (iter != HighsStringOptions::intOptions.end()) {
        *ret = *(T*)((void*)&iter->second);
      } else {
        *ret = NULL;
      }
    } else if (typeid(T) == typeid(bool)) {
      std::map<char*, bool, char_cmp>::iterator iter =
          HighsStringOptions::boolOptions.find(key);
      if (iter != HighsStringOptions::boolOptions.end()) {
        *ret = *(T*)((void*)&iter->second);
      } else {
        *ret = NULL;
      }
    } else if (typeid(T) == typeid(FILE*)) {
      std::map<char*, FILE*, char_cmp>::iterator iter =
          HighsStringOptions::fileOptions.find(key);
      if (iter != HighsStringOptions::fileOptions.end()) {
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
          HighsStringOptions::doubleOptions.find(key);
      if (it == HighsStringOptions::doubleOptions.end()) {
        HighsStringOptions::doubleOptions.insert(
            std::map<char*, double, char_cmp>::value_type(key, value));
      } else {
        it->second = value;
      }
    } else if (typeid(value) == typeid(int)) {
      std::map<char*, int, char_cmp>::iterator it =
          HighsStringOptions::intOptions.find(key);
      if (it == HighsStringOptions::intOptions.end()) {
        HighsStringOptions::intOptions.insert(
            std::map<char*, int, char_cmp>::value_type(key, value));
      } else {
        it->second = value;
      }
    } else if (typeid(value) == typeid(bool)) {
      std::map<char*, bool, char_cmp>::iterator it =
          HighsStringOptions::boolOptions.find(key);
      if (it == HighsStringOptions::boolOptions.end()) {
        HighsStringOptions::boolOptions.insert(
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
          HighsStringOptions::fileOptions.find(key);
      if (it == HighsStringOptions::fileOptions.end()) {
        HighsStringOptions::fileOptions.insert(
            std::map<char*, FILE*, char_cmp>::value_type(key, value));
      } else {
        it->second = value;
      }
    }
  }
};

std::map<char*, double, char_cmp> HighsStringOptions::doubleOptions;
std::map<char*, int, char_cmp> HighsStringOptions::intOptions;
std::map<char*, bool, char_cmp> HighsStringOptions::boolOptions;
std::map<char*, FILE*, char_cmp> HighsStringOptions::fileOptions;

#endif