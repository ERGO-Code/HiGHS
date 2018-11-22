#ifndef HIGHSOPTIONS_H
#define HIGHSOPTIONS_H

#include <string.h>
#include <map>

struct char_cmp {
  bool operator()(const char* a, const char* b) const {
    return strcmp(a, b) < 0;
  }
};

std::map<char*, double, char_cmp> doubleOptions;
std::map<char*, int, char_cmp> intOptions;
std::map<char*, bool, char_cmp> boolOptions;

class HighsStringOptions {
 public:
  template <class T>
  static T getValue(char* key) {
    if (key == NULL) {
      T t;
      return t;
    }
    if (typeid(T) == typeid(double)) {
      std::map<char*, double, char_cmp>::iterator iter =
          doubleOptions.find(key);
      return iter->second;
    } else if (typeid(T) == typeid(int)) {
      std::map<char*, int, char_cmp>::iterator iter = intOptions.find(key);
      return iter->second;
    } else if (typeid(T) == typeid(bool)) {
      std::map<char*, bool, char_cmp>::iterator iter = boolOptions.find(key);
      return iter->second;
    }
  }

  template <class T>
  static void setValue(char* key, T value) {
    // does nothing if key already exists
    if (typeid(value) == typeid(double)) {
      doubleOptions.insert(
          std::map<char*, double, char_cmp>::value_type(key, value));
    } else if (typeid(value) == typeid(int)) {
      intOptions.insert(std::map<char*, int, char_cmp>::value_type(key, value));
    } else if (typeid(value) == typeid(bool)) {
      boolOptions.insert(
          std::map<char*, bool, char_cmp>::value_type(key, value));
    }
  }
};
#endif