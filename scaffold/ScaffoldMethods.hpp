#ifndef SCAFFOLD_METHODS_HPP_
#define SCAFFOLD_METHODS_HPP_

#include <iostream>

namespace scaffold {

class ScaffoldUtils {
 public:
  static void printHelloHighssScaffold() {
    std::cout << "HiGHS scaffold linked." << std::endl;
    return;
  }
};

}  // namespace scaffold

#endif