#ifndef SCAFFOLD_METHODS_HPP_
#define SCAFFOLD_METHODS_HPP_

#include <iostream>

namespace scaffold {

class ScaffoldUtils {
 public:
  static void scaffoldHello() {
    std::cout << "HiGHS scaffold linked." << std::endl << std::endl;
    return;
  }
};

}  // namespace scaffold

#endif