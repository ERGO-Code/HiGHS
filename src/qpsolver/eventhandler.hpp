#ifndef __SRC_LIB_EVENTHANDLER_HPP__
#define __SRC_LIB_EVENTHANDLER_HPP__

#include <vector>

template<typename T> // T: void (*fncptr)(int, double)
class Eventhandler {

   std::vector<void(*)(T)> subscribers;

public:
   void subscribe(void(*subscriber)(T)) {
      subscribers.push_back(subscriber);
   }

   void fire(T args) {
      for (void(*fun)(T) : subscribers) {
         fun(args);
      }
   }
};

#endif
