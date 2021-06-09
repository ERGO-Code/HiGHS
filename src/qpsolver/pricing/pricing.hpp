#ifndef __SRC_LIB_PRICING_HPP__
#define __SRC_LIB_PRICING_HPP__

#include "../vector.hpp"
#include "../matrix.hpp"

class Pricing {  
public:
   virtual int price(const Vector& x, const Vector& gradient) = 0;
   virtual void update_weights(const Vector& aq, const Vector& ep, int p, int q) = 0;
};

#endif
