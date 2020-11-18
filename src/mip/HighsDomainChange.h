
#ifndef HIGHS_DOMAIN_CHANGE_H_
#define HIGHS_DOMAIN_CHANGE_H_

enum class HighsBoundType { Lower, Upper };

struct HighsDomainChange {
  HighsBoundType boundtype;
  int column;
  double boundval;
};

#endif