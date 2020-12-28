
#ifndef HIGHS_DOMAIN_CHANGE_H_
#define HIGHS_DOMAIN_CHANGE_H_

enum class HighsBoundType { Lower, Upper };

struct HighsDomainChange {
  HighsBoundType boundtype;
  int column;
  double boundval;

  bool operator<(const HighsDomainChange& other) const {
    if (column < other.column) return true;
    if (other.column < column) return false;
    if ((int)boundtype < (int)other.boundtype) return true;
    return false;
  }
};

struct HighsSubstitution {
  int substcol;
  int staycol;
  double scale;
  double offset;
};

#endif