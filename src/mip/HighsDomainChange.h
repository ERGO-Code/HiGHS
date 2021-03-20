/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

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

  bool operator==(const HighsDomainChange& other) const {
    return boundtype == other.boundtype && column == other.column &&
           boundval == other.boundval;
  }
};

struct HighsSubstitution {
  int substcol;
  int staycol;
  double scale;
  double offset;
};

#endif
