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

#include "util/HighsInt.h"

enum class HighsBoundType { Lower, Upper };

struct HighsDomainChange {
  HighsBoundType boundtype;
  HighsInt column;
  double boundval;

  bool operator<(const HighsDomainChange& other) const {
    if (column < other.column) return true;
    if (other.column < column) return false;
    if ((HighsInt)boundtype < (HighsInt)other.boundtype) return true;
    return false;
  }

  bool operator==(const HighsDomainChange& other) const {
    return boundtype == other.boundtype && column == other.column &&
           boundval == other.boundval;
  }
};

struct HighsSubstitution {
  HighsInt substcol;
  HighsInt staycol;
  double scale;
  double offset;
};

#endif
