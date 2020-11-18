#ifndef HIGHS_IMPLICATIONS_H_
#define HIGHS_IMPLICATIONS_H_

#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>

#include "mip/HighsDomain.h"
#include "mip/HighsDomainChange.h"

class HighsCliqueTable;

class HighsImplications {
  std::vector<HighsDomainChange> implications;

  struct Implics {
    int start;
    int num;
  };
  std::vector<Implics> implicationmap;

  bool computeImplications(int col, bool val);

 public:
  HighsDomain& globaldomain;
  HighsCliqueTable& cliquetable;
  HighsImplications(HighsDomain& globaldom, HighsCliqueTable& cliquetable)
      : globaldomain(globaldom), cliquetable(cliquetable) {
    implicationmap.resize(2 * globaldom.colLower_.size(), {-1, 0});
  }

  int getImplications(int col, bool val,
                      const HighsDomainChange*& implicationsstart,
                      bool& infeasible) {
    int loc = 2 * col + val;
    if (implicationmap[loc].start == -1) {
      infeasible = computeImplications(col, val);

      if (infeasible) return 0;
    } else
      infeasible = false;

    assert(implicationmap[loc].start != -1);

    implicationsstart = &implications[implicationmap[loc].start];

    return implicationmap[loc].num;
  }

  bool implicationsCached(int col, bool val) {
    int loc = 2 * col + val;
    return implicationmap[loc].start != -1;
  }
};

#endif