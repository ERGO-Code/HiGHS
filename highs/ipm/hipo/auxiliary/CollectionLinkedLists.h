#ifndef HIGHS_COLLECTION_LINKED_LISTS_H
#define HIGHS_COLLECTION_LINKED_LISTS_H

#include <vector>

#include "util/HighsType.h"

namespace highs {

namespace folding {

// Collection of linked lists.
// See highs/ipm/basiclu/lu_list.h for an explanation.

class CollectionLinkedLists {
  std::vector<HighsInt> forward_;
  std::vector<HighsInt> backward_;
  std::vector<HighsInt> length_;

  HighsInt n_elem_{};
  HighsInt n_lists_{};

  void print() const;

 public:
  // initialise and clear
  void init(HighsInt n_elem, HighsInt n_lists);
  void clear(HighsInt list);
  void clear();

  // modify the lists
  void append(HighsInt elem, HighsInt list);
  void remove(HighsInt elem, HighsInt list);

  // read the lists
  HighsInt head(HighsInt list) const { return forward_[n_elem_ + list]; }
  HighsInt next(HighsInt elem) const { return forward_[elem]; }
  HighsInt tail(HighsInt list) const { return backward_[n_elem_ + list]; }
  HighsInt prev(HighsInt elem) const { return backward_[elem]; }
  HighsInt length(HighsInt list) const { return length_[list]; }
  bool cont(HighsInt v) const { return v < n_elem_; }
};

/*
To go through list i:

HighsInt v = head(i);
while (cont(v)){
  ...
  v = next(v);
}

or reverse

HighsInt v = tail(i);
while (cont(v)){
  ...
  v = prev(v);
}
*/

}  // namespace folding

}  // namespace highs

#endif