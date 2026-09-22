#ifndef HIPO_COLLECTION_LINKED_LISTS_H
#define HIPO_COLLECTION_LINKED_LISTS_H

#include <vector>

#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

// Collection of linked lists.
// See highs/ipm/basiclu/lu_list.h for an explanation.

class CollectionLinkedLists {
  std::vector<Int> forward_;
  std::vector<Int> backward_;
  std::vector<Int> length_;

  Int n_elem_{};
  Int n_lists_{};

  void print() const;

 public:
  // initialise and clear
  void init(Int n_elem, Int n_lists);
  void clear(Int list);
  void clear();

  // modify the lists
  void append(Int elem, Int list);
  void remove(Int elem, Int list);

  // read the lists
  Int head(Int list) const { return forward_[n_elem_ + list]; }
  Int next(Int elem) const { return forward_[elem]; }
  Int tail(Int list) const { return backward_[n_elem_ + list]; }
  Int prev(Int elem) const { return backward_[elem]; }
  Int length(Int list) const { return length_[list]; }
  bool cont(Int v) const { return v < n_elem_; }

  void test();
};

/*
To go through list i:

Int v = head(i);
while (cont(v)){
  ...
  v = next(v);
}

or reverse

Int v = tail(i);
while (cont(v)){
  ...
  v = prev(v);
}
*/

}  // namespace hipo

#endif