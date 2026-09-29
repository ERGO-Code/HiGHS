#ifndef HIGHS_COLLECTION_LINKED_LISTS_H
#define HIGHS_COLLECTION_LINKED_LISTS_H

#include <vector>

#include "util/HighsType.h"

namespace highs {

namespace folding {

// Collection of linked lists.
// See highs/ipm/basiclu/lu_list.h for an explanation.

class LinkedLists {
  std::vector<HighsInt> forward_;
  std::vector<HighsInt> backward_;
  std::vector<HighsInt> length_;

  HighsInt n_elem_{};
  HighsInt n_lists_{};

  void print() const;

 public:
  void init(HighsInt n_elem, HighsInt n_lists);
  void clear(HighsInt list);

  HighsInt n() const { return n_elem_; }
  void append(HighsInt elem, HighsInt list);
  void remove(HighsInt elem, HighsInt list);
  HighsInt length(HighsInt list) const { return length_[list]; }

  const HighsInt& head(HighsInt list) const { return forward_[n_elem_ + list]; }
  const HighsInt& next(HighsInt elem) const { return forward_[elem]; }

  // Define iterator for range-based loop:
  //  for (HighsInt v : list(i))
  //
  struct List {
    const LinkedLists* owner;
    const HighsInt list;

    struct Iterator {
      const LinkedLists* owner;
      HighsInt current;

      HighsInt operator*() const { return current; }
      Iterator& operator++() {
        current = owner->next(current);
        return *this;
      }
      bool operator!=(const Iterator& o) const { return current != o.current; }
      bool operator==(const Iterator& o) const { return current == o.current; }
    };

    Iterator begin() const { return {owner, owner->head(list)}; }
    Iterator end() const { return {owner, owner->n() + list}; }
  };

  List list(HighsInt l) const { return {this, l}; }
};

}  // namespace folding

}  // namespace highs

#endif