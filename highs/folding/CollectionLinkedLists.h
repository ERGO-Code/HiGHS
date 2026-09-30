#ifndef HIGHS_COLLECTION_LINKED_LISTS_H
#define HIGHS_COLLECTION_LINKED_LISTS_H

#include <vector>

#include "util/HighsType.h"

namespace highs {

namespace folding {

// Collection of linked lists.
// See highs/ipm/basiclu/lu_list.h for an explanation.

class LinkedLists {
  HighsInt n_elem_{};
  HighsInt n_lists_{};

  std::vector<HighsInt> forward_;
  std::vector<HighsInt> backward_;
  std::vector<HighsInt> length_;

 public:
  void init(HighsInt n_elem, HighsInt n_lists);
  void clear(HighsInt list);

  void append(HighsInt elem, HighsInt list);
  void remove(HighsInt elem, HighsInt list);
  HighsInt length(HighsInt list) const { return length_[list]; }

  const HighsInt& head(HighsInt list) const { return forward_[n_elem_ + list]; }
  const HighsInt& next(HighsInt elem) const { return forward_[elem]; }
  const HighsInt& tail(HighsInt list) const {
    return backward_[n_elem_ + list];
  }
  const HighsInt& prev(HighsInt elem) const { return backward_[elem]; }

  // Define iterator for range-based loop:
  //  for (HighsInt v : list(i))
  //
  // and reverse:
  //  for (HighsInt v : listReverse(i))
  //
  template <bool reverse>
  struct List {
    const LinkedLists* owner;
    const HighsInt list;

    struct Iterator {
      const LinkedLists* owner;
      HighsInt current;

      HighsInt operator*() const { return current; }
      Iterator& operator++() {
        current = reverse ? owner->prev(current) : owner->next(current);
        return *this;
      }
      bool operator!=(const Iterator& o) const { return current != o.current; }
      bool operator==(const Iterator& o) const { return current == o.current; }
    };

    Iterator begin() const {
      return {owner, reverse ? owner->tail(list) : owner->head(list)};
    }
    Iterator end() const { return {owner, owner->n_elem_ + list}; }
  };

  List<false> list(HighsInt l) const { return List<false>{this, l}; }
  List<true> listReverse(HighsInt l) const { return List<true>{this, l}; }
};

}  // namespace folding

}  // namespace highs

#endif