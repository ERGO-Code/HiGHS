#include "CollectionLinkedLists.h"

#include <cassert>
#include <cstdio>

namespace highs {

namespace folding {

void CollectionLinkedLists::init(HighsInt n_elem, HighsInt n_lists) {
  n_elem_ = n_elem;
  n_lists_ = n_lists;
  forward_.resize(n_elem + n_lists);
  backward_.resize(n_elem + n_lists);
  for (HighsInt i = 0; i < n_elem + n_lists; ++i) {
    forward_[i] = i;
    backward_[i] = i;
  }
  length_.assign(n_lists, 0);
}

void CollectionLinkedLists::clear(HighsInt list) {
  /*HighsInt current = forward_[n_elem_ + list];
  while (current < n_elem_) {
    const HighsInt temp = forward_[current];
    forward_[current] = current;
    backward_[current] = current;
    current = temp;
  }
  */

  forward_[n_elem_ + list] = n_elem_ + list;
  backward_[n_elem_ + list] = n_elem_ + list;
  length_[list] = 0;
}

void CollectionLinkedLists::append(HighsInt elem, HighsInt list) {
  const HighsInt temp = backward_[n_elem_ + list];
  backward_[n_elem_ + list] = elem;
  backward_[elem] = temp;
  forward_[temp] = elem;
  forward_[elem] = n_elem_ + list;
  length_[list]++;
}

void CollectionLinkedLists::remove(HighsInt elem, HighsInt list) {
  forward_[backward_[elem]] = forward_[elem];
  backward_[forward_[elem]] = backward_[elem];
  forward_[elem] = elem;
  backward_[elem] = elem;
  length_[list]--;
}

void CollectionLinkedLists::print() const {
  printf("H:  ");
  for (HighsInt i = 0; i < n_lists_; ++i) printf("%3d", forward_[n_elem_ + i]);
  printf("\n");

  printf("T:  ");
  for (HighsInt i = 0; i < n_lists_; ++i) printf("%3d", backward_[n_elem_ + i]);
  printf("\n");

  printf("L:  ");
  for (HighsInt i = 0; i < n_lists_; ++i) printf("%3d", length_[i]);
  printf("\n");

  printf("N:  ");
  for (HighsInt i = 0; i < n_elem_; ++i) printf("%3d", forward_[i]);
  printf("\n");

  printf("P:  ");
  for (HighsInt i = 0; i < n_elem_; ++i) printf("%3d", backward_[i]);
  printf("\n");

  printf("\n\n");
}

}  // namespace folding

}  // namespace highs