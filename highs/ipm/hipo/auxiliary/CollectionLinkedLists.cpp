#include "CollectionLinkedLists.h"

#include <cassert>

namespace hipo {

void CollectionLinkedLists::init(Int n_elem, Int n_lists) {
  n_elem_ = n_elem;
  n_lists_ = n_lists;
  forward_.resize(n_elem + n_lists);
  backward_.resize(n_elem + n_lists);
  for (Int i = 0; i < n_elem + n_lists; ++i) {
    forward_[i] = i;
    backward_[i] = i;
  }
  length_.assign(n_lists, 0);
}

void CollectionLinkedLists::clear(Int list) {
  Int current = forward_[n_elem_ + list];
  while (current < n_elem_) {
    const Int temp = forward_[current];
    forward_[current] = current;
    backward_[current] = current;
    current = temp;
  }
  forward_[n_elem_ + list] = n_elem_ + list;
  backward_[n_elem_ + list] = n_elem_ + list;
  length_[list] = 0;
}

void CollectionLinkedLists::clear() {
  for (Int i = 0; i < n_lists_; ++i) clear(i);
}

void CollectionLinkedLists::append(Int elem, Int list) {
  const Int temp = backward_[n_elem_ + list];
  backward_[n_elem_ + list] = elem;
  backward_[elem] = temp;
  forward_[temp] = elem;
  forward_[elem] = n_elem_ + list;
  length_[list]++;
}

void CollectionLinkedLists::remove(Int elem, Int list) {
  forward_[backward_[elem]] = forward_[elem];
  backward_[forward_[elem]] = backward_[elem];
  forward_[elem] = elem;
  backward_[elem] = elem;
  length_[list]--;
}

void CollectionLinkedLists::test() {
  init(16, 4);
  print();

  append(0, 0);
  append(7, 2);
  append(1, 0);
  append(4, 0);
  append(2, 1);
  append(3, 3);
  append(6, 3);
  append(5, 1);
  append(10, 0);
  append(8, 2);
  append(9, 2);
  append(11, 3);

  append(12, 1);
  append(13, 3);
  remove(12, 1);
  append(12, 2);

  clear(1);
  print();

  exit(1);
}

void CollectionLinkedLists::print() const {
  printf("H:  ");
  for (Int i = 0; i < n_lists_; ++i) printf("%3d", forward_[n_elem_ + i]);
  printf("\n");

  printf("T:  ");
  for (Int i = 0; i < n_lists_; ++i) printf("%3d", backward_[n_elem_ + i]);
  printf("\n");

  printf("L:  ");
  for (Int i = 0; i < n_lists_; ++i) printf("%3d", length_[i]);
  printf("\n");

  printf("N:  ");
  for (Int i = 0; i < n_elem_; ++i) printf("%3d", forward_[i]);
  printf("\n");

  printf("P:  ");
  for (Int i = 0; i < n_elem_; ++i) printf("%3d", backward_[i]);
  printf("\n");

  printf("\n\n");
}

}  // namespace hipo