#include "HighsLinkedLists.h"

#include <cassert>
#include <cstdio>

void HighsLinkedLists::init(HighsInt n_elem, HighsInt n_lists) {
  n_elem_ = n_elem;
  n_lists_ = n_lists;
  forward_.resize(n_elem + n_lists);
  backward_.resize(n_elem + n_lists);
  length_.assign(n_lists, 0);
  for (HighsInt i = 0; i < n_lists; ++i) {
    forward_[n_elem_ + i] = n_elem_ + i;
    backward_[n_elem_ + i] = n_elem_ + i;
  }
}

void HighsLinkedLists::clear(HighsInt list) {
  forward_[n_elem_ + list] = n_elem_ + list;
  backward_[n_elem_ + list] = n_elem_ + list;
  length_[list] = 0;
}

void HighsLinkedLists::append(HighsInt elem, HighsInt list) {
  const HighsInt temp = backward_[n_elem_ + list];
  backward_[n_elem_ + list] = elem;
  backward_[elem] = temp;
  forward_[temp] = elem;
  forward_[elem] = n_elem_ + list;
  length_[list]++;
}

void HighsLinkedLists::remove(HighsInt elem, HighsInt list) {
  forward_[backward_[elem]] = forward_[elem];
  backward_[forward_[elem]] = backward_[elem];
  forward_[elem] = elem;
  backward_[elem] = elem;
  length_[list]--;
}
