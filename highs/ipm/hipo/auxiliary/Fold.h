#ifndef HIPO_FOLDING_H
#define HIPO_FOLDING_H

#include "CollectionLinkedLists.h"

namespace hipo {

std::vector<Int> colourRefinement(const std::vector<Int> ptr, const std::vector<Int> adj);

void test_folding();

}  // namespace hipo

#endif