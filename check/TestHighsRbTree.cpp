
#include <numeric>

#include "catch.hpp"
#include "util/HighsRandom.h"
#include "util/HighsRbTree.h"

using namespace highs;

struct Node {
  HighsInt key;
  RbTreeLinks<HighsInt> links;
};

class MyRbTree;

namespace highs {
template <>
struct RbTreeTraits<MyRbTree> {
  using KeyType = HighsInt;
  using LinkType = HighsInt;
};
}  // namespace highs

class MyRbTree : public RbTree<MyRbTree> {
 public:
  HighsInt root = -1;
  std::vector<Node> nodes;

  MyRbTree() : RbTree<MyRbTree>(root) {}

  RbTreeLinks<HighsInt>& getRbTreeLinks(HighsInt node) {
    return nodes[node].links;
  }
  const RbTreeLinks<HighsInt>& getRbTreeLinks(HighsInt node) const {
    return nodes[node].links;
  }
  HighsInt getKey(HighsInt node) const { return nodes[node].key; }

  bool insert(HighsInt x) {
    std::pair<HighsInt, bool> p = find(x);

    if (p.second) return false;

    HighsInt newNode = nodes.size();
    nodes.emplace_back();
    nodes.back().key = x;
    link(newNode, p.first);
    return true;
  }

  void erase(HighsInt node) { unlink(node); }

  bool contains(HighsInt x) { return find(x).second; }
};

static void checkRbTree(MyRbTree& tree, HighsInt* expectedKeys,
                        HighsInt numExpectedKeys) {
  std::vector<HighsInt> keys;
  keys.reserve(numExpectedKeys);
  if (tree.root != -1) REQUIRE(tree.nodes[tree.root].links.isBlack());

  HighsInt x = tree.first();
  while (x != -1) {
    keys.push_back(tree.nodes[x].key);
    if (tree.nodes[x].links.isRed()) {
      HighsInt lChild = tree.nodes[x].links.child[0];
      HighsInt rChild = tree.nodes[x].links.child[1];
      if (lChild != -1) REQUIRE(tree.nodes[lChild].links.isBlack());
      if (rChild != -1) REQUIRE(tree.nodes[rChild].links.isBlack());
    }
    x = tree.successor(x);
    REQUIRE((HighsInt)keys.size() <= numExpectedKeys);
  }

  REQUIRE((HighsInt)keys.size() == numExpectedKeys);
  std::sort(expectedKeys, expectedKeys + numExpectedKeys);
  bool isOk = std::equal(keys.begin(), keys.end(), expectedKeys);
  REQUIRE(isOk);
}

TEST_CASE("HighsRbTree", "[util]") {
  std::vector<HighsInt> keys;
  keys.resize(1000);
  std::iota(keys.begin(), keys.end(), 0);

  HighsRandom rand;
  rand.shuffle(keys.data(), keys.size());
  MyRbTree rbTree;

  for (size_t i = 0; i < keys.size(); ++i) {
    HighsInt x = keys[i];
    bool inserted = rbTree.insert(x);
    REQUIRE(inserted);
    checkRbTree(rbTree, keys.data(), i + 1);
  }

  // randomly delete half of the elements and check the tree after each deletion

  for (size_t i = keys.size() - 1; i > keys.size() / 2; --i) {
    HighsInt k = rand.integer(i + 1);
    std::swap(keys[k], keys[i]);
    HighsInt x = keys[i];
    std::pair<HighsInt, bool> node = rbTree.find(x);
    REQUIRE(node.second);
    rbTree.erase(node.first);
    checkRbTree(rbTree, keys.data(), i);
  }
}
