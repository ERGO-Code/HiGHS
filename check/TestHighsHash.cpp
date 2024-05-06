#include <numeric>

#include "HCheckConfig.h"
#include "catch.hpp"
#include "util/HighsHash.h"
#include "util/HighsHashTree.h"

TEST_CASE("Highs_log2i", "[util]") {
  // test 32 bit and 64 bit values whoes log2 value should be floored
  uint32_t x = 12345;
  uint64_t y = 1234567891011;
  REQUIRE(HighsHashHelpers::log2i(x) == 13);
  REQUIRE(HighsHashHelpers::log2i(y) == 40);

  // test 2^13 and 2^40 which should yield the same result
  x = 8192;
  y = 1099511627776;
  REQUIRE(HighsHashHelpers::log2i(x) == 13);
  REQUIRE(HighsHashHelpers::log2i(y) == 40);
}

TEST_CASE("Highs_HashTree", "[util]") {
  // test 32 bit and 64 bit values whoes log2 value should be floored
  HighsHashTree<int> htree;

  htree.insert(1);
  htree.insert(2);
  htree.insert(3);

  REQUIRE(!htree.contains(4));
  REQUIRE(htree.contains(1));
  REQUIRE(htree.contains(2));
  REQUIRE(htree.contains(3));

  htree.erase(2);
  REQUIRE(htree.contains(1));
  REQUIRE(!htree.contains(2));
  REQUIRE(htree.contains(3));

  REQUIRE(!htree.insert(1));
  REQUIRE(!htree.insert(3));

  constexpr int NUM_CHECK = 10000;
  for (int i = 0; i < NUM_CHECK; ++i) {
    htree.insert(i);
  }

  std::vector<int> test;
  std::vector<int> testReference;
  testReference.resize(NUM_CHECK);
  std::iota(testReference.begin(), testReference.end(), 0);

  htree.for_each([&](const HighsHashTableEntry<int>& entry) {
    test.push_back(entry.key());
    return false;
  });

  REQUIRE(test.size() == NUM_CHECK);
  std::sort(test.begin(), test.end());
  REQUIRE(std::equal(test.begin(), test.end(), testReference.begin()));

  for (int i = 0; i < NUM_CHECK; ++i) {
    REQUIRE(htree.contains(i));
  }

  testReference.clear();

  for (int i = 0; i < NUM_CHECK; ++i) {
    if ((i % 7) == 0) {
      REQUIRE(htree.contains(i));
      htree.erase(i);
      REQUIRE(!htree.contains(i));
    } else {
      testReference.push_back(i);
    }
  }

  test.clear();
  htree.for_each([&](const HighsHashTableEntry<int>& entry) {
    test.push_back(entry.key());
    return false;
  });

  REQUIRE(test.size() == testReference.size());
  std::sort(test.begin(), test.end());
  REQUIRE(std::equal(test.begin(), test.end(), testReference.begin()));

  HighsHashTree<int> htree2;
  for (int i = 0; i < NUM_CHECK; ++i) {
    if (htree.contains(i)) {
      REQUIRE((i % 7) != 0);
    } else {
      REQUIRE((i % 7) == 0);
      htree2.insert(i);
    }
  }

  const auto* commonElement = htree.find_common(htree2);
  REQUIRE(commonElement == nullptr);

  for (int i = 0; i < NUM_CHECK; ++i) {
    if (i % 7 != 0) {
      htree2.insert(i);
      commonElement = htree.find_common(htree2);
      REQUIRE(commonElement != nullptr);
      REQUIRE(commonElement->key() == i);
      htree.erase(i);

      commonElement = htree.find_common(htree2);
      REQUIRE(commonElement == nullptr);
    }
  }

  HighsHashTree<int> htree3 = htree2;

  for (int i = 0; i < NUM_CHECK; ++i) {
    const int* a = htree2.find(i);
    const int* b = htree3.find(i);

    if (a == b) {
      REQUIRE(a == nullptr);
    } else {
      REQUIRE(*a == *b);
    }
  }
}
