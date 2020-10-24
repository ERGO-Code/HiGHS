#include "catch.hpp"
#include "util/HSet.h"

const bool dev_run = true;

TEST_CASE("HSet", "[highs_test_hset]") {
  const int size = 5;
  const int max_value = 10;
  HSet set;
  REQUIRE(!set.setup(-1, max_value));
  REQUIRE(!set.setup(0, max_value));
  REQUIRE(!set.setup(size, -1));
  REQUIRE(set.setup(size, max_value));
  set.print();
  REQUIRE(!set.add(-1));
  REQUIRE(set.add(0));
  REQUIRE(set.add(1));
  REQUIRE(!set.add(1));
  set.print();
  REQUIRE(set.add(3));
  REQUIRE(set.add(8));
  REQUIRE(set.add(7));
  REQUIRE(set.add(5));
  set.print();
  REQUIRE(set.add(11));
  set.print();
  REQUIRE(!set.remove(4));
  REQUIRE(set.remove(3));
  set.print();
  REQUIRE(set.remove(11));
  REQUIRE(set.remove(0));
  REQUIRE(set.remove(8));
  REQUIRE(set.remove(5));
  REQUIRE(set.remove(1));
  REQUIRE(set.remove(7));
  set.print();
  REQUIRE(set.add(3));
  REQUIRE(set.add(8));
  REQUIRE(set.add(7));
  REQUIRE(set.add(5));
  set.print();
  set.clear();
  set.print();
}
