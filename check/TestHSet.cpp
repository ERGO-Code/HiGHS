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
  bool debug = true;
  // Don't allow asserts so that debugging can be tested without
  // throwing an assert!
  bool allow_assert = false;
  FILE* output = NULL;
  if (dev_run) output = stdout;
  REQUIRE(set.setup(size, max_value, debug, allow_assert, output));
  //  set.print();
  REQUIRE(!set.add(-1));
  REQUIRE(set.add(0));
  REQUIRE(set.add(1));
  REQUIRE(!set.add(1));
  //  set.print();
  REQUIRE(set.add(3));
  REQUIRE(set.add(8));
  REQUIRE(set.add(7));
  REQUIRE(set.add(5));
  //  set.print();
  REQUIRE(set.add(11));
  //  set.print();
  REQUIRE(set.add(99));
  //  set.print();
  REQUIRE(!set.remove(4));
  REQUIRE(set.remove(3));
  //  set.print();
  REQUIRE(set.remove(11));
  REQUIRE(set.remove(0));
  REQUIRE(set.remove(8));
  REQUIRE(set.remove(5));
  REQUIRE(set.remove(99));
  REQUIRE(set.remove(1));
  REQUIRE(set.remove(7));
  //  set.print();
  set.clear();
  REQUIRE(set.add(3));
  REQUIRE(set.add(8));
  REQUIRE(set.add(7));
  REQUIRE(set.add(5));
  //  set.print();
  const int& count = set.count();
  const vector<int>& value = set.value();
  int value_size = value.size();
  REQUIRE(count==4);
  REQUIRE(value_size==8);
  REQUIRE(value[0]==3);
  REQUIRE(value[1]==8);
  REQUIRE(value[2]==7);
  REQUIRE(value[3]==5);
  if (dev_run) {
    printf("Set(%d, %d)\nValues: ", value_size, count);
    for (int ix=0; ix < count; ix++) printf(" %d", value[ix]);
    printf("\n");
  }

}
