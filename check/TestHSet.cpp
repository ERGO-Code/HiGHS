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
  REQUIRE(set.add(3));
  REQUIRE(set.add(8));
  REQUIRE(set.add(7));
  REQUIRE(set.add(5));
  //  set.print();
  set.clear();
  //  set.print();

  if (dev_run) printf("\nChange count illegally\n");
  REQUIRE(set.add(3));
  REQUIRE(set.add(8));
  REQUIRE(set.add(7));
  REQUIRE(set.add(5));
  set.count_ = 1;
  //  set.print();
  REQUIRE(!set.debug());
  set.clear();
  
  REQUIRE(set.add(3));
  REQUIRE(set.add(8));
  REQUIRE(set.add(7));
  REQUIRE(set.add(5));
  //  set.print();
  REQUIRE(set.debug());

  if (dev_run) printf("\nChange entry illegally\n");
  set.value_[0] = 1;  
  //  set.print();
  REQUIRE(!set.debug());
  set.clear();

  REQUIRE(set.add(3));
  REQUIRE(set.add(8));
  REQUIRE(set.add(7));
  REQUIRE(set.add(5));
  //  set.print();

  if (dev_run) printf("\nCreate duplicate\n");
  set.value_[1] = 3;  
  //  set.print();
  REQUIRE(!set.debug());
  set.clear();



}
