#include "catch.hpp"
#include "util/HSet.h"

const bool dev_run = false;

TEST_CASE("HSet", "[highs_test_hset]") {
  const HighsInt size = 5;
  const HighsInt max_entry = 10;
  HSet set;
  REQUIRE(!set.setup(-1, max_entry));
  REQUIRE(!set.setup(0, max_entry));
  REQUIRE(!set.setup(size, -1));
  bool debug = true;
  // Don't allow asserts so that debugging can be tested without
  // triggering an assert!
  bool allow_assert = false;
  FILE* log_file = NULL;
  if (dev_run) log_file = stdout;
  bool output_flag = false;
  if (dev_run) output_flag = true;
  REQUIRE(
      set.setup(size, max_entry, output_flag, log_file, debug, allow_assert));
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
  // Confirm that 8 is in the set
  REQUIRE(set.in(8));
  REQUIRE(set.remove(8));
  // Confirm that 8 is not in the set
  REQUIRE(!set.in(8));
  // Confirm that entry too small is not in the set
  REQUIRE(!set.in(-1));
  // Confirm that entry too big is not in the set
  REQUIRE(!set.in(999));
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
  set.print();
  const HighsInt& count = set.count();
  const vector<HighsInt>& entry = set.entry();
  HighsInt entry_size = entry.size();
  REQUIRE(count == 4);
  REQUIRE(entry_size == 8);
  REQUIRE(entry[0] == 3);
  REQUIRE(entry[1] == 8);
  REQUIRE(entry[2] == 7);
  REQUIRE(entry[3] == 5);
  if (dev_run) {
    printf("Set(%" HIGHSINT_FORMAT ", %" HIGHSINT_FORMAT ")\nEntries: ",
           entry_size, count);
    for (HighsInt ix = 0; ix < count; ix++)
      printf(" %" HIGHSINT_FORMAT "", entry[ix]);
    printf("\n");
  }
}
