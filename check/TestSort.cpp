#include "catch.hpp"
#include "util/HighsSort.h"
#include "util/HighsRandom.h"

#include <vector>
#include <algorithm>

// No commas in test case name.
TEST_CASE("HiGHS_sort", "[highs_data]") {

  int num_values = 10;
  std::vector<int> indices;
  std::vector<int> int_values;
  std::vector<double> double_values;
  std::vector<double> original_double_values;
  indices.resize(1+num_values);
  int_values.resize(num_values);
  double_values.resize(1+num_values);
  original_double_values.resize(1+num_values);

  // Set up a vector of random number and their corresponding indices
  HighsRandom random;
  for (int ix = 0; ix < num_values; ix++) {
    double_values[1+ix] = random.fraction();
    original_double_values[1+ix] = double_values[1+ix];
    indices[1+ix] = ix;
  }
  // Sort the vector of random number and their corresponding indices
  maxheapsort(&double_values[0], &indices[0], num_values);

  // Check that the random numbers are ascending and that the indices
  // point from the original values to their new positions
  bool error0 = false;
  bool error1 = false;
  double previous_double = -1e200;
  for (int ix = 0; ix < num_values; ix++) {
    //    printf("%2d: %2d %12g %12g\n", ix, indices[1+ix], double_values[1+ix], original_double_values[1+ix]);
    error0 = error0 || double_values[1+ix] < previous_double;
    previous_double = double_values[1+ix];
    error1 = error1 || double_values[1+ix] == original_double_values[indices[1+ix]];
  }
  
  REQUIRE(error0 == false);  
  REQUIRE(error1 == false);  

  // Use the indices of the previous sort as a vector of integers to sort
  for (int ix = 0; ix < num_values; ix++) {
    double_values[ix] = double_values[ix+1];
    int_values[ix] = indices[1+ix];
  }
  std::make_heap(int_values.begin(), int_values.end());
  std::sort_heap(int_values.begin(), int_values.end());
  //  maxheapsort(&int_values[0], num_values);

  bool ok;
  // Check that the values in the vector of doubles are ascending - can do strict test
  ok = increasing_set_ok(&double_values[0], num_values, 0, 1, true);
  REQUIRE(ok == true);  

  // Check that the values in the vector of integers are ascending - maybe can't do strict test
  ok = increasing_set_ok(&int_values[0], num_values, 0, num_values, false);
  REQUIRE(ok == true);  
}
