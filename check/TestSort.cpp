#include <algorithm>
#include <vector>

#include "catch.hpp"
#include "lp_data/HConst.h"
#include "util/HighsRandom.h"
#include "util/HighsSort.h"

using std::vector;

const bool dev_run = false;
const bool use_ekk = true;

// No commas in test case name.
void getRandomValues(const int num_values, vector<double>& values,
                     vector<int>& indices) {
  // Set up a vector of random number and their corresponding indices
  HighsRandom random;
  for (int ix = 0; ix < num_values; ix++) {
    values[1 + ix] = random.fraction();
    indices[1 + ix] = ix;
  }
}

void doFullSort(const int num_values, vector<double>& values,
                vector<int>& indices) {
  // Sort the vector of random number and their corresponding indices
  maxheapsort(&values[0], &indices[0], num_values);
}

void doAddSort(int& num_values_sorted, const int& max_num_values_sorted,
               vector<double>& best_double_values, vector<int>& best_indices,
               const int num_values, const vector<double>& values,
               const vector<int>& indices) {
  num_values_sorted = 0;
  for (int ix = 1; ix <= num_values; ix++) {
    addToDecreasingHeap(num_values_sorted, max_num_values_sorted,
                        best_double_values, best_indices, values[ix],
                        indices[ix]);
  }
  sortDecreasingHeap(num_values_sorted, best_double_values, best_indices);
}

void reportValuesIndices(const int num_values, const vector<double>& values,
                         const vector<int>& indices) {
  if (!dev_run) return;
  printf("\n  Ix      Value Index\n");
  for (int ix = 1; ix <= num_values; ix++) {
    printf("%4d %10.8f  %4d\n", ix, values[ix], indices[ix]);
  }
}

void checkIncreasingSort(const int num_sorted, const vector<double>& values,
                         const vector<int>& indices,
                         const vector<double>& original_values) {
  // Check that the random numbers are ascending and that the indices
  // point from the original values to their new positions
  bool error0 = false;
  bool error1 = false;
  double previous = -HIGHS_CONST_INF;
  for (int ix = 0; ix < num_sorted; ix++) {
    if (values[1 + ix] < previous) {
      printf("Values[%2d] = %f5.4 < %f5.4 = previous\n", 1 + ix, values[1 + ix],
             previous);
      error0 = true;
    }
    previous = values[1 + ix];
    if (values[1 + ix] != original_values[1 + indices[1 + ix]]) {
      printf("Values[%2d] = %f5.4 != %f5.4 = original_values[indices[%2d]]\n",
             1 + ix, values[1 + ix], original_values[indices[1 + ix]], 1 + ix);
      error1 = true;
    }
  }

  REQUIRE(error0 == false);
  REQUIRE(error1 == false);
}

void checkDecreasingSort(const int num_sorted, const vector<double>& values,
                         const vector<int>& indices,
                         const vector<double>& original_values) {
  // Check that the random numbers are ascending and that the indices
  // point from the original values to their new positions
  bool error0 = false;
  bool error1 = false;
  double previous = HIGHS_CONST_INF;
  for (int ix = 0; ix < num_sorted; ix++) {
    if (values[1 + ix] > previous) {
      printf("Values[%2d] = %f5.4 < %f5.4 = previous\n", 1 + ix, values[1 + ix],
             previous);
      error0 = true;
    }
    previous = values[1 + ix];
    if (values[1 + ix] != original_values[1 + indices[1 + ix]]) {
      printf("Values[%2d] = %f5.4 != %f5.4 = original_values[indices[%2d]]\n",
             1 + ix, values[1 + ix], original_values[indices[1 + ix]], 1 + ix);
      error1 = true;
    }
  }

  REQUIRE(error0 == false);
  REQUIRE(error1 == false);
}

TEST_CASE("HiGHS_sort", "[highs_data]") {
  int num_values = 10;
  vector<int> indices;
  vector<int> int_values;
  vector<double> double_values;
  vector<double> original_double_values;
  indices.resize(1 + num_values);
  int_values.resize(num_values);
  double_values.resize(1 + num_values);
  original_double_values.resize(1 + num_values);

  getRandomValues(num_values, double_values, indices);
  reportValuesIndices(num_values, double_values, indices);
  original_double_values = double_values;
  doFullSort(num_values, double_values, indices);
  reportValuesIndices(num_values, double_values, indices);
  checkIncreasingSort(num_values, double_values, indices,
                      original_double_values);

  // Use the indices of the previous sort as a vector of integers to sort
  for (int ix = 0; ix < num_values; ix++) {
    double_values[ix] = double_values[ix + 1];
    int_values[ix] = indices[1 + ix];
  }
  std::make_heap(int_values.begin(), int_values.end());
  std::sort_heap(int_values.begin(), int_values.end());
  //  maxheapsort(&int_values[0], num_values);

  bool ok;
  // Check that the values in the vector of doubles are ascending - can do
  // strict test
  ok = increasingSetOk(&double_values[0], num_values, 0, 1, true);
  REQUIRE(ok == true);

  // Check that the values in the vector of integers are ascending - maybe can't
  // do strict test
  ok = increasingSetOk(&int_values[0], num_values, 0, num_values, false);
  REQUIRE(ok == true);

  num_values = 14;
  vector<int> set;
  vector<double> lb;
  vector<double> ub;
  set.resize(num_values);
  lb.resize(num_values);
  ub.resize(num_values);

  set[0] = 6;
  lb[0] = 60;
  ub[0] = 60;
  set[1] = 7;
  lb[1] = 6;
  ub[1] = 6;
  set[2] = 8;
  lb[2] = 60;
  ub[2] = 60;
  set[3] = 13;
  lb[3] = 600;
  ub[3] = 1200;
  set[4] = 4;
  lb[4] = 70;
  ub[4] = 70;
  set[5] = 5;
  lb[5] = 16;
  ub[5] = 16;
  set[6] = 2;
  lb[6] = 70;
  ub[6] = 70;
  set[7] = 3;
  lb[7] = 7;
  ub[7] = 7;
  set[8] = 11;
  lb[8] = 200;
  ub[8] = 1400;
  set[9] = 0;
  lb[9] = 75;
  ub[9] = 75;
  set[10] = 1;
  lb[10] = 12;
  ub[10] = 12;
  set[11] = 14;
  lb[11] = 0;
  ub[11] = 1400;
  set[12] = 9;
  lb[12] = 6;
  ub[12] = 6;
  set[13] = 15;
  lb[13] = 600;
  ub[13] = 1200;

  vector<int> sorted_set = set;
  vector<double> sorted_lb;
  vector<double> sorted_ub;
  sorted_lb.resize(num_values);
  sorted_ub.resize(num_values);

  sortSetData(num_values, &sorted_set[0], &lb[0], &ub[0], NULL, &sorted_lb[0],
              &sorted_ub[0], NULL);

  int prev_ix = -HIGHS_CONST_I_INF;
  for (int k0 = 0; k0 < num_values; k0++) {
    int ix = sorted_set[k0];
    REQUIRE(ix >= prev_ix);
    int k1 = -HIGHS_CONST_I_INF;
    for (int k_1 = 0; k_1 < num_values; k_1++) {
      if (set[k_1] == ix) {
        k1 = k_1;
        break;
      }
    }
    REQUIRE(k1 > -HIGHS_CONST_I_INF);
    REQUIRE(sorted_lb[k0] == lb[k1]);
    REQUIRE(sorted_ub[k0] == ub[k1]);
  }

  num_values = 10;
  int max_num_values_sorted = 5;
  vector<int> best_indices;
  vector<double> best_double_values;
  indices.assign(1 + num_values, 0);
  double_values.assign(1 + num_values, 0);
  best_indices.assign(1 + max_num_values_sorted, 0);
  best_double_values.assign(1 + max_num_values_sorted, 0);
  int num_values_sorted;
  getRandomValues(num_values, double_values, indices);
  reportValuesIndices(num_values, double_values, indices);
  doAddSort(num_values_sorted, max_num_values_sorted, best_double_values,
            best_indices, num_values, double_values, indices);
  reportValuesIndices(num_values_sorted, best_double_values, best_indices);
  checkDecreasingSort(num_values_sorted, best_double_values, best_indices,
                      double_values);
}
