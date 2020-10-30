/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsSort.cpp
 * @brief Sorting routines for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include <cstddef>
#include "lp_data/HConst.h"
#include "util/HighsSort.h"

using std::vector;

void addToDecreasingHeap(int& n, int mx_n,
			 vector<double>& heap_v,
			 vector<int>& heap_ix,
			 const double v,
			 const int ix) {
  int cd_p, pa_p;
  if (n < mx_n) {
    // The heap is not full so put the new value at the bottom of the
    // heap and let it rise up to its correct level.
    n++;
    cd_p = n;
    pa_p = cd_p/2;
    // 10
    for (;;) {
      if (pa_p > 0) {
	if (v < heap_v[pa_p]) {
	  heap_v[cd_p] = heap_v[pa_p];
	  heap_ix[cd_p] = heap_ix[pa_p];
	  cd_p = pa_p;
	  pa_p = pa_p/2;
	  continue;
	}
      }
      break;
    }
    heap_v[cd_p] = v;
    heap_ix[cd_p] = ix;
  } else if (v > heap_v[1]) {
    // The heap is full so replace the least value with the new value
    // and let it sink down to its correct level.
    pa_p = 1;
    cd_p = pa_p + pa_p;
    // 20
    for (;;) {
      if (cd_p <= n) {
	if (cd_p < n) {
	  if (heap_v[cd_p] > heap_v[cd_p+1]) cd_p++;
	}
	if (v > heap_v[cd_p]) {
	  heap_v[pa_p] = heap_v[cd_p];
	  heap_ix[pa_p] = heap_ix[cd_p];
	  pa_p = cd_p;
	  cd_p = cd_p + cd_p;
	  continue;
	}
      }
      break;
    }
    heap_v[pa_p] = v;
    heap_ix[pa_p] = ix;
    // Set heap_ix(0)=1 to indicate that the values form a heap.
  }
  heap_ix[0] = 1;
  return;
}

void sortDecreasingHeap(const int n, vector<double>& heap_v, vector<int>& heap_ix) {
  int fo_p, srt_p;
  int cd_p, pa_p;
  int ix;
  double v;
  if (n <= 1) return;
  if (heap_ix[0] != 1) {
    // The data are assumed to be completely unordered. A heap will be formed and { sorted.
    fo_p = n/2 + 1;
    srt_p = n;
  } else {
    // The data are assumed to form a heap which is to be sorted.
    fo_p = 1;
    srt_p = n;
  }
  // 10   continue
  for (;;) {
    if (fo_p > 1) {
      fo_p = fo_p - 1;
      v = heap_v[fo_p];
      ix = heap_ix[fo_p];
    } else {
      v = heap_v[srt_p];
      ix = heap_ix[srt_p];
      heap_v[srt_p] = heap_v[1];
      heap_ix[srt_p] = heap_ix[1];
      srt_p--;
      if (srt_p == 1) {
	heap_v[1] = v;
	heap_ix[1] = ix;
	return;
      }
    }
    pa_p = fo_p;
    cd_p = fo_p + fo_p;
    for (;;) {
      if (cd_p <= srt_p) {
	if (cd_p < srt_p) {
	  if (heap_v[cd_p] > heap_v[cd_p+1]) cd_p = cd_p + 1;
	}
	if (v > heap_v[cd_p]) {
	  heap_v[pa_p] = heap_v[cd_p];
	  heap_ix[pa_p] = heap_ix[cd_p];
	  pa_p = cd_p;
	  cd_p = cd_p + cd_p;
	  continue;
	}
      }
      break;
    }
      heap_v[pa_p] = v;
      heap_ix[pa_p] = ix;
  }
    return;
}

void maxheapsort(int* heap_v, int n) {
  buildMaxheap(heap_v, n);
  maxHeapsort(heap_v, n);
}

void maxheapsort(int* heap_v, int* heap_i, int n) {
  buildMaxheap(heap_v, heap_i, n);
  maxHeapsort(heap_v, heap_i, n);
}

void maxheapsort(double* heap_v, int* heap_i, int n) {
  buildMaxheap(heap_v, heap_i, n);
  maxHeapsort(heap_v, heap_i, n);
}

void buildMaxheap(int* heap_v, int n) {
  int i;
  for (i = n / 2; i >= 1; i--) {
    maxHeapify(heap_v, i, n);
  }
}

void buildMaxheap(int* heap_v, int* heap_i, int n) {
  int i;
  for (i = n / 2; i >= 1; i--) {
    maxHeapify(heap_v, heap_i, i, n);
  }
}

void buildMaxheap(double* heap_v, int* heap_i, int n) {
  int i;
  for (i = n / 2; i >= 1; i--) {
    maxHeapify(heap_v, heap_i, i, n);
  }
}

void maxHeapsort(int* heap_v, int n) {
  int temp_v;
  int i;
  for (i = n; i >= 2; i--) {
    temp_v = heap_v[i];
    heap_v[i] = heap_v[1];
    heap_v[1] = temp_v;
    maxHeapify(heap_v, 1, i - 1);
  }
}

void maxHeapsort(int* heap_v, int* heap_i, int n) {
  int temp_v;
  int i, temp_i;
  for (i = n; i >= 2; i--) {
    temp_v = heap_v[i];
    heap_v[i] = heap_v[1];
    heap_v[1] = temp_v;
    temp_i = heap_i[i];
    heap_i[i] = heap_i[1];
    heap_i[1] = temp_i;
    maxHeapify(heap_v, heap_i, 1, i - 1);
  }
}

void maxHeapsort(double* heap_v, int* heap_i, int n) {
  double temp_v;
  int i, temp_i;
  for (i = n; i >= 2; i--) {
    temp_v = heap_v[i];
    heap_v[i] = heap_v[1];
    heap_v[1] = temp_v;
    temp_i = heap_i[i];
    heap_i[i] = heap_i[1];
    heap_i[1] = temp_i;
    maxHeapify(heap_v, heap_i, 1, i - 1);
  }
}

void maxHeapify(int* heap_v, int i, int n) {
  int temp_v;
  int j;
  temp_v = heap_v[i];
  j = 2 * i;
  while (j <= n) {
    if (j < n && heap_v[j + 1] > heap_v[j]) j = j + 1;
    if (temp_v > heap_v[j])
      break;
    else if (temp_v <= heap_v[j]) {
      heap_v[j / 2] = heap_v[j];
      j = 2 * j;
    }
  }
  heap_v[j / 2] = temp_v;
  return;
}

void maxHeapify(int* heap_v, int* heap_i, int i, int n) {
  int temp_v;
  int j, temp_i;
  temp_v = heap_v[i];
  temp_i = heap_i[i];
  j = 2 * i;
  while (j <= n) {
    if (j < n && heap_v[j + 1] > heap_v[j]) j = j + 1;
    if (temp_v > heap_v[j])
      break;
    else if (temp_v <= heap_v[j]) {
      heap_v[j / 2] = heap_v[j];
      heap_i[j / 2] = heap_i[j];
      j = 2 * j;
    }
  }
  heap_v[j / 2] = temp_v;
  heap_i[j / 2] = temp_i;
  return;
}

void maxHeapify(double* heap_v, int* heap_i, int i, int n) {
  double temp_v;
  int j, temp_i;
  temp_v = heap_v[i];
  temp_i = heap_i[i];
  j = 2 * i;
  while (j <= n) {
    if (j < n && heap_v[j + 1] > heap_v[j]) j = j + 1;
    if (temp_v > heap_v[j])
      break;
    else if (temp_v <= heap_v[j]) {
      heap_v[j / 2] = heap_v[j];
      heap_i[j / 2] = heap_i[j];
      j = 2 * j;
    }
  }
  heap_v[j / 2] = temp_v;
  heap_i[j / 2] = temp_i;
  return;
}

bool increasingSetOk(const int* set, const int set_num_entries,
                     const int set_entry_lower, const int set_entry_upper,
                     bool strict) {
  if (set_num_entries < 0) return false;
  if (set == NULL) return false;
  bool check_bounds = set_entry_lower <= set_entry_upper;
  int previous_entry;
  if (check_bounds) {
    if (strict) {
      previous_entry = set_entry_lower - 1;
    } else {
      previous_entry = set_entry_lower;
    }
  } else {
    previous_entry = -HIGHS_CONST_I_INF;
  }
  for (int k = 0; k < set_num_entries; k++) {
    int entry = set[k];
    if (strict) {
      if (entry <= previous_entry) return false;
    } else {
      if (entry < previous_entry) return false;
    }
    if (check_bounds && entry > set_entry_upper) return false;
    previous_entry = entry;
  }
  return true;
}

bool increasingSetOk(const double* set, const int set_num_entries,
                     const double set_entry_lower, const double set_entry_upper,
                     bool strict) {
  if (set_num_entries < 0) return false;
  if (set == NULL) return false;
  bool check_bounds = set_entry_lower <= set_entry_upper;
  double previous_entry;
  if (check_bounds) {
    if (strict) {
      if (set_entry_lower < 0) {
        previous_entry = (1 + HIGHS_CONST_TINY) * set_entry_lower;
      } else if (set_entry_lower > 0) {
        previous_entry = (1 - HIGHS_CONST_TINY) * set_entry_lower;
      } else {
        previous_entry = -HIGHS_CONST_TINY;
      }
    } else {
      previous_entry = set_entry_lower;
    }
  } else {
    previous_entry = -HIGHS_CONST_INF;
  }
  for (int k = 0; k < set_num_entries; k++) {
    double entry = set[k];
    if (strict) {
      if (entry <= previous_entry) return false;
    } else {
      if (entry < previous_entry) return false;
    }
    if (check_bounds && entry > set_entry_upper) return false;
    previous_entry = entry;
  }
  return true;
}

void sortSetData(const int num_entries, int* set, const double* data0,
                 const double* data1, const double* data2, double* sorted_data0,
                 double* sorted_data1, double* sorted_data2) {
  vector<int> sort_set_vec(1 + num_entries);
  vector<int> perm_vec(1 + num_entries);

  int* sort_set = &sort_set_vec[0];
  int* perm = &perm_vec[0];

  for (int ix = 0; ix < num_entries; ix++) {
    sort_set[1 + ix] = set[ix];
    perm[1 + ix] = ix;
  }
  maxheapsort(sort_set, perm, num_entries);
  for (int ix = 0; ix < num_entries; ix++) {
    set[ix] = sort_set[1 + ix];
    if (data0 != NULL) sorted_data0[ix] = data0[perm[1 + ix]];
    if (data1 != NULL) sorted_data1[ix] = data1[perm[1 + ix]];
    if (data2 != NULL) sorted_data2[ix] = data2[perm[1 + ix]];
  }
}
