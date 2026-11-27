// Property tests for the sparse graph routines used by IPX and BASICLU:
//   ipx::DepthFirstSearch, ipx::AugmentingPath and lu_dfs.
//
// Each routine is checked against its documented contract on random graphs,
// using small independent reference computations (breadth-first reachability
// and a simple maximum bipartite matching), not against another
// implementation of the same routine.

#include <algorithm>
#include <functional>
#include <random>
#include <vector>

#include "catch.hpp"
#include "ipm/basiclu/basiclu.h"
#include "ipm/ipx/sparse_utils.h"

// lu_internal.h cannot be included from C++, so declare lu_dfs here.
extern "C" lu_int lu_dfs(lu_int i, const lu_int* begin, const lu_int* end,
                         const lu_int* index, lu_int top, lu_int* xi,
                         lu_int* pstack, lu_int* marked, const lu_int M);

namespace {

using Adjacency = std::vector<std::vector<HighsInt>>;

// Random adjacency lists for num_lists lists over targets 0..num_targets-1.
// If acyclic, list j only contains targets > j, so that the graph on
// num_targets == num_lists nodes has no cycles.
Adjacency randomAdjacency(HighsInt num_lists, HighsInt num_targets,
                          double density, bool acyclic, std::mt19937& rng) {
  std::uniform_real_distribution<double> coin(0.0, 1.0);
  Adjacency adj(num_lists);
  for (HighsInt j = 0; j < num_lists; j++) {
    for (HighsInt i = acyclic ? j + 1 : 0; i < num_targets; i++)
      if (coin(rng) < density) adj[j].push_back(i);
    // Neighbour order matters for the search order; do not keep it sorted.
    std::shuffle(adj[j].begin(), adj[j].end(), rng);
  }
  return adj;
}

// Column pointers and indices (CSC) of an adjacency structure.
template <typename I>
void toCsc(const Adjacency& adj, std::vector<I>& ptr, std::vector<I>& idx) {
  ptr.assign(1, 0);
  idx.clear();
  for (const auto& list : adj) {
    for (HighsInt i : list) idx.push_back(i);
    ptr.push_back(idx.size());
  }
}

// Nodes reachable from start without passing through a blocked node.
std::vector<HighsBool> reachableFrom(const Adjacency& adj, HighsInt start,
                                const std::vector<HighsBool>& blocked) {
  std::vector<HighsBool> seen(adj.size(), false);
  std::vector<HighsInt> queue{start};
  seen[start] = true;
  for (size_t q = 0; q < queue.size(); q++) {
    for (HighsInt next : adj[queue[q]]) {
      if (!seen[next] && !blocked[next]) {
        seen[next] = true;
        queue.push_back(next);
      }
    }
  }
  return seen;
}

// Checks one call of a reach routine: the output segment out[newtop..top-1]
// is exactly reach(start) in the graph minus blocked nodes, starts with
// start, is topologically ordered if the graph is acyclic, and the marks
// were set for exactly those nodes.
template <typename I>
void checkReach(const Adjacency& adj, HighsInt start,
                const std::vector<HighsBool>& blocked_before, bool acyclic,
                const std::vector<I>& out, I newtop, I top,
                const std::vector<I>& marked, I marker) {
  const HighsInt n = adj.size();
  const std::vector<HighsBool> expected = reachableFrom(adj, start, blocked_before);

  REQUIRE(newtop <= top);
  REQUIRE(newtop >= 0);
  REQUIRE(out[newtop] == start);

  std::vector<HighsInt> position(n, -1);
  for (I k = newtop; k < top; k++) {
    const HighsInt node = out[k];
    REQUIRE(node >= 0);
    REQUIRE(node < n);
    REQUIRE(position[node] == -1);  // no duplicates
    position[node] = k;
  }
  for (HighsInt node = 0; node < n; node++) {
    const bool in_output = position[node] >= 0;
    REQUIRE(in_output == (expected[node] && !blocked_before[node]));
    REQUIRE((marked[node] == marker) == (blocked_before[node] || in_output));
  }
  if (acyclic) {
    for (HighsInt u = 0; u < n; u++) {
      if (position[u] < 0) continue;
      for (HighsInt v : adj[u])
        if (position[v] >= 0) REQUIRE(position[u] < position[v]);
    }
  }
}

// Size of a maximum matching between columns and non-excluded rows, by
// repeated simple augmentation (Kuhn's method).
HighsInt maximumMatchingSize(const Adjacency& cols, HighsInt num_rows,
                             const std::vector<HighsBool>& excluded) {
  std::vector<HighsInt> row_to_col(num_rows, -1);
  std::vector<HighsInt> visit_stamp(num_rows, -1);
  // Recursive augmentation; depth is bounded by the small test sizes.
  std::function<bool(HighsInt, HighsInt)> augment = [&](HighsInt col,
                                                        HighsInt stamp) {
    for (HighsInt row : cols[col]) {
      if (excluded[row] || visit_stamp[row] == stamp) continue;
      visit_stamp[row] = stamp;
      if (row_to_col[row] < 0 || augment(row_to_col[row], stamp)) {
        row_to_col[row] = col;
        return true;
      }
    }
    return false;
  };
  HighsInt size = 0;
  for (HighsInt col = 0; col < (HighsInt)cols.size(); col++)
    if (augment(col, col)) size++;
  return size;
}

}  // namespace

TEST_CASE("ipx-depth-first-search", "[highs_sparse_graph]") {
  using ipx::Int;
  std::mt19937 rng(20261126);
  for (int trial = 0; trial < 400; trial++) {
    const Int m = 1 + rng() % 50;
    const bool acyclic = trial % 2 == 0;
    const bool use_colmap = trial % 3 == 0;
    const double density = (trial % 4 == 0) ? 0.02 : 0.15;

    // With a colmap, node i uses column colmap[i] (or none if negative);
    // otherwise node i uses column i.
    Adjacency columns = randomAdjacency(m, m, density, acyclic, rng);
    std::vector<Int> colmap(m);
    Adjacency adj(m);
    for (Int i = 0; i < m; i++) {
      if (use_colmap) {
        // Keep the graph acyclic by mapping i to itself or to no column.
        colmap[i] = (rng() % 4 == 0) ? -1 : (acyclic ? i : (Int)(rng() % m));
        if (colmap[i] >= 0) adj[i] = columns[colmap[i]];
      } else {
        adj[i] = columns[i];
      }
    }
    std::vector<Int> Ap, Ai;
    toCsc(columns, Ap, Ai);

    // Some nodes are marked before the first search; they must be skipped.
    const Int marker = 3;
    std::vector<Int> marked(m, -1);
    for (Int i = 0; i < m; i++)
      if (rng() % 5 == 0) marked[i] = marker;

    std::vector<Int> istack(m, -7), work(m, -7);
    Int top = m;
    std::vector<Int> order(m);
    for (Int i = 0; i < m; i++) order[i] = i;
    std::shuffle(order.begin(), order.end(), rng);

    for (Int start : order) {
      if (marked[start] == marker) continue;
      std::vector<HighsBool> blocked(m);
      for (Int i = 0; i < m; i++) blocked[i] = marked[i] == marker;
      const std::vector<Int> output_before(istack.begin() + top, istack.end());

      const Int newtop = ipx::DepthFirstSearch(
          start, Ap.data(), Ai.data(), use_colmap ? colmap.data() : nullptr,
          top, istack.data(), marked.data(), marker, work.data());

      checkReach<Int>(adj, start, blocked, acyclic, istack, newtop, top,
                      marked, marker);
      // Output from earlier calls is left in place.
      REQUIRE(std::equal(output_before.begin(), output_before.end(),
                         istack.begin() + top));
      top = newtop;
    }
  }
}

TEST_CASE("ipx-augmenting-path", "[highs_sparse_graph]") {
  using ipx::Int;
  std::mt19937 rng(17);
  for (int trial = 0; trial < 400; trial++) {
    const Int m = 1 + rng() % 40;
    const Int n = 1 + rng() % 50;
    const double density = (trial % 3 == 0) ? 0.03 : 0.12;
    const Adjacency cols = randomAdjacency(n, m, density, false, rng);
    std::vector<Int> Ap, Ai;
    toCsc(cols, Ap, Ai);

    // Rows with jmatch < -1 are excluded from matching.
    std::vector<Int> jmatch(m, -1);
    std::vector<HighsBool> excluded(m, false);
    for (Int i = 0; i < m; i++) {
      if (rng() % 6 == 0) {
        jmatch[i] = -2;
        excluded[i] = true;
      }
    }
    std::vector<Int> cheap(Ap.begin(), Ap.end() - 1);
    std::vector<Int> marked(n, -1);
    std::vector<Int> work(m), work2(m + 1), work3(m + 1);

    std::vector<Int> col_order(n);
    for (Int j = 0; j < n; j++) col_order[j] = j;
    std::shuffle(col_order.begin(), col_order.end(), rng);

    Int num_matched = 0;
    for (Int jstart : col_order) {
      const std::vector<Int> jmatch_before = jmatch;
      const bool matched =
          ipx::AugmentingPath(jstart, Ap.data(), Ai.data(), jmatch.data(),
                              cheap.data(), marked.data(), work.data(),
                              work2.data(), work3.data());

      // The matching is valid: each matched row lies in its column, no
      // column is matched twice, and excluded rows are untouched.
      std::vector<Int> rows_of_col(n, 0);
      for (Int i = 0; i < m; i++) {
        if (excluded[i]) {
          REQUIRE(jmatch[i] == -2);
          continue;
        }
        REQUIRE(jmatch[i] >= -1);
        if (jmatch[i] >= 0) {
          const Int j = jmatch[i];
          REQUIRE(std::find(cols[j].begin(), cols[j].end(), i) !=
                  cols[j].end());
          rows_of_col[j]++;
        }
      }
      for (Int j = 0; j < n; j++) REQUIRE(rows_of_col[j] <= 1);

      // Previously matched columns stay matched; jstart is matched iff the
      // call returned true; on failure nothing changes.
      for (Int i = 0; i < m; i++) {
        if (jmatch_before[i] >= 0) {
          const Int j = jmatch_before[i];
          REQUIRE(std::count(jmatch.begin(), jmatch.end(), j) == 1);
        }
      }
      const bool jstart_matched =
          std::count(jmatch.begin(), jmatch.end(), jstart) == 1;
      REQUIRE(matched == jstart_matched);
      if (!matched) REQUIRE(jmatch == jmatch_before);
      num_matched += matched;

      // cheap[j] only passes over rows that are no longer free.
      for (Int j = 0; j < n; j++) {
        REQUIRE(cheap[j] >= Ap[j]);
        REQUIRE(cheap[j] <= Ap[j + 1]);
        for (Int p = Ap[j]; p < cheap[j]; p++) REQUIRE(jmatch[Ai[p]] != -1);
      }
    }

    // Trying every column once gives a maximum matching.
    REQUIRE(num_matched == maximumMatchingSize(cols, m, excluded));
  }
}

TEST_CASE("basiclu-lu-dfs", "[highs_sparse_graph]") {
  std::mt19937 rng(4242);
  for (int trial = 0; trial < 400; trial++) {
    const lu_int m = 1 + rng() % 50;
    const bool acyclic = trial % 2 == 0;
    // lu_dfs accepts either explicit list ends or lists terminated by a
    // negative index (end == NULL).
    const bool use_end = trial % 4 < 2;
    const double density = (trial % 3 == 0) ? 0.02 : 0.15;
    const Adjacency adj = randomAdjacency(m, m, density, acyclic, rng);

    // Store the lists in a shuffled order in index[] so begin[] is not
    // monotone, as happens in BASICLU's row and column files.
    std::vector<lu_int> list_order(m);
    for (lu_int j = 0; j < m; j++) list_order[j] = j;
    std::shuffle(list_order.begin(), list_order.end(), rng);
    std::vector<lu_int> begin(m), end(m), index;
    for (lu_int j : list_order) {
      begin[j] = index.size();
      for (HighsInt i : adj[j]) index.push_back(i);
      end[j] = index.size();
      if (!use_end) index.push_back(-1);
    }
    if (index.empty()) index.push_back(-1);

    const lu_int M = 5;
    std::vector<lu_int> marked(m, -1);
    for (lu_int i = 0; i < m; i++)
      if (rng() % 5 == 0) marked[i] = M;

    std::vector<lu_int> xi(m, -7), pstack(m);
    lu_int top = m;
    for (lu_int start = 0; start < m; start++) {
      std::vector<HighsBool> blocked(m);
      for (lu_int i = 0; i < m; i++) blocked[i] = marked[i] == M;
      const std::vector<lu_int> marked_before = marked;
      const std::vector<lu_int> output_before(xi.begin() + top, xi.end());

      const lu_int newtop =
          lu_dfs(start, begin.data(), use_end ? end.data() : nullptr,
                 index.data(), top, xi.data(), pstack.data(), marked.data(), M);

      if (blocked[start]) {
        // A marked start node is a no-op.
        REQUIRE(newtop == top);
        REQUIRE(marked == marked_before);
      } else {
        checkReach<lu_int>(adj, start, blocked, acyclic, xi, newtop, top,
                           marked, M);
      }
      REQUIRE(std::equal(output_before.begin(), output_before.end(),
                         xi.begin() + top));
      top = newtop;
    }
  }
}
