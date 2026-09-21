#include "Fold.h"

#include <algorithm>
#include <stack>

namespace hipo {

std::vector<Int> colourRefinement(const std::vector<Int> ptr,
                                  const std::vector<Int> adj) {
  // Taken from "Tight Lower and Upper Bounds for the Complexity
  // of Canonical Colour Refinement", Berkholz, Bonsma, Grohe

  // initial uniform colour 0

  const Int n = ptr.size() - 1;

  std::vector<Int> maxcdeg(n, 0), mincdeg(n, 0), cdeg(n, 0), colour(n, 0);

  std::vector<Int> colours_adj(n, 0);
  std::vector<HighsBool> colours_adj_belong(n, 0);
  Int colours_adj_top = 0;

  std::vector<Int> colours_split(n, 0);
  Int colours_split_top = 0;

  CollectionLinkedLists C, A;
  C.init(n, n);
  A.init(n, n);
  for (Int i = 0; i < n; ++i) C.append(i, 0);

  Int k = 0;

  std::stack<Int> S_refine;
  std::vector<HighsBool> in_stack(n, 0);
  S_refine.push(0);
  in_stack[0] = 1;

  auto split_up_colour = [&](Int s) {
    const Int maxcdegs = maxcdeg[s];
    std::vector<Int> numcdeg(maxcdegs + 1, 0);
    numcdeg[0] = C.length(s) - A.length(s);

    Int v = A.head(s);
    while (v < n) {
      numcdeg[cdeg[v]]++;
      v = A.next(v);
    }

    Int b = 0;
    for (Int i = 1; i <= maxcdegs; ++i) {
      if (numcdeg[i] > numcdeg[b]) b = i;
    }

    std::vector<Int> f(maxcdegs + 1, 0);
    for (Int i = 0; i <= maxcdegs; ++i) {
      if (numcdeg[i] >= 1) {
        if (i == mincdeg[s]) {
          f[i] = s;
          if (!in_stack[s] && b != i) {
            S_refine.push(f[i]);
            in_stack[f[i]] = 1;
          }

        } else {
          k++;
          f[i] = k;
          if (in_stack[s] || i != b) {
            S_refine.push(f[i]);
            in_stack[f[i]] = 1;
          }
        }
      }
    }

    v = A.head(s);
    while (v < n) {
      if (f[cdeg[v]] != s) {
        C.remove(v, s);
        C.append(v, f[cdeg[v]]);
        colour[v] = f[cdeg[v]];
      }
      v = A.next(v);
    }
  };

  while (!S_refine.empty()) {
    printf("\n");
    for (Int c : colour) printf("%d", c);
    printf("\n");

    const Int r = S_refine.top();
    S_refine.pop();
    in_stack[r] = 0;

    printf("Refine with r = %d\n", r);

    Int v = C.head(r);
    while (v < n) {
      for (Int el = ptr[v]; el < ptr[v + 1]; ++el) {
        const Int w = adj[el];
        cdeg[w]++;
        if (cdeg[w] == 1) A.append(w, colour[w]);

        if (!colours_adj_belong[colour[w]]) {
          colours_adj[colours_adj_top] = colour[w];
          colours_adj_top++;
          colours_adj_belong[colour[w]] = 1;
        }

        if (cdeg[w] > maxcdeg[colour[w]]) maxcdeg[colour[w]] = cdeg[w];
      }
      v = C.next(v);
    }

    for (Int c_ind = 0; c_ind < colours_adj_top; ++c_ind) {
      const Int c = colours_adj[c_ind];
      if (C.length(c) != A.length(c))
        mincdeg[c] = 0;
      else {
        mincdeg[c] = maxcdeg[c];
        Int v = A.head(c);
        while (v < n) {
          if (cdeg[v] < mincdeg[c]) mincdeg[c] = cdeg[v];
          v = A.next(v);
        }
      }
    }

    colours_split_top = 0;
    for (Int c_ind = 0; c_ind < colours_adj_top; ++c_ind) {
      const Int c = colours_adj[c_ind];
      if (mincdeg[c] < maxcdeg[c]) {
        colours_split[colours_split_top] = c;
        colours_split_top++;
      }
    }
    std::sort(colours_split.begin(), colours_split.begin() + colours_split_top);
    for (Int s_ind = 0; s_ind < colours_split_top; ++s_ind) {
      const Int s = colours_split[s_ind];
      printf("\tSplit %d\n", s);
      split_up_colour(s);
    }

    for (Int c_ind = 0; c_ind < colours_adj_top; ++c_ind) {
      const Int c = colours_adj[c_ind];
      Int v = A.head(c);
      while (v < n) {
        cdeg[v] = 0;
        v = A.next(v);
      }
      maxcdeg[c] = 0;
      A.clear(c);
      colours_adj_belong[c] = 0;
    }
    colours_adj_top = 0;
  }

  return colour;
}

void test_folding() {
  const std::vector<Int> ptr = {0, 3, 5, 8, 10, 14, 16, 19, 21, 24};
  const std::vector<Int> adj = {1, 2, 3, 0, 4, 0, 3, 4, 0, 2, 1, 2,
                                5, 8, 4, 6, 5, 7, 8, 6, 8, 4, 6, 7};

  const std::vector<Int> colour = colourRefinement(ptr, adj);

  printf("\n\n");
  for (Int c : colour) printf("%d", c);
  printf("\n");
}

}  // namespace hipo