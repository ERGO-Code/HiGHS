#include "Fold.h"

#include <algorithm>
#include <stack>

namespace hipo {

// pick a colour
// compute colour degrees
// compute colours that are split
// split colours

std::vector<Int> colourRefinement(const std::vector<Int> ptr,
                                  const std::vector<Int> adj) {
  // Taken from "Tight Lower and Upper Bounds for the Complexity
  // of Canonical Colour Refinement", Berkholz, Bonsma, Grohe

  // initial uniform colour 0

  const Int n = ptr.size() - 1;

  std::vector<Int> colour(n, 0);
  std::vector<Int> colour_degree(n, 0);
  std::vector<Int> max_colour_degree(n, 0);
  std::vector<Int> min_colour_degree(n, 0);

  std::vector<Int> colours_touched(n, 0);
  std::vector<HighsBool> in_colours_touched(n, 0);
  Int top_touched = 0;

  std::vector<Int> colours_split(n, 0);
  Int top_split = 0;

  CollectionLinkedLists colour_classes;
  colour_classes.init(n, n);
  for (Int i = 0; i < n; ++i) colour_classes.append(i, 0);

  CollectionLinkedLists colour_classes_touched;
  colour_classes_touched.init(n, n);

  Int latest_colour = 0;

  std::stack<Int> stack_refine;
  std::vector<HighsBool> in_stack(n, 0);
  stack_refine.push(0);
  in_stack[0] = 1;

  auto split_up_colour = [&](Int s) {
    const Int maxcdegs = max_colour_degree[s];
    std::vector<Int> numcdeg(maxcdegs + 1, 0);
    numcdeg[0] = colour_classes.length(s) - colour_classes_touched.length(s);

    Int v = colour_classes_touched.head(s);
    while (colour_classes_touched.cont(v)) {
      numcdeg[colour_degree[v]]++;
      v = colour_classes_touched.next(v);
    }

    Int b = 0;
    for (Int i = 1; i <= maxcdegs; ++i) {
      if (numcdeg[i] > numcdeg[b]) b = i;
    }

    std::vector<Int> f(maxcdegs + 1, 0);
    for (Int i = 0; i <= maxcdegs; ++i) {
      if (numcdeg[i] >= 1) {
        if (i == min_colour_degree[s]) {
          f[i] = s;
          if (!in_stack[s] && b != i) {
            stack_refine.push(f[i]);
            in_stack[f[i]] = 1;
          }

        } else {
          latest_colour++;
          f[i] = latest_colour;
          if (in_stack[s] || i != b) {
            stack_refine.push(f[i]);
            in_stack[f[i]] = 1;
          }
        }
      }
    }

    v = colour_classes_touched.head(s);
    while (colour_classes_touched.cont(v)) {
      if (f[colour_degree[v]] != s) {
        colour_classes.remove(v, s);
        colour_classes.append(v, f[colour_degree[v]]);
        colour[v] = f[colour_degree[v]];
      }
      v = colour_classes_touched.next(v);
    }
  };

  while (!stack_refine.empty()) {
    printf("\n");
    for (Int c : colour) printf("%d", c);
    printf("\n");

    const Int refining_colour = stack_refine.top();
    stack_refine.pop();
    in_stack[refining_colour] = 0;

    printf("Refine with r = %d\n", refining_colour);

    Int v = colour_classes.head(refining_colour);
    while (colour_classes.cont(v)) {
      for (Int el = ptr[v]; el < ptr[v + 1]; ++el) {
        const Int w = adj[el];
        colour_degree[w]++;
        if (colour_degree[w] == 1) colour_classes_touched.append(w, colour[w]);

        if (!in_colours_touched[colour[w]]) {
          colours_touched[top_touched] = colour[w];
          top_touched++;
          in_colours_touched[colour[w]] = 1;
        }

        if (colour_degree[w] > max_colour_degree[colour[w]])
          max_colour_degree[colour[w]] = colour_degree[w];
      }
      v = colour_classes.next(v);
    }

    for (Int el = 0; el < top_touched; ++el) {
      const Int c = colours_touched[el];
      if (colour_classes.length(c) != colour_classes_touched.length(c))
        min_colour_degree[c] = 0;
      else {
        min_colour_degree[c] = max_colour_degree[c];
        Int v = colour_classes_touched.head(c);
        while (colour_classes_touched.cont(v)) {
          if (colour_degree[v] < min_colour_degree[c])
            min_colour_degree[c] = colour_degree[v];
          v = colour_classes_touched.next(v);
        }
      }
    }

    top_split = 0;
    for (Int el = 0; el < top_touched; ++el) {
      const Int c = colours_touched[el];
      if (min_colour_degree[c] < max_colour_degree[c]) {
        colours_split[top_split] = c;
        top_split++;
      }
    }
    std::sort(colours_split.begin(), colours_split.begin() + top_split);
    for (Int el = 0; el < top_split; ++el) {
      const Int s = colours_split[el];
      printf("\tSplit %d\n", s);
      split_up_colour(s);
    }

    for (Int el = 0; el < top_touched; ++el) {
      const Int c = colours_touched[el];
      Int v = colour_classes_touched.head(c);
      while (colour_classes_touched.cont(v)) {
        colour_degree[v] = 0;
        v = colour_classes_touched.next(v);
      }
      max_colour_degree[c] = 0;
      colour_classes_touched.clear(c);
      in_colours_touched[c] = 0;
    }
    top_touched = 0;
  }

  return colour;
}

void test_folding() {
  const std::vector<Int> ptr = {0,  3,  6,  9,  19, 22, 25,
                                28, 31, 34, 37, 40, 41, 42};
  const std::vector<Int> adj = {1, 2, 3, 0, 2,  3, 0, 1, 3, 0, 1,  2,  4,  5,
                                6, 7, 8, 9, 10, 3, 7, 9, 3, 6, 8,  3,  5,  8,
                                3, 4, 9, 3, 5,  6, 3, 4, 7, 3, 11, 12, 10, 10};

  const std::vector<Int> colour = colourRefinement(ptr, adj);

  printf("\n\n");
  for (Int c : colour) printf("%d", c);
  printf("\n");
}

}  // namespace hipo