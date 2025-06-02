#include "Auxiliary.h"

#include <stack>

namespace highspm {

void counts2Ptr(std::vector<Int>& ptr, std::vector<Int>& w) {
  // Given the column counts in the vector w (of size n),
  // compute the column pointers in the vector ptr (of size n+1),
  // and copy the first n pointers back into w.

  Int temp_nz{};
  Int n = w.size();
  for (Int j = 0; j < n; ++j) {
    ptr[j] = temp_nz;
    temp_nz += w[j];
    w[j] = ptr[j];
  }
  ptr[n] = temp_nz;
}

void inversePerm(const std::vector<Int>& perm, std::vector<Int>& iperm) {
  // Given the permutation perm, produce the inverse permutation iperm.
  // perm[i] : i-th entry to use in the new order.
  // iperm[i]: where entry i is located in the new order.

  for (Int i = 0; i < perm.size(); ++i) {
    iperm[perm[i]] = i;
  }
}

void subtreeSize(const std::vector<Int>& parent, std::vector<Int>& sizes) {
  // Compute sizes of subtrees of the tree given by parent

  Int n = parent.size();
  sizes.assign(n, 1);

  for (Int i = 0; i < n; ++i) {
    Int k = parent[i];
    if (k != -1) sizes[k] += sizes[i];
  }
}

void transpose(const std::vector<Int>& ptr, const std::vector<Int>& rows,
               std::vector<Int>& ptrT, std::vector<Int>& rowsT) {
  // Compute the transpose of the matrix and return it in rowsT and ptrT

  Int n = ptr.size() - 1;

  std::vector<Int> work(n);

  // count the entries in each row into work
  for (Int i = 0; i < ptr.back(); ++i) {
    ++work[rows[i]];
  }

  // sum row sums to obtain pointers
  counts2Ptr(ptrT, work);

  for (Int j = 0; j < n; ++j) {
    for (Int el = ptr[j]; el < ptr[j + 1]; ++el) {
      Int i = rows[el];

      // entry (i,j) becomes entry (j,i)
      Int pos = work[i]++;
      rowsT[pos] = j;
    }
  }
}

void transpose(const std::vector<Int>& ptr, const std::vector<Int>& rows,
               const std::vector<double>& val, std::vector<Int>& ptrT,
               std::vector<Int>& rowsT, std::vector<double>& valT) {
  // Compute the transpose of the matrix and return it in rowsT, ptrT and valT

  Int n = ptr.size() - 1;

  std::vector<Int> work(n);

  // count the entries in each row into work
  for (Int i = 0; i < ptr.back(); ++i) {
    ++work[rows[i]];
  }

  // sum row sums to obtain pointers
  counts2Ptr(ptrT, work);

  for (Int j = 0; j < n; ++j) {
    for (Int el = ptr[j]; el < ptr[j + 1]; ++el) {
      Int i = rows[el];

      // entry (i,j) becomes entry (j,i)
      Int pos = work[i]++;
      rowsT[pos] = j;
      valT[pos] = val[el];
    }
  }
}

void symProduct(const std::vector<Int>& ptr, const std::vector<Int>& rows,
                const std::vector<double>& vals, const std::vector<double>& x,
                std::vector<double>& y, double alpha) {
  // Matrix-vector product in CSC format, for symmetric matrix which stores only
  // the lower triangle.
  // Compute y = y + alpha * M * x

  const Int n = ptr.size() - 1;

  for (Int col = 0; col < n; ++col) {
    for (Int el = ptr[col]; el < ptr[col + 1]; ++el) {
      Int row = rows[el];
      double val = vals[el];

      y[row] += alpha * val * x[col];
      if (row != col) y[col] += alpha * val * x[row];
    }
  }
}

void symProductQuad(const std::vector<Int>& ptr, const std::vector<Int>& rows,
                    const std::vector<double>& vals,
                    const std::vector<double>& x, std::vector<HighsCDouble>& y,
                    double alpha) {
  // Matrix-vector product in CSC format, for symmetric matrix which stores only
  // the lower triangle.
  // Compute y = y + alpha * M * x

  const Int n = ptr.size() - 1;

  for (Int col = 0; col < n; ++col) {
    for (Int el = ptr[col]; el < ptr[col + 1]; ++el) {
      Int row = rows[el];
      HighsCDouble val = vals[el];

      y[row] += val * (HighsCDouble)x[col] * (HighsCDouble)alpha;
      if (row != col)
        y[col] += val * (HighsCDouble)x[row] * (HighsCDouble)alpha;
    }
  }
}

void childrenLinkedList(const std::vector<Int>& parent, std::vector<Int>& head,
                        std::vector<Int>& next) {
  // Create linked lists of children in elimination tree.
  // parent gives the dependencies of the tree,
  // head[node] is the first child of node,
  // next[head[node]] is the second child,
  // next[next[head[node]]] is the third child...
  // until -1 is reached.

  Int n = parent.size();
  head.assign(n, -1);
  next.assign(n, -1);
  for (Int node = n - 1; node >= 0; --node) {
    if (parent[node] == -1) continue;
    next[node] = head[parent[node]];
    head[parent[node]] = node;
  }
}

void reverseLinkedList(std::vector<Int>& head, std::vector<Int>& next) {
  // Reverse the linked list of children of each node.
  // If a node has children (a -> b -> c -> -1), the reverse list contains
  // children (c -> b -> a -> -1).

  const Int n = head.size();

  for (Int node = 0; node < n; ++node) {
    Int prev_node = -1;
    Int curr_node = head[node];
    Int next_node = -1;

    while (curr_node != -1) {
      next_node = next[curr_node];
      next[curr_node] = prev_node;
      prev_node = curr_node;
      curr_node = next_node;
    }

    head[node] = prev_node;
  }
}

void dfsPostorder(Int node, Int& start, std::vector<Int>& head,
                  const std::vector<Int>& next, std::vector<Int>& order) {
  // Perform depth first search starting from root node and order the nodes
  // starting from the value start. head and next contain the linked list of
  // children.

  std::stack<Int> stack;
  stack.push(node);

  while (!stack.empty()) {
    const Int current = stack.top();
    const Int child = head[current];

    if (child == -1) {
      // no children left to order,
      // remove from the stack and order
      stack.pop();
      order[start++] = current;
    } else {
      // at least one child left to order,
      // add it to the stack and remove it from the list of children
      stack.push(child);
      head[current] = next[child];
    }
  }
}

void processEdge(Int j, Int i, const std::vector<Int>& first,
                 std::vector<Int>& maxfirst, std::vector<Int>& delta,
                 std::vector<Int>& prevleaf, std::vector<Int>& ancestor) {
  // Process edge of skeleton matrix.
  // Taken from Tim Davis "Direct Methods for Sparse Linear Systems".

  // j not a leaf of ith row subtree
  if (i <= j || first[j] <= maxfirst[i]) {
    return;
  }

  // max first[j] so far
  maxfirst[i] = first[j];

  // previous leaf of ith row subtree
  Int jprev = prevleaf[i];

  // A(i,j) is in the skeleton matrix
  delta[j]++;

  if (jprev != -1) {
    // find least common ancestor of jprev and j
    Int q = jprev;
    while (q != ancestor[q]) {
      q = ancestor[q];
    }

    // path compression
    Int sparent;
    for (Int s = jprev; s != q; s = sparent) {
      sparent = ancestor[s];
      ancestor[s] = q;
    }

    // consider overlap
    delta[q]--;
  }

  // previous leaf of ith subtree set to j
  prevleaf[i] = j;
}

double getDiagStart(Int n, Int k, Int nb, Int n_blocks, std::vector<Int>& start,
                    bool triang) {
  // start position of diagonal blocks for blocked dense formats
  start.assign(n_blocks, 0);
  for (Int i = 1; i < n_blocks; ++i) {
    start[i] = start[i - 1] + nb * (n - (i - 1) * nb);
    if (triang) start[i] -= nb * (nb - 1) / 2;
  }

  Int jb = std::min(nb, k - (n_blocks - 1) * nb);
  double result = (double)start.back() + (double)(n - (n_blocks - 1) * nb) * jb;
  if (triang) result -= (double)jb * (jb - 1) / 2;
  return result;
}

Clock::Clock() { start(); }
void Clock::start() { t0 = std::chrono::high_resolution_clock::now(); }
double Clock::stop() const {
  auto t1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> d = t1 - t0;
  return d.count();
}

}  // namespace highspm
