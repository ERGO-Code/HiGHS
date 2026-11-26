#include "ipm/ipx/sparse_utils.h"
#include <cassert>

namespace ipx {

// Iterative depth-first search. istack[0..depth] holds the nodes on the path
// from istart to the current node, and work[d] is the position in the
// neighbour list of istack[d] where the scan resumes. A node is written to
// the output (istack below top, growing downward) once all its neighbours
// have been visited, so the output is in reverse postorder.
//
// The path and the output share istack. They do not overlap, because each
// node in either part was unmarked on entry and the caller guarantees that
// istack[0..top-1] has room for all such nodes.
Int DepthFirstSearch(Int istart, const Int* Ap, const Int* Ai,
                     const Int* colmap, Int top, Int* istack, Int* marked,
                     Int marker, Int* work) {
    assert(marked[istart] != marker);
    auto neighbours_begin = [&](Int i) -> Int {
        const Int j = colmap ? colmap[i] : i;
        return j >= 0 ? Ap[j] : 0;
    };
    auto neighbours_end = [&](Int i) -> Int {
        const Int j = colmap ? colmap[i] : i;
        return j >= 0 ? Ap[j+1] : 0;
    };

    Int depth = 0;
    istack[0] = istart;
    marked[istart] = marker;
    work[0] = neighbours_begin(istart);
    while (depth >= 0) {
        const Int i = istack[depth];
        const Int pend = neighbours_end(i);
        Int p = work[depth];
        while (p < pend && marked[Ai[p]] == marker)
            p++;
        if (p < pend) {
            // descend to the first unmarked neighbour
            const Int inext = Ai[p];
            work[depth] = p+1;
            depth++;
            istack[depth] = inext;
            marked[inext] = marker;
            work[depth] = neighbours_begin(inext);
        } else {
            // all neighbours of i visited
            depth--;
            istack[--top] = i;
        }
    }
    return top;
}

// Depth-first search for an alternating path that starts at column jstart,
// alternates between unmatched and matched entries, and ends at a free row.
//
// Each level d of the path uses two entries of @work: work[2*d] is the column
// on the path, and work[2*d+1] is the position in that column where the scan
// resumes. When the path leaves the column at level d, it does so through row
// Ai[work[2*d+1]-1], so the rows on the path do not need their own stack.
// The path has at most m+1 columns, because each column after jstart is
// matched to a different row; hence @work needs 2*(m+1) entries.
//
// When a column is entered for the first time, its entries from cheap[j]
// onward are checked for a free row ("cheap assignment"). Matched rows never
// become free, so cheap[j] only moves forward and the total cost of these
// checks over all calls is O(nnz(A)).
bool AugmentingPath(Int jstart, const Int* Ap, const Int* Ai, Int* jmatch,
                    Int* cheap, Int* marked, Int* work) {
    Int depth = 0;
    Int j = jstart;
    bool entered = true;

    for (;;) {
        if (entered) {
            entered = false;
            work[2*depth] = j;
            marked[j] = jstart;
            while (cheap[j] < Ap[j+1]) {
                const Int i = Ai[cheap[j]++];
                if (jmatch[i] == -1) {
                    // flip the path: match i to the last column and the
                    // row leaving each earlier column to that column
                    jmatch[i] = j;
                    for (Int d = depth-1; d >= 0; d--)
                        jmatch[Ai[work[2*d+1]-1]] = work[2*d];
                    return true;
                }
            }
            work[2*depth+1] = Ap[j];
        }

        // look for a matched row in column j whose column is not on the path
        // and has not been visited in this search
        const Int pend = Ap[j+1];
        Int p = work[2*depth+1];
        Int jnext = -1;
        for (; p < pend; p++) {
            const Int jm = jmatch[Ai[p]];
            assert(jm != -1);  // no free rows remain in column j
            if (jm >= 0 && marked[jm] != jstart) {
                jnext = jm;
                break;
            }
        }
        if (jnext >= 0) {
            work[2*depth+1] = p+1;
            depth++;
            j = jnext;
            entered = true;
        } else {
            // no alternating path through column j
            if (depth == 0)
                return false;
            depth--;
            j = work[2*depth];
        }
    }
}

}  // namespace ipx
