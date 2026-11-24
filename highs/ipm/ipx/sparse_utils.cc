#include "ipm/ipx/sparse_utils.h"
#include <cassert>

namespace ipx {

Int DepthFirstSearch(Int istart, const Int* Ap, const Int* Ai,
                     const Int* colmap, Int top, Int* istack, Int* marked,
                     Int marker, Int* work) {
    // The current node i and the position p in its list of neighbours are
    // kept in local variables. When moving to a neighbour, (i,p) is saved in
    // istack[depth], work[depth]. Nodes are marked when reached and moved to
    // the output when all their neighbours have been processed.
    assert(marked[istart] != marker);
    Int depth = 0;
    Int i = istart;
    marked[i] = marker;
    Int j = colmap ? colmap[i] : i;
    Int p = j >= 0 ? Ap[j] : 0;
    Int pend = j >= 0 ? Ap[j+1] : 0;
    for (;;) {
        if (p < pend) {
            const Int inext = Ai[p++];
            if (marked[inext] != marker) {
                // save current node and go to neighbour
                marked[inext] = marker;
                istack[depth] = i;
                work[depth] = p;
                depth++;
                i = inext;
                j = colmap ? colmap[i] : i;
                p = j >= 0 ? Ap[j] : 0;
                pend = j >= 0 ? Ap[j+1] : 0;
            }
        } else {
            // node finished, go back to the last saved node
            istack[--top] = i;
            if (depth == 0)
                break;
            depth--;
            i = istack[depth];
            p = work[depth];
            j = colmap ? colmap[i] : i;
            pend = j >= 0 ? Ap[j+1] : 0;
        }
    }
    return top;
}

// Returns a free row in column j, or -1 if there is none. Rows of column j
// before position cheap[j] are known to be matched.
static Int FindFreeRow(Int j, const Int* Ap, const Int* Ai, const Int* jmatch,
                       Int* cheap) {
    for (Int p = cheap[j]; p < Ap[j+1]; p++) {
        if (jmatch[Ai[p]] == -1) {
            cheap[j] = p+1;
            return Ai[p];
        }
    }
    cheap[j] = Ap[j+1];
    return -1;
}

bool AugmentingPath(Int jstart, const Int* Ap, const Int* Ai, Int* jmatch,
                    Int* cheap, Int* marked, Int* work, Int* work2, Int* work3){
    // Depth first search for an alternating path from column jstart to a free
    // row. The current column j and the position p in it are kept in local
    // variables. When the path moves along row i to column jmatch[i],
    // (i,j,p) is saved in work[depth], work2[depth], work3[depth].
    Int depth = 0;
    Int j = jstart;
    marked[j] = jstart;
    Int ifree = FindFreeRow(j, Ap, Ai, jmatch, cheap);
    Int p = Ap[j];
    while (ifree < 0) {
        if (p < Ap[j+1]) {
            const Int i = Ai[p++];
            const Int jnext = jmatch[i];
            assert(jnext != -1);
            // follow row i unless it is to be ignored or its column is
            // already in the path
            if (jnext >= 0 && marked[jnext] != jstart) {
                work[depth] = i;
                work2[depth] = j;
                work3[depth] = p;
                depth++;
                j = jnext;
                marked[j] = jstart;
                ifree = FindFreeRow(j, Ap, Ai, jmatch, cheap);
                p = Ap[j];
            }
        } else {
            // no path through column j, go back
            if (depth == 0)
                return false;
            depth--;
            j = work2[depth];
            p = work3[depth];
        }
    }
    // match ifree to the last column and each saved row to its saved column
    jmatch[ifree] = j;
    while (depth > 0) {
        depth--;
        jmatch[work[depth]] = work2[depth];
    }
    return true;
}

}  // namespace ipx
