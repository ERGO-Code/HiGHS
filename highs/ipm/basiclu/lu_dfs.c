/*
 * lu_dfs.c
 *
 * Copyright (C) 2016-2018  ERGO-Code
 *
 * Depth first search in a graph
 *
 */

#include "ipm/basiclu/lu_internal.h"

static lu_int dfs_end(
    lu_int i, const lu_int *begin, const lu_int *end, const lu_int *index,
    lu_int top, lu_int *xi, lu_int *pstack, lu_int *marked, const lu_int M);

static lu_int dfs(
    lu_int i, const lu_int *begin, const lu_int *index,
    lu_int top, lu_int *xi, lu_int *pstack, lu_int *marked, const lu_int M);


/*
 * lu_dfs() - compute reach(i) in a graph by depth first search
 *
 * @begin, @end, @index define the graph. When end is not NULL, then node j
 * has neighbours index[begin[j]..end[j]-1]. When end is NULL, then the
 * neighbour list is terminated by a negative index.
 *
 * On return xi[newtop..top-1] hold reach(i) in topological order; newtop is
 * the function return value. Nodes that were already marked are excluded from
 * the reach.
 *
 * @pstack is size m workspace (m the number of nodes in the graph); the
 * contents of pstack is undefined on entry/return.
 *
 * @marked is size m array. Node j is marked iff marked[j] == M.
 * On return nodes xi[newtop..top-1] are marked.
 *
 * If node i is marked on entry, the function does nothing.
 *
 */

lu_int lu_dfs(
    lu_int i, const lu_int *begin, const lu_int *end, const lu_int *index,
    lu_int top, lu_int *xi, lu_int *pstack, lu_int *marked, const lu_int M)
{
    if (marked[i] == M)
        return top;

    return end ?
        dfs_end(i, begin, end, index, top, xi, pstack, marked, M) :
        dfs(i, begin, index, top, xi, pstack, marked, M);
}


/* ==========================================================================
 * dfs_end
 *
 * The current node i and the position p in its list of neighbours are kept in
 * local variables. When moving to a neighbour, (i,p) is saved in xi[depth],
 * pstack[depth]. Nodes are marked when reached and moved to the output when
 * all their neighbours have been processed.
 * ========================================================================== */

static lu_int dfs_end(
    lu_int i, const lu_int *begin, const lu_int *end, const lu_int *index,
    lu_int top, lu_int *xi, lu_int *pstack, lu_int *marked, const lu_int M)
{
    lu_int inext, p, depth = 0;
    assert(marked[i] != M);

    marked[i] = M;
    p = begin[i];
    for (;;)
    {
        if (p < end[i])
        {
            inext = index[p++];
            if (marked[inext] != M) /* save node i and go to neighbour */
            {
                marked[inext] = M;
                xi[depth] = i;
                pstack[depth] = p;
                depth++;
                i = inext;
                p = begin[i];
            }
        }
        else                    /* node i finished, go back */
        {
            xi[--top] = i;
            if (depth == 0)
                break;
            depth--;
            i = xi[depth];
            p = pstack[depth];
        }
    }

    return top;
}


/* ==========================================================================
 * dfs
 *
 * Same as dfs_end, but the neighbour lists are terminated by a negative index.
 * ========================================================================== */

static lu_int dfs(
    lu_int i, const lu_int *begin, const lu_int *index,
    lu_int top, lu_int *xi, lu_int *pstack, lu_int *marked, const lu_int M)
{
    lu_int inext, p, depth = 0;
    assert(marked[i] != M);

    marked[i] = M;
    p = begin[i];
    for (;;)
    {
        if ((inext = index[p]) >= 0)
        {
            p++;
            if (marked[inext] != M) /* save node i and go to neighbour */
            {
                marked[inext] = M;
                xi[depth] = i;
                pstack[depth] = p;
                depth++;
                i = inext;
                p = begin[i];
            }
        }
        else                    /* node i finished, go back */
        {
            xi[--top] = i;
            if (depth == 0)
                break;
            depth--;
            i = xi[depth];
            p = pstack[depth];
        }
    }

    return top;
}
