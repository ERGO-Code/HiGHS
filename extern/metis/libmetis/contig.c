/*!
\file 
\brief Functions that deal with eliminating disconnected partitions

\date Started 7/15/98
\author George
\author Copyright 1997-2009, Regents of the University of Minnesota 
\version $Id: contig.c 10513 2011-07-07 22:06:03Z karypis $
*/

#include "metislib.h"

/*************************************************************************/
/*! This function identifies the number of connected components in a graph
    that result after removing the vertices that belong to the vertex 
    separator (i.e., graph->where[i] == 2).
    The connected component memberships are returned in the CSR-style 
    pair of arrays cptr, cind.
*/
/**************************************************************************/
idx_t FindSepInducedComponents(ctrl_t *ctrl, graph_t *graph, idx_t *cptr, 
          idx_t *cind)
{
  idx_t i, j, k, nvtxs, first, last, nleft, ncmps, wgt;
  const idx_t *xadj, *adjncy;
  idx_t *where, *touched, *queue;

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj_orig ? graph->xadj_orig : graph->xadj;
  adjncy = graph->adjncy_orig ? graph->adjncy_orig : graph->adjncy;
  where  = graph->where;

  touched = ismalloc(nvtxs, 0);

  for (i=0; i<graph->nbnd; i++)
    touched[graph->bndind[i]] = 1;

  queue = cind;

  nleft = 0;
  for (i=0; i<nvtxs; i++) {
    if (where[i] != 2) 
      nleft++;
  }

  for (i=0; i<nvtxs; i++) {
    if (where[i] != 2)
      break;
  }

  touched[i] = 1;
  queue[0]   = i;
  first      = 0; 
  last       = 1;
  cptr[0]    = 0;  /* This actually points to queue */
  ncmps      = 0;

  while (first != nleft) {
    if (first == last) { /* Find another starting vertex */
      cptr[++ncmps] = first;
      for (i=0; i<nvtxs; i++) {
        if (!touched[i])
          break;
      }
      queue[last++] = i;
      touched[i] = 1;
    }

    i = queue[first++];
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];
      if (!touched[k]) {
        queue[last++] = k;
        touched[k] = 1;
      }
    }
  }
  cptr[++ncmps] = first;

  gk_free((void **)&touched);

  return ncmps;
}
