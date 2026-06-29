/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * separator.c
 *
 * This file contains code for separator extraction
 *
 * Started 8/1/97
 * George
 *
 * $Id: separator.c 10481 2011-07-05 18:01:23Z karypis $
 *
 */

#include "metislib.h"

/*************************************************************************
* This function takes a bisection and constructs a minimum weight vertex 
* separator out of it. It uses the node-based separator refinement for it.
**************************************************************************/
void ConstructSeparator(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, j, k, nvtxs, nbnd;
  const idx_t *xadj;
  idx_t *where, *bndind;

  WCOREPUSH;

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj_orig ? graph->xadj_orig : graph->xadj;
  nbnd   = graph->nbnd;
  bndind = graph->bndind;

  where = icopy(nvtxs, graph->where, iwspacemalloc(ctrl, nvtxs));

  /* Put the nodes in the boundary into the separator */
  for (i=0; i<nbnd; i++) {
    j = bndind[i];
    if (xadj[j+1]-xadj[j] > 0)  /* Ignore islands */
      where[j] = 2;
  }

  FreeRData(graph);

  Allocate2WayNodePartitionMemory(ctrl, graph);
  icopy(nvtxs, where, graph->where);

  WCOREPOP;


  Compute2WayNodePartitionParams(ctrl, graph);

  FM_2WayNodeRefine2Sided(ctrl, graph, 1); 
  FM_2WayNodeRefine1Sided(ctrl, graph, 4); 


}
