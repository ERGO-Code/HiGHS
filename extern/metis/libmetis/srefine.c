/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * srefine.c
 *
 * This file contains code for the separator refinement algorithms
 *
 * Started 8/1/97
 * George
 *
 * $Id: srefine.c 14362 2013-05-21 21:35:23Z karypis $
 *
 */

#include "metislib.h"


/*************************************************************************/
/*! This function is the entry point of the separator refinement. 
    It does not perform any refinement on graph, but it starts by first
    projecting it to the next level finer graph and proceeds from there. */
/*************************************************************************/
void Refine2WayNode(ctrl_t *ctrl, graph_t *orggraph, graph_t *graph)
{


  if (graph == orggraph) {
    Compute2WayNodePartitionParams(ctrl, graph);
  }
  else {
    do {
      graph = graph->finer;

      Project2WayNodePartition(ctrl, graph);

      FM_2WayNodeBalance(ctrl, graph); 


      switch (ctrl->rtype) {
        case METIS_RTYPE_SEP2SIDED:
          FM_2WayNodeRefine2Sided(ctrl, graph, ctrl->niter); 
          break;
        case METIS_RTYPE_SEP1SIDED:
          FM_2WayNodeRefine1Sided(ctrl, graph, ctrl->niter); 
          break;
        default:
          gk_errexit("Unknown rtype of %d\n", ctrl->rtype);
      }

    } while (graph != orggraph);
  }

}


/*************************************************************************/
/*! This function allocates memory for 2-way node-based refinement */
/**************************************************************************/
void Allocate2WayNodePartitionMemory(ctrl_t *ctrl, graph_t *graph)
{
  idx_t nvtxs;

  nvtxs = graph->nvtxs;

  graph->pwgts  = imalloc(3);
  graph->where  = imalloc(nvtxs);
  graph->bndptr = imalloc(nvtxs);
  graph->bndind = imalloc(nvtxs);
  graph->nrinfo = (nrinfo_t *)malloc(nvtxs*sizeof(nrinfo_t));
}


/*************************************************************************/
/*! This function computes the edegrees[] to the left & right sides */
/*************************************************************************/
void Compute2WayNodePartitionParams(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, j, nvtxs, nbnd;
  const idx_t *xadj, *adjncy;
  idx_t *vwgt;
  idx_t *where, *pwgts, *bndind, *bndptr, *edegrees;
  nrinfo_t *rinfo;
  idx_t me, other;

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj_orig ? graph->xadj_orig : graph->xadj;
  vwgt   = graph->vwgt;
  adjncy = graph->adjncy_orig ? graph->adjncy_orig : graph->adjncy;

  where  = graph->where;
  rinfo  = graph->nrinfo;
  pwgts  = iset(3, 0, graph->pwgts);
  bndind = graph->bndind;
  bndptr = iset(nvtxs, -1, graph->bndptr);


  /*------------------------------------------------------------
  / Compute now the separator external degrees
  /------------------------------------------------------------*/
  nbnd = 0;
  for (i=0; i<nvtxs; i++) {
    me = where[i];
    pwgts[me] += vwgt[i];

    if (me == 2) { /* If it is on the separator do some computations */
      BNDInsert(nbnd, bndind, bndptr, i);

      edegrees = rinfo[i].edegrees;
      edegrees[0] = edegrees[1] = 0;

      for (j=xadj[i]; j<xadj[i+1]; j++) {
        other = where[adjncy[j]];
        if (other != 2)
          edegrees[other] += vwgt[adjncy[j]];
      }
    }
  }


  graph->mincut = pwgts[2];
  graph->nbnd   = nbnd;
}


/*************************************************************************/
/*! This function projects the node-based bisection */
/*************************************************************************/
void Project2WayNodePartition(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, j, nvtxs;
  idx_t *cmap, *where, *cwhere;
  graph_t *cgraph;

  cgraph = graph->coarser;
  cwhere = cgraph->where;

  nvtxs = graph->nvtxs;
  cmap  = graph->cmap;

  Allocate2WayNodePartitionMemory(ctrl, graph);
  where = graph->where;
  
  /* Project the partition */
  for (i=0; i<nvtxs; i++) {
    where[i] = cwhere[cmap[i]];
  }

  FreeGraph(&graph->coarser);
  graph->coarser = NULL;

  Compute2WayNodePartitionParams(ctrl, graph);
}
