/*
\file
\brief This file contains the driving routines for multilevel refinement

\date   Started 7/24/1997
\author George  
\author Copyright 1997-2009, Regents of the University of Minnesota 
\version\verbatim $Id: refine.c 14362 2013-05-21 21:35:23Z karypis $ \endverbatim
*/

#include "metislib.h"


/*************************************************************************/
/*! This function allocates memory for 2-way edge refinement */
/*************************************************************************/
void Allocate2WayPartitionMemory(ctrl_t *ctrl, graph_t *graph)
{
  idx_t nvtxs, ncon;

  nvtxs = graph->nvtxs;
  ncon  = graph->ncon;

  graph->pwgts  = imalloc(2*ncon);
  graph->where  = imalloc(nvtxs);
  graph->bndptr = imalloc(nvtxs);
  graph->bndind = imalloc(nvtxs);
  graph->id     = imalloc(nvtxs);
  graph->ed     = imalloc(nvtxs);
}


/*************************************************************************/
/*! This function computes the initial id/ed */
/*************************************************************************/
void Compute2WayPartitionParams(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, j, nvtxs, ncon, nbnd, mincut, istart, iend, tid, ted, me;
  const idx_t *xadj, *adjncy;
  idx_t *vwgt, *adjwgt, *pwgts;
  idx_t *where, *bndptr, *bndind, *id, *ed;

  nvtxs  = graph->nvtxs;
  ncon   = graph->ncon;
  xadj   = graph->xadj_orig ? graph->xadj_orig : graph->xadj;
  vwgt   = graph->vwgt;
  adjncy = graph->adjncy_orig ? graph->adjncy_orig : graph->adjncy;
  adjwgt = graph->adjwgt;

  where  = graph->where;
  id     = graph->id;
  ed     = graph->ed;

  pwgts  = iset(2*ncon, 0, graph->pwgts);
  bndptr = iset(nvtxs, -1, graph->bndptr);
  bndind = graph->bndind;

  /* Compute pwgts */
  if (ncon == 1) {
    for (i=0; i<nvtxs; i++) {
      pwgts[where[i]] += vwgt[i];
    }
  }
  else {
    for (i=0; i<nvtxs; i++) {
      me = where[i];
      for (j=0; j<ncon; j++)
        pwgts[me*ncon+j] += vwgt[i*ncon+j];
    }
  }


  /* Compute the required info for refinement  */
  for (nbnd=0, mincut=0, i=0; i<nvtxs; i++) {
    istart = xadj[i];
    iend   = xadj[i+1];

    me = where[i];
    tid = ted = 0;

    for (j=istart; j<iend; j++) {
      if (me == where[adjncy[j]])
        tid += adjwgt[j];
      else
        ted += adjwgt[j];
    }
    id[i] = tid;
    ed[i] = ted;
  
    if (ted > 0 || istart == iend) {
      BNDInsert(nbnd, bndind, bndptr, i);
      mincut += ted;
    }
  }

  graph->mincut = mincut/2;
  graph->nbnd   = nbnd;

}
