/**
\file
\brief Functions that deal with setting up the graphs for METIS.

\date   Started 7/25/1997
\author George  
\author Copyright 1997-2009, Regents of the University of Minnesota 
\version\verbatim $Id: graph.c 15817 2013-11-25 14:58:41Z karypis $ \endverbatim
*/

#include "metislib.h"


/*************************************************************************/
/*! This function sets up the graph from the user input */
/*************************************************************************/
graph_t *SetupGraph(ctrl_t *ctrl, idx_t nvtxs, idx_t ncon, const idx_t *xadj, 
             const idx_t *adjncy, idx_t *vwgt, idx_t *vsize, idx_t *adjwgt) 
{
  idx_t i, j, k, sum;
  real_t *nvwgt;
  graph_t *graph;

  /* allocate the graph and fill in the fields */
  graph = CreateGraph();

  graph->nvtxs  = nvtxs;
  graph->nedges = xadj[nvtxs];
  graph->ncon   = ncon;

  graph->xadj_orig = xadj;
  graph->xadj      = NULL;
  graph->free_xadj = 0;

  graph->adjncy_orig = adjncy;
  graph->adjncy      = NULL;
  graph->free_adjncy = 0;

  graph->droppedewgt = 0;

  /* setup the vertex weights */
  if (vwgt) {
    graph->vwgt      = vwgt;
    graph->free_vwgt = 0;
  }
  else {
    vwgt = graph->vwgt = ismalloc(ncon*nvtxs, 1);
  }

  graph->tvwgt    = imalloc(ncon);
  graph->invtvwgt = rmalloc(ncon);
  for (i=0; i<ncon; i++) {
    graph->tvwgt[i]    = isum(nvtxs, vwgt+i, ncon);
    graph->invtvwgt[i] = 1.0/(graph->tvwgt[i] > 0 ? graph->tvwgt[i] : 1);
  }


  /* setup the edge weights */
  if (adjwgt) {
    graph->adjwgt      = adjwgt;
    graph->free_adjwgt = 0;
  }
  else {
    adjwgt = graph->adjwgt = ismalloc(graph->nedges, 1);
  }
  


  /* setup various derived info */
  SetupGraph_tvwgt(graph);

  if (ctrl->optype == METIS_OP_OMETIS) 
    SetupGraph_label(graph);


  return graph;
}


/*************************************************************************/
/*! Set's up the tvwgt/invtvwgt info */
/*************************************************************************/
void SetupGraph_tvwgt(graph_t *graph)
{
  idx_t i;

  if (graph->tvwgt == NULL) 
    graph->tvwgt  = imalloc(graph->ncon);
  if (graph->invtvwgt == NULL) 
    graph->invtvwgt = rmalloc(graph->ncon);

  for (i=0; i<graph->ncon; i++) {
    graph->tvwgt[i]    = isum(graph->nvtxs, graph->vwgt+i, graph->ncon);
    graph->invtvwgt[i] = 1.0/(graph->tvwgt[i] > 0 ? graph->tvwgt[i] : 1);
  }
}


/*************************************************************************/
/*! Set's up the label info */
/*************************************************************************/
void SetupGraph_label(graph_t *graph)
{
  idx_t i;

  if (graph->label == NULL)
    graph->label = imalloc(graph->nvtxs);

  for (i=0; i<graph->nvtxs; i++)
    graph->label[i] = i;
}


/*************************************************************************/
/*! Setup the various arrays for the split graph */
/*************************************************************************/
graph_t *SetupSplitGraph(graph_t *graph, idx_t snvtxs, idx_t snedges)
{
  graph_t *sgraph;

  sgraph = CreateGraph();

  sgraph->nvtxs  = snvtxs;
  sgraph->nedges = snedges;
  sgraph->ncon   = graph->ncon;

  /* Allocate memory for the split graph */
  sgraph->xadj        = imalloc(snvtxs+1);
  sgraph->vwgt        = imalloc(sgraph->ncon*snvtxs);
  sgraph->adjncy      = imalloc(snedges);
  sgraph->adjwgt      = imalloc(snedges);
  sgraph->label	      = imalloc(snvtxs);
  sgraph->tvwgt       = imalloc(sgraph->ncon);
  sgraph->invtvwgt    = rmalloc(sgraph->ncon);

  if (graph->vsize)
    sgraph->vsize     = imalloc(snvtxs);

  return sgraph;
}


/*************************************************************************/
/*! This function creates and initializes a graph_t data structure */
/*************************************************************************/
graph_t *CreateGraph(void)
{
  graph_t *graph;

  graph = (graph_t *)malloc(sizeof(graph_t));

  InitGraph(graph);

  return graph;
}


/*************************************************************************/
/*! This function initializes a graph_t data structure */
/*************************************************************************/
void InitGraph(graph_t *graph) 
{
  memset((void *)graph, 0, sizeof(graph_t));

  /* graph size constants */
  graph->nvtxs     = -1;
  graph->nedges    = -1;
  graph->ncon      = -1;
  graph->mincut    = -1;
  graph->minvol    = -1;
  graph->nbnd      = -1;

  /* memory for the graph structure */
  graph->xadj      = NULL;
  graph->vwgt      = NULL;
  graph->vsize     = NULL;
  graph->adjncy    = NULL;
  graph->adjwgt    = NULL;
  graph->label     = NULL;
  graph->cmap      = NULL;
  graph->tvwgt     = NULL;
  graph->invtvwgt  = NULL;

  graph->xadj_orig   = NULL;
  graph->adjncy_orig = NULL;

  /* by default these are set to true, but the can be explicitly changed afterwards */
  graph->free_xadj   = 1;
  graph->free_vwgt   = 1;
  graph->free_vsize  = 1;
  graph->free_adjncy = 1;
  graph->free_adjwgt = 1;


  /* memory for the partition/refinement structure */
  graph->where     = NULL;
  graph->pwgts     = NULL;
  graph->id        = NULL;
  graph->ed        = NULL;
  graph->bndptr    = NULL;
  graph->bndind    = NULL;
  graph->nrinfo    = NULL;

  /* linked-list structure */
  graph->coarser   = NULL;
  graph->finer     = NULL;

}


/*************************************************************************/
/*! This function frees the memory storing the structure of the graph */
/*************************************************************************/
void FreeSData(graph_t *graph) 
{
  /* free graph structure */
  if (graph->free_xadj)
    gk_free((void **)&graph->xadj);
  if (graph->free_vwgt)
    gk_free((void **)&graph->vwgt);
  if (graph->free_vsize)
    gk_free((void **)&graph->vsize);
  if (graph->free_adjncy)
    gk_free((void **)&graph->adjncy);
  if (graph->free_adjwgt)
    gk_free((void **)&graph->adjwgt);
}


/*************************************************************************/
/*! This function frees the refinement/partition memory stored in a graph */
/*************************************************************************/
void FreeRData(graph_t *graph) 
{
  /* free partition/refinement structure */
  gk_free((void **)&graph->where);
  gk_free((void**)&graph->pwgts);
  gk_free((void**)&graph->id);
  gk_free((void**)&graph->ed);
  gk_free((void**)&graph->bndptr);
  gk_free((void**)&graph->bndind);
  gk_free((void**)&graph->nrinfo);
}


/*************************************************************************/
/*! This function deallocates any memory stored in a graph */
/*************************************************************************/
void FreeGraph(graph_t **r_graph) 
{
  graph_t *graph;

  graph = *r_graph;

  /* free the graph structure's fields */
  FreeSData(graph);

  /* free the partition/refinement fields */
  FreeRData(graph);

  gk_free((void **)&graph->tvwgt);
  gk_free((void**)&graph->invtvwgt);
  gk_free((void**)&graph->label);
  gk_free((void**)&graph->cmap);
  gk_free((void**)&graph);

  *r_graph = NULL;
}
