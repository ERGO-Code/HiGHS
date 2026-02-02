/*!
\file  
\brief Functions for computing matchings during graph coarsening

\date Started 7/23/97
\author George  
\author Copyright 1997-2011, Regents of the University of Minnesota 
\version\verbatim $Id: coarsen.c 20398 2016-11-22 17:17:12Z karypis $ \endverbatim
*/


#include "metislib.h"

#define UNMATCHEDFOR2HOP  0.10  /* The fraction of unmatched vertices that triggers 2-hop */
                                  

/*************************************************************************/
/*! This function takes a graph and creates a sequence of coarser graphs.
    It implements the coarsening phase of the multilevel paradigm. 
 */
/*************************************************************************/
graph_t *CoarsenGraph(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, eqewgts, level=0;


  /* determine if the weights on the edges are all the same */
  for (eqewgts=1, i=1; i<graph->nedges; i++) {
    if (graph->adjwgt[0] != graph->adjwgt[i]) {
      eqewgts = 0;
      break;
    }
  }

  /* set the maximum allowed coarsest vertex weight */
  for (i=0; i<graph->ncon; i++)
    ctrl->maxvwgt[i] = 1.5*graph->tvwgt[i]/ctrl->CoarsenTo;

  do {
    IFSET(ctrl->dbglvl, METIS_DBG_COARSEN, PrintCGraphStats(ctrl, graph));

    /* allocate memory for cmap, if it has not already been done due to
       multiple cuts */
    if (graph->cmap == NULL)
      graph->cmap = imalloc(graph->nvtxs);

    /* determine which matching scheme you will use */
    switch (ctrl->ctype) {
      case METIS_CTYPE_RM:
        Match_RM(ctrl, graph);
        break;
      case METIS_CTYPE_SHEM:
        if (eqewgts || graph->nedges == 0)
          Match_RM(ctrl, graph);
        else
          Match_SHEM(ctrl, graph);
        break;
      default:
        gk_errexit("Unknown ctype: %d\n", ctrl->ctype);
    }

    graph = graph->coarser;
    eqewgts = 0;
    level++;


  } while (graph->nvtxs > ctrl->CoarsenTo && 
           graph->nvtxs < COARSEN_FRACTION*graph->finer->nvtxs && 
           graph->nedges > graph->nvtxs/2);

  IFSET(ctrl->dbglvl, METIS_DBG_COARSEN, PrintCGraphStats(ctrl, graph));

  return graph;
}


/*************************************************************************/
/*! This function takes a graph and creates a sequence of nlevels coarser 
    graphs, where nlevels is an input parameter.
 */
/*************************************************************************/
graph_t *CoarsenGraphNlevels(ctrl_t *ctrl, graph_t *graph, idx_t nlevels)
{
  idx_t i, eqewgts, level;


  /* determine if the weights on the edges are all the same */
  for (eqewgts=1, i=1; i<graph->nedges; i++) {
    if (graph->adjwgt[0] != graph->adjwgt[i]) {
      eqewgts = 0;
      break;
    }
  }

  /* set the maximum allowed coarsest vertex weight */
  for (i=0; i<graph->ncon; i++)
    ctrl->maxvwgt[i] = 1.5*graph->tvwgt[i]/ctrl->CoarsenTo;

  for (level=0; level<nlevels; level++) {
    IFSET(ctrl->dbglvl, METIS_DBG_COARSEN, PrintCGraphStats(ctrl, graph));

    /* allocate memory for cmap, if it has not already been done due to
       multiple cuts */
    if (graph->cmap == NULL)
      graph->cmap = imalloc(graph->nvtxs);

    /* determine which matching scheme you will use */
    switch (ctrl->ctype) {
      case METIS_CTYPE_RM:
        Match_RM(ctrl, graph);
        break;
      case METIS_CTYPE_SHEM:
        if (eqewgts || graph->nedges == 0)
          Match_RM(ctrl, graph);
        else
          Match_SHEM(ctrl, graph);
        break;
      default:
        gk_errexit("Unknown ctype: %d\n", ctrl->ctype);
    }

    graph = graph->coarser;
    eqewgts = 0;


    if (graph->nvtxs < ctrl->CoarsenTo || 
        graph->nvtxs > COARSEN_FRACTION*graph->finer->nvtxs || 
        graph->nedges < graph->nvtxs/2)
      break; 
  } 

  IFSET(ctrl->dbglvl, METIS_DBG_COARSEN, PrintCGraphStats(ctrl, graph));

  return graph;
}


/*************************************************************************/
/*! This function finds a matching by randomly selecting one of the 
    unmatched adjacent vertices. 
 */
/**************************************************************************/
idx_t Match_RM(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, pi, ii, j, jj, jjinc, k, nvtxs, ncon, cnvtxs, maxidx, 
        last_unmatched, avgdegree, bnum;
  const idx_t *xadj, *adjncy;
  idx_t *vwgt, *adjwgt, *maxvwgt;
  idx_t *match, *cmap, *degrees, *perm, *tperm;
  size_t nunmatched=0;

  WCOREPUSH;


  nvtxs  = graph->nvtxs;
  ncon   = graph->ncon;
  xadj   = graph->xadj_orig ? graph->xadj_orig : graph->xadj;
  vwgt   = graph->vwgt;
  adjncy = graph->adjncy_orig ? graph->adjncy_orig : graph->adjncy;
  adjwgt = graph->adjwgt;
  cmap   = graph->cmap;

  maxvwgt  = ctrl->maxvwgt;

  match   = iset(nvtxs, UNMATCHED, iwspacemalloc(ctrl, nvtxs));
  perm    = iwspacemalloc(ctrl, nvtxs);
  tperm   = iwspacemalloc(ctrl, nvtxs);
  degrees = iwspacemalloc(ctrl, nvtxs);

  /* Determine a "random" traversal order that is biased towards 
     low-degree vertices */
  irandArrayPermute(nvtxs, tperm, nvtxs/8, 1, &ctrl->rng_state);

  avgdegree = 4.0*(xadj[nvtxs]/nvtxs);
  for (i=0; i<nvtxs; i++) {
    bnum = sqrt(1+xadj[i+1]-xadj[i]);
    degrees[i] = (bnum > avgdegree ? avgdegree : bnum);
  }
  BucketSortKeysInc(ctrl, nvtxs, avgdegree, degrees, tperm, perm);


  /* Traverse the vertices and compute the matching */
  for (cnvtxs=0, last_unmatched=0, pi=0; pi<nvtxs; pi++) {
    i = perm[pi];

    if (match[i] == UNMATCHED) {  /* Unmatched */
      maxidx = i;

      if ((ncon == 1 ? vwgt[i] < maxvwgt[0] : ivecle(ncon, vwgt+i*ncon, maxvwgt))) {
        /* Deal with island vertices. Find a non-island and match it with. 
           The matching ignores ctrl->maxvwgt requirements */
        if (xadj[i] == xadj[i+1]) {
          last_unmatched = gk_max(pi, last_unmatched)+1;
          for (; last_unmatched<nvtxs; last_unmatched++) {
            j = perm[last_unmatched];
            if (match[j] == UNMATCHED) {
              maxidx = j;
              break;
            }
          }
        }
        else {
          /* Find a random matching, subject to maxvwgt constraints */
          if (ncon == 1) {
            /* single constraint version */
            for (j=xadj[i]; j<xadj[i+1]; j++) {
              k = adjncy[j];
              if (match[k] == UNMATCHED && vwgt[i]+vwgt[k] <= maxvwgt[0]) {
                maxidx = k;
                break;
              }
            }

            /* If it did not match, record for a 2-hop matching. */
            if (maxidx == i && 2*vwgt[i] < maxvwgt[0]) {
              nunmatched++;
              maxidx = UNMATCHED;
            }
          }
          else {
            /* multi-constraint version */
            for (j=xadj[i]; j<xadj[i+1]; j++) {
              k = adjncy[j];
              if (match[k] == UNMATCHED && 
                  ivecaxpylez(ncon, 1, vwgt+i*ncon, vwgt+k*ncon, maxvwgt)) {
                maxidx = k;
                break;
              }
            }

            /* If it did not match, record for a 2-hop matching. */
            if (maxidx == i && ivecaxpylez(ncon, 2, vwgt+i*ncon, vwgt+i*ncon, maxvwgt)) {
              nunmatched++;
              maxidx = UNMATCHED;
            }
          }
        }
      }

      if (maxidx != UNMATCHED) {
        cmap[i]  = cmap[maxidx] = cnvtxs++;
        match[i] = maxidx;
        match[maxidx] = i;
      }
    }
  }

  //printf("nunmatched: %zu\n", nunmatched);

  /* see if a 2-hop matching is required/allowed */
  if (!ctrl->no2hop && nunmatched > UNMATCHEDFOR2HOP*nvtxs) 
    cnvtxs = Match_2Hop(ctrl, graph, perm, match, cnvtxs, nunmatched);


  /* match the final unmatched vertices with themselves and reorder the vertices 
     of the coarse graph for memory-friendly contraction */
  for (cnvtxs=0, i=0; i<nvtxs; i++) {
    if (match[i] == UNMATCHED) {
      match[i] = i;
      cmap[i]  = cnvtxs++;
    }
    else {
      if (i <= match[i]) 
        cmap[i] = cmap[match[i]] = cnvtxs++;
    }
  }


  CreateCoarseGraph(ctrl, graph, cnvtxs, match);

  WCOREPOP;

  return cnvtxs;
}


/**************************************************************************/
/*! This function finds a matching using the HEM heuristic. The vertices 
    are visited based on increasing degree to ensure that all vertices are 
    given a chance to match with something. 
 */
/**************************************************************************/
idx_t Match_SHEM(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, pi, ii, j, jj, jjinc, k, nvtxs, ncon, cnvtxs, maxidx, maxwgt, 
        last_unmatched, avgdegree, bnum;
  const idx_t *xadj, *adjncy;
  idx_t *vwgt, *adjwgt, *maxvwgt;
  idx_t *match, *cmap, *degrees, *perm, *tperm;
  size_t nunmatched=0;

  WCOREPUSH;


  nvtxs  = graph->nvtxs;
  ncon   = graph->ncon;
  xadj   = graph->xadj_orig ? graph->xadj_orig : graph->xadj;
  vwgt   = graph->vwgt;
  adjncy = graph->adjncy_orig ? graph->adjncy_orig : graph->adjncy;
  adjwgt = graph->adjwgt;
  cmap   = graph->cmap;

  maxvwgt  = ctrl->maxvwgt;

  match   = iset(nvtxs, UNMATCHED, iwspacemalloc(ctrl, nvtxs));
  perm    = iwspacemalloc(ctrl, nvtxs);
  tperm   = iwspacemalloc(ctrl, nvtxs);
  degrees = iwspacemalloc(ctrl, nvtxs);

  /* Determine a "random" traversal order that is biased towards low-degree vertices */
  irandArrayPermute(nvtxs, tperm, nvtxs/8, 1, &ctrl->rng_state);

  avgdegree = 4.0*(xadj[nvtxs]/nvtxs);
  for (i=0; i<nvtxs; i++) {
    bnum = sqrt(1+xadj[i+1]-xadj[i]);
    degrees[i] = (bnum > avgdegree ? avgdegree : bnum);
  }
  BucketSortKeysInc(ctrl, nvtxs, avgdegree, degrees, tperm, perm);


  /* Traverse the vertices and compute the matching */
  for (cnvtxs=0, last_unmatched=0, pi=0; pi<nvtxs; pi++) {
    i = perm[pi];

    if (match[i] == UNMATCHED) {  /* Unmatched */
      maxidx = i;
      maxwgt = -1;

      if ((ncon == 1 ? vwgt[i] < maxvwgt[0] : ivecle(ncon, vwgt+i*ncon, maxvwgt))) {
        /* Deal with island vertices. Find a non-island and match it with. 
           The matching ignores ctrl->maxvwgt requirements */
        if (xadj[i] == xadj[i+1]) { 
          last_unmatched = gk_max(pi, last_unmatched)+1;
          for (; last_unmatched<nvtxs; last_unmatched++) {
            j = perm[last_unmatched];
            if (match[j] == UNMATCHED) {
              maxidx = j;
              break;
            }
          }
        }
        else {
          /* Find a heavy-edge matching, subject to maxvwgt constraints */
          if (ncon == 1) {
            /* single constraint version */
            for (j=xadj[i]; j<xadj[i+1]; j++) {
              k = adjncy[j];
              if (match[k] == UNMATCHED && 
                  maxwgt < adjwgt[j] && vwgt[i]+vwgt[k] <= maxvwgt[0]) {
                maxidx = k;
                maxwgt = adjwgt[j];
              }
            }

            /* If it did not match, record for a 2-hop matching. */
            if (maxidx == i && 2*vwgt[i] < maxvwgt[0]) {
              nunmatched++;
              maxidx = UNMATCHED;
            }
          }
          else {
            /* multi-constraint version */
            for (j=xadj[i]; j<xadj[i+1]; j++) {
              k = adjncy[j];
              if (match[k] == UNMATCHED && 
                  ivecaxpylez(ncon, 1, vwgt+i*ncon, vwgt+k*ncon, maxvwgt) &&
                  (maxwgt < adjwgt[j] || 
                   (maxwgt == adjwgt[j] && 
                    BetterVBalance(ncon, graph->invtvwgt, vwgt+i*ncon, 
                        vwgt+maxidx*ncon, vwgt+k*ncon)))) {
                maxidx = k;
                maxwgt = adjwgt[j];
              }
            }

            /* If it did not match, record for a 2-hop matching. */
            if (maxidx == i && ivecaxpylez(ncon, 2, vwgt+i*ncon, vwgt+i*ncon, maxvwgt)) {
              nunmatched++;
              maxidx = UNMATCHED;
            }
          }
        }
      }

      if (maxidx != UNMATCHED) {
        cmap[i]  = cmap[maxidx] = cnvtxs++;
        match[i] = maxidx;
        match[maxidx] = i;
      }
    }
  }

  //printf("nunmatched: %zu\n", nunmatched);

  /* see if a 2-hop matching is required/allowed */
  if (!ctrl->no2hop && nunmatched > UNMATCHEDFOR2HOP*nvtxs) 
    cnvtxs = Match_2Hop(ctrl, graph, perm, match, cnvtxs, nunmatched);


  /* match the final unmatched vertices with themselves and reorder the vertices 
     of the coarse graph for memory-friendly contraction */
  for (cnvtxs=0, i=0; i<nvtxs; i++) {
    if (match[i] == UNMATCHED) {
      match[i] = i;
      cmap[i] = cnvtxs++;
    }
    else {
      if (i <= match[i]) 
        cmap[i] = cmap[match[i]] = cnvtxs++;
    }
  }


  CreateCoarseGraph(ctrl, graph, cnvtxs, match);

  WCOREPOP;

  return cnvtxs;
}


/*************************************************************************/
/*! This function matches the unmatched vertices using a 2-hop matching 
    that involves vertices that are two hops away from each other. */
/**************************************************************************/
idx_t Match_2Hop(ctrl_t *ctrl, graph_t *graph, idx_t *perm, idx_t *match, 
          idx_t cnvtxs, size_t nunmatched)
{

  cnvtxs = Match_2HopAny(ctrl, graph, perm, match, cnvtxs, &nunmatched, 2);
  cnvtxs = Match_2HopAll(ctrl, graph, perm, match, cnvtxs, &nunmatched, 64);
  if (nunmatched > 1.5*UNMATCHEDFOR2HOP*graph->nvtxs) 
    cnvtxs = Match_2HopAny(ctrl, graph, perm, match, cnvtxs, &nunmatched, 3);
  if (nunmatched > 2.0*UNMATCHEDFOR2HOP*graph->nvtxs) 
    cnvtxs = Match_2HopAny(ctrl, graph, perm, match, cnvtxs, &nunmatched, graph->nvtxs);

  return cnvtxs;
}


/*************************************************************************/
/*! This function matches the unmatched vertices whose degree is less than
    maxdegree using a 2-hop matching that involves vertices that are two 
    hops away from each other. 
    The requirement of the 2-hop matching is a simple non-empty overlap
    between the adjacency lists of the vertices. */
/**************************************************************************/
idx_t Match_2HopAny(ctrl_t *ctrl, graph_t *graph, idx_t *perm, idx_t *match, 
          idx_t cnvtxs, size_t *r_nunmatched, size_t maxdegree)
{
  idx_t i, pi, ii, j, jj, k, nvtxs;
  const idx_t *xadj, *adjncy;
  idx_t *colptr, *rowind;
  idx_t *cmap;
  size_t nunmatched;


  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj_orig ? graph->xadj_orig : graph->xadj;
  adjncy = graph->adjncy_orig ? graph->adjncy_orig : graph->adjncy;
  cmap   = graph->cmap;

  nunmatched = *r_nunmatched;

  /*IFSET(ctrl->dbglvl, METIS_DBG_COARSEN, printf("IN: nunmatched: %zu\t", nunmatched)); */

  /* create the inverted index */
  WCOREPUSH;
  colptr = iset(nvtxs, 0, iwspacemalloc(ctrl, nvtxs+1));
  for (i=0; i<nvtxs; i++) {
    if (match[i] == UNMATCHED && xadj[i+1]-xadj[i] < maxdegree) {
      for (j=xadj[i]; j<xadj[i+1]; j++)
        colptr[adjncy[j]]++;
    }
  }
  MAKECSR(i, nvtxs, colptr);

  rowind = iwspacemalloc(ctrl, colptr[nvtxs]);
  for (pi=0; pi<nvtxs; pi++) {
    i = perm[pi];
    if (match[i] == UNMATCHED && xadj[i+1]-xadj[i] < maxdegree) {
      for (j=xadj[i]; j<xadj[i+1]; j++)
        rowind[colptr[adjncy[j]]++] = i;
    }
  }
  SHIFTCSR(i, nvtxs, colptr);

  /* compute matchings by going down the inverted index */
  for (pi=0; pi<nvtxs; pi++) {
    i = perm[pi];
    if (colptr[i+1]-colptr[i] < 2)
      continue;

    for (jj=colptr[i+1], j=colptr[i]; j<jj; j++) {
      if (match[rowind[j]] == UNMATCHED) {
        for (jj--; jj>j; jj--) {
          if (match[rowind[jj]] == UNMATCHED) {
            cmap[rowind[j]] = cmap[rowind[jj]] = cnvtxs++;
            match[rowind[j]]  = rowind[jj];
            match[rowind[jj]] = rowind[j];
            nunmatched -= 2;
            break;
          }
        }
      }
    }
  }
  WCOREPOP;

  *r_nunmatched = nunmatched;
  return cnvtxs;
}


/*************************************************************************/
/*! This function matches the unmatched vertices whose degree is less than
    maxdegree using a 2-hop matching that involves vertices that are two 
    hops away from each other. 
    The requirement of the 2-hop matching is that of identical adjacency
    lists.
 */
/**************************************************************************/
idx_t Match_2HopAll(ctrl_t *ctrl, graph_t *graph, idx_t *perm, idx_t *match, 
          idx_t cnvtxs, size_t *r_nunmatched, size_t maxdegree)
{
  idx_t i, pi, pk, ii, j, jj, k, nvtxs, mask, idegree;
  const idx_t *xadj, *adjncy;
  idx_t *cmap, *mark;
  ikv_t *keys;
  size_t nunmatched, ncand;

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj_orig ? graph->xadj_orig : graph->xadj;
  adjncy = graph->adjncy_orig ? graph->adjncy_orig : graph->adjncy;
  cmap   = graph->cmap;

  nunmatched = *r_nunmatched;
  mask = IDX_MAX/maxdegree;

  /*IFSET(ctrl->dbglvl, METIS_DBG_COARSEN, printf("IN: nunmatched: %zu\t", nunmatched)); */

  WCOREPUSH;

  /* collapse vertices with identical adjacency lists */
  keys = ikvwspacemalloc(ctrl, nunmatched);
  for (ncand=0, pi=0; pi<nvtxs; pi++) {
    i = perm[pi];
    idegree = xadj[i+1]-xadj[i];
    if (match[i] == UNMATCHED && idegree > 1 && idegree < maxdegree) {
      for (k=0, j=xadj[i]; j<xadj[i+1]; j++) 
        k += adjncy[j]%mask;
      keys[ncand].val = i;
      keys[ncand].key = (k%mask)*maxdegree + idegree;
      ncand++;
    }
  }
  ikvsorti(ncand, keys);

  mark = iset(nvtxs, 0, iwspacemalloc(ctrl, nvtxs));
  for (pi=0; pi<ncand; pi++) {
    i = keys[pi].val;
    if (match[i] != UNMATCHED)
      continue;

    for (j=xadj[i]; j<xadj[i+1]; j++)
      mark[adjncy[j]] = i;

    for (pk=pi+1; pk<ncand; pk++) {
      k = keys[pk].val;
      if (match[k] != UNMATCHED)
        continue;

      if (keys[pi].key != keys[pk].key)
        break;
      if (xadj[i+1]-xadj[i] != xadj[k+1]-xadj[k])
        break;

      for (jj=xadj[k]; jj<xadj[k+1]; jj++) {
        if (mark[adjncy[jj]] != i)
          break;
      }
      if (jj == xadj[k+1]) {
        cmap[i] = cmap[k] = cnvtxs++;
        match[i] = k;
        match[k] = i;
        nunmatched -= 2;
        break;
      }
    }
  }
  WCOREPOP;

  *r_nunmatched = nunmatched;
  return cnvtxs;
}

/*************************************************************************/
/*! This function prints various stats for each graph during coarsening 
 */
/*************************************************************************/
void PrintCGraphStats(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i;

  printf("%10"PRIDX" %10"PRIDX" %10"PRIDX" [%"PRIDX"] [", 
      graph->nvtxs, graph->nedges, isum(graph->nedges, graph->adjwgt, 1), ctrl->CoarsenTo);

  for (i=0; i<graph->ncon; i++)
    printf(" %8"PRIDX":%8"PRIDX, ctrl->maxvwgt[i], graph->tvwgt[i]);
  printf(" ]\n");
}


/*************************************************************************/
/*! This function creates the coarser graph. Depending on the size of the
    candidate adjacency lists it either uses a hash table or an array
    to do duplicate detection.
 */
/*************************************************************************/
void CreateCoarseGraph(ctrl_t *ctrl, graph_t *graph, idx_t cnvtxs, 
         idx_t *match)
{
  idx_t j, jj, k, kk, l, m, istart, iend, nvtxs, nedges, ncon, 
        cnedges, v, u, mask;
  const idx_t *xadj, *adjncy;
  idx_t *vwgt, *vsize, *adjwgt;
  idx_t *cmap, *htable, *dtable;
  idx_t *cxadj, *cvwgt, *cvsize, *cadjncy, *cadjwgt;
  graph_t *cgraph;
  int dovsize, dropedges;
  idx_t cv, nkeys, droppedewgt;
  idx_t *keys=NULL, *medianewgts=NULL, *noise=NULL;

  WCOREPUSH;

  dovsize   = 0;
  dropedges = ctrl->dropedges;

  mask = HTLENGTH;

  nvtxs   = graph->nvtxs;
  ncon    = graph->ncon;
  xadj    = graph->xadj_orig ? graph->xadj_orig : graph->xadj;
  vwgt    = graph->vwgt;
  vsize   = graph->vsize;
  adjncy  = graph->adjncy_orig ? graph->adjncy_orig : graph->adjncy;
  adjwgt  = graph->adjwgt;
  cmap    = graph->cmap;

  /* Setup structures for dropedges */
  if (dropedges) {
    for (nkeys=0, v=0; v<nvtxs; v++) 
      nkeys = gk_max(nkeys, xadj[v+1]-xadj[v]);
    nkeys = 2*nkeys+1;

    keys        = iwspacemalloc(ctrl, nkeys);
    noise       = iwspacemalloc(ctrl, cnvtxs);
    medianewgts = iset(cnvtxs, -1, iwspacemalloc(ctrl, cnvtxs));

    for (v=0; v<cnvtxs; v++) 
      noise[v] = irandInRange(128, &ctrl->rng_state);
  }

  /* Initialize the coarser graph */
  cgraph   = SetupCoarseGraph(graph, cnvtxs, dovsize);
  cxadj    = cgraph->xadj;
  cvwgt    = cgraph->vwgt;
  cvsize   = cgraph->vsize;
  cadjncy  = cgraph->adjncy;
  cadjwgt  = cgraph->adjwgt;

  htable = iset(mask+1, -1, iwspacemalloc(ctrl, mask+1));   /* hash table */
  dtable = iset(cnvtxs, -1, iwspacemalloc(ctrl, cnvtxs));   /* direct table */

  cxadj[0] = cnvtxs = cnedges = 0;
  for (v=0; v<nvtxs; v++) {
    if ((u = match[v]) < v)
      continue;


    /* take care of the vertices */
    if (ncon == 1)
      cvwgt[cnvtxs] = vwgt[v];
    else
      icopy(ncon, vwgt+v*ncon, cvwgt+cnvtxs*ncon);

    if (dovsize)
      cvsize[cnvtxs] = vsize[v];

    if (v != u) { 
      if (ncon == 1)
        cvwgt[cnvtxs] += vwgt[u];
      else
        iaxpy(ncon, 1, vwgt+u*ncon, 1, cvwgt+cnvtxs*ncon, 1);

      if (dovsize)
        cvsize[cnvtxs] += vsize[u];
    }


    /* take care of the edges */ 
    if ((xadj[v+1]-xadj[v] + xadj[u+1]-xadj[u]) < (mask>>2)) { /* use mask */
      /* put the ID of the contracted node itself at the start, so that it can be 
       * removed easily */
      htable[cnvtxs&mask] = 0;
      cadjncy[0] = cnvtxs;
      nedges = 1;

      istart = xadj[v];
      iend   = xadj[v+1];
      for (j=istart; j<iend; j++) {
        k = cmap[adjncy[j]];
        for (kk=k&mask; htable[kk]!=-1 && cadjncy[htable[kk]]!=k; kk=((kk+1)&mask));
        if ((m = htable[kk]) == -1) {
          cadjncy[nedges] = k;
          cadjwgt[nedges] = adjwgt[j];
          htable[kk] = nedges++;
        }
        else {
          cadjwgt[m] += adjwgt[j];
        }
      }
  
      if (v != u) { 
        istart = xadj[u];
        iend   = xadj[u+1];
        for (j=istart; j<iend; j++) {
          k = cmap[adjncy[j]];
          for (kk=k&mask; htable[kk]!=-1 && cadjncy[htable[kk]]!=k; kk=((kk+1)&mask));
          if ((m = htable[kk]) == -1) {
            cadjncy[nedges] = k;
            cadjwgt[nedges] = adjwgt[j];
            htable[kk]      = nedges++;
          }
          else {
            cadjwgt[m] += adjwgt[j];
          }
        }
      }

      /* reset the htable -- reverse order (LIFO) is critical to prevent cadjncy[-1]
       * indexing due to a remove of an earlier entry */
      for (j=nedges-1; j>=0; j--) {
        k = cadjncy[j];
        for (kk=k&mask; cadjncy[htable[kk]]!=k; kk=((kk+1)&mask));
        htable[kk] = -1;  
      }

      /* remove the contracted vertex from the list */
      cadjncy[0] = cadjncy[--nedges];
      cadjwgt[0] = cadjwgt[nedges];
    }
    else {
      nedges = 0;
      istart = xadj[v];
      iend   = xadj[v+1];
      for (j=istart; j<iend; j++) {
        k = cmap[adjncy[j]];
        if ((m = dtable[k]) == -1) {
          cadjncy[nedges] = k;
          cadjwgt[nedges] = adjwgt[j];
          dtable[k] = nedges++;
        }
        else {
          cadjwgt[m] += adjwgt[j];
        }
      }

      if (v != u) { 
        istart = xadj[u];
        iend   = xadj[u+1];
        for (j=istart; j<iend; j++) {
          k = cmap[adjncy[j]];
          if ((m = dtable[k]) == -1) {
            cadjncy[nedges] = k;
            cadjwgt[nedges] = adjwgt[j];
            dtable[k] = nedges++;
          }
          else {
            cadjwgt[m] += adjwgt[j];
          }
        }

        /* Remove the contracted self-loop, when present */
        if ((j = dtable[cnvtxs]) != -1) {
          cadjncy[j]        = cadjncy[--nedges];
          cadjwgt[j]        = cadjwgt[nedges];
          dtable[cnvtxs] = -1;
        }
      }

      /* Zero out the dtable */
      for (j=0; j<nedges; j++)
        dtable[cadjncy[j]] = -1;  
    }


    /* Determine the median weight of the incident edges, which will be used
       to keep an edge (u, v) iff wgt(u, v) >= min(medianewgts[u], medianewgts[v]) */
    if (dropedges) {
      medianewgts[cnvtxs] = 8;  /* default for island nodes */ 
      if (nedges > 0) {
        for (j=0; j<nedges; j++) 
          keys[j] = (cadjwgt[j]<<8) + noise[cnvtxs] + noise[cadjncy[j]];
        isortd(nedges, keys);
        medianewgts[cnvtxs] = keys[gk_min(nedges-1, ((xadj[v+1]-xadj[v] + xadj[u+1]-xadj[u])>>1))];
      }
    }

    cadjncy         += nedges;
    cadjwgt         += nedges;
    cnedges         += nedges;
    cxadj[++cnvtxs]  = cnedges;
  }


  /* compact the adjacency structure of the coarser graph to keep only +ve edges */
  if (dropedges) { 
    droppedewgt = 0;

    cadjncy  = cgraph->adjncy;
    cadjwgt  = cgraph->adjwgt;

    cnedges = 0;
    for (u=0; u<cnvtxs; u++) {
      istart = cxadj[u];
      iend   = cxadj[u+1];
      for (j=istart; j<iend; j++) {
        v = cadjncy[j];
        if ((cadjwgt[j]<<8) + noise[u] + noise[v] >= gk_min(medianewgts[u], medianewgts[v])) {
          cadjncy[cnedges]   = cadjncy[j];
          cadjwgt[cnedges++] = cadjwgt[j];
        }
        else 
          droppedewgt += cadjwgt[j];
      }
      cxadj[u] = cnedges;
    }
    SHIFTCSR(j, cnvtxs, cxadj);

    cgraph->droppedewgt = droppedewgt;
  }

  cgraph->nedges = cnedges;

  for (j=0; j<ncon; j++) {
    cgraph->tvwgt[j]    = isum(cgraph->nvtxs, cgraph->vwgt+j, ncon);
    cgraph->invtvwgt[j] = 1.0/(cgraph->tvwgt[j] > 0 ? cgraph->tvwgt[j] : 1);
  }

  ReAdjustMemory(ctrl, graph, cgraph);

  WCOREPOP;
}


/*************************************************************************/
/*! Setup the various arrays for the coarse graph 
 */
/*************************************************************************/
graph_t *SetupCoarseGraph(graph_t *graph, idx_t cnvtxs, int dovsize)
{
  graph_t *cgraph;

  cgraph = CreateGraph();

  cgraph->nvtxs = cnvtxs;
  cgraph->ncon  = graph->ncon;

  cgraph->finer  = graph;
  graph->coarser = cgraph;

  /* Allocate memory for the coarser graph.
     NOTE: The +1 in the adjwgt/adjncy is to allow the optimization of self-loop
           detection by adding ahead of time the self-loop. That optimization
           requires a +1 adjncy/adjwgt array for the limit case where the 
           coarser graph is of the same size of the previous graph. */
  cgraph->xadj     = imalloc(cnvtxs+1);
  cgraph->adjncy   = imalloc(graph->nedges+1);
  cgraph->adjwgt   = imalloc(graph->nedges+1);
  cgraph->vwgt     = imalloc(cgraph->ncon*cnvtxs);
  cgraph->tvwgt    = imalloc(cgraph->ncon);
  cgraph->invtvwgt = rmalloc(cgraph->ncon);

  if (dovsize)
    cgraph->vsize = imalloc(cnvtxs);

  return cgraph;
}


/*************************************************************************/
/*! This function re-adjusts the amount of memory that was allocated if
    it will lead to significant savings 
 */
/*************************************************************************/
void ReAdjustMemory(ctrl_t *ctrl, graph_t *graph, graph_t *cgraph) 
{
  if (cgraph->nedges > 10000 && cgraph->nedges < 0.9*graph->nedges) {
    cgraph->adjncy = irealloc(cgraph->adjncy, cgraph->nedges);
    cgraph->adjwgt = irealloc(cgraph->adjwgt, cgraph->nedges);
  }
}
