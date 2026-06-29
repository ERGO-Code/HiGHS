/*
 * mutil.c 
 *
 * This file contains various utility functions for the MOC portion of the
 * code
 *
 * Started 2/15/98
 * George
 *
 * $Id: mcutil.c 13901 2013-03-24 16:17:03Z karypis $
 *
 */

#include "metislib.h"


/*************************************************************************/
/*! This function returns true if \forall i, x[i] <= z[i]. */
/**************************************************************************/
int ivecle(idx_t n, idx_t *x, idx_t *z)
{
  for (n--; n>=0; n--) {
    if (x[n] > z[n]) 
      return 0;
  }

  return  1;
}

/*************************************************************************/
/*! This function returns true if \forall i, a*x[i]+y[i] <= z[i]. */
/**************************************************************************/
int ivecaxpylez(idx_t n, idx_t a, idx_t *x, idx_t *y, idx_t *z)
{
  for (n--; n>=0; n--) {
    if (a*x[n]+y[n] > z[n]) 
      return 0;
  }

  return  1;
}

/*************************************************************************/
/*! This function checks if v+u2 provides a better balance in the weight 
     vector that v+u1 */
/*************************************************************************/
int BetterVBalance(idx_t ncon, real_t *invtvwgt, idx_t *v_vwgt, idx_t *u1_vwgt, 
        idx_t *u2_vwgt)
{
  idx_t i;
  real_t sum1=0.0, sum2=0.0, diff1=0.0, diff2=0.0;

  for (i=0; i<ncon; i++) {
    sum1 += (v_vwgt[i]+u1_vwgt[i])*invtvwgt[i];
    sum2 += (v_vwgt[i]+u2_vwgt[i])*invtvwgt[i];
  }
  sum1 = sum1/ncon;
  sum2 = sum2/ncon;

  for (i=0; i<ncon; i++) {
    diff1 += rabs(sum1 - (v_vwgt[i]+u1_vwgt[i])*invtvwgt[i]);
    diff2 += rabs(sum2 - (v_vwgt[i]+u2_vwgt[i])*invtvwgt[i]);
  }

  return (diff1 - diff2 >= 0);
}

/*************************************************************************/
/*! This function takes two ubfactor-centered load imbalance vectors x & y, 
    and returns true if y is better balanced than x. */
/*************************************************************************/ 
int BetterBalance2Way(idx_t n, real_t *x, real_t *y)
{
  real_t nrm1=0.0, nrm2=0.0;

  for (--n; n>=0; n--) {
    if (x[n] > 0) nrm1 += x[n]*x[n];
    if (y[n] > 0) nrm2 += y[n]*y[n];
  }
  return nrm2 < nrm1;
}

/*************************************************************************/
/*! Computes the maximum load imbalance of a partitioning solution over 
    all the constraints. */
/**************************************************************************/ 
real_t ComputeLoadImbalance(graph_t *graph, idx_t nparts, real_t *pijbm)
{
  idx_t i, j, ncon, *pwgts;
  real_t max, cur;

  ncon  = graph->ncon;
  pwgts = graph->pwgts;

  max = 1.0;
  for (i=0; i<ncon; i++) {
    for (j=0; j<nparts; j++) {
      cur = pwgts[j*ncon+i]*pijbm[j*ncon+i];
      if (cur > max)
        max = cur;
    }
  }

  return max;
}

/*************************************************************************/
/*! Computes the maximum load imbalance difference of a partitioning 
    solution over all the constraints. 
    The difference is defined with respect to the allowed maximum 
    unbalance for the respective constraint. 
 */
/**************************************************************************/ 
real_t ComputeLoadImbalanceDiff(graph_t *graph, idx_t nparts, real_t *pijbm,
           real_t *ubvec)
{
  idx_t i, j, ncon, *pwgts;
  real_t max, cur;

  ncon  = graph->ncon;
  pwgts = graph->pwgts;

  max = -1.0;
  for (i=0; i<ncon; i++) {
    for (j=0; j<nparts; j++) {
      cur = pwgts[j*ncon+i]*pijbm[j*ncon+i] - ubvec[i];
      if (cur > max)
        max = cur;
    }
  }

  return max;
}

/*************************************************************************/
/*! Computes the difference between load imbalance of each constraint across 
    the partitions minus the desired upper bound on the load imabalnce.
    It also returns the maximum load imbalance across the partitions &
    constraints. */
/**************************************************************************/ 
real_t ComputeLoadImbalanceDiffVec(graph_t *graph, idx_t nparts, real_t *pijbm, 
         real_t *ubfactors, real_t *diffvec)
{
  idx_t i, j, ncon, *pwgts;
  real_t cur, max;

  ncon  = graph->ncon;
  pwgts = graph->pwgts;

  for (max=-1.0, i=0; i<ncon; i++) {
    diffvec[i] = pwgts[i]*pijbm[i] - ubfactors[i];
    for (j=1; j<nparts; j++) {
      cur = pwgts[j*ncon+i]*pijbm[j*ncon+i] - ubfactors[i];
      if (cur > diffvec[i])
        diffvec[i] = cur;
    }
    if (max < diffvec[i])
      max = diffvec[i];
  }

  return max;
}
