/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * util.c
 *
 * This function contains various utility routines
 *
 * Started 9/28/95
 * George
 *
 * $Id: util.c 10495 2011-07-06 16:04:45Z karypis $
 */

#include "metislib.h"

/*************************************************************************/
/*! Returns the highest weight index of x[i]*y[i] 
 */
/*************************************************************************/
idx_t iargmax_nrm(size_t n, idx_t *x, real_t *y)
{
  idx_t i, max=0;
      
  for (i=1; i<n; i++)
     max = (x[i]*y[i] > x[max]*y[max] ? i : max);
                
  return max;
}

/*************************************************************************/
/*! These functions return the index of the second largest elements in the
    vector formed by x.y where '.' is element-wise multiplication */
/*************************************************************************/
idx_t iargmax2_nrm(size_t n, idx_t *x, real_t *y)
{
  size_t i, max1, max2;

  if (x[0]*y[0] > x[1]*y[1]) {
    max1 = 0;
    max2 = 1;
  }
  else {
    max1 = 1;
    max2 = 0;
  }

  for (i=2; i<n; i++) {
    if (x[i]*y[i] > x[max1]*y[max1]) {
      max2 = max1;
      max1 = i;
    }
    else if (x[i]*y[i] > x[max2]*y[max2])
      max2 = i;
  }

  return max2;
}


