/*!
\file  gklib.c
\brief Various helper routines generated using GKlib's templates

\date   Started 4/12/2007
\author George  
\author Copyright 1997-2009, Regents of the University of Minnesota 
\version\verbatim $Id: gklib.c 10395 2011-06-23 23:28:06Z karypis $ \endverbatim
*/


#include "metislib.h"


/*************************************************************************/
/*! BLAS routines */
/*************************************************************************/
GK_MKBLAS(i,  idx_t,  idx_t)

/*************************************************************************/
/*! Memory allocation routines */
/*************************************************************************/
GK_MKALLOC(i,    idx_t)
GK_MKALLOC(r,    real_t)
GK_MKALLOC(ikv,  ikv_t)
GK_MKALLOC(rkv,  rkv_t)

/*************************************************************************/
/*! Priority queues routines */
/*************************************************************************/
#define key_gt(a, b) ((a) > (b))
GK_MKPQUEUE(rpq, rpq_t, rkv_t, real_t, idx_t, rkvmalloc, REAL_MAX, key_gt)
#undef key_gt

/*************************************************************************/
/*! Random number generation routines */
/*************************************************************************/
GK_MKRANDOM(i, idx_t, idx_t)

/*************************************************************************/
/*! Sorting routines */
/*************************************************************************/
void isortd(size_t n, idx_t *base)
{
#define i_gt(a, b) ((*a) > (*b))
  GK_MKQSORT(idx_t, base, n, i_gt);
#undef i_gt
}

void ikvsorti(size_t n, ikv_t *base)
{
#define ikey_lt(a, b) ((a)->key < (b)->key)
  GK_MKQSORT(ikv_t, base, n, ikey_lt);
#undef ikey_lt
}
