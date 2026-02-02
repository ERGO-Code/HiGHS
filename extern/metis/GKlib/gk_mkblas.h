/*!
\file  gk_mkblas.h
\brief Templates for BLAS-like routines

\date   Started 3/28/07
\author George
\version\verbatim $Id: gk_mkblas.h 16304 2014-02-25 14:27:19Z karypis $ \endverbatim
*/

#ifndef _GK_MKBLAS_H_
#define _GK_MKBLAS_H_


#define GK_MKBLAS(PRFX, TYPE, OUTTYPE) \
/*************************************************************************/\
/*! The macro for gk_?sum()-class of routines */\
/**************************************************************************/\
OUTTYPE PRFX ## sum(size_t n, TYPE *x, size_t incx)\
{\
  size_t i;\
  OUTTYPE sum = 0;\
\
  for (i=0; i<n; i++, x+=incx)\
    sum += (*x);\
\
  return sum;\
}\
\
\
/*************************************************************************/\
/*! The macro for gk_?axpy()-class of routines */\
/**************************************************************************/\
TYPE *PRFX ## axpy(size_t n, TYPE alpha, TYPE *x, size_t incx, TYPE *y, size_t incy)\
{\
  size_t i;\
  TYPE *y_in = y;\
\
  for (i=0; i<n; i++, x+=incx, y+=incy)\
    *y += alpha*(*x);\
\
  return y_in;\
}\



#define GK_MKBLAS_PROTO(PRFX, TYPE, OUTTYPE) \
  OUTTYPE  PRFX ## sum(size_t n, TYPE *x, size_t incx);\
  TYPE    *PRFX ## axpy(size_t n, TYPE alpha, TYPE *x, size_t incx, TYPE *y, size_t incy);\


#endif
