/*!
\file  gk_mkmemory.h
\brief Templates for memory allocation routines

\date   Started 3/29/07
\author George
\version\verbatim $Id: gk_mkmemory.h 10711 2011-08-31 22:23:04Z karypis $ \endverbatim
*/

#ifndef _GK_MKMEMORY_H_
#define _GK_MKMEMORY_H_


#define GK_MKALLOC(PRFX, TYPE)\
/*************************************************************************/\
/*! The macro for gk_?malloc()-class of routines */\
/**************************************************************************/\
TYPE *PRFX ## malloc(size_t n)\
{\
  return (TYPE *)malloc(sizeof(TYPE)*n);\
}\
\
\
/*************************************************************************/\
/*! The macro for gk_?realloc()-class of routines */\
/**************************************************************************/\
TYPE *PRFX ## realloc(TYPE *ptr, size_t n)\
{\
  return (TYPE *)realloc((void *)ptr, sizeof(TYPE)*n);\
}\
\
\
/*************************************************************************/\
/*! The macro for gk_?smalloc()-class of routines */\
/**************************************************************************/\
TYPE *PRFX ## smalloc(size_t n, TYPE ival)\
{\
  TYPE *ptr;\
\
  ptr = (TYPE *)malloc(sizeof(TYPE)*n);\
  if (ptr == NULL) \
    return NULL; \
\
  return PRFX ## set(n, ival, ptr); \
}\
\
\
/*************************************************************************/\
/*! The macro for gk_?set()-class of routines */\
/*************************************************************************/\
TYPE *PRFX ## set(size_t n, TYPE val, TYPE *x)\
{\
  size_t i;\
\
  for (i=0; i<n; i++)\
    x[i] = val;\
\
  return x;\
}\
\
\
/*************************************************************************/\
/*! The macro for gk_?set()-class of routines */\
/*************************************************************************/\
TYPE *PRFX ## copy(size_t n, TYPE *a, TYPE *b)\
{\
  return (TYPE *)memmove((void *)b, (void *)a, sizeof(TYPE)*n);\
}\


#define GK_MKALLOC_PROTO(PRFX, TYPE)\
  TYPE  *PRFX ## malloc(size_t n);\
  TYPE  *PRFX ## realloc(TYPE *ptr, size_t n);\
  TYPE  *PRFX ## smalloc(size_t n, TYPE ival);\
  TYPE  *PRFX ## set(size_t n, TYPE val, TYPE *x);\
  TYPE  *PRFX ## copy(size_t n, TYPE *a, TYPE *b);\



#endif
