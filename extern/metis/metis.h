/*!
\file metis.h
\brief This file contains function prototypes and constant definitions for METIS
 *
\author George
\date   Started 8/9/02
\version\verbatim $Id$\endverbatim
*/

#ifndef HIGHS_METIS_ORDERING_H
#define HIGHS_METIS_ORDERING_H

/****************************************************************************
* A set of defines that can be modified by the user
*****************************************************************************/

/*--------------------------------------------------------------------------
 Specifies the width of the elementary data type that will hold information
 about vertices and their adjacency lists.

 Possible values:
   32 : Use 32 bit signed integers
   64 : Use 64 bit signed integers

 A width of 64 should be specified if the number of vertices or the total
 number of edges in the graph exceed the limits of a 32 bit signed integer
 i.e., 2^31-1.
 Proper use of 64 bit integers requires that the c99 standard datatypes
 int32_t and int64_t are supported by the compiler.
 GCC does provides these definitions in stdint.h, but it may require some
 modifications on other architectures.
--------------------------------------------------------------------------*/
#include "HighsMetis.h"

/*--------------------------------------------------------------------------
 Specifies the data type that will hold floating-point style information.

 Possible values:
   32 : single precision floating point (float)
   64 : double precision floating point (double)
--------------------------------------------------------------------------*/
#define REALTYPEWIDTH 32

/****************************************************************************
* In principle, nothing needs to be changed beyond this point, unless the
* int32_t and int64_t cannot be found in the normal places.
*****************************************************************************/

/* Uniform definitions for various compilers */
#if defined(_MSC_VER)
  #define COMPILER_MSC
#endif
#if defined(__ICC)
  #define COMPILER_ICC
#endif
#if defined(__GNUC__)
  #define COMPILER_GCC
#endif

/* Include c99 int definitions and needed constants. When building the library,
 * these are already defined by GKlib; hence the test for _GKLIB_H_ */
#ifndef _GKLIB_H_
#ifdef __cplusplus
#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#endif
#include <stdint.h>
#include <inttypes.h>
#endif


/*------------------------------------------------------------------------
* Setup the basic datatypes
*-------------------------------------------------------------------------*/
#if IDXTYPEWIDTH == 32
  #define IDX_MAX   INT32_MAX
  #define IDX_MIN   INT32_MIN

  #define SCIDX  SCNd32
  #define PRIDX  PRId32

  #define strtoidx      strtol
  #define iabs          abs
#elif IDXTYPEWIDTH == 64
  #define IDX_MAX   INT64_MAX
  #define IDX_MIN   INT64_MIN

  #define SCIDX  SCNd64
  #define PRIDX  PRId64

#ifdef COMPILER_MSC
  #define strtoidx      _strtoi64
#else
  #define strtoidx      strtoll
#endif
  #define iabs          labs
#else
  #error "Incorrect user-supplied value of IDXTYPEWIDTH"
#endif


#if REALTYPEWIDTH == 32
  typedef float real_t;

  #define SCREAL         "f"
  #define PRREAL         "f"
  #define REAL_MAX       FLT_MAX
  #define REAL_MIN       FLT_MIN
  #define REAL_EPSILON   FLT_EPSILON

  #define rabs          fabsf
  #define REALEQ(x,y) ((rabs((x)-(y)) <= FLT_EPSILON))

#ifdef COMPILER_MSC
  #define strtoreal     (float)strtod
#else
  #define strtoreal     strtof
#endif
#elif REALTYPEWIDTH == 64
  typedef double real_t;

  #define SCREAL         "lf"
  #define PRREAL         "lf"
  #define REAL_MAX       DBL_MAX
  #define REAL_MIN       DBL_MIN
  #define REAL_EPSILON   DBL_EPSILON

  #define rabs          fabs
  #define REALEQ(x,y) ((rabs((x)-(y)) <= DBL_EPSILON))

  #define strtoreal     strtod
#else
  #error "Incorrect user-supplied value for REALTYPEWIDTH"
#endif


/*------------------------------------------------------------------------
* Constant definitions 
*-------------------------------------------------------------------------*/
/* Metis's version number */
#define METIS_VER_MAJOR         5
#define METIS_VER_MINOR         2
#define METIS_VER_SUBMINOR      1


/*------------------------------------------------------------------------
* Enum type definitions 
*-------------------------------------------------------------------------*/


/*! Operation type codes */
typedef enum {
  METIS_OP_OMETIS
} moptype_et;


/*! Coarsening Schemes */
typedef enum {
  METIS_CTYPE_RM,
  METIS_CTYPE_SHEM
} mctype_et;

/*! Initial partitioning schemes */
typedef enum {
  METIS_IPTYPE_EDGE,
  METIS_IPTYPE_NODE,
} miptype_et;


/*! Refinement schemes */
typedef enum {
  METIS_RTYPE_SEP2SIDED,
  METIS_RTYPE_SEP1SIDED
} mrtype_et;


/*! Debug Levels */
typedef enum {
  METIS_DBG_INFO       = 1,       /*!< Shows various diagnostic messages */
  METIS_DBG_COARSEN    = 2,	      /*!< Show the coarsening progress */
  METIS_DBG_REFINE     = 4,	      /*!< Show the refinement progress */
  METIS_DBG_IPART      = 8, 	    /*!< Show info on initial partitioning */
  METIS_DBG_MOVEINFO   = 16, 	    /*!< Show info on vertex moves during refinement */
  METIS_DBG_SEPINFO    = 32, 	    /*!< Show info on vertex moves during sep refinement */
} mdbglvl_et;


#endif  /* _METIS_H_ */
