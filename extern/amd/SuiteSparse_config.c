//------------------------------------------------------------------------------
// SuiteSparse_config/SuiteSparse_config.c: common utilites for SuiteSparse
//------------------------------------------------------------------------------

// SuiteSparse_config, Copyright (c) 2012-2023, Timothy A. Davis.
// All Rights Reserved.
// SPDX-License-Identifier: BSD-3-clause

//------------------------------------------------------------------------------

/* SuiteSparse configuration : memory manager and printf functions.
 */

#include "SuiteSparse_config.h"

void *SuiteSparse_malloc    /* pointer to allocated block of memory */
(
    size_t nitems,          /* number of items to malloc */
    size_t size_of_item     /* sizeof each item */
)
{
    void *p ;
    size_t size ;
    if (nitems < 1) nitems = 1 ;
    if (size_of_item < 1) size_of_item = 1 ;
    size = nitems * size_of_item  ;

    if (size != ((double) nitems) * size_of_item)
    {
        /* size_t overflow */
        p = NULL ;
    }
    else
    {
        p = (void *) malloc (size) ;
    }
    return (p) ;
}

/* -------------------------------------------------------------------------- */
/* SuiteSparse_free: free wrapper */
/* -------------------------------------------------------------------------- */

void *SuiteSparse_free      /* always returns NULL */
(
    void *p                 /* block to free */
)
{
    if (p)
    {
        free(p) ;
    }
    return (NULL) ;
}
