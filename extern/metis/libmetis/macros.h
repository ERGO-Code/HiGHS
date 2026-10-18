/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * macros.h
 *
 * This file contains macros used in multilevel
 *
 * Started 9/25/94
 * George
 *
 * $Id: macros.h 10060 2011-06-02 18:56:30Z karypis $
 *
 */

#ifndef _LIBMETIS_MACROS_H_
#define _LIBMETIS_MACROS_H_

#define SWAP gk_SWAP

/* gets the appropriate option value */
#define GETOPTION(options, idx, defval) \
            ((options) == NULL || (options)[idx] == -1 ? defval : (options)[idx]) 

/* converts a user provided ufactor into a real ubfactor */
#define I2RUBFACTOR(ufactor) (1.0+0.001*(ufactor))

/* set/reset the current workspace core */
#define WCOREPUSH    wspacepush(ctrl)
#define WCOREPOP     wspacepop(ctrl)



/*************************************************************************
* These macros insert and remove nodes from a Direct Access list 
**************************************************************************/
#define ListInsert(n, lind, lptr, i) \
   do { \
     lind[n] = i; \
     lptr[i] = (n)++;\
   } while(0) 

#define ListDelete(n, lind, lptr, i) \
   do { \
     lind[lptr[i]] = lind[--(n)]; \
     lptr[lind[n]] = lptr[i]; \
     lptr[i] = -1; \
   } while(0) 


/*************************************************************************
* These macros insert and remove nodes from the boundary list
**************************************************************************/
#define BNDInsert(nbnd, bndind, bndptr, vtx) \
  ListInsert(nbnd, bndind, bndptr, vtx)

#define BNDDelete(nbnd, bndind, bndptr, vtx) \
  ListDelete(nbnd, bndind, bndptr, vtx)

#endif
