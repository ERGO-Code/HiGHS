/*!
\file 
\brief Functions dealing with creating and allocating mcores

\date Started 5/30/11
\author George
\author Copyright 1997-2011, Regents of the University of Minnesota 
\version $Id: mcore.c 13953 2013-03-30 16:20:07Z karypis $
*/

#include "GKlib.h"


/*************************************************************************/
/*! This function creates an mcore 
 */
/*************************************************************************/
gk_mcore_t *gk_mcoreCreate(size_t coresize)
{
  gk_mcore_t *mcore;

  mcore = (gk_mcore_t *)malloc(sizeof(gk_mcore_t));
  memset(mcore, 0, sizeof(gk_mcore_t));

  mcore->coresize = coresize;
  mcore->corecpos = 0;

  mcore->core = (coresize == 0 ? NULL : malloc(mcore->coresize));

  /* allocate the memory for keeping track of malloc ops */
  mcore->nmops = 2048;
  mcore->cmop  = 0;
  mcore->mops  = (gk_mop_t *)malloc(mcore->nmops*sizeof(gk_mop_t));

  return mcore;
}

/*************************************************************************/
/*! This function destroys an mcore.
 */
/*************************************************************************/
void gk_mcoreDestroy(gk_mcore_t **r_mcore, int showstats)
{
  gk_mcore_t *mcore = *r_mcore;

  if (mcore == NULL)
    return;

  gk_free((void **)&mcore->core);
  gk_free((void**)&mcore->mops);
  gk_free((void**)&mcore);

  *r_mcore = NULL;
}


/*************************************************************************/
/*! This function allocate space from the core/heap 
 */
/*************************************************************************/
void *gk_mcoreMalloc(gk_mcore_t *mcore, size_t nbytes)
{
  void *ptr;

  /* pad to make pointers 8-byte aligned */
  nbytes += (nbytes%8 == 0 ? 0 : 8 - nbytes%8);

  if (mcore->corecpos + nbytes < mcore->coresize) {
    /* service this request from the core */
    ptr = ((char *)mcore->core)+mcore->corecpos;
    mcore->corecpos += nbytes;

    gk_mcoreAdd(mcore, GK_MOPT_CORE, nbytes, ptr);
  }
  else {
    /* service this request from the heap */
    ptr = malloc(nbytes);

    gk_mcoreAdd(mcore, GK_MOPT_HEAP, nbytes, ptr);
  }

  /*
  printf("MCMALLOC: %zu %d %8zu\n", mcore->cmop-1, 
      mcore->mops[mcore->cmop-1].type, mcore->mops[mcore->cmop-1].nbytes);
  */

  return ptr;
}


/*************************************************************************/
/*! This function sets a marker in the stack of malloc ops to be used
    subsequently for freeing purposes 
 */
/*************************************************************************/
void gk_mcorePush(gk_mcore_t *mcore)
{
  gk_mcoreAdd(mcore, GK_MOPT_MARK, 0, NULL);
  /* printf("MCPPUSH:   %zu\n", mcore->cmop-1); */
}


/*************************************************************************/
/*! This function frees all mops since the last push 
 */
/*************************************************************************/
void gk_mcorePop(gk_mcore_t *mcore)
{
  while (mcore->cmop > 0) {
    mcore->cmop--;
    switch (mcore->mops[mcore->cmop].type) {
      case GK_MOPT_MARK: /* push marker */
        goto DONE;
        break; 

      case GK_MOPT_CORE: /* core free */
        if (mcore->corecpos < mcore->mops[mcore->cmop].nbytes)
          gk_errexit("Internal Error: wspace's core is about to be over-freed [%zu, %zu, %zd]\n",
              mcore->coresize, mcore->corecpos, mcore->mops[mcore->cmop].nbytes);

        mcore->corecpos    -= mcore->mops[mcore->cmop].nbytes;
        break;

      case GK_MOPT_HEAP: /* heap free */
        gk_free((void **)&mcore->mops[mcore->cmop].ptr);
        break;

      default:
        gk_errexit("Unknown mop type of %d\n", mcore->mops[mcore->cmop].type);
    }
  }

DONE:
  ;
  /*printf("MCPPOP:    %zu\n", mcore->cmop); */
}


/*************************************************************************/
/*! Adds a memory allocation at the end of the list.
 */
/*************************************************************************/
void gk_mcoreAdd(gk_mcore_t *mcore, int type, size_t nbytes, void *ptr)
{
  if (mcore->cmop == mcore->nmops) {
    mcore->nmops *= 2;
    mcore->mops = realloc(mcore->mops, mcore->nmops*sizeof(gk_mop_t));
    if (mcore->mops == NULL) 
      gk_errexit("***Memory allocation for gkmcore failed.\n");
  }

  mcore->mops[mcore->cmop].type   = type;
  mcore->mops[mcore->cmop].nbytes = nbytes;
  mcore->mops[mcore->cmop].ptr    = ptr;
  mcore->cmop++;
}


/*************************************************************************/
/*! This function deletes the mop associated with the supplied pointer.
    The mop has to be a heap allocation, otherwise it fails violently.
 */
/*************************************************************************/
void gk_mcoreDel(gk_mcore_t *mcore, void *ptr)
{
  int i;

  for (i=mcore->cmop-1; i>=0; i--) {
    if (mcore->mops[i].type == GK_MOPT_MARK)
      gk_errexit("Could not find pointer %p in mcore\n", ptr);

    if (mcore->mops[i].ptr == ptr) {
      if (mcore->mops[i].type != GK_MOPT_HEAP)
        gk_errexit("Trying to delete a non-HEAP mop.\n");

      mcore->mops[i] = mcore->mops[--mcore->cmop];
      return;
    }
  }

  gk_errexit("mcoreDel should never have been here!\n");
}

