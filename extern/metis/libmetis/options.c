/**
  \file
  \brief This file contains various routines for dealing with options and ctrl_t.

  \date   Started 5/12/2011
  \author George  
  \author Copyright 1997-2011, Regents of the University of Minnesota 
  \version\verbatim $Id: options.c 17717 2014-10-03 19:09:31Z dominique $ \endverbatim
  */

#include "metislib.h"


/*************************************************************************/
/*! This function creates and sets the run parameters (ctrl_t) */
/*************************************************************************/
ctrl_t *SetupCtrl(moptype_et optype, idx_t *options, idx_t ncon, idx_t nparts, 
            real_t *tpwgts, real_t *ubvec)
{
  idx_t i, j;
  ctrl_t *ctrl;

  ctrl = (ctrl_t *)malloc(sizeof(ctrl_t));
  
  memset((void *)ctrl, 0, sizeof(ctrl_t));

  switch (optype) {
    case METIS_OP_OMETIS:
      ctrl->rtype    = GETOPTION(options, METIS_OPTION_RTYPE,    METIS_RTYPE_SEP1SIDED);
      ctrl->iptype   = GETOPTION(options, METIS_OPTION_IPTYPE,   METIS_IPTYPE_EDGE);
      ctrl->nseps    = GETOPTION(options, METIS_OPTION_NSEPS,    1);
      ctrl->niter    = GETOPTION(options, METIS_OPTION_NITER,    10);
      ctrl->ufactor  = GETOPTION(options, METIS_OPTION_UFACTOR,  OMETIS_DEFAULT_UFACTOR);
      ctrl->compress = GETOPTION(options, METIS_OPTION_COMPRESS, 1);
      ctrl->ccorder  = GETOPTION(options, METIS_OPTION_CCORDER,  0);
      ctrl->pfactor  = 0.1*GETOPTION(options, METIS_OPTION_PFACTOR,  0);

      ctrl->CoarsenTo = 100;
      break;

    default:
      gk_errexit("Unknown optype of %d\n", optype);
  }

  /* common options */
  ctrl->ctype     = GETOPTION(options, METIS_OPTION_CTYPE, METIS_CTYPE_SHEM);
  ctrl->no2hop    = GETOPTION(options, METIS_OPTION_NO2HOP, 0);
  ctrl->seed      = GETOPTION(options, METIS_OPTION_SEED, -1);
  ctrl->dbglvl    = GETOPTION(options, METIS_OPTION_DBGLVL, 0);
  ctrl->dropedges = GETOPTION(options, METIS_OPTION_DROPEDGES, 0);

  /* set non-option information */
  ctrl->optype  = optype;
  ctrl->ncon    = ncon;
  ctrl->nparts  = nparts;
  ctrl->maxvwgt = ismalloc(ncon, 0);

  ctrl->rng_state = (ctrl->seed == -1 ? 4321 : ctrl->seed);


  /* setup the ubfactors */
  ctrl->ubfactors = rsmalloc(ctrl->ncon, I2RUBFACTOR(ctrl->ufactor));
  if (ubvec)
    rcopy(ctrl->ncon, ubvec, ctrl->ubfactors);
  for (i=0; i<ctrl->ncon; i++)
    ctrl->ubfactors[i] += 0.0000499;

  /* Allocate memory for balance multipliers. 
     Note that for PMETIS/OMETIS routines the memory allocated is more 
     than required as balance multipliers for 2 parts is sufficient. */
  ctrl->pijbm = rmalloc(nparts*ncon);

  IFSET(ctrl->dbglvl, METIS_DBG_INFO, PrintCtrl(ctrl));

  if (!CheckParams(ctrl)) {
    FreeCtrl(&ctrl);
    return NULL;
  }
  else {
    return ctrl;
  }
}

/*************************************************************************/
/*! Computes the per-partition/constraint balance multipliers */
/*************************************************************************/
void Setup2WayBalMultipliers(ctrl_t *ctrl, graph_t *graph, real_t *tpwgts)
{
  idx_t i, j;

  for (i=0; i<2; i++) {
    for (j=0; j<graph->ncon; j++)
      ctrl->pijbm[i*graph->ncon+j] = graph->invtvwgt[j]/tpwgts[i*graph->ncon+j];
  }
}


/*************************************************************************/
/*! This function prints the various control fields */
/*************************************************************************/
void PrintCtrl(ctrl_t *ctrl)
{
  idx_t i, j, modnum;

  printf(" Runtime parameters:\n");

  printf("   Coarsening type: ");
  switch (ctrl->ctype) {
    case METIS_CTYPE_RM:
      printf("METIS_CTYPE_RM\n");
      break;
    case METIS_CTYPE_SHEM:
      printf("METIS_CTYPE_SHEM\n");
      break;
    default:
      printf("Unknown!\n");
  }

  printf("   Initial partitioning type: ");
  switch (ctrl->iptype) {
    case METIS_IPTYPE_EDGE:
      printf("METIS_IPTYPE_EDGE\n");
      break;
    case METIS_IPTYPE_NODE:
      printf("METIS_IPTYPE_NODE\n");
      break;
    default:
      printf("Unknown!\n");
  }

  printf("   Refinement type: ");
  switch (ctrl->rtype) {
    case METIS_RTYPE_SEP2SIDED:
      printf("METIS_RTYPE_SEP2SIDED\n");
      break;
    case METIS_RTYPE_SEP1SIDED:
      printf("METIS_RTYPE_SEP1SIDED\n");
      break;
    default:
      printf("Unknown!\n");
  }

  printf("   Perform a 2-hop matching: %s\n", (ctrl->no2hop ? "No" : "Yes"));

  printf("   Drop edges: %s\n", (ctrl->dropedges ? "Yes" : "No"));

  printf("   Number of balancing constraints: %"PRIDX"\n", ctrl->ncon);
  printf("   Number of refinement iterations: %"PRIDX"\n", ctrl->niter);
  printf("   Random number seed: %"PRIDX"\n", ctrl->seed);

  if (ctrl->optype == METIS_OP_OMETIS) {
    printf("   Number of separators: %"PRIDX"\n", ctrl->nseps);
    printf("   Compress graph prior to ordering: %s\n", (ctrl->compress ? "Yes" : "No"));
    printf("   Detect & order connected components separately: %s\n", (ctrl->ccorder ? "Yes" : "No"));
    printf("   Prunning factor for high degree vertices: %"PRREAL"\n", ctrl->pfactor);
  }

  printf("   Allowed maximum load imbalance: ");
  for (i=0; i<ctrl->ncon; i++) 
    printf("%.3"PRREAL" ", ctrl->ubfactors[i]);
  printf("\n");

  printf("\n");
}


/*************************************************************************/
/*! This function checks the validity of user-supplied parameters */
/*************************************************************************/
int CheckParams(ctrl_t *ctrl)
{
  idx_t i, j;
  real_t sum;
  mdbglvl_et  dbglvl=METIS_DBG_INFO;

  switch (ctrl->optype) {
    case METIS_OP_OMETIS:
      if (ctrl->ctype != METIS_CTYPE_RM && ctrl->ctype != METIS_CTYPE_SHEM) {
        IFSET(dbglvl, METIS_DBG_INFO, printf("Input Error: Incorrect coarsening scheme.\n"));
        return 0;
      }
      if (ctrl->iptype != METIS_IPTYPE_EDGE && ctrl->iptype != METIS_IPTYPE_NODE) {
        IFSET(dbglvl, METIS_DBG_INFO, printf("Input Error: Incorrect initial partitioning scheme.\n"));
        return 0;
      }
      if (ctrl->rtype != METIS_RTYPE_SEP1SIDED && ctrl->rtype != METIS_RTYPE_SEP2SIDED) {
        IFSET(dbglvl, METIS_DBG_INFO, printf("Input Error: Incorrect refinement scheme.\n"));
        return 0;
      }
      if (ctrl->nseps <= 0) {
        IFSET(dbglvl, METIS_DBG_INFO, printf("Input Error: Incorrect nseps.\n"));
        return 0;
      }
      if (ctrl->niter <= 0) {
        IFSET(dbglvl, METIS_DBG_INFO, printf("Input Error: Incorrect niter.\n"));
        return 0;
      }
      if (ctrl->ufactor <= 0) {
        IFSET(dbglvl, METIS_DBG_INFO, printf("Input Error: Incorrect ufactor.\n"));
        return 0;
      }
      if (ctrl->nparts != 3) {
        IFSET(dbglvl, METIS_DBG_INFO, printf("Input Error: Incorrect nparts.\n"));
        return 0;
      }
      if (ctrl->ncon != 1) {
        IFSET(dbglvl, METIS_DBG_INFO, printf("Input Error: Incorrect ncon.\n"));
        return 0;
      }
      if (ctrl->compress != 0 && ctrl->compress != 1) {
        IFSET(dbglvl, METIS_DBG_INFO, printf("Input Error: Incorrect compress.\n"));
        return 0;
      }
      if (ctrl->ccorder != 0 && ctrl->ccorder != 1) {
        IFSET(dbglvl, METIS_DBG_INFO, printf("Input Error: Incorrect ccorder.\n"));
        return 0;
      }
      if (ctrl->pfactor < 0.0 ) {
        IFSET(dbglvl, METIS_DBG_INFO, printf("Input Error: Incorrect pfactor.\n"));
        return 0;
      }

      for (i=0; i<ctrl->ncon; i++) {
        if (ctrl->ubfactors[i] <= 1.0) {
          IFSET(dbglvl, METIS_DBG_INFO, 
              printf("Input Error: Incorrect ubfactor for constraint %"PRIDX".\n", i));
          return 0;
        }
      }

      break;

    default:
      IFSET(dbglvl, METIS_DBG_INFO, printf("Input Error: Incorrect optype\n"));
      return 0;
  }

  return 1;
}

  
/*************************************************************************/
/*! This function frees the memory associated with a ctrl_t */
/*************************************************************************/
void FreeCtrl(ctrl_t **r_ctrl)
{
  ctrl_t *ctrl = *r_ctrl;

  FreeWorkSpace(ctrl);

  gk_free((void **)&ctrl->tpwgts);
  gk_free((void**)&ctrl->pijbm);
  gk_free((void**)&ctrl->ubfactors);
  gk_free((void**)&ctrl->maxvwgt);
  gk_free((void**)&ctrl);

  *r_ctrl = NULL;
}


