/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * defs.h
 *
 * This file contains constant definitions
 *
 * Started 8/27/94
 * George
 *
 * $Id: defs.h 20398 2016-11-22 17:17:12Z karypis $
 *
 */

#ifndef _LIBMETIS_DEFS_H_
#define _LIBMETIS_DEFS_H_

#define METISTITLE              "METIS 5.2.1 Copyright 1998-22, Regents of the University of Minnesota\n"

#define HTLENGTH		((1<<13)-1)

#define INIT_MAXNAD             200     /* Initial number of maximum number of 
                                           adjacent domains. This number will be
                                           adjusted as required. */

#define UNMATCHED		-1

#define LARGENIPARTS		7	/* Number of random initial partitions */
#define SMALLNIPARTS		5	/* Number of random initial partitions */

#define COARSEN_FRACTION	0.85	/* Node reduction between successive coarsening levels */

#define COMPRESSION_FRACTION		0.85

#define MMDSWITCH		        120

#define OMETIS_DEFAULT_UFACTOR          200

#endif
