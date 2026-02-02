//------------------------------------------------------------------------------
// AMD/Source/amd_control: print control parameters for AMD
//------------------------------------------------------------------------------

// AMD, Copyright (c) 1996-2022, Timothy A. Davis, Patrick R. Amestoy, and
// Iain S. Duff.  All Rights Reserved.
// SPDX-License-Identifier: BSD-3-clause

//------------------------------------------------------------------------------

/* User-callable.  Prints the control parameters for AMD.  See amd.h
 * for details.  If the Control array is not present, the defaults are
 * printed instead.
 */

#include "amd_internal.h"

void amd_control
(
    double Control [ ]
)
{
    double alpha ;
    amd_int aggressive ;

    if (Control != (double *) NULL)
    {
	alpha = Control [AMD_DENSE] ;
	aggressive = Control [AMD_AGGRESSIVE] != 0 ;
    }
    else
    {
	alpha = AMD_DEFAULT_DENSE ;
	aggressive = AMD_DEFAULT_AGGRESSIVE ;
    }

    SUITESPARSE_PRINTF ((
        "\nAMD version %d.%d.%d, %s: approximate minimum degree ordering\n"
	"    dense row parameter: %g\n", AMD_MAIN_VERSION, AMD_SUB_VERSION,
	AMD_SUBSUB_VERSION, AMD_DATE, alpha)) ;

    if (alpha < 0)
    {
	SUITESPARSE_PRINTF (("    no rows treated as dense\n")) ;
    }
    else
    {
	SUITESPARSE_PRINTF ((
	"    (rows with more than max (%g * sqrt (n), 16) entries are\n"
	"    considered \"dense\", and placed last in output permutation)\n",
	alpha)) ;
    }

    if (aggressive)
    {
	SUITESPARSE_PRINTF (("    aggressive absorption:  yes\n")) ;
    }
    else
    {
	SUITESPARSE_PRINTF (("    aggressive absorption:  no\n")) ;
    }

    SUITESPARSE_PRINTF (("    size of AMD integer: %lu\n\n", sizeof (amd_int))) ;
}
