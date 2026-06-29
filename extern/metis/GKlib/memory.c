/*!
\file  memory.c
\brief This file contains various allocation routines 

The allocation routines included are for 1D and 2D arrays of the 
most datatypes that GKlib support. Many of these routines are 
defined with the help of the macros in gk_memory.h. These macros 
can be used to define other memory allocation routines.

\date   Started 4/3/2007
\author George
\version\verbatim $Id: memory.c 21050 2017-05-25 03:53:58Z karypis $ \endverbatim
*/


#include "GKlib.h"

void gk_free(void **ptr1)
{
  free(*ptr1);
  *ptr1 = NULL;
}          

