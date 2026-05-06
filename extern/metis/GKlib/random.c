/*!
\file  
\brief Various routines for providing portable 32 and 64 bit random number
       generators.

\date   Started 5/17/2007
\author George
\version\verbatim $Id: random.c 18796 2015-06-02 11:39:45Z karypis $ \endverbatim
*/

#include "GKlib.h"

/*************************************************************************/
/*! Define function rand_r, which may not exist on certain machines */
/*************************************************************************/
int my_rand_r(unsigned *rng_state) {
  // Linear congruential generator, with values from wikipedia
  int result = ((*rng_state * 1103515245) + 12345) & 0x7fffffff;
  *rng_state = result;
  return result;
}

/* generates a random number on [0, 2^64-1]-interval */
uint64_t gk_randint64(unsigned *rng_state)
{
uint64_t piece_1 = ((uint64_t) my_rand_r(rng_state)) << 32;
uint64_t piece_2 = ((uint64_t) my_rand_r(rng_state));
  return (uint64_t)(piece_1 | piece_2);
}

/* generates a random number on [0, 2^32-1]-interval */
uint32_t gk_randint32(unsigned *rng_state)
{
  return (uint32_t)my_rand_r(rng_state);
}


