/*!
\file gk_macros.h
\brief This file contains various macros

\date   Started 3/27/2007
\author George
\version\verbatim $Id: gk_macros.h 15048 2013-08-31 19:38:14Z karypis $ \endverbatim
*/

#ifndef _GK_MACROS_H_
#define _GK_MACROS_H_

/*-------------------------------------------------------------
 * Usefull commands 
 *-------------------------------------------------------------*/
#define gk_max(a, b) ((a) >= (b) ? (a) : (b))
#define gk_min(a, b) ((a) >= (b) ? (b) : (a))
#define gk_SWAP(a, b, tmp) do {(tmp) = (a); (a) = (b); (b) = (tmp);} while(0) 
#define INC_DEC(a, b, val) do {(a) += (val); (b) -= (val);} while(0)

#define ONEOVERRANDMAX (1.0/(RAND_MAX+1.0))
#define RandomInRange(u, rng_state) ((int) (ONEOVERRANDMAX*(u)*my_rand_r(rng_state)))
#define RandomInRange_r(s, u) ((int) (ONEOVERRANDMAX*(u)*rand_r(s)))

/*-------------------------------------------------------------
 * dbglvl handling macros
 *-------------------------------------------------------------*/
#define IFSET(a, flag, cmd) if ((a)&(flag)) (cmd);

/*-------------------------------------------------------------
 * CSR conversion macros
 *-------------------------------------------------------------*/
#define MAKECSR(i, n, a) \
   do { \
     for (i=1; i<n; i++) a[i] += a[i-1]; \
     for (i=n; i>0; i--) a[i] = a[i-1]; \
     a[0] = 0; \
   } while(0) 

#define SHIFTCSR(i, n, a) \
   do { \
     for (i=n; i>0; i--) a[i] = a[i-1]; \
     a[0] = 0; \
   } while(0) 

#endif
