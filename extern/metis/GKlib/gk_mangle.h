#ifndef HIGHS_GKLIB_MANGLE_H
#define HIGHS_GKLIB_MANGLE_H


/* memory.c */
#define gk_free             Highs_gk_free

/* error.c */
#define gk_errexit          Highs_gk_errexit

/* random.c */
#define gk_randint64        Highs_gk_randint64
#define gk_randint32        Highs_gk_randint32
#define gk_rand_r           Highs_gk_rand_r


/* mcore.c */
#define gk_mcoreCreate      Highs_gk_mcoreCreate
#define gk_mcoreDestroy     Highs_gk_mcoreDestroy
#define gk_mcoreMalloc      Highs_gk_mcoreMalloc
#define gk_mcorePush        Highs_gk_mcorePush
#define gk_mcorePop         Highs_gk_mcorePop
#define gk_mcoreAdd         Highs_gk_mcoreAdd
#define gk_mcoreDel         Highs_gk_mcoreDel


#endif