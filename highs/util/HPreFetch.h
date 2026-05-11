/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HPreFetch.h
 * @brief Platform-independent prefetch wrappers for cache optimization
 */
#ifndef HIGHS_UTIL_HPREFETCH_H_
#define HIGHS_UTIL_HPREFETCH_H_

#include "util/HighsInt.h"

#ifndef HPC_PREFETCH_DIST
#define HPC_PREFETCH_DIST 8
#endif

#if defined(__GNUC__) || defined(__clang__)

#define HPC_PREFETCH_RD(addr) __builtin_prefetch((const void*)(addr), 0, 1)
#define HPC_PREFETCH_WR(addr) __builtin_prefetch((const void*)(addr), 1, 1)

#elif defined(_MSC_VER)

#include <intrin.h>
#define HPC_PREFETCH_RD(addr) _mm_prefetch((const char*)(addr), _MM_HINT_T1)
#define HPC_PREFETCH_WR(addr) _mm_prefetch((const char*)(addr), _MM_HINT_T1)

#else

#define HPC_PREFETCH_RD(addr) ((void)0)
#define HPC_PREFETCH_WR(addr) ((void)0)

#endif

// Branch prediction hints for hot paths
#if defined(__GNUC__) || defined(__clang__)
#define HIGHS_LIKELY(x)   __builtin_expect(!!(x), 1)
#define HIGHS_UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
#define HIGHS_LIKELY(x)   (x)
#define HIGHS_UNLIKELY(x) (x)
#endif

#endif /* HIGHS_UTIL_HPREFETCH_H_ */
