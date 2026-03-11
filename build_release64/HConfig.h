#ifndef HCONFIG_H_
#define HCONFIG_H_

/* #undef FAST_BUILD */
#define ZLIB_FOUND
#define CUPDLP_CPU
/* #undef CUPDLP_GPU */
/* #undef HIPO */
#define CMAKE_BUILD_TYPE "Release"
#define HIGHSINT64
/* #undef HIGHS_NO_DEFAULT_THREADS */
#define HIGHS_HAVE_MM_PAUSE
#define HIGHS_HAVE_BUILTIN_CLZ
/* #undef HIGHS_HAVE_BITSCAN_REVERSE */
/* #undef BLAS_LIBRARIES */

#define HIGHS_GITHASH "63ab434b2"
#define HIGHS_VERSION_MAJOR 1
#define HIGHS_VERSION_MINOR 13
#define HIGHS_VERSION_PATCH 1

#endif /* HCONFIG_H_ */
