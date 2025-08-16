#ifndef HCONFIG_H_
#define HCONFIG_H_

// Minimal configuration for testing
#define FAST_BUILD
/* #undef ZLIB_FOUND */
/* #undef CUPDLP_CPU */
/* #undef CUPDLP_GPU */
#define CMAKE_BUILD_TYPE "Debug"
/* #undef HIGHSINT64 */
/* #undef HIGHS_NO_DEFAULT_THREADS */
/* #undef HIGHS_HAVE_MM_PAUSE */
/* #undef HIGHS_HAVE_BUILTIN_CLZ */
/* #undef HIGHS_HAVE_BITSCAN_REVERSE */

#define HIGHS_GITHASH "test-integration"
#define HIGHS_VERSION_MAJOR 1
#define HIGHS_VERSION_MINOR 7
#define HIGHS_VERSION_PATCH 2

#endif /* HCONFIG_H_ */