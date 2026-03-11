
#ifndef HIGHS_EXPORT_H
#define HIGHS_EXPORT_H

#ifdef HIGHS_STATIC_DEFINE
#  define HIGHS_EXPORT
#  define HIGHS_NO_EXPORT
#else
#  ifndef HIGHS_EXPORT
#    ifdef highs_EXPORTS
        /* We are building this library */
#      define HIGHS_EXPORT __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define HIGHS_EXPORT __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef HIGHS_NO_EXPORT
#    define HIGHS_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef HIGHS_DEPRECATED
#  define HIGHS_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef HIGHS_DEPRECATED_EXPORT
#  define HIGHS_DEPRECATED_EXPORT HIGHS_EXPORT HIGHS_DEPRECATED
#endif

#ifndef HIGHS_DEPRECATED_NO_EXPORT
#  define HIGHS_DEPRECATED_NO_EXPORT HIGHS_NO_EXPORT HIGHS_DEPRECATED
#endif

/* NOLINTNEXTLINE(readability-avoid-unconditional-preprocessor-if) */
#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef HIGHS_NO_DEPRECATED
#    define HIGHS_NO_DEPRECATED
#  endif
#endif

#endif /* HIGHS_EXPORT_H */
