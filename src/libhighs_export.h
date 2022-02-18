#ifndef LIBHIGHS_EXPORT_H
#define LIBHIGHS_EXPORT_H

#if defined _WIN32 || defined __CYGWIN__
  #ifdef BUILDING_LIBHIGHS
    #define LIBHIGHS_EXPORT __declspec(dllexport)
  #else
    #define LIBHIGHS_EXPORT __declspec(dllimport)
  #endif
  #define LIBHIGHS_NO_EXPORT
#else
  #define LIBHIGHS_EXPORT __attribute__ ((visibility ("default")))
  #define LIBHIGHS_NO_EXPORT  __attribute__ ((visibility ("hidden")))
#endif

#endif /* LIBHIGHS_EXPORT_H */
