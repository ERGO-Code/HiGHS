## HiGHS CMake configuration file

set(_VERSION 1.13.1)



include(CMakeFindDependencyMacro)
find_dependency(Threads)

# Only look for BLAS if highs was built with it
set(_HIGHS_HAVE_BLAS OFF)

if(_HIGHS_HAVE_BLAS)
    find_dependency(BLAS)
endif()

# Let users know about optional features
set(HIGHS_HAVE_BLAS ${_HIGHS_HAVE_BLAS})



include("${CMAKE_CURRENT_LIST_DIR}/highs-targets.cmake")
