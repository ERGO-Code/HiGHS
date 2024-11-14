

set(HIGHS_DIR "")

if(NOT TARGET highs)
  include("${CMAKE_CURRENT_LIST_DIR}/highs-targets.cmake")
endif()

set(HIGHS_LIBRARIES highs)

set(HIGHS_INCLUDE_DIRS "/home/ivet/code/HiGHS/src;/home/ivet/code/HiGHS/build-cpu")

set(HIGHS_FOUND TRUE)

include(CMakeFindDependencyMacro)
find_dependency(Threads)
find_dependency(ZLIB)
