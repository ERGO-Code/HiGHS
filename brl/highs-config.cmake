if(NOT TARGET libhighs)
  include("${CMAKE_CURRENT_LIST_DIR}/highs-targets.cmake")
endif()

set(HIGHS_LIBRARIES libhighs)
set(HIGHS_INCLUDE_DIRS "/Users/mac/projects/HiGHS/src;/Users/mac/projects/HiGHS/brl")
set(HIGHS_FOUND TRUE)
