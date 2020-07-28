if(NOT TARGET libhighs)
  include("${CMAKE_CURRENT_LIST_DIR}/highs-targets.cmake")
endif()

set(HIGHS_LIBRARIES libhighs)
set(HIGHS_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/../../../include")
set(HIGHS_FOUND TRUE)
