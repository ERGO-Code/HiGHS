if(NOT TARGET libh2gmp)
  include("${CMAKE_CURRENT_LIST_DIR}/h2gmp-targets.cmake")
endif()

set(H2GMP_LIBRARIES libh2gmp)
set(H2GMP_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/../../../include")
set(H2GMP_FOUND TRUE)
