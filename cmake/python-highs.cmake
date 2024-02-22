set(CMAKE_VERBOSE_MAKEFILE ON)

include(sources)

# Find Python 3

find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)

include(FetchContent)

message(CHECK_START "Fetching pybind11")
list(APPEND CMAKE_MESSAGE_INDENT "  ")
set(PYBIND11_INSTALL ON)
set(PYBIND11_TEST OFF)
FetchContent_Declare(
  pybind11
  GIT_REPOSITORY "https://github.com/pybind/pybind11.git"
  GIT_TAG "v2.11.1"
)
FetchContent_MakeAvailable(pybind11)
list(POP_BACK CMAKE_MESSAGE_INDENT)
message(CHECK_PASS "fetched")

function(search_python_module)
  set(options NO_VERSION)
  set(oneValueArgs NAME PACKAGE)
  set(multiValueArgs "")
  cmake_parse_arguments(MODULE
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
  )
  message(STATUS "Searching python module: \"${MODULE_NAME}\"")
  if(${MODULE_NO_VERSION})
    execute_process(
      COMMAND ${Python3_EXECUTABLE} -c "import ${MODULE_NAME}"
      RESULT_VARIABLE _RESULT
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    set(MODULE_VERSION "unknown")
  else()
    execute_process(
      COMMAND ${Python3_EXECUTABLE} -c "import ${MODULE_NAME}; print(${MODULE_NAME}.__version__)"
      RESULT_VARIABLE _RESULT
      OUTPUT_VARIABLE MODULE_VERSION
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  endif()
  if(${_RESULT} STREQUAL "0")
    message(STATUS "Found python module: \"${MODULE_NAME}\" (found version \"${MODULE_VERSION}\")")
  else()
      message(FATAL_ERROR "Can't find python module: \"${MODULE_NAME}\", please install it using your system package manager.")
  endif()
endfunction()

search_python_module(
  NAME setuptools
  PACKAGE setuptools)
search_python_module(
  NAME wheel
  PACKAGE wheel)

set(PYTHON_PROJECT "highspy")
message(STATUS "Python project: ${PYTHON_PROJECT}")
set(PYTHON_PROJECT_DIR ${PROJECT_BINARY_DIR}/${PYTHON_PROJECT})
message(STATUS "Python project build path: ${PYTHON_PROJECT_DIR}")

pybind11_add_module(highspy highspy/highs_bindings.cpp)

target_compile_definitions(highspy 
                PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})

# set_target_properties(highspy PROPERTIES
#   LIBRARY_OUTPUT_NAME "highspy")

target_include_directories(highspy PUBLIC ${include_dirs})

target_sources(highspy PUBLIC
    ${cupdlp_sources}
    ${ipx_sources}
    ${basiclu_sources}
    ${highs_sources}
)

# target_include_directories(highs_bindings PUBLIC src)
#  target_include_directories(highs_bindings PUBLIC ${CMAKE_SOURCE_DIR}/src)

  # target_include_directories(highs_bindings PUBLIC
  #   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  #   $<BUILD_INTERFACE:${HIGHS_BINARY_DIR}>
  #   $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/highs>
  # )


  # if(APPLE)
#   set_target_properties(highs_bindings PROPERTIES
#     SUFFIX ".so"
#     INSTALL_RPATH "@loader_path;@loader_path/../../../${PYTHON_PROJECT_DIR}/.libs"
#     )
# elseif(UNIX)
#   set_target_properties(highs_bindings PROPERTIES
#     INSTALL_RPATH "$ORIGIN:$ORIGIN/../../../${PYTHON_PROJECT_DIR}/.libs"
#     )
# endif()

add_library(${PROJECT_NAMESPACE}::highspy ALIAS highspy)

target_compile_definitions(highspy
                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})

# target_link_libraries(highs_bindings PRIVATE
#   ${PROJECT_NAMESPACE}::highs
# )

# file(GENERATE OUTPUT ${PYTHON_PROJECT_DIR}/__init__.py CONTENT "")

# file(COPY
#   highspy/setup.py 
#   highspy/pyproject.toml
#   highspy/README.md
#   DESTINATION ${PYTHON_PROJECT_DIR})

# file(COPY
#   highspy/__init__.py
#   highspy/highs.py
#   highspy/highs_bindings.cpp
#   DESTINATION ${PYTHON_PROJECT_DIR}/highspy)

# add_custom_command(
#   OUTPUT highspy/dist/timestamp
#   COMMAND ${CMAKE_COMMAND} -E remove_directory dist
#   COMMAND ${CMAKE_COMMAND} -E make_directory ${PYTHON_PROJECT_DIR}/.libs
#   # # Don't need to copy static lib on Windows.
#   COMMAND ${CMAKE_COMMAND} -E $<IF:$<STREQUAL:$<TARGET_PROPERTY:highs,TYPE>,SHARED_LIBRARY>,copy,true>
#   $<$<STREQUAL:$<TARGET_PROPERTY:highs,TYPE>,SHARED_LIBRARY>:$<TARGET_SONAME_FILE:highs>>
#   ${PYTHON_PROJECT_DIR}/.libs
  
#   #COMMAND ${Python3_EXECUTABLE} setup.py bdist_egg bdist_wheel
#   # COMMAND ${Python3_EXECUTABLE} setup.py bdist_wheel
#   # COMMAND ${CMAKE_COMMAND} -E touch ${PROJECT_BINARY_DIR}/highspy/dist/timestamp

#   # BYPRODUCTS
#   #   highspy/${PYTHON_PROJECT}.egg-info
#   #   highspy/build
#   #   highspy/distoutput_output_flagflag
#   DEPENDS
#     ${PROJECT_NAMESPACE}::highs
#   WORKING_DIRECTORY ${PYTHON_PROJECT_DIR}
#   COMMAND_EXPAND_LISTS)

# # main target
# add_custom_target(python_package all
#   DEPENDS depends
#     python/dist/timestamp
#   WORKING_DIRECTORY ${PYTHON_PROJECT_DIR})
