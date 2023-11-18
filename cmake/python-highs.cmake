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
# search_python_module(
#   NAME pybind11
#   PACKAGE pybind11)

set(PYTHON_PROJECT "highspy")
message(STATUS "Python project: ${PYTHON_PROJECT}")
set(PYTHON_PROJECT_DIR ${PROJECT_BINARY_DIR}/python/${PYTHON_PROJECT})
message(STATUS "Python project build path: ${PYTHON_PROJECT_DIR}")



pybind11_add_module(highs_bindings highspy/highs_bindings.cpp)
set_target_properties(highs_bindings PROPERTIES
  LIBRARY_OUTPUT_NAME "highs_bindings")


if(APPLE)
  set_target_properties(highs_bindings PROPERTIES
    SUFFIX ".so"
    INSTALL_RPATH "@loader_path;@loader_path/../../../${PYTHON_PROJECT}/.libs"
    )
elseif(UNIX)
  set_target_properties(highs_bindings PROPERTIES
    INSTALL_RPATH "$ORIGIN:$ORIGIN/../../../${PYTHON_PROJECT}/.libs"
    )
endif()

add_library(${PROJECT_NAMESPACE}::highs_bindings ALIAS highs_bindings)

target_link_libraries(highs_bindings PRIVATE
  ${PROJECT_NAMESPACE}::highs
)

file(GENERATE OUTPUT ${PYTHON_PROJECT_DIR}/__init__.py CONTENT "")

file(COPY
  highspy/highs_bindings.cpp
  setup.py 
  pyproject.toml
  DESTINATION ${PYTHON_PROJECT_DIR})

# add_custom_command(
#   OUTPUT python/dist/timestamp
#   # Don't need to copy static lib on Windows.
#   COMMAND ${CMAKE_COMMAND} -E $<IF:$<STREQUAL:$<TARGET_PROPERTY:highs,TYPE>,SHARED_LIBRARY>,copy,true>
#   $<$<STREQUAL:$<TARGET_PROPERTY:highs,TYPE>,SHARED_LIBRARY>:$<TARGET_SONAME_FILE:highs>>
#   libs

#   COMMAND ${CMAKE_COMMAND} -E copy $<${TARGET_FILE}::highs> python/
#   WORKING_DIRECTORY ${PYTHON_PROJECT_DIR}
# )

add_custom_command(
  OUTPUT python/dist/timestamp
  COMMAND ${CMAKE_COMMAND} -E remove_directory dist
  COMMAND ${CMAKE_COMMAND} -E make_directory ${PYTHON_PROJECT}/.libs
  # Don't need to copy static lib on Windows.
  COMMAND ${CMAKE_COMMAND} -E $<IF:$<STREQUAL:$<TARGET_PROPERTY:highs,TYPE>,SHARED_LIBRARY>,copy,true>
  $<$<STREQUAL:$<TARGET_PROPERTY:highs,TYPE>,SHARED_LIBRARY>:$<TARGET_SONAME_FILE:highs>>
  ${PYTHON_PROJECT}/.libs
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:highs> ${PYTHON_PROJECT}/.libs
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:highs_bindings> ${PYTHON_PROJECT}/.libs

  #COMMAND ${Python3_EXECUTABLE} setup.py bdist_egg bdist_wheel
  COMMAND ${Python3_EXECUTABLE} setup.py bdist_wheel

  COMMAND ${CMAKE_COMMAND} -E touch ${PROJECT_BINARY_DIR}/python/dist/timestamp

  DEPENDS
    ${PROJECT_NAMESPACE}::highs
  BYPRODUCTS
    python/${PYTHON_PROJECT}
    python/${PYTHON_PROJECT}.egg-info
    python/build
    python/dist
  WORKING_DIRECTORY python/highspy
  COMMAND_EXPAND_LISTS)

# main target
add_custom_target(python_package all
  depends
    python/dist/timestamp
  working_directory python)

  set(VENV_EXECUTABLE ${Python3_EXECUTABLE} -m virtualenv)


set(INSTALL_PYTHON ON)

if(INSTALL_PYTHON)
  # make a virtualenv to install our python package in it
  add_custom_command(TARGET python_package POST_BUILD
    # Clean previous install otherwise pip install may do nothing

    # Must NOT call it in a folder containing the setup.py otherwise pip call it
    # (i.e. "python setup.py bdist") while we want to consume the wheel package
    COMMAND ${VENV_Python3_EXECUTABLE} -m pip install ${CMAKE_CURRENT_BINARY_DIR}/python
      
    # install modules only required to run examples
    COMMAND ${VENV_Python3_EXECUTABLE} -m pip install pytest

    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Install ${PYTHON_PROJECT}"
    VERBATIM)
endif()
