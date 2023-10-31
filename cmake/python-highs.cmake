find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 REQUIRED)

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
    if(FETCH_PYTHON_DEPS)
      message(WARNING "Can't find python module: \"${MODULE_NAME}\", install it using pip...")
      execute_process(
        COMMAND ${Python3_EXECUTABLE} -m pip install --user ${MODULE_PACKAGE}
        OUTPUT_STRIP_TRAILING_WHITESPACE
        COMMAND_ERROR_IS_FATAL ANY
      )
    else()
      message(FATAL_ERROR "Can't find python module: \"${MODULE_NAME}\", please install it using your system package manager.")
    endif()
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
    )
endif()

# if(APPLE)
#   set_target_properties(highs_bindings PROPERTIES
#     SUFFIX ".so"
#     INSTALL_RPATH "@loader_path;@loader_path/../../${PYTHON_PROJECT}/libs"
#     )
# elseif(UNIX)
#   set_target_properties(highs_bindings PROPERTIES
#     INSTALL_RPATH "$ORIGIN:$ORIGIN/../../${PYTHON_PROJECT}/libs"
#     )
# endif()

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



add_custom_command(
  OUTPUT python/dist/timestamp
  # Don't need to copy static lib on Windows.
  COMMAND ${CMAKE_COMMAND} -E $<IF:$<STREQUAL:$<TARGET_PROPERTY:highs,TYPE>,SHARED_LIBRARY>,copy,true>
  $<$<STREQUAL:$<TARGET_PROPERTY:highs,TYPE>,SHARED_LIBRARY>:$<TARGET_SONAME_FILE:highs>>
  libs

  COMMAND ${CMAKE_COMMAND} -E copy $<${TARGET_FILE}::highs> python/
  WORKING_DIRECTORY ${PYTHON_PROJECT_DIR}
)

  # add_custom_command(
  # OUTPUT python/dist/timestamp
  # COMMAND ${CMAKE_COMMAND} -E remove_directory dist
  # COMMAND ${CMAKE_COMMAND} -E make_directory ${PYTHON_PROJECT}/.libs
  # # Don't need to copy static lib on Windows.
  # COMMAND ${CMAKE_COMMAND} -E $<IF:$<STREQUAL:$<TARGET_PROPERTY:highs,TYPE>,SHARED_LIBRARY>,copy,true>
  # $<$<STREQUAL:$<TARGET_PROPERTY:highs,TYPE>,SHARED_LIBRARY>:$<TARGET_SONAME_FILE:highs>>
  # ${PYTHON_PROJECT}/.libs
  # COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:highs> ${PYTHON_PROJECT}/libs
  # COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:highs_bindings> ${PYTHON_PROJECT}/libs

  # #COMMAND ${Python3_EXECUTABLE} setup.py bdist_egg bdist_wheel
  # COMMAND ${Python3_EXECUTABLE} setup.py bdist_wheel

  # COMMAND ${CMAKE_COMMAND} -E touch ${PROJECT_BINARY_DIR}/python/dist/timestamp

  # DEPENDS
  #   python/setup.py
  #   ${PROJECT_NAMESPACE}::highs
  # BYPRODUCTS
  #   python/${PYTHON_PROJECT}
  #   python/${PYTHON_PROJECT}.egg-info
  #   python/build
  #   python/dist
  # WORKING_DIRECTORY python
  # COMMAND_EXPAND_LISTS)

# Main Target
# add_custom_target(python_package ALL
#   DEPENDS
#     python/dist/timestamp
#   WORKING_DIRECTORY python)

# Install rules
# configure_file(
#   ${PROJECT_SOURCE_DIR}/cmake/python-install.cmake.in
#   ${PROJECT_BINARY_DIR}/python/python-install.cmake
#   @ONLY)
# install(SCRIPT ${PROJECT_BINARY_DIR}/python/python-install.cmake)

# if(BUILD_VENV)
#   # make a virtualenv to install our python package in it
#   add_custom_command(TARGET python_package POST_BUILD
#     # Clean previous install otherwise pip install may do nothing
#     COMMAND ${CMAKE_COMMAND} -E remove_directory ${VENV_DIR}
#     COMMAND ${VENV_EXECUTABLE} -p ${Python3_EXECUTABLE}
#     $<IF:$<BOOL:${VENV_USE_SYSTEM_SITE_PACKAGES}>,--system-site-packages,-q>
#       ${VENV_DIR}
#     #COMMAND ${VENV_EXECUTABLE} ${VENV_DIR}
#     # Must NOT call it in a folder containing the setup.py otherwise pip call it
#     # (i.e. "python setup.py bdist") while we want to consume the wheel package
#     COMMAND ${VENV_Python3_EXECUTABLE} -m pip install
#       --find-links=${CMAKE_CURRENT_BINARY_DIR}/python/dist ${PYTHON_PROJECT}==${PROJECT_VERSION}
#     # install modules only required to run examples
#     COMMAND ${VENV_Python3_EXECUTABLE} -m pip install pandas matplotlib pytest scipy
#     BYPRODUCTS ${VENV_DIR}
#     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#     COMMENT "Create venv and install ${PYTHON_PROJECT}"
#     VERBATIM)
# endif()
