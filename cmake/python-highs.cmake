if(NOT PYTHON)
  return()
endif()

if(NOT TARGET ${PROJECT_NAMESPACE}::highs)
  message(FATAL_ERROR "Python: missing ortools TARGET")
endif()

# Find Python 3
find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)

# Find if the python module is available,
# otherwise install it (PACKAGE_NAME) to the Python3 user install directory.
# If CMake option FETCH_PYTHON_DEPS is OFF then issue a fatal error instead.
# e.g
# search_python_module(
#   NAME
#     mypy_protobuf
#   PACKAGE
#     mypy-protobuf
#   NO_VERSION
# )
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
        )
    else()
      message(FATAL_ERROR "Can't find python module: \"${MODULE_NAME}\", please install it using your system package manager.")
    endif()
  endif()
endfunction()

# Find if a python builtin module is available.
# e.g
# search_python_internal_module(
#   NAME
#     mypy_protobuf
# )
function(search_python_internal_module)
  set(options "")
  set(oneValueArgs NAME)
  set(multiValueArgs "")
  cmake_parse_arguments(MODULE
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
  )
  message(STATUS "Searching python module: \"${MODULE_NAME}\"")
  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import ${MODULE_NAME}"
    RESULT_VARIABLE _RESULT
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  if(${_RESULT} STREQUAL "0")
    message(STATUS "Found python internal module: \"${MODULE_NAME}\"")
  else()
    message(FATAL_ERROR "Can't find python internal module \"${MODULE_NAME}\", please install it using your system package manager.")
  endif()
endfunction()

if(BUILD_TESTING)
  # add_python_test()
  # CMake function to generate and build python test.
  # Parameters:
  #  the python filename
  # e.g.:
  # add_python_test(foo.py)
  function(add_python_test FILE_NAME)
    message(STATUS "Configuring test ${FILE_NAME} ...")
    get_filename_component(TEST_NAME ${FILE_NAME} NAME_WE)
    get_filename_component(WRAPPER_DIR ${FILE_NAME} DIRECTORY)
    get_filename_component(COMPONENT_DIR ${WRAPPER_DIR} DIRECTORY)
    get_filename_component(COMPONENT_NAME ${COMPONENT_DIR} NAME)

    add_test(
      NAME python_${COMPONENT_NAME}_${TEST_NAME}
      COMMAND ${VENV_Python3_EXECUTABLE} -m pytest ${FILE_NAME}
      WORKING_DIRECTORY ${VENV_DIR})
    message(STATUS "Configuring test ${FILE_NAME} done")
  endfunction()
endif()

#######################
##  PYTHON WRAPPERS  ##
#######################
list(APPEND CMAKE_SWIG_FLAGS "-I${PROJECT_SOURCE_DIR}")

set(PYTHON_PROJECT ${PROJECT_NAME})
message(STATUS "Python project: ${PYTHON_PROJECT}")
set(PYTHON_PROJECT_DIR ${PROJECT_BINARY_DIR}/python/${PYTHON_PROJECT})
message(STATUS "Python project build path: ${PYTHON_PROJECT_DIR}")

# #######################
# ## Python Packaging  ##
# #######################
# #file(MAKE_DIRECTORY python/${PYTHON_PROJECT})
# file(GENERATE OUTPUT ${PYTHON_PROJECT_DIR}/__init__.py CONTENT "__version__ = \"${PROJECT_VERSION}\"\n")
# file(GENERATE OUTPUT ${PYTHON_PROJECT_DIR}/algorithms/__init__.py CONTENT "")
# file(GENERATE OUTPUT ${PYTHON_PROJECT_DIR}/bop/__init__.py CONTENT "")
# file(GENERATE OUTPUT ${PYTHON_PROJECT_DIR}/constraint_solver/__init__.py CONTENT "")
# file(GENERATE OUTPUT ${PYTHON_PROJECT_DIR}/glop/__init__.py CONTENT "")
# file(GENERATE OUTPUT ${PYTHON_PROJECT_DIR}/graph/__init__.py CONTENT "")
# file(GENERATE OUTPUT ${PYTHON_PROJECT_DIR}/graph/python/__init__.py CONTENT "")
# file(GENERATE OUTPUT ${PYTHON_PROJECT_DIR}/init/__init__.py CONTENT "")
# file(GENERATE OUTPUT ${PYTHON_PROJECT_DIR}/linear_solver/__init__.py CONTENT "")
# file(GENERATE OUTPUT ${PYTHON_PROJECT_DIR}/linear_solver/python/__init__.py CONTENT "")
# file(GENERATE OUTPUT ${PYTHON_PROJECT_DIR}/packing/__init__.py CONTENT "")
# file(GENERATE OUTPUT ${PYTHON_PROJECT_DIR}/pdlp/__init__.py CONTENT "")
# file(GENERATE OUTPUT ${PYTHON_PROJECT_DIR}/sat/__init__.py CONTENT "")
# file(GENERATE OUTPUT ${PYTHON_PROJECT_DIR}/sat/python/__init__.py CONTENT "")
# file(GENERATE OUTPUT ${PYTHON_PROJECT_DIR}/scheduling/__init__.py CONTENT "")
# file(GENERATE OUTPUT ${PYTHON_PROJECT_DIR}/util/__init__.py CONTENT "")
# file(GENERATE OUTPUT ${PYTHON_PROJECT_DIR}/util/python/__init__.py CONTENT "")

# file(COPY
#   ortools/linear_solver/linear_solver_natural_api.py
#   DESTINATION ${PYTHON_PROJECT_DIR}/linear_solver)
# file(COPY
#   ortools/linear_solver/python/model_builder.py
#   ortools/linear_solver/python/model_builder_helper.py
#   DESTINATION ${PYTHON_PROJECT_DIR}/linear_solver/python)
# file(COPY
#   ortools/sat/python/cp_model.py
#   ortools/sat/python/cp_model_helper.py
#   ortools/sat/python/visualization.py
#   DESTINATION ${PYTHON_PROJECT_DIR}/sat/python)

# # setup.py.in contains cmake variable e.g. @PYTHON_PROJECT@ and
# # generator expression e.g. $<TARGET_FILE_NAME:pyFoo>
# configure_file(
#   ${PROJECT_SOURCE_DIR}/ortools/python/setup.py.in
#   ${PROJECT_BINARY_DIR}/python/setup.py.in
#   @ONLY)
# file(GENERATE
#   OUTPUT ${PROJECT_BINARY_DIR}/python/setup.py
#   INPUT ${PROJECT_BINARY_DIR}/python/setup.py.in)

# #add_custom_command(
# #  OUTPUT python/setup.py
# #  DEPENDS ${PROJECT_BINARY_DIR}/python/setup.py
# #  COMMAND ${CMAKE_COMMAND} -E copy setup.py setup.py
# #  WORKING_DIRECTORY python)

# configure_file(
#   ${PROJECT_SOURCE_DIR}/ortools/python/README.pypi.txt
#   ${PROJECT_BINARY_DIR}/python/README.txt
#   COPYONLY)

# # Look for required python modules
# search_python_module(
#   NAME setuptools
#   PACKAGE setuptools)
# search_python_module(
#   NAME wheel
#   PACKAGE wheel)

# add_custom_command(
#   OUTPUT python/dist/timestamp
#   COMMAND ${CMAKE_COMMAND} -E remove_directory dist
#   COMMAND ${CMAKE_COMMAND} -E make_directory ${PYTHON_PROJECT}/.libs
#   # Don't need to copy static lib on Windows.
#   COMMAND ${CMAKE_COMMAND} -E $<IF:$<STREQUAL:$<TARGET_PROPERTY:ortools,TYPE>,SHARED_LIBRARY>,copy,true>
#   $<$<STREQUAL:$<TARGET_PROPERTY:ortools,TYPE>,SHARED_LIBRARY>:$<TARGET_SONAME_FILE:ortools>>
#   ${PYTHON_PROJECT}/.libs
#   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:pywrapinit> ${PYTHON_PROJECT}/init
#   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:pywrapknapsack_solver> ${PYTHON_PROJECT}/algorithms
#   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:linear_sum_assignment_pybind11> ${PYTHON_PROJECT}/graph/python
#   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:max_flow_pybind11> ${PYTHON_PROJECT}/graph/python
#   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:min_cost_flow_pybind11> ${PYTHON_PROJECT}/graph/python
#   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:pywrapcp> ${PYTHON_PROJECT}/constraint_solver
#   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:pywraplp> ${PYTHON_PROJECT}/linear_solver
#   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:pywrap_model_builder_helper> ${PYTHON_PROJECT}/linear_solver/python
#   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:swig_helper> ${PYTHON_PROJECT}/sat/python
#   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:pywraprcpsp> ${PYTHON_PROJECT}/scheduling
#   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:sorted_interval_list> ${PYTHON_PROJECT}/util/python
#   #COMMAND ${Python3_EXECUTABLE} setup.py bdist_egg bdist_wheel
#   COMMAND ${Python3_EXECUTABLE} setup.py bdist_wheel
#   COMMAND ${CMAKE_COMMAND} -E touch ${PROJECT_BINARY_DIR}/python/dist/timestamp
#   MAIN_DEPENDENCY
#     ortools/python/setup.py.in
#   DEPENDS
#     python/setup.py
#     Py${PROJECT_NAME}_proto
#     ${PROJECT_NAMESPACE}::ortools
#     pywrapinit
#     pywrapknapsack_solver
#     linear_sum_assignment_pybind11
#     max_flow_pybind11
#     min_cost_flow_pybind11
#     pywrapcp
#     pywraplp
#     pywrap_model_builder_helper
#     swig_helper
#     pywraprcpsp
#     sorted_interval_list
#   BYPRODUCTS
#     python/${PYTHON_PROJECT}
#     python/${PYTHON_PROJECT}.egg-info
#     python/build
#     python/dist
#   WORKING_DIRECTORY python
#   COMMAND_EXPAND_LISTS)

# # Main Target
# add_custom_target(python_package ALL
#   DEPENDS
#     python/dist/timestamp
#   WORKING_DIRECTORY python)

# # Install rules
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
#     COMMAND ${VENV_Python3_EXECUTABLE} -m pip install pandas matplotlib pytest
#     BYPRODUCTS ${VENV_DIR}
#     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#     COMMENT "Create venv and install ${PYTHON_PROJECT}"
#     VERBATIM)
# endif()

# if(BUILD_TESTING)
#   configure_file(
#     ${PROJECT_SOURCE_DIR}/ortools/init/python/version_test.py.in
#     ${PROJECT_BINARY_DIR}/python/version_test.py
#     @ONLY)

#   # run the tests within the virtualenv
#   add_test(NAME python_init_version_test
#     COMMAND ${VENV_Python3_EXECUTABLE} ${PROJECT_BINARY_DIR}/python/version_test.py)
# endif()

# #####################
# ##  Python Sample  ##
# #####################
# # add_python_sample()
# # CMake function to generate and build python sample.
# # Parameters:
# #  the python filename
# # e.g.:
# # add_python_sample(foo.py)
# function(add_python_sample FILE_NAME)
#   message(STATUS "Configuring sample ${FILE_NAME} ...")
#   get_filename_component(SAMPLE_NAME ${FILE_NAME} NAME_WE)
#   get_filename_component(SAMPLE_DIR ${FILE_NAME} DIRECTORY)
#   get_filename_component(COMPONENT_DIR ${SAMPLE_DIR} DIRECTORY)
#   get_filename_component(COMPONENT_NAME ${COMPONENT_DIR} NAME)

#   if(BUILD_TESTING)
#     add_test(
#       NAME python_${COMPONENT_NAME}_${SAMPLE_NAME}
#       COMMAND ${VENV_Python3_EXECUTABLE} ${FILE_NAME}
#       WORKING_DIRECTORY ${VENV_DIR})
#   endif()
#   message(STATUS "Configuring sample ${FILE_NAME} done")
# endfunction()

######################
##  Python Example  ##
######################
# add_python_example()
# CMake function to generate and build python example.
# Parameters:
#  the python filename
# e.g.:
# add_python_example(foo.py)
function(add_python_example FILE_NAME)
  message(STATUS "Configuring example ${FILE_NAME} ...")
  get_filename_component(EXAMPLE_NAME ${FILE_NAME} NAME_WE)
  get_filename_component(COMPONENT_DIR ${FILE_NAME} DIRECTORY)
  get_filename_component(COMPONENT_NAME ${COMPONENT_DIR} NAME)

  if(BUILD_TESTING)
    add_test(
      NAME python_${COMPONENT_NAME}_${EXAMPLE_NAME}
      COMMAND ${VENV_Python3_EXECUTABLE} ${FILE_NAME}
      WORKING_DIRECTORY ${VENV_DIR})
  endif()
  message(STATUS "Configuring example ${FILE_NAME} done")
endfunction()
