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



pybind11_add_module(highs_bindings highspy/highs_bindings.cpp)


# note: macOS is APPLE and also UNIX !
# if(APPLE)
#   set_target_properties(model_builder_helper_pybind11 PROPERTIES
#     SUFFIX ".so"
#     INSTALL_RPATH "@loader_path;@loader_path/../../../${PYTHON_PROJECT}/.libs"
#     )
#   set_property(TARGET model_builder_helper_pybind11 APPEND PROPERTY
#     LINK_FLAGS "-flat_namespace -undefined suppress"
#     )
# elseif(UNIX)
#   set_target_properties(model_builder_helper_pybind11 PROPERTIES
#     INSTALL_RPATH "$ORIGIN:$ORIGIN/../../../${PYTHON_PROJECT}/.libs"
#     )
# endif()

target_link_libraries(highs_bindings PRIVATE
  ${PROJECT_NAMESPACE}::highs
)

# add_library(${PROJECT_NAMESPACE}::model_builder_helper_pybind11 ALIAS model_builder_helper_pybind11)

# if(BUILD_TESTING)
#   file(GLOB PYTHON_SRCS "*_test.py")
#   foreach(FILE_NAME IN LISTS PYTHON_SRCS)
#     add_python_test(${FILE_NAME})
#   endforeach()
# endif()
