if (NOT PYTHON_BUILD_SETUP)
  return()
endif()

set(CMAKE_VERBOSE_MAKEFILE ON)

include(sources-python)

set(sources_python ${highs_sources_python} 
                   ${cupdlp_sources_python} 
                   ${ipx_sources_python} 
                   ${basiclu_sources_python})

set(headers_python ${highs_headers_python} 
                   ${cupdlp_headers_python} 
                   ${ipx_headers_python} 
                   ${basiclu_headers_python})

# Find Python 3
find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)

# Pybind11
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

# add module
pybind11_add_module(highspy highspy/highs_bindings.cpp)

# todo is this version required?
# target_compile_definitions(highspy 
#                 PRIVATE VERSION_INFO=${VERSION_INFO})

# sources for python 
target_sources(highspy PUBLIC ${sources_python} ${headers_python})

# include directories for python 
target_include_directories(highspy PUBLIC ${include_dirs_python})
