# CMakeLists.txt for highspy wheel
# Builds python bindings for HiGHS library
#
# We use static linking for easier distribution, but using a shared library can simplify 
# debugging and development, i.e., can simply update libhighs without rebuilding the wheel. 
# NOTE: HiGHS library can statically link to highs_extras if built with BUILD_SHARED_EXTRAS_LIB=OFF
cmake_minimum_required(VERSION 3.15...3.27)

project(_core LANGUAGES CXX)

# Find Python 3
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG)

# static link to libhighs
set(BUILD_SHARED_EXTRAS_LIB ON)
set(BUILD_SHARED_LIBS OFF)
set(BUILD_CXX ON)
set(BUILD_TESTING OFF)
set(BUILD_SAMPLES OFF)
set(BUILD_EXAMPLES OFF)
set(BUILD_CXX_EXE OFF)

# build highs library with output into build/highs directory
# add_subdirectory(
#     "${CMAKE_CURRENT_LIST_DIR}/../"
#     "${CMAKE_CURRENT_BINARY_DIR}/highs")

python_add_library(_core MODULE highs_bindings.cpp WITH_SOABI)
target_link_libraries(_core PRIVATE pybind11::headers highs ${CMAKE_DL_LIBS})

if (APPLE)
  set_target_properties(highs PROPERTIES INSTALL_RPATH "@loader_path")
else()
  set_target_properties(_core PROPERTIES INSTALL_RPATH "$ORIGIN")
endif()

set_target_properties(_core PROPERTIES
    BUILD_WITH_INSTALL_RPATH TRUE
)

target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})

if (${BUILD_SHARED_EXTRAS_LIB})
  target_compile_definitions(_core PRIVATE HIGHS_SHARED_EXTRAS_LIBRARY)
else()
  target_link_libraries(_core PRIVATE highs_extras)
endif()

target_include_directories(_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../highs>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../extern>)

if(MSVC)
  target_compile_options(_core PRIVATE "/bigobj")
endif()

if (NOT MSVC)
  target_compile_options(_core PRIVATE "-ftemplate-depth=2048")
endif()

# The install directory is the output (wheel) directory
install(TARGETS _core
        RUNTIME DESTINATION highspy
        LIBRARY DESTINATION highspy)

# Only needed if the wheel is built with a shared library, but can be left in for static linking as well.
# Include highs library (but ignore .lib file)
install(TARGETS highs
        RUNTIME DESTINATION highspy
        LIBRARY DESTINATION highspy
        NAMELINK_SKIP)
