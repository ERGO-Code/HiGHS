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
  set_target_properties(_core PROPERTIES INSTALL_RPATH "@loader_path")
else()
  set_target_properties(_core PROPERTIES INSTALL_RPATH "$ORIGIN")
endif()

set_target_properties(_core PROPERTIES
    BUILD_WITH_INSTALL_RPATH TRUE
)

target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})


target_include_directories(_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../highs>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../extern>)

if(MSVC)
  target_compile_options(_core PRIVATE "/bigobj")
endif()

if (NOT MSVC)
  target_compile_options(_core PRIVATE "-ftemplate-depth=2048")
endif()

if(CUPDLP_GPU)
  target_link_libraries(_core PRIVATE cudalin ${CUDA_LIBRARY})
endif()

# The install directory is the output (wheel) directory
install(TARGETS _core
        RUNTIME DESTINATION highspy
        LIBRARY DESTINATION highspy
        COMPONENT python)

# install(TARGETS highs/cudalin ...) is silently ignored in CMake 4.x when the
# target is owned by a different directory, and also when BUILD_SHARED_LIBS is
# OFF in this scope (reset above). Use install(FILES $<TARGET_FILE:...>) which
# resolves the actual DLL path at install time with no directory restrictions.
if (CUPDLP_GPU)
  # Use install(FILES) to avoid CMake 4.x cross-directory install(TARGETS)
  # restrictions (silently ignored for targets owned by other directories).
  # $<TARGET_SONAME_FILE:highs> gives the SONAME file on Linux (libhighs.so.1.14)
  # or the DLL on Windows — the exact file the extension's linker records.
  install(FILES $<TARGET_SONAME_FILE:highs>
          DESTINATION highspy
          COMPONENT python)
  if (NOT HIGHS_GPU_LIB)
    install(FILES $<TARGET_FILE:cudalin>
            DESTINATION highspy
            COMPONENT python)
  endif()
else()
  # Non-GPU: on Unix highs is shared by default; LIBRARY+NAMELINK_SKIP installs
  # the SONAME file that the extension's RPATH resolves.
  install(TARGETS highs
          LIBRARY DESTINATION highspy
          NAMELINK_SKIP
          COMPONENT python)
endif()
