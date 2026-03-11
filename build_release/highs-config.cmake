## HiGHS CMake configuration file

set(HIGHS_VERSION 1.13.1)


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was highs-config.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

####################################################################################

include(CMakeFindDependencyMacro)
find_dependency(Threads)

# Only look for BLAS if highs was built with it
set(_HIGHS_HAVE_BLAS OFF)

if(_HIGHS_HAVE_BLAS)
    find_dependency(BLAS)
endif()

# Let users know about optional features
set(HIGHS_HAVE_BLAS ${_HIGHS_HAVE_BLAS})

find_dependency(ZLIB)

include("${CMAKE_CURRENT_LIST_DIR}/highs-targets.cmake")
