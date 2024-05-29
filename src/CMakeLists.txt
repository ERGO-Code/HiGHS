if (NOT BUILD_CXX)
  return()
endif()
  
# Define library.
include(sources)
set(sources ${highs_sources} ${cupdlp_sources} ${ipx_sources} ${basiclu_sources})
set(headers ${highs_headers} ${cupdlp_headers} ${ipx_headers} ${basiclu_headers})

# Configure the config windows version file
if(MSVC)
  string(REPLACE "." ","  PROJECT_RC_VERSION "${PROJECT_VERSION}")
  configure_file(${HIGHS_SOURCE_DIR}/version.rc.in
    "${HIGHS_BINARY_DIR}/version.rc" @ONLY)
  set(win_version_file ${HIGHS_BINARY_DIR}/version.rc)
else()
  set(win_version_file)
endif()

# Outdated CMake approach: update in progress
if(NOT FAST_BUILD)

  add_library(libhighs ${sources} ${headers} ${win_version_file})
  target_include_directories(libhighs PRIVATE ${include_dirs})

  if(${BUILD_SHARED_LIBS})
    # put version information into shared library file
    set_target_properties(libhighs PROPERTIES
      VERSION
      ${HIGHS_VERSION_MAJOR}.${HIGHS_VERSION_MINOR}.${HIGHS_VERSION_PATCH}
      SOVERSION ${HIGHS_VERSION_MAJOR}.${HIGHS_VERSION_MINOR})
    if(MINGW)
      target_compile_definitions(libhighs PUBLIC LIBHIGHS_STATIC_DEFINE)
    endif()
  else()
    # create static highs library with pic
    set_target_properties(libhighs PROPERTIES
      POSITION_INDEPENDENT_CODE on)
    target_compile_definitions(libhighs PUBLIC LIBHIGHS_STATIC_DEFINE)
  endif()

  # on UNIX system the 'lib' prefix is automatically added
  set_target_properties(libhighs PROPERTIES
    OUTPUT_NAME "highs"
    PDB_NAME "libhighs"
    MACOSX_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

  if(ZLIB AND ZLIB_FOUND)
    target_link_libraries(libhighs ZLIB::ZLIB)
    set(CONF_DEPENDENCIES "include(CMakeFindDependencyMacro)\nfind_dependency(ZLIB)")
  endif()

  # set the install rpath to the installed destination
  set_target_properties(libhighs PROPERTIES INSTALL_RPATH
    "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

  # install the header files of highs
  foreach(file ${headers})
    get_filename_component(dir ${file} DIRECTORY)

    if(NOT dir STREQUAL "")
      string(REPLACE ../extern/ "" dir ${dir})
    endif()

    install(FILES ${file} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs/${dir})
  endforeach()
  install(FILES ${HIGHS_BINARY_DIR}/HConfig.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs)

  if(UNIX)
    # target_compile_options(libhighs PRIVATE "-Wno-defaulted-function-deleted")
    # target_compile_options(libhighs PRIVATE "-Wno-return-type-c-linkage")
    target_compile_options(libhighs PRIVATE "-Wno-return-type" "-Wno-switch")

    target_compile_options(libhighs PRIVATE "-Wno-unused-variable")
    target_compile_options(libhighs PRIVATE "-Wno-unused-const-variable")

    # target_compile_options(libhighs PRIVATE "-Wno-sign-compare")
    # target_compile_options(libhighs PRIVATE "-Wno-logical-op-parentheses")

    # target_compile_options(libipx PRIVATE "-Wno-defaulted-function-deleted")
    # target_compile_options(libipx PRIVATE "-Wno-return-type-c-linkage")
    # target_compile_options(libipx PRIVATE "-Wno-return-type" "-Wno-switch")

    # target_compile_options(libipx PRIVATE "-Wno-sign-compare")
    # target_compile_options(libipx PRIVATE "-Wno-logical-op-parentheses")
  endif()
  
  install(TARGETS libhighs EXPORT highs-targets
    LIBRARY
    ARCHIVE
    RUNTIME
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs)

  # Add library targets to the build-tree export set
  export(TARGETS libhighs
    FILE "${HIGHS_BINARY_DIR}/highs-targets.cmake")

  # Configure the config file for the build tree:
  # Either list all the src/* directories here, or put explicit paths in all the
  # include statements.
  # M reckons that the latter is more transparent, and I'm inclined to agree.
  set(CONF_INCLUDE_DIRS "${HIGHS_SOURCE_DIR}/src" "${HIGHS_BINARY_DIR}")
  configure_file(${HIGHS_SOURCE_DIR}/highs-config.cmake.in
    "${HIGHS_BINARY_DIR}/highs-config.cmake" @ONLY)

  # Configure the config file for the install
  set(CONF_INCLUDE_DIRS "\${CMAKE_CURRENT_LIST_DIR}/../../../${CMAKE_INSTALL_INCLUDEDIR}/highs")
  configure_file(${HIGHS_SOURCE_DIR}/highs-config.cmake.in
    "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs-config.cmake" @ONLY)

  # Configure the pkg-config file for the install
  configure_file(${HIGHS_SOURCE_DIR}/highs.pc.in
    "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs.pc" @ONLY)

  # Install the targets of the highs export group, the config file so that other
  # cmake-projects can link easily against highs, and the pkg-config flie so that
  # other projects can easily build against highs
  install(EXPORT highs-targets FILE highs-targets.cmake DESTINATION
    ${CMAKE_INSTALL_LIBDIR}/cmake/highs)
  install(FILES "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs-config.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/highs)
  install(FILES "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs.pc"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

else()
  # FAST_BUILD is set to on.
  # At the moment used only for gradually updating the CMake targets build and
  # install / export.
  # Define library in modern CMake using target_*()
  # No interfaces (apart from c); No ipx; New (short) ctest instances.
  add_library(highs)
  add_library(${PROJECT_NAMESPACE}::highs ALIAS highs)

  if(${BUILD_SHARED_LIBS})
    # put version information into shared library file
    set_target_properties(highs PROPERTIES
      VERSION
      ${HIGHS_VERSION_MAJOR}.${HIGHS_VERSION_MINOR}.${HIGHS_VERSION_PATCH}
      SOVERSION ${HIGHS_VERSION_MAJOR}.${HIGHS_VERSION_MINOR})
  endif()

  set_target_properties(highs PROPERTIES POSITION_INDEPENDENT_CODE ON)

  if(APPLE)
      set_target_properties(highs PROPERTIES INSTALL_RPATH "@loader_path")
  elseif (UNIX)
      set_target_properties(highs PROPERTIES INSTALL_RPATH "$ORIGIN")
  endif()

  target_include_directories(highs PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${HIGHS_BINARY_DIR}>
    # $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/highs>
  )

  # target_include_directories(highs PRIVATE
  #   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interfaces>
  #   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/io>
  #   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/pdlp>
  #   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/pdlp/cupdlp>
  #   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ipm>
  #   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ipm/ipx>
  #   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ipm/basiclu>
  #   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/lp_data>
  #   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mip>
  #   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/model>
  #   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/parallel>
  #   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/presolve>
  #   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/qpsolver>
  #   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/simplex>
  #   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/util>
  #   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/test>)

  # target_include_directories(highs PRIVATE
  #   $<BUILD_INTERFACE:${HIGHS_SOURCE_DIR}/extern/>
  #   $<BUILD_INTERFACE:${HIGHS_SOURCE_DIR}/extern/filereader>
  #   $<BUILD_INTERFACE:${HIGHS_SOURCE_DIR}/extern/pdqsort>)

  target_sources(highs PRIVATE ${sources} ${headers} ${win_version_file})
  target_include_directories(highs PRIVATE ${include_dirs})

  if(MSVC)
    list(APPEND highs_compile_opts
      "/bigobj" # Allow big object
      "/DNOMINMAX"
      "/DWIN32_LEAN_AND_MEAN=1"
      "/D_CRT_SECURE_NO_WARNINGS"
      "/D_CRT_SECURE_NO_DEPRECATE"
      "/MP" # Build with multiple processes
      "/Zc:preprocessor" # Enable preprocessor conformance mode
      "/fp:precise"
      )
    if (CMAKE_BUILD_TYPE STREQUAL Release)
      list(APPEND highs_compile_opts
        "/DNDEBUG")
    endif()
    # MSVC warning suppressions
    list(APPEND highs_compile_opts
      "/wd4005" # 'macro-redefinition'
      "/wd4018" # 'expression' : signed/unsigned mismatch
      "/wd4065" # switch statement contains 'default' but no 'case' labels
      "/wd4068" # 'unknown pragma'
      "/wd4101" # 'identifier' : unreferenced local variable
      "/wd4146" # unary minus operator applied to unsigned type, result still unsigned
      "/wd4200" # nonstandard extension used : zero-sized array in struct/union
      "/wd4244" # 'conversion' conversion from 'type1' to 'type2', possible loss of data
      "/wd4251" # 'identifier' : class 'type' needs to have dll-interface to be used by clients of class 'type2'
      "/wd4267" # 'var' : conversion from 'size_t' to 'type', possible loss of data
      "/wd4305" # 'identifier' : truncation from 'type1' to 'type2'
      "/wd4307" # 'operator' : integral constant overflow
      "/wd4309" # 'conversion' : truncation of constant value
      "/wd4334" # 'operator' : result of 32-bit shift implicitly converted to 64 bits (was 64-bit shift intended?)
      "/wd4355" # 'this' : used in base member initializer list
      "/wd4477" # 'fwprintf' : format string '%s' requires an argument of type 'wchar_t *'
      "/wd4506" # no definition for inline function 'function'
      "/wd4715" # function' : not all control paths return a value
      "/wd4800" # 'type' : forcing value to bool 'true' or 'false' (performance warning)
      "/wd4996" # The compiler encountered a deprecated declaration.
      )
    target_compile_options(highs PRIVATE ${highs_compile_opts})
  # else()
  #   list(APPEND highs_compile_opts "-fwrapv")
  endif()


  if(ZLIB AND ZLIB_FOUND)
    target_include_directories(highs PRIVATE
      $<BUILD_INTERFACE:${HIGHS_SOURCE_DIR}/extern/zstr>
    )
    target_link_libraries(highs ZLIB::ZLIB)
    set(CONF_DEPS 
      "include(CMakeFindDependencyMacro)\nfind_dependency(Threads)\nfind_dependency(ZLIB)")
    set(CONF_DEPENDENCIES ${CONF_DEPS})
  else() 
    set(CONF_DEPENDENCIES "include(CMakeFindDependencyMacro)\nfind_dependency(Threads)")
  endif()
    

  # # on UNIX system the 'lib' prefix is automatically added
  # set_target_properties(highs PROPERTIES
  # OUTPUT_NAME "highs"
  # MACOSX_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

  # if (UNIX)
  # set_target_properties(highs PROPERTIES
  # LIBRARY_OUTPUT_DIRECTORY "${HIGHS_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
  # endif()
  
  # set_target_properties(highs PROPERTIES PUBLIC_HEADER "src/Highs.h;src/lp_data/HighsLp.h;src/lp_data/HighsLpSolverObject.h")

  # install the header files of highs
  foreach(file ${headers})
    get_filename_component(dir ${file} DIRECTORY)

    if(NOT dir STREQUAL "")
      string(REPLACE ../extern/ "" dir ${dir})
    endif()

    install(FILES ${file} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs/${dir})
  endforeach()
  install(FILES ${HIGHS_BINARY_DIR}/HConfig.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs)

  # target_compile_options(highs PRIVATE "-Wall")
  # target_compile_options(highs PRIVATE "-Wunused")
  
  if (UNIX)
    target_compile_options(highs PRIVATE "-Wno-unused-variable")
    target_compile_options(highs PRIVATE "-Wno-unused-const-variable")
  endif()

  # Configure the config file for the build tree:
  # Either list all the src/* directories here, or put explicit paths in all the
  # include statements.
  # M reckons that the latter is more transparent, and I'm inclined to agree.
  set(CONF_INCLUDE_DIRS "${HIGHS_SOURCE_DIR}/src" "${HIGHS_BINARY_DIR}")
  configure_file(${HIGHS_SOURCE_DIR}/highs-config.cmake.in
    "${HIGHS_BINARY_DIR}/highs-config.cmake" @ONLY)

  # Configure the config file for the install
  set(CONF_INCLUDE_DIRS "\${CMAKE_CURRENT_LIST_DIR}/../../../${CMAKE_INSTALL_INCLUDEDIR}")
  configure_file(${HIGHS_SOURCE_DIR}/highs-config.cmake.in
    "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs-config.cmake" @ONLY)

  # Configure the pkg-config file for the install
  configure_file(${HIGHS_SOURCE_DIR}/highs.pc.in
    "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs.pc" @ONLY)

  # Install the targets of the highs export group, the config file so that other
  # cmake-projects can link easily against highs, and the pkg-config flie so that
  # other projects can easily build against highs
  install(FILES "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs-config.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/highs)
  install(FILES "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs.pc"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)


  if (BUILD_DOTNET)
  
      # see: https://docs.microsoft.com/en-us/dotnet/core/rid-catalog
      if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64)")
        set(DOTNET_PLATFORM arm64)
      else()
        set(DOTNET_PLATFORM x64)
      endif()

      if(APPLE)
        set(DOTNET_RID osx-${DOTNET_PLATFORM})
      elseif(UNIX)
        set(DOTNET_RID linux-${DOTNET_PLATFORM})
      elseif(WIN32)
        set(DOTNET_RID win-${DOTNET_PLATFORM})
      else()
        message(FATAL_ERROR "Unsupported system !")
      endif()
      message(STATUS ".Net RID: ${DOTNET_RID}")

      set(DOTNET_PROJECT Highs.Native)
      set(DOTNET_PROJECT_DIR ${PROJECT_BINARY_DIR}/dotnet/${DOTNET_PROJECT})
      file(MAKE_DIRECTORY ${DOTNET_PROJECT_DIR}/runtimes/${DOTNET_RID}/native)

      if (APPLE)
        set(TARGET_FILE_NAME "highs.dylib")
      elseif(UNIX)
        set(TARGET_FILE_NAME "highs.so")
      elseif(WIN32)
        set(TARGET_FILE_NAME "highs.dll")
      endif()

      add_custom_command(TARGET highs POST_BUILD 
        COMMAND "${CMAKE_COMMAND}" -E copy 
          "$<TARGET_FILE:highs>"
        ${DOTNET_PROJECT_DIR}/runtimes/${DOTNET_RID}/native/${TARGET_FILE_NAME}
        COMMENT "Copying to output directory")
  endif()

endif()

if(FORTRAN_FOUND)
  set(fortransources interfaces/highs_fortran_api.f90)
  set(CMAKE_Fortran_MODULE_DIRECTORY ${HIGHS_BINARY_DIR}/modules)
  add_library(FortranHighs interfaces/highs_fortran_api.f90)

  if(NOT FAST_BUILD)
    target_link_libraries(FortranHighs PUBLIC libhighs)
  else()
    target_link_libraries(FortranHighs PUBLIC highs)
  endif()

  install(TARGETS FortranHighs 
    LIBRARY
    ARCHIVE
    RUNTIME
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs
    MODULES DESTINATION modules)
  if(NOT MSVC)
    install(FILES ${HIGHS_BINARY_DIR}/modules/highs_fortran_api.mod DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs/fortran)
  else()
    install(FILES ${HIGHS_BINARY_DIR}/modules/${CMAKE_BUILD_TYPE}/highs_fortran_api.mod DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs/fortran)
  endif()
  # use link rpath
  # set_target_properties(FortranHighs PROPERTIES INSTALL_RPATH
  #   "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
endif(FORTRAN_FOUND)

if(CSHARP_FOUND)
  message(STATUS "CSharp supported")
  set(csharpsources
    interfaces/highs_csharp_api.cs)
  add_library(HighsCsharp SHARED interfaces/highs_csharp_api.cs)
  add_library(${PROJECT_NAMESPACE}::HighsCsharp ALIAS HighsCsharp)
  target_compile_options(HighsCsharp PUBLIC "/unsafe")
  add_executable(csharpexample ../examples/call_highs_from_csharp.cs)
  target_compile_options(csharpexample PUBLIC "/unsafe")
  target_link_libraries(csharpexample PUBLIC HighsCsharp)
else()
  message(STATUS "No CSharp support")
endif()

find_package(Threads)
if(Threads_FOUND)
  include(CheckAtomic)
  if(HAVE_CXX_ATOMICS64_WITH_LIB)
    if(FAST_BUILD)
      target_link_libraries(highs atomic)
    else()
      target_link_libraries(libhighs atomic)
    endif()
  endif()
endif()

if(FAST_BUILD)
  target_link_libraries(highs Threads::Threads)
else()
  target_link_libraries(libhighs Threads::Threads)
endif()
