if(FAST_BUILD)
  if (BUILD_CXX)
    # create highs binary using library without pic
    add_executable(highs-bin)

    # Configure the config windows version file
    if(MSVC)
      string(REPLACE "." ","  PROJECT_RC_VERSION "${PROJECT_VERSION}")
      configure_file(${HIGHS_SOURCE_DIR}/version.rc.in
        "${HIGHS_BINARY_DIR}/version.rc" @ONLY)
      set(win_version_file ${HIGHS_BINARY_DIR}/version.rc)
    else()
      set(win_version_file)
    endif()

    target_sources(highs-bin PRIVATE RunHighs.cpp HighsRuntimeOptions.h ${win_version_file})

    target_include_directories(highs-bin PRIVATE
      $<BUILD_INTERFACE:${HIGHS_SOURCE_DIR}/app>
    )

    if(UNIX)
      target_compile_options(highs-bin PUBLIC "-Wno-unused-variable")
      target_compile_options(highs-bin PUBLIC "-Wno-unused-const-variable")
    endif()

    set_target_properties(highs-bin
      PROPERTIES OUTPUT_NAME highs)

    target_link_libraries(highs-bin highs)

    if(APPLE)
        set_target_properties(highs-bin PROPERTIES INSTALL_RPATH 
        "@loader_path/../${CMAKE_INSTALL_LIBDIR};@loader_path")
    elseif (UNIX)
        set_target_properties(highs-bin PROPERTIES INSTALL_RPATH 
        "$ORIGIN:$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
    endif()

    # install the binary
    install(TARGETS highs-bin EXPORT highs-targets
      RUNTIME)
  endif()
else()
  # create highs binary using library without pic
  add_executable(highs)

  target_sources(highs PRIVATE RunHighs.cpp)

  if(UNIX)
    target_compile_options(highs PUBLIC "-Wno-unused-variable")
    target_compile_options(highs PUBLIC "-Wno-unused-const-variable")
  endif()

  target_link_libraries(highs libhighs)


  target_include_directories(highs PRIVATE
    $<BUILD_INTERFACE:${HIGHS_SOURCE_DIR}/app>
  )

  # install the binary
  install(TARGETS highs EXPORT highs-targets
    RUNTIME)
endif()

# Add demo to FAST_BUILD as well. 
if(EMSCRIPTEN AND EMSCRIPTEN_HTML)
  set(CMAKE_EXECUTABLE_SUFFIX ".html")
  set_target_properties(highs PROPERTIES LINK_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/highs_webdemo_shell.html)
endif()