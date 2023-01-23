# create highs binary using library without pic
add_executable(highs-bin)

target_sources(highs-bin PRIVATE RunHighs.cpp)

target_include_directories(highs-bin PRIVATE
        $<BUILD_INTERFACE:${HIGHS_SOURCE_DIR}/app>  
        )

if (UNIX)
        target_compile_options(highs-bin PUBLIC "-Wno-unused-variable")
        target_compile_options(highs-bin PUBLIC "-Wno-unused-const-variable")
endif()

set_target_properties(highs-bin
        PROPERTIES OUTPUT_NAME highs)

target_link_libraries(highs-bin highs)


# install the binary
install(TARGETS highs-bin EXPORT highs-targets
        RUNTIME)
