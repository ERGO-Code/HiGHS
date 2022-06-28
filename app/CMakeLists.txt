# create highs binary using library without pic
add_executable(highs)

target_sources(highs PRIVATE RunHighs.cpp)

if (UNIX)
        target_compile_options(highs PUBLIC "-Wno-unused-variable")
        target_compile_options(highs PUBLIC "-Wno-unused-const-variable")
        target_compile_options(highs PUBLIC "-Wno-sign-compare")
endif()

target_link_libraries(highs libhighs)

target_include_directories(highs PRIVATE
        $<BUILD_INTERFACE:${HIGHS_SOURCE_DIR}/app>  
        )

# install the binary
install(TARGETS highs EXPORT highs-targets
        RUNTIME)
