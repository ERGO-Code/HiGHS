# create highs binary using library without pic
add_executable(highs)

target_sources(highs PRIVATE RunHighs.cpp)

target_compile_options(highs PUBLIC "-Wno-defaulted-function-deleted")
target_compile_options(highs PUBLIC "-Wno-unused-variable")
target_compile_options(highs PUBLIC "-Wno-unused-const-variable")
target_compile_options(highs PUBLIC "-Wno-sign-compare")

target_include_directories(highs PRIVATE .)
target_link_libraries(highs libhighs)

# install the binary
install(TARGETS highs EXPORT highs-targets
        RUNTIME DESTINATION bin)
