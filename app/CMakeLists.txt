# create highs binary using library without pic
add_executable(highs)

target_sources(highs PRIVATE RunHighs.cpp)
target_include_directories(highs PRIVATE .)
target_link_libraries(highs libhighs)

set_target_properties(highs PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# install the binary
install(TARGETS highs EXPORT highs-targets
        RUNTIME DESTINATION bin)
