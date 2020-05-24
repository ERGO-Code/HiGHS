# create highs binary using library without pic
add_executable(highs)

target_sources(highs PRIVATE RunHighs.cpp)
target_include_directories(highs PRIVATE .)
target_link_libraries(highs libhighs)

# install the binary
# so it is generated in build/bin/ rather than build/app
install(TARGETS highs EXPORT highs-targets
        RUNTIME DESTINATION bin)
