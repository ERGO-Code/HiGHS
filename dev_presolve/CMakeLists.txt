# Create the interface cmake library for dev-presolve to be used in scaffold/
# and HiGHS/.
add_library(dev-presolve-lib INTERFACE)

target_sources(
  dev-presolve-lib
  INTERFACE ${HIGHS_SOURCE_DIR}/dev-presolve/DevPresolveMethods.hpp)

#target_link_libraries(libhighs PRIVATE dev-presolve)
target_include_directories(test-presolve
                           INTERFACE ${HIGHS_SOURCE_DIR}/component_test)

# Creates the executable for dev-presolve
add_executable(dev-presolve)

find_package(Boost REQUIRED OPTIONAL_COMPONENTS iostreams serialization program_options)

target_include_directories(
  dev-presolve
  PUBLIC ${HIGHS_SOURCE_DIR}/dev_presolve ${HIGHS_SOURCE_DIR}/component_test
         ${HIGHS_SOURCE_DIR}/scaffold
         ${Boost_INCLUDE_DIRS})

target_sources(dev-presolve
               PUBLIC ${HIGHS_SOURCE_DIR}/dev_presolve/DevPresolveMain.cpp)

target_link_libraries(dev-presolve libhighs ${Boost_LIBRARIES}) 

add_test(NAME dev-presolve COMMAND  dev-presolve)
        set_tests_properties (dev-presolve PROPERTIES
                PASS_REGULAR_EXPRESSION
                "dev-presolve: pass")