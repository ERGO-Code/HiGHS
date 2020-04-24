# Define executable(s) to build. The scaffold executable will run tests we would
# like to run on each update of master. For test or dev it can also define and
# link with dev- and test- components.

# In this first set up the executable dev-presolve builds and runs the target
# defined in dev-presolve/

add_executable(scaffold)
target_include_directories(scaffold PUBLIC ${HIGHS_SOURCE_DIR}/dev-presolve
                                           ${HIGHS_SOURCE_DIR}/scaffold)

target_sources(scaffold PUBLIC ${HIGHS_SOURCE_DIR}/scaffold/ScaffoldMain.cpp)