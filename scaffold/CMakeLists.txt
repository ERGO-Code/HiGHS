# Define executable(s) to build
# In this first set up the executable dev-presolve builds and runs 
# the target defined in dev-presolve/ 

# Will define scaffold exe once CMake is More Modern.

# add_executable(dev-presolve)
# target_include_directories(dev-presolve
#                            PUBLIC 
#                            ${HIGHS_SOURCE_DIR}/dev-presolve
#                            ${HIGHS_SOURCE_DIR}/scaffold)
# target_sources(dev-presolve PUBLIC
#                ${HIGHS_SOURCE_DIR}/dev_presolve/DevPresolveMain.cpp)

######
# Define scaffold exe target here since it has to be linked against HiGHS and for the 
# moment the include issues are still present due to the outdated CMAKE

add_executable(scaffold)
target_include_directories(scaffold
                           PUBLIC 
                           ${HIGHS_SOURCE_DIR}/dev-presolve
                           ${HIGHS_SOURCE_DIR}/scaffold)

 target_sources(scaffold PUBLIC
                ${HIGHS_SOURCE_DIR}/scaffold/ScaffoldMain.cpp)