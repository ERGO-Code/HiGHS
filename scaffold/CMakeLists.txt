# Define executable(s) to build
# In this first set up the executable dev-presolve builds and runs 
# the target defined in dev-presolve/ 

add_executable(dev-presolve)
target_include_directories(dev-presolve
                           PUBLIC 
                           ${PROJECT_SOURCE_DIR}/dev-presolve
                           ${PROJECT_SOURCE_DIR}/scaffold)
target_sources(dev-presolve PUBLIC
               ${PROJECT_SOURCE_DIR}/dev-presolve/DevPresolveMain.cpp)

target_link_libraries(dev-presolve libhighs)