add_executable(dev-presolve)

target_sources(dev-presolve PUBLIC
               ${PROJECT_SOURCE_DIR}/dev-presolve/DevPresolveMain.cpp)

target_include_directories(dev-presolve
                           PUBLIC 
                           ${PROJECT_SOURCE_DIR}/dev-presolve
                           ${PROJECT_SOURCE_DIR}/scaffold)


target_link_libraries(dev-presolve libhighs)