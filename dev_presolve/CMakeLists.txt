add_executable(dev-presolve)

target_sources(dev-presolve PUBLIC
               ${HIGHS_SOURCE_DIR}/dev-presolve/DevPresolveMain.cpp)

target_include_directories(dev-presolve
                           PUBLIC 
                           ${HIGHS_SOURCE_DIR}/dev_presolve
                           ${HIGHS_SOURCE_DIR}/scaffold)