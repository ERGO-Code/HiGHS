add_executable(dev-presolve)

target_include_directories(dev-presolve
                           PUBLIC 
                           ${HIGHS_SOURCE_DIR}/dev_presolve
                           ${HIGHS_SOURCE_DIR}/component_test
                           ${HIGHS_SOURCE_DIR}/scaffold)

target_sources(dev-presolve PUBLIC
               ${HIGHS_SOURCE_DIR}/dev_presolve/DevPresolveMain.cpp)