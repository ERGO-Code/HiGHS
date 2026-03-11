# CMake generated Testfile for 
# Source directory: /home/yzhou/Github/HiGHS/check
# Build directory: /home/yzhou/Github/HiGHS/build_release64/check
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[capi_unit_tests]=] "/home/yzhou/Github/HiGHS/build_release64/bin/capi_unit_tests")
set_tests_properties([=[capi_unit_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;152;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([=[unit-test-build]=] "/snap/cmake/1515/bin/cmake" "--build" "/home/yzhou/Github/HiGHS/build_release64" "--target" "unit_tests")
set_tests_properties([=[unit-test-build]=] PROPERTIES  RESOURCE_LOCK "unittestbin" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;158;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([=[unit_tests_all]=] "/home/yzhou/Github/HiGHS/build_release64/bin/unit_tests" "--success")
set_tests_properties([=[unit_tests_all]=] PROPERTIES  DEPENDS "unit-test-build" TIMEOUT "10000" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;171;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[25fv47--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "/home/yzhou/Github/HiGHS/check/instances/25fv47.mps")
set_tests_properties([==[25fv47--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Simplex   iterations: 3149
Objective value     :  5.5018458883" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[25fv47--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "/home/yzhou/Github/HiGHS/check/instances/25fv47.mps")
set_tests_properties([==[25fv47--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  5.5018458883" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[25fv47--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "/home/yzhou/Github/HiGHS/check/instances/25fv47.mps")
set_tests_properties([==[25fv47--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  5.5018458883" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[25fv47--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "/home/yzhou/Github/HiGHS/check/instances/25fv47.mps")
set_tests_properties([==[25fv47--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  5.5018458883" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[25fv47--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "/home/yzhou/Github/HiGHS/check/instances/25fv47.mps")
set_tests_properties([==[25fv47--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  5.5018458883" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[80bau3b--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "/home/yzhou/Github/HiGHS/check/instances/80bau3b.mps")
set_tests_properties([==[80bau3b--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Simplex   iterations: 3686
Objective value     :  9.8722419241" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[80bau3b--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "/home/yzhou/Github/HiGHS/check/instances/80bau3b.mps")
set_tests_properties([==[80bau3b--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  9.8722419241" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[80bau3b--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "/home/yzhou/Github/HiGHS/check/instances/80bau3b.mps")
set_tests_properties([==[80bau3b--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  9.8722419241" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[80bau3b--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "/home/yzhou/Github/HiGHS/check/instances/80bau3b.mps")
set_tests_properties([==[80bau3b--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  9.8722419241" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[80bau3b--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "/home/yzhou/Github/HiGHS/check/instances/80bau3b.mps")
set_tests_properties([==[80bau3b--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  9.8722419241" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[adlittle--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "/home/yzhou/Github/HiGHS/check/instances/adlittle.mps")
set_tests_properties([==[adlittle--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Simplex   iterations: 74
Objective value     :  2.2549496316" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[adlittle--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "/home/yzhou/Github/HiGHS/check/instances/adlittle.mps")
set_tests_properties([==[adlittle--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  2.2549496316" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[adlittle--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "/home/yzhou/Github/HiGHS/check/instances/adlittle.mps")
set_tests_properties([==[adlittle--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  2.2549496316" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[adlittle--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "/home/yzhou/Github/HiGHS/check/instances/adlittle.mps")
set_tests_properties([==[adlittle--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  2.2549496316" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[adlittle--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "/home/yzhou/Github/HiGHS/check/instances/adlittle.mps")
set_tests_properties([==[adlittle--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  2.2549496316" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[afiro--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "/home/yzhou/Github/HiGHS/check/instances/afiro.mps")
set_tests_properties([==[afiro--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Simplex   iterations: 22
Objective value     : -4.6475314286" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[afiro--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "/home/yzhou/Github/HiGHS/check/instances/afiro.mps")
set_tests_properties([==[afiro--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     : -4.6475314286" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[afiro--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "/home/yzhou/Github/HiGHS/check/instances/afiro.mps")
set_tests_properties([==[afiro--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     : -4.6475314286" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[afiro--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "/home/yzhou/Github/HiGHS/check/instances/afiro.mps")
set_tests_properties([==[afiro--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     : -4.6475314286" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[afiro--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "/home/yzhou/Github/HiGHS/check/instances/afiro.mps")
set_tests_properties([==[afiro--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     : -4.6475314286" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[etamacro--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "/home/yzhou/Github/HiGHS/check/instances/etamacro.mps")
set_tests_properties([==[etamacro--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Simplex   iterations: 532
Objective value     : -7.5571523330" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[etamacro--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "/home/yzhou/Github/HiGHS/check/instances/etamacro.mps")
set_tests_properties([==[etamacro--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     : -7.5571523330" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[etamacro--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "/home/yzhou/Github/HiGHS/check/instances/etamacro.mps")
set_tests_properties([==[etamacro--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     : -7.5571523330" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[etamacro--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "/home/yzhou/Github/HiGHS/check/instances/etamacro.mps")
set_tests_properties([==[etamacro--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     : -7.5571523330" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[etamacro--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "/home/yzhou/Github/HiGHS/check/instances/etamacro.mps")
set_tests_properties([==[etamacro--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     : -7.5571523330" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[greenbea--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "/home/yzhou/Github/HiGHS/check/instances/greenbea.mps")
set_tests_properties([==[greenbea--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Simplex   iterations: 5109
Objective value     : -7.2555248130" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[greenbea--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "/home/yzhou/Github/HiGHS/check/instances/greenbea.mps")
set_tests_properties([==[greenbea--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     : -7.2555248130" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[greenbea--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "/home/yzhou/Github/HiGHS/check/instances/greenbea.mps")
set_tests_properties([==[greenbea--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     : -7.2555248130" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[greenbea--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "/home/yzhou/Github/HiGHS/check/instances/greenbea.mps")
set_tests_properties([==[greenbea--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     : -7.2555248130" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[greenbea--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "/home/yzhou/Github/HiGHS/check/instances/greenbea.mps")
set_tests_properties([==[greenbea--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     : -7.2555248130" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[shell--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "/home/yzhou/Github/HiGHS/check/instances/shell.mps")
set_tests_properties([==[shell--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Simplex   iterations: 623
Objective value     :  1.2088253460" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[shell--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "/home/yzhou/Github/HiGHS/check/instances/shell.mps")
set_tests_properties([==[shell--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  1.2088253460" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[shell--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "/home/yzhou/Github/HiGHS/check/instances/shell.mps")
set_tests_properties([==[shell--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  1.2088253460" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[shell--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "/home/yzhou/Github/HiGHS/check/instances/shell.mps")
set_tests_properties([==[shell--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  1.2088253460" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[shell--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "/home/yzhou/Github/HiGHS/check/instances/shell.mps")
set_tests_properties([==[shell--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  1.2088253460" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[stair--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "/home/yzhou/Github/HiGHS/check/instances/stair.mps")
set_tests_properties([==[stair--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Simplex   iterations: 529
Objective value     : -2.5126695119" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[stair--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "/home/yzhou/Github/HiGHS/check/instances/stair.mps")
set_tests_properties([==[stair--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     : -2.5126695119" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[stair--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "/home/yzhou/Github/HiGHS/check/instances/stair.mps")
set_tests_properties([==[stair--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     : -2.5126695119" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[stair--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "/home/yzhou/Github/HiGHS/check/instances/stair.mps")
set_tests_properties([==[stair--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     : -2.5126695119" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[stair--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "/home/yzhou/Github/HiGHS/check/instances/stair.mps")
set_tests_properties([==[stair--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     : -2.5126695119" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[standata--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "/home/yzhou/Github/HiGHS/check/instances/standata.mps")
set_tests_properties([==[standata--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Simplex   iterations: 72
Objective value     :  1.2576995000" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[standata--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "/home/yzhou/Github/HiGHS/check/instances/standata.mps")
set_tests_properties([==[standata--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  1.2576995000" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[standata--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "/home/yzhou/Github/HiGHS/check/instances/standata.mps")
set_tests_properties([==[standata--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  1.2576995000" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[standata--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "/home/yzhou/Github/HiGHS/check/instances/standata.mps")
set_tests_properties([==[standata--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  1.2576995000" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[standata--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "/home/yzhou/Github/HiGHS/check/instances/standata.mps")
set_tests_properties([==[standata--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  1.2576995000" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[standgub--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "/home/yzhou/Github/HiGHS/check/instances/standgub.mps")
set_tests_properties([==[standgub--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Simplex   iterations: 68
Objective value     :  1.2576995000" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[standgub--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "/home/yzhou/Github/HiGHS/check/instances/standgub.mps")
set_tests_properties([==[standgub--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  1.2576995000" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[standgub--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "/home/yzhou/Github/HiGHS/check/instances/standgub.mps")
set_tests_properties([==[standgub--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  1.2576995000" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[standgub--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "/home/yzhou/Github/HiGHS/check/instances/standgub.mps")
set_tests_properties([==[standgub--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  1.2576995000" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[standgub--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "/home/yzhou/Github/HiGHS/check/instances/standgub.mps")
set_tests_properties([==[standgub--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  1.2576995000" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[standmps--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "/home/yzhou/Github/HiGHS/check/instances/standmps.mps")
set_tests_properties([==[standmps--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Simplex   iterations: 218
Objective value     :  1.4060175000" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[standmps--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "/home/yzhou/Github/HiGHS/check/instances/standmps.mps")
set_tests_properties([==[standmps--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  1.4060175000" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[standmps--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "/home/yzhou/Github/HiGHS/check/instances/standmps.mps")
set_tests_properties([==[standmps--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  1.4060175000" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[standmps--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "/home/yzhou/Github/HiGHS/check/instances/standmps.mps")
set_tests_properties([==[standmps--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  1.4060175000" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[standmps--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "/home/yzhou/Github/HiGHS/check/instances/standmps.mps")
set_tests_properties([==[standmps--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  1.4060175000" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;447;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[bgetam--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "/home/yzhou/Github/HiGHS/check/instances/bgetam.mps")
set_tests_properties([==[bgetam--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[bgetam--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "/home/yzhou/Github/HiGHS/check/instances/bgetam.mps")
set_tests_properties([==[bgetam--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[bgetam--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "/home/yzhou/Github/HiGHS/check/instances/bgetam.mps")
set_tests_properties([==[bgetam--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[bgetam--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "/home/yzhou/Github/HiGHS/check/instances/bgetam.mps")
set_tests_properties([==[bgetam--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[bgetam--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "/home/yzhou/Github/HiGHS/check/instances/bgetam.mps")
set_tests_properties([==[bgetam--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[box1--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "/home/yzhou/Github/HiGHS/check/instances/box1.mps")
set_tests_properties([==[box1--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[box1--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "/home/yzhou/Github/HiGHS/check/instances/box1.mps")
set_tests_properties([==[box1--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[box1--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "/home/yzhou/Github/HiGHS/check/instances/box1.mps")
set_tests_properties([==[box1--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[box1--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "/home/yzhou/Github/HiGHS/check/instances/box1.mps")
set_tests_properties([==[box1--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[box1--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "/home/yzhou/Github/HiGHS/check/instances/box1.mps")
set_tests_properties([==[box1--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[ex72a--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "/home/yzhou/Github/HiGHS/check/instances/ex72a.mps")
set_tests_properties([==[ex72a--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[ex72a--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "/home/yzhou/Github/HiGHS/check/instances/ex72a.mps")
set_tests_properties([==[ex72a--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[ex72a--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "/home/yzhou/Github/HiGHS/check/instances/ex72a.mps")
set_tests_properties([==[ex72a--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[ex72a--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "/home/yzhou/Github/HiGHS/check/instances/ex72a.mps")
set_tests_properties([==[ex72a--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[ex72a--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "/home/yzhou/Github/HiGHS/check/instances/ex72a.mps")
set_tests_properties([==[ex72a--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[forest6--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "/home/yzhou/Github/HiGHS/check/instances/forest6.mps")
set_tests_properties([==[forest6--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[forest6--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "/home/yzhou/Github/HiGHS/check/instances/forest6.mps")
set_tests_properties([==[forest6--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[forest6--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "/home/yzhou/Github/HiGHS/check/instances/forest6.mps")
set_tests_properties([==[forest6--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[forest6--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "/home/yzhou/Github/HiGHS/check/instances/forest6.mps")
set_tests_properties([==[forest6--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[forest6--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "/home/yzhou/Github/HiGHS/check/instances/forest6.mps")
set_tests_properties([==[forest6--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[galenet--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "/home/yzhou/Github/HiGHS/check/instances/galenet.mps")
set_tests_properties([==[galenet--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[galenet--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "/home/yzhou/Github/HiGHS/check/instances/galenet.mps")
set_tests_properties([==[galenet--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[galenet--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "/home/yzhou/Github/HiGHS/check/instances/galenet.mps")
set_tests_properties([==[galenet--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[galenet--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "/home/yzhou/Github/HiGHS/check/instances/galenet.mps")
set_tests_properties([==[galenet--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[galenet--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "/home/yzhou/Github/HiGHS/check/instances/galenet.mps")
set_tests_properties([==[galenet--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[gams10am--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "/home/yzhou/Github/HiGHS/check/instances/gams10am.mps")
set_tests_properties([==[gams10am--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[gams10am--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "/home/yzhou/Github/HiGHS/check/instances/gams10am.mps")
set_tests_properties([==[gams10am--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[gams10am--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "/home/yzhou/Github/HiGHS/check/instances/gams10am.mps")
set_tests_properties([==[gams10am--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[gams10am--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "/home/yzhou/Github/HiGHS/check/instances/gams10am.mps")
set_tests_properties([==[gams10am--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[gams10am--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "/home/yzhou/Github/HiGHS/check/instances/gams10am.mps")
set_tests_properties([==[gams10am--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[refinery--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "/home/yzhou/Github/HiGHS/check/instances/refinery.mps")
set_tests_properties([==[refinery--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[refinery--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "/home/yzhou/Github/HiGHS/check/instances/refinery.mps")
set_tests_properties([==[refinery--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[refinery--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "/home/yzhou/Github/HiGHS/check/instances/refinery.mps")
set_tests_properties([==[refinery--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[refinery--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "/home/yzhou/Github/HiGHS/check/instances/refinery.mps")
set_tests_properties([==[refinery--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[refinery--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "/home/yzhou/Github/HiGHS/check/instances/refinery.mps")
set_tests_properties([==[refinery--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[woodinfe--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "/home/yzhou/Github/HiGHS/check/instances/woodinfe.mps")
set_tests_properties([==[woodinfe--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[woodinfe--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "/home/yzhou/Github/HiGHS/check/instances/woodinfe.mps")
set_tests_properties([==[woodinfe--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[woodinfe--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "/home/yzhou/Github/HiGHS/check/instances/woodinfe.mps")
set_tests_properties([==[woodinfe--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[woodinfe--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "/home/yzhou/Github/HiGHS/check/instances/woodinfe.mps")
set_tests_properties([==[woodinfe--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[woodinfe--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "/home/yzhou/Github/HiGHS/check/instances/woodinfe.mps")
set_tests_properties([==[woodinfe--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Infeasible" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;451;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[gas11--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "/home/yzhou/Github/HiGHS/check/instances/gas11.mps")
set_tests_properties([==[gas11--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Unbounded" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;452;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[gas11--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "/home/yzhou/Github/HiGHS/check/instances/gas11.mps")
set_tests_properties([==[gas11--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Unbounded" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;452;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[gas11--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "/home/yzhou/Github/HiGHS/check/instances/gas11.mps")
set_tests_properties([==[gas11--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Unbounded" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;452;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[gas11--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "/home/yzhou/Github/HiGHS/check/instances/gas11.mps")
set_tests_properties([==[gas11--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Unbounded" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;452;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[gas11--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "/home/yzhou/Github/HiGHS/check/instances/gas11.mps")
set_tests_properties([==[gas11--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model status        : Unbounded" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;410;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;452;add_instancetests;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[small_mip--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/small_mip.mps")
set_tests_properties([==[small_mip--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      3.2368421.*
  Dual bound        3.2368421.*
  Solution status   feasible
                    3.2368421.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[small_mip--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/small_mip.mps")
set_tests_properties([==[small_mip--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      3.2368421.*
  Dual bound        3.2368421.*
  Solution status   feasible
                    3.2368421.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[small_mip--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/small_mip.mps")
set_tests_properties([==[small_mip--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      3.2368421.*
  Dual bound        3.2368421.*
  Solution status   feasible
                    3.2368421.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[small_mip--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/small_mip.mps")
set_tests_properties([==[small_mip--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      3.2368421.*
  Dual bound        3.2368421.*
  Solution status   feasible
                    3.2368421.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[small_mip--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/small_mip.mps")
set_tests_properties([==[small_mip--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      3.2368421.*
  Dual bound        3.2368421.*
  Solution status   feasible
                    3.2368421.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[flugpl--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/flugpl.mps")
set_tests_properties([==[flugpl--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      1201500.*
  Dual bound        1201500.*
  Solution status   feasible
                    1201500.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[flugpl--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/flugpl.mps")
set_tests_properties([==[flugpl--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      1201500.*
  Dual bound        1201500.*
  Solution status   feasible
                    1201500.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[flugpl--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/flugpl.mps")
set_tests_properties([==[flugpl--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      1201500.*
  Dual bound        1201500.*
  Solution status   feasible
                    1201500.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[flugpl--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/flugpl.mps")
set_tests_properties([==[flugpl--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      1201500.*
  Dual bound        1201500.*
  Solution status   feasible
                    1201500.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[flugpl--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/flugpl.mps")
set_tests_properties([==[flugpl--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      1201500.*
  Dual bound        1201500.*
  Solution status   feasible
                    1201500.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[lseu--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/lseu.mps")
set_tests_properties([==[lseu--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      1120|1119.9999999.*
  Dual bound        1120|1119.9999999.*
  Solution status   feasible
                    1120|1119.9999999.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[lseu--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/lseu.mps")
set_tests_properties([==[lseu--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      1120|1119.9999999.*
  Dual bound        1120|1119.9999999.*
  Solution status   feasible
                    1120|1119.9999999.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[lseu--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/lseu.mps")
set_tests_properties([==[lseu--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      1120|1119.9999999.*
  Dual bound        1120|1119.9999999.*
  Solution status   feasible
                    1120|1119.9999999.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[lseu--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/lseu.mps")
set_tests_properties([==[lseu--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      1120|1119.9999999.*
  Dual bound        1120|1119.9999999.*
  Solution status   feasible
                    1120|1119.9999999.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[lseu--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/lseu.mps")
set_tests_properties([==[lseu--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      1120|1119.9999999.*
  Dual bound        1120|1119.9999999.*
  Solution status   feasible
                    1120|1119.9999999.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[egout--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/egout.mps")
set_tests_properties([==[egout--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      (568.1007|568.1006999).*
  Dual bound        (568.1007|568.1006999).*
  Solution status   feasible
                    (568.1007|568.1006999).* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[egout--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/egout.mps")
set_tests_properties([==[egout--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      (568.1007|568.1006999).*
  Dual bound        (568.1007|568.1006999).*
  Solution status   feasible
                    (568.1007|568.1006999).* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[egout--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/egout.mps")
set_tests_properties([==[egout--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      (568.1007|568.1006999).*
  Dual bound        (568.1007|568.1006999).*
  Solution status   feasible
                    (568.1007|568.1006999).* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[egout--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/egout.mps")
set_tests_properties([==[egout--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      (568.1007|568.1006999).*
  Dual bound        (568.1007|568.1006999).*
  Solution status   feasible
                    (568.1007|568.1006999).* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[egout--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/egout.mps")
set_tests_properties([==[egout--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      (568.1007|568.1006999).*
  Dual bound        (568.1007|568.1006999).*
  Solution status   feasible
                    (568.1007|568.1006999).* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[gt2--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/gt2.mps")
set_tests_properties([==[gt2--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      21166.*
  Dual bound        21166.*
  Solution status   feasible
                    21166.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[gt2--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/gt2.mps")
set_tests_properties([==[gt2--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      21166.*
  Dual bound        21166.*
  Solution status   feasible
                    21166.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[gt2--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/gt2.mps")
set_tests_properties([==[gt2--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      21166.*
  Dual bound        21166.*
  Solution status   feasible
                    21166.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[gt2--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/gt2.mps")
set_tests_properties([==[gt2--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      21166.*
  Dual bound        21166.*
  Solution status   feasible
                    21166.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[gt2--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/gt2.mps")
set_tests_properties([==[gt2--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      21166.*
  Dual bound        21166.*
  Solution status   feasible
                    21166.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[rgn--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/rgn.mps")
set_tests_properties([==[rgn--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      82.19999.*
  Dual bound        82.19999.*
  Solution status   feasible
                    82.19999.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[rgn--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/rgn.mps")
set_tests_properties([==[rgn--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      82.19999.*
  Dual bound        82.19999.*
  Solution status   feasible
                    82.19999.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[rgn--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/rgn.mps")
set_tests_properties([==[rgn--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      82.19999.*
  Dual bound        82.19999.*
  Solution status   feasible
                    82.19999.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[rgn--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/rgn.mps")
set_tests_properties([==[rgn--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      82.19999.*
  Dual bound        82.19999.*
  Solution status   feasible
                    82.19999.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[rgn--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/rgn.mps")
set_tests_properties([==[rgn--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      82.19999.*
  Dual bound        82.19999.*
  Solution status   feasible
                    82.19999.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[bell5--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/bell5.mps")
set_tests_properties([==[bell5--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      (8966406.49152|8966406.491519|8966406.49151).*
  Dual bound        (8966406.49152|8966406.491519|8966406.49151).*
  Solution status   feasible
                    (8966406.49152|8966406.491519|8966406.49151).* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[bell5--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/bell5.mps")
set_tests_properties([==[bell5--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      (8966406.49152|8966406.491519|8966406.49151).*
  Dual bound        (8966406.49152|8966406.491519|8966406.49151).*
  Solution status   feasible
                    (8966406.49152|8966406.491519|8966406.49151).* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[bell5--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/bell5.mps")
set_tests_properties([==[bell5--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      (8966406.49152|8966406.491519|8966406.49151).*
  Dual bound        (8966406.49152|8966406.491519|8966406.49151).*
  Solution status   feasible
                    (8966406.49152|8966406.491519|8966406.49151).* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[bell5--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/bell5.mps")
set_tests_properties([==[bell5--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      (8966406.49152|8966406.491519|8966406.49151).*
  Dual bound        (8966406.49152|8966406.491519|8966406.49151).*
  Solution status   feasible
                    (8966406.49152|8966406.491519|8966406.49151).* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[bell5--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/bell5.mps")
set_tests_properties([==[bell5--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      (8966406.49152|8966406.491519|8966406.49151).*
  Dual bound        (8966406.49152|8966406.491519|8966406.49151).*
  Solution status   feasible
                    (8966406.49152|8966406.491519|8966406.49151).* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[sp150x300d--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/sp150x300d.mps")
set_tests_properties([==[sp150x300d--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      (69|68.9999999).*
  Dual bound        (69|68.9999999).*
  Solution status   feasible
                    (69|68.9999999).* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[sp150x300d--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/sp150x300d.mps")
set_tests_properties([==[sp150x300d--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      (69|68.9999999).*
  Dual bound        (69|68.9999999).*
  Solution status   feasible
                    (69|68.9999999).* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[sp150x300d--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/sp150x300d.mps")
set_tests_properties([==[sp150x300d--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      (69|68.9999999).*
  Dual bound        (69|68.9999999).*
  Solution status   feasible
                    (69|68.9999999).* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[sp150x300d--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/sp150x300d.mps")
set_tests_properties([==[sp150x300d--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      (69|68.9999999).*
  Dual bound        (69|68.9999999).*
  Solution status   feasible
                    (69|68.9999999).* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[sp150x300d--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/sp150x300d.mps")
set_tests_properties([==[sp150x300d--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      (69|68.9999999).*
  Dual bound        (69|68.9999999).*
  Solution status   feasible
                    (69|68.9999999).* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[p0548--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/p0548.mps")
set_tests_properties([==[p0548--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      (8691|8690.9999999).*
  Dual bound        (8691|8690.9999999).*
  Solution status   feasible
                    (8691|8690.9999999).* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[p0548--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/p0548.mps")
set_tests_properties([==[p0548--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      (8691|8690.9999999).*
  Dual bound        (8691|8690.9999999).*
  Solution status   feasible
                    (8691|8690.9999999).* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[p0548--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/p0548.mps")
set_tests_properties([==[p0548--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      (8691|8690.9999999).*
  Dual bound        (8691|8690.9999999).*
  Solution status   feasible
                    (8691|8690.9999999).* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[p0548--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/p0548.mps")
set_tests_properties([==[p0548--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      (8691|8690.9999999).*
  Dual bound        (8691|8690.9999999).*
  Solution status   feasible
                    (8691|8690.9999999).* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[p0548--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/p0548.mps")
set_tests_properties([==[p0548--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      (8691|8690.9999999).*
  Dual bound        (8691|8690.9999999).*
  Solution status   feasible
                    (8691|8690.9999999).* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[dcmulti--presolve=off]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=off" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/dcmulti.mps")
set_tests_properties([==[dcmulti--presolve=off]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      188182.*
  Dual bound        188182.*
  Solution status   feasible
                    188182.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[dcmulti--presolve=on]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--presolve=on" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/dcmulti.mps")
set_tests_properties([==[dcmulti--presolve=on]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      188182.*
  Dual bound        188182.*
  Solution status   feasible
                    188182.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[dcmulti--random_seed=1]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=1" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/dcmulti.mps")
set_tests_properties([==[dcmulti--random_seed=1]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      188182.*
  Dual bound        188182.*
  Solution status   feasible
                    188182.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[dcmulti--random_seed=2]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=2" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/dcmulti.mps")
set_tests_properties([==[dcmulti--random_seed=2]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      188182.*
  Dual bound        188182.*
  Solution status   feasible
                    188182.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
add_test([==[dcmulti--random_seed=3]==] "/home/yzhou/Github/HiGHS/build_release64/bin/highs" "--random_seed=3" "--options_file" "/home/yzhou/Github/HiGHS/build_release64/testoptions.txt" "/home/yzhou/Github/HiGHS/check/instances/dcmulti.mps")
set_tests_properties([==[dcmulti--random_seed=3]==] PROPERTIES  DEPENDS "unit_tests_all" FAIL_REGULAR_EXPRESSION "Solution status   infeasible" PASS_REGULAR_EXPRESSION "Status            Optimal
  Primal bound      188182.*
  Dual bound        188182.*
  Solution status   feasible
                    188182.* \\(objective\\)" _BACKTRACE_TRIPLES "/home/yzhou/Github/HiGHS/check/CMakeLists.txt;466;add_test;/home/yzhou/Github/HiGHS/check/CMakeLists.txt;0;")
