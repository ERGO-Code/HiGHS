# CMake generated Testfile for 
# Source directory: /home/ivet/code/HiGHS/examples
# Build directory: /home/ivet/code/HiGHS/build-cpu/examples
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[cxx_examples_call_highs_from_cpp]=] "/home/ivet/code/HiGHS/build-cpu/bin/call_highs_from_cpp")
set_tests_properties([=[cxx_examples_call_highs_from_cpp]=] PROPERTIES  _BACKTRACE_TRIPLES "/home/ivet/code/HiGHS/cmake/cpp-highs.cmake;126;add_test;/home/ivet/code/HiGHS/examples/CMakeLists.txt;8;highs_cxx_test;/home/ivet/code/HiGHS/examples/CMakeLists.txt;0;")
add_test([=[c_examples_call_highs_from_c_minimal]=] "/home/ivet/code/HiGHS/build-cpu/bin/call_highs_from_c_minimal")
set_tests_properties([=[c_examples_call_highs_from_c_minimal]=] PROPERTIES  _BACKTRACE_TRIPLES "/home/ivet/code/HiGHS/cmake/cpp-highs.cmake;159;add_test;/home/ivet/code/HiGHS/examples/CMakeLists.txt;24;highs_c_test;/home/ivet/code/HiGHS/examples/CMakeLists.txt;0;")
