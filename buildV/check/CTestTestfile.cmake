# CMake generated Testfile for 
# Source directory: /home/jajhall/h2gmp/check
# Build directory: /home/jajhall/h2gmp/buildV/check
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(brandy  "/home/jajhall/h2gmp/buildV/bin/h2gmp" " " "-f" "/home/jajhall/h2gmp/check/brandy.mps")
set_tests_properties(brandy  PROPERTIES  PASS_REGULAR_EXPRESSION "Solve plain: OPTIMAL")
