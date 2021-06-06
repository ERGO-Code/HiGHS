---
title: Link with HiGHS Library
permalink: /link-library/
---
HiGHS is implemented in C++11.
It is build with CMake and can be easily linked with other CMake projects.

Call HiGHS library from another CMake project:

Example.cpp
Makefile
- build/
- install/
# Linking

## From another CMake project
todo:
To use the library from a CMake project use

`find_package(HiGHS)` 

and add the correct path to HIGHS_DIR.


#### build/ folder

#### cmake/ folder

#### install/ folder

## From another C++ code using g++

Call HiGHS library from another C++ source:

The default library type is dynamic. For static libraries, see [install.md].
Comments in CMake file?

CMake project
Example.cpp
- cmake/
- build/
- install/

PATH=install_folder/lib/ ./use_highs` 
Make sure the library is in `build/lib/libhighs.so`. todo: versions, several libs generated currently.

todo:
An executable defined in the file `use_highs.cpp` is linked with the HiGHS library as follows. After running the code above, compile and run with

`g++ -o use_highs use_highs.cpp -I install_folder/include/ -L install_folder/lib/ -lhighs` 

`LD_LIBRARY_