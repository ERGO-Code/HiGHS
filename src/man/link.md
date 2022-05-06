There are several ways the HiGHS library can be used within another C++ project. 

## Use HiGHS from another CMake Project
To use the library from a CMake project use

`find_package(HiGHS)`
`find_package(Threads)`

and add the correct path to HIGHS_DIR.

### Compiling and linking without CMake

An executable defined in the file `use_highs.cpp` is linked with the HiGHS library as follows. Make sure HiGHS is installed in `install folder` following todo: add link to install. Afterwards, compile and run with

`g++ -o use_highs use_highs.cpp -I install_folder/include/ -L install_folder/lib/ -lhighs`

`LD_LIBRARY_PATH=install_folder/lib/ ./use_highs`