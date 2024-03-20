#!/bin/bash
sudo apt update
sudo apt upgrade -y
sudo apt install build-essential -y
sudo apt install cmake -y

mkdir build_linux
cd build_linux
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DCMAKE_CXX_FLAGS="-O3 -march=native" -DCMAKE_C_FLAGS="-O3 -march=native"
make -j8