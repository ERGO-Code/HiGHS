#!/bin/bash

sudo apt update
sudo apt upgrade -y
sudo apt install build-essential -y
sudo apt install cmake -y
sudo apt install g++-aarch64-linux-gnu gcc-aarch64-linux-gnu -y

mkdir build_arm
cd build_arm
cmake ../.. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON  -DCMAKE_TOOLCHAIN_FILE=../arm-toolchain.cmake
make
make -j8