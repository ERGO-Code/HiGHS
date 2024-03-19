How To for nuget
===

This document describes one way of how to create a nuget package for .NET containing the highs solver for different runtimes.

Steps
---

### Build shared libraries
For all desired targets, build a **shared library** of the c++ highs solver, e.g.
  - win-x64
  - linux-x64
  - linux-arm64

### Move the fils to the nuget folder
Copy the resulting libraries into the **runtimes** sub folder
  - nuget/runtimes/win-x64/highs.dll
  - nuget/runtimes/linux-x64/libhighs.so
  - nuget/runtimes/linux-arm64/libhighs.so

### Create the package
Run the dotet command to build and pack the package in the nuget directory. The required csproj. file is already present in the nuget folder
- `dotnet pack -c Release /p:Version=$version`

### Checking the package
In order to check if the runtimes are contained in the nuget package, one can open the nupkg file in a tool like 7zip and check the runtimes folder.

## nuget structure for native libraries
The nuget package is required to look like this for the native libraries
```
package/
|-- lib/
|   |-- netstandard2.0/
|       |-- highs_csharp_api.dll
|-- runtimes/
|   |-- linux-x64/
|   |   |-- native/
|   |       |-- [linux-x64 native libraries]
|   |-- linux-arm64/
|   |   |-- native/
|   |       |-- [linux-arm64 native libraries]
|   |-- win-x64/
|   |   |-- native/
|   |       |-- [win-x64 native libraries]
```

## Examples for the builds
These are examples, how one might run the builds on an ubuntu system, e.g. as a wsl
Compiler flags should be adjusted to provide the best performance
### linux-x64
This should run on a linux system
```shell
mkdir build_linux
cd build_linux
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DCMAKE_CXX_FLAGS="-O3 -march=native" -DCMAKE_C_FLAGS="-O3 -march=native"
make
```

### linux-arm64
This might run on a linux-x64 system (or use a linux-arm64 system and skipp the toolchain part)
```shell
mkdir build_arm
cd build_arm
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DCMAKE_TOOLCHAIN_FILE=../arm-toolchain.cmake
make

```

## win-x64
This should run on a windows system
```shell
mkdir build_windows
cd build_windows
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DCMAKE_SYSTEM_NAME=Windows -A x64
make

```
