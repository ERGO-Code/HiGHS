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

### Move the files to the runtimes subfolder of the nuget folder
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


