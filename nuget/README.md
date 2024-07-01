This is the documentation page for the .NET wrapper of HiGHS.

## NuGet 

The nuget package Highs.Native is on https://www.nuget.org, at https://www.nuget.org/packages/HiGHS/. 

It can be added to your C# project with `dotnet`

```shell
dotnet add package Highs.Native --version 1.7.2
```

The nuget package contains runtime libraries for 

* `win-x64`
* `win-x32`
* `linux-x64`
* `linux-arm64`
* `macos-x64`
* `macos-arm64`

#### Local build

To build the wrapper locally, you would need `cmake` and `dotnet`. CMake can be configured to generate the files required for the dotnet package locally, wtih the `BUILD_DOTNET` cmake variable. Assuming the build directory is called `build`, the package is generated in `build/dotnet/Highs.Native`, with a single runtime library, depending on the platform. From the HiGHS root directory, run 

``` bash
cmake -S. -Bbuild -DCSHARP=ON -DBUILD_DOTNET=ON
```

Then, from `build/dotnet/Highs.Native`, run 

```shell
dotnet pack -c Release /p:Version=$version
```

At the moment version is set manually.
