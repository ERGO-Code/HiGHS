This is the documentation page for the .NET wrapper of HiGHS.

## NuGet

The nuget package Highs.Native is on https://www.nuget.org, at https://www.nuget.org/packages/HiGHS/.

It can be added to your C# project with `dotnet`

```shell
dotnet add package Highs.Native --version 1.16.0
```

The nuget package contains runtime libraries for

* `win-x64`
* `win-x86`
* `win-arm64`
* `linux-x64`
* `linux-arm64`
* `macos-x64`
* `macos-arm64`

#### HiPO solver

The HiPO interior point solver needs the optional `Highs.Native.Extras` package, which contains third-party dependencies under the Apache-2.0 and BSD-3-Clause licenses. It is the counterpart of `highspy-extras` for Python. Add it with the same version as `Highs.Native`:

```shell
dotnet add package Highs.Native.Extras --version 1.16.0
```

Then select HiPO with `solver.setStringOptionValue("solver", "hipo")`. Without the package, that call returns `HighsStatus.kError`.

#### Local build

To build the wrapper locally, you would need `cmake` and `dotnet`. CMake can be configured to generate the files required for the dotnet package locally, with the `BUILD_DOTNET` cmake variable. Assuming the build directory is called `build`, the package is generated in `build/dotnet/Highs.Native`, with a single runtime library, depending on the platform. From the HiGHS root directory, run

``` bash
cmake -S. -Bbuild -DCSHARP=ON -DBUILD_DOTNET=ON
```

To also generate `build/dotnet/Highs.Native.Extras`, add `-DHIPO=ON -DBUILD_OPENBLAS=ON`.

Then, from `build/dotnet/Highs.Native` (and `build/dotnet/Highs.Native.Extras`), run

```shell
dotnet pack -c Release /p:Version=$version
```

At the moment version is set manually.

Many thanks to Markus-Ludwig Wermer and the team at MENTZ GmbH for the NuGet packaging support!
