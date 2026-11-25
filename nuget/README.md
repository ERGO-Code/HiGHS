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
* `linux-arm`
* `macos-x64`
* `macos-arm64`

#### HiPO solver

The HiPO interior point solver needs the optional `Highs.Native.Extras` package, which contains third-party dependencies under the Apache-2.0 and BSD-3-Clause licenses. It is the counterpart of `highspy-extras` for Python. It depends on the `Highs.Native` package of the same version, so adding it is enough:

```shell
dotnet add package Highs.Native.Extras
```

Then select HiPO with `solver.setStringOptionValue("solver", "hipo")`. Without the package, that call returns `HighsStatus.kError`. `Highs.Native.Extras` has no `linux-arm` runtime, so HiPO is unavailable there.

#### Local build

To build the wrapper locally, you would need `cmake` and `dotnet`. CMake can be configured to generate the files required for the dotnet package locally, with the `BUILD_DOTNET` cmake variable. Assuming the build directory is called `build`, the package is generated in `build/dotnet/Highs.Native`, with a single runtime library, depending on the platform. From the HiGHS root directory, run

``` bash
cmake -S. -Bbuild -DCSHARP=ON -DBUILD_DOTNET=ON
```

To also generate `build/dotnet/Highs.Native.Extras`, add `-DHIPO=ON -DBUILD_OPENBLAS=ON`.

Then, from `build/dotnet/Highs.Native`, run

```shell
dotnet pack -c Release /p:Version=$version -o $nupkgs
```

`Highs.Native.Extras` depends on that `Highs.Native` package, so pack it from `build/dotnet/Highs.Native.Extras` with the same output folder as the package source:

```shell
dotnet pack -c Release /p:Version=$version -p:RestoreSources=$nupkgs -o $nupkgs
```

At the moment version is set manually.

Many thanks to Markus-Ludwig Wermer and the team at MENTZ GmbH for the NuGet packaging support!
