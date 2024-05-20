# [Getting started](@id csharp-getting-started)

## Install 

### Build from source

There is a C# example code in `examples/call_highs_from_csharp.cs`. From the HiGHS root directory, run 

``` bash
cmake -S. -Bbuild -DCSHARP=ON
```

If a CSharp compiler is available, this builds the example using cmake and generates a binary in the build directory (`build/bin/csharpexample`).

### NuGet

The nuget package Highs.Native is on https://www.nuget.org, at https://www.nuget.org/packages/Highs.Native/. 

It can be added to your C# project with `dotnet`

```bash
dotnet add package Highs.Native --version 1.7.0
```

The nuget package contains runtime libraries for 

* `win-x64`
* `win-x32`
* `linux-x64`
* `linux-arm64`
* `macos-x64`
* `macos-arm64`

Details for building locally can be found in `nuget/README.md`.

