#### [CSharp](@id csharp)

#### Build from source

There is a C# example code in `examples/call_highs_from_csharp.cs`. From the HiGHS root directory, run 

``` bash
cmake -S. -Bbuild -DCSHARP=ON
```

If a CSharp compiler is available, this builds the example using cmake and generates a binary in the build directory (`build/bin/csharpexample`).

#### [NuGet](@id nuget)

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

#### C# API

The C# API can be called directly. Here are observations on calling the HiGHS C# API from C#:

 * The file `HiGHS/src/interfaces/highs_csharp_api.cs` contains all the PInvoke you need. 
 * Make sure, that the native HiGHS library (`highs.dll`, `libhighs.dll`,
   `libhighs.so`, ... depending on your platform) can be found at runtime. How
   to do this is platform dependent, copying it next to your C# executable
   should work in most cases. You can use msbuild for that. On linux, installing
   HiGHS system wide should work.
 * Make sure that all dependencies of the HiGHS library can be found, too. For
   example, if HiGHS was build using `Visual C++` make sure that the
   `MSVCRuntime` is installed on the machine you want to run your application
   on.
 * Depending on the name of your HiGHS library, it might be necessary to change
   the constant "highslibname". See [document](https://learn.microsoft.com/en-us/dotnet/standard/native-interop/cross-platform)
   on writing cross platform P/Invoke code if necessary.
 * Call the Methods in `highs_csharp_api.cs` and have fun with HiGHS.

This is the normal way to call plain old C from C# with the great simplification
that you don't have to write the PInvoke declarations yourself.

