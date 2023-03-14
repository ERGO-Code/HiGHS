Precompiled binaries are available for a variety of
platforms at the [JuliaBinaryWrappers HiGHS
repository](https://github.com/JuliaBinaryWrappers/HiGHS_jll.jl/releases). Each
includes library files for linking to external projects, and a
stand-alone executable.

**Installation instructions**

To install, download the appropriate file and extract the executable located at `/bin/highs`.

* For Windows users: if in doubt, choose the file ending in `x86_64-w64-mingw32-cxx11.tar.gz`
* For M1 macOS users: choose the file ending in `aarch64-apple-darwin.tar.gz`
* For Intel macOS users: choose the file ending in `x86_64-apple-darwin.tar.gz`

 * These files link against `libstdc++`. If you do not have one installed, download the platform-specific libraries from the [JuliaBinaryWrappers CompilerSupportLibraries repository](https://github.com/JuliaBinaryWrappers/CompilerSupportLibraries_jll.jl/releases/tag/CompilerSupportLibraries-v0.5.1%2B0) and copy all the libraries into the same folder as the `highs` executable.
