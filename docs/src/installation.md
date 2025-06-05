# Install HiGHS

## Compile from source 

HiGHS uses CMake as build system, and requires at least version
3.15. Details about building from source using CMake can be found in `HiGHS/cmake/README.md`.

## Install via a package manager

HiGHS can be installed using a package manager in the cases of
[`Julia`](@ref HiGHS.jl), [`Python`](@ref python-getting-started), [`CSharp`](@ref nuget) and [`Rust`](@ref Rust).

## Precompiled Binaries

_These binaries are provided by the Julia community and are not officially
supported by the HiGHS development team. If you have trouble using these
libraries, please open a GitHub issue and tag `@odow` in your question._

Precompiled static executables are available for a variety of platforms at

 * [https://github.com/JuliaBinaryWrappers/HiGHSstatic_jll.jl/releases](https://github.com/JuliaBinaryWrappers/HiGHSstatic_jll.jl/releases)

Multiple versions are available. Each version has the form `vX.Y.Z`. In
general, you should choose the most recent version.

To install a precompiled binary, download the appropriate `HiGHSstatic.vX.Y.Z.[platform-string].tar.gz`
file and extract the executable located at `/bin/highs`.

Do not download the file starting with `HiGHSstatic-logs`. These files contain
information from the automated compilation system. Click "Show all N assets"
to see more files.

### Platform strings

The GitHub releases contain precompiled binaries for a number of different
platforms. These are indicated by the platform-specific string in each
filename.

 * For Windows users: choose the file ending in `x86_64-w64-mingw32-cxx11.tar.gz`
 * For M1 macOS users: choose the file ending in `aarch64-apple-darwin.tar.gz`
 * For Intel macOS users: choose the file ending in `x86_64-apple-darwin.tar.gz`
