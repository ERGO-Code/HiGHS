# !!! note
#     This file is a version of the file we use to package HiGHS for the Julia
#     ecosystem. If you make changes to this file during the development of
#     HiGHS, please tag `@odow` so we can make the correponding changes to:
#     https://github.com/JuliaPackaging/Yggdrasil/blob/master/H/HiGHS

using BinaryBuilder, Pkg

name = "HiGHS"
version = VersionNumber(ENV["HIGHS_RELEASE"])

sources = [GitSource(ENV["HIGHS_URL"], ENV["HIGHS_COMMIT"])]

script = raw"""
export BUILD_SHARED="ON"
export BUILD_STATIC="OFF"

cd $WORKSPACE/srcdir/HiGHS

# Remove system CMake to use the jll version
apk del cmake

mkdir -p build
cd build

# Do fully static build only on Windows
if [[ "${BUILD_SHARED}" == "OFF" ]] && [[ "${target}" == *-mingw* ]]; then
    export CXXFLAGS="-static"
fi

cmake -DCMAKE_INSTALL_PREFIX=${prefix} \
    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TARGET_TOOLCHAIN} \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_SHARED_LIBS=${BUILD_SHARED} \
    -DZLIB_USE_STATIC_LIBS=${BUILD_STATIC} \
    -DFAST_BUILD=ON ..

if [[ "${target}" == *-linux-* ]]; then
        make -j ${nproc}
else
    if [[ "${target}" == *-mingw* ]]; then
        cmake --build . --config Release
    else
        cmake --build . --config Release --parallel
    fi
fi
make install

install_license ../LICENSE.txt
"""

products = [
    LibraryProduct("libhighs", :libhighs),
    ExecutableProduct("highs", :highs),
]

platforms = supported_platforms()
platforms = expand_cxxstring_abis(platforms)

dependencies = [
    Dependency("CompilerSupportLibraries_jll"),
    Dependency("Zlib_jll"),
    HostBuildDependency(PackageSpec(; name="CMake_jll")),
]

build_tarballs(
    ARGS,
    name,
    version,
    sources,
    script,
    platforms,
    products,
    dependencies;
    preferred_gcc_version = v"6",
    julia_compat = "1.6",
)
