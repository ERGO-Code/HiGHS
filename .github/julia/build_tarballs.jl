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
cd $WORKSPACE/srcdir/HiGHS

# Remove system CMake to use the jll version
apk del cmake

rm -rf build
mkdir build

if [[ "${target}" == *-mingw* ]]; then
    LBT=blastrampoline-5
else
    LBT=blastrampoline
fi

cmake -S . -B build \
    -DCMAKE_INSTALL_PREFIX=${prefix} \
    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TARGET_TOOLCHAIN} \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_SHARED_LIBS=ON \
    -DHIPO=ON \
    -DBUILD_SHARED_EXTRAS_LIB=OFF \
    -DBLA_VENDOR=libblastrampoline

if [[ "${target}" == *-linux-* ]]; then
    make -C build -j ${nproc}
else
    if [[ "${target}" == *-mingw* ]]; then
        cmake --build build --config Release
    else
        cmake --build build --config Release --parallel
    fi
fi
cmake --install build

install_license LICENSE.txt
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
    Dependency("libblastrampoline_jll"),
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
    preferred_gcc_version = v"11",
    julia_compat = "1.10",
)
