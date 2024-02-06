import os
import re
import sys
import sysconfig
import platform
import subprocess
from pathlib import Path

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from setuptools.command.test import test as TestCommand

# Convert distutils Windows platform specifiers to CMake -A arguments
PLAT_TO_CMAKE = {
    "win32": "Win32",
    "win-amd64": "x64",
    "win-arm32": "ARM",
    "win-arm64": "ARM64",
}

# A CMakeExtension needs a sourcedir instead of a file list.
# The name must be the _single_ output extension from the CMake build.
# If you need multiple extensions, see scikit-build.


class CMakeExtension(Extension):
    def __init__(self, name: str, sourcedir: str = "") -> None:
        super().__init__(name, sources=[])
        self.sourcedir = os.fspath(Path(sourcedir).resolve())

# class CMakeExtension(Extension):
#     def __init__(self, name):
#         Extension.__init__(self, name, sources=[])


class CMakeBuild(build_ext):
    def build_extension(self, ext: CMakeExtension) -> None:
        # Must be in this form due to bug in .resolve() only fixed in Python 3.10+
        ext_fullpath = Path.cwd() / self.get_ext_fullpath(ext.name)
        extdir = ext_fullpath.parent.resolve()

        # Using this requires trailing slash for auto-detection & inclusion of
        # auxiliary "native" libs

        debug = int(os.environ.get("DEBUG", 0)) if self.debug is None else self.debug
        cfg = "Debug" if debug else "Release"

        # CMake lets you override the generator - we need to check this.
        # Can be set with Conda-Build, for example.
        cmake_generator = os.environ.get("CMAKE_GENERATOR", "")

        # Set Python_EXECUTABLE instead if you use PYBIND11_FINDPYTHON
        # EXAMPLE_VERSION_INFO shows you how to pass a value into the C++ code
        # from Python.
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}{os.sep}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE={cfg}",  # not used on MSVC, but no harm
            "-DPYTHON=ON"
        ]
        build_args = []
        # Adding CMake arguments set as environment variable
        # (needed e.g. to build for ARM OSx on conda-forge)
        if "CMAKE_ARGS" in os.environ:
            cmake_args += [item for item in os.environ["CMAKE_ARGS"].split(" ") if item]

        # In this example, we pass in the version to C++. You might not need to.
        cmake_args += [f"-DEXAMPLE_VERSION_INFO={self.distribution.get_version()}"]

        if self.compiler.compiler_type != "msvc":
            # Using Ninja-build since it a) is available as a wheel and b)
            # multithreads automatically. MSVC would require all variables be
            # exported for Ninja to pick it up, which is a little tricky to do.
            # Users can override the generator with CMAKE_GENERATOR in CMake
            # 3.15+.
            if not cmake_generator or cmake_generator == "Ninja":
                try:
                    import ninja

                    ninja_executable_path = Path(ninja.BIN_DIR) / "ninja"
                    cmake_args += [
                        "-GNinja",
                        f"-DCMAKE_MAKE_PROGRAM:FILEPATH={ninja_executable_path}",
                    ]
                except ImportError:
                    pass

        else:
            # Single config generators are handled "normally"
            single_config = any(x in cmake_generator for x in {"NMake", "Ninja"})

            # CMake allows an arch-in-generator style for backward compatibility
            contains_arch = any(x in cmake_generator for x in {"ARM", "Win64"})

            # Specify the arch if using MSVC generator, but only if it doesn't
            # contain a backward-compatibility arch spec already in the
            # generator name.
            if not single_config and not contains_arch:
                cmake_args += ["-A", PLAT_TO_CMAKE[self.plat_name]]

            # Multi-config generators have a different way to specify configs
            if not single_config:
                cmake_args += [
                    f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}"
                ]
                build_args += ["--config", cfg]

        if sys.platform.startswith("darwin"):
            # Cross-compile support for macOS - respect ARCHFLAGS if set
            archs = re.findall(r"-arch (\S+)", os.environ.get("ARCHFLAGS", ""))
            if archs:
                cmake_args += ["-DCMAKE_OSX_ARCHITECTURES={}".format(";".join(archs))]

        # Set CMAKE_BUILD_PARALLEL_LEVEL to control the parallel build level
        # across all generators.
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            # self.parallel is a Python 3 only way to set parallel jobs by hand
            # using -j in the build_ext call, not supported by pip or PyPA-build.
            if hasattr(self, "parallel") and self.parallel:
                # CMake 3.12+ only.
                build_args += [f"-j{self.parallel}"]

        build_temp = Path(self.build_temp) / ext.name
        if not build_temp.exists():
            build_temp.mkdir(parents=True)

        subprocess.run(
            ["cmake", ext.sourcedir, *cmake_args], cwd=build_temp, check=True
        )
        subprocess.run(
            ["cmake", "--build", ".", *build_args], cwd=build_temp, check=True
        )

# The information here can also be placed in setup.cfg - better separation of
# logic and declaration, and simpler if you include description/version in a file.
setup(
    name="highspy",
    version="1.6.0.dev8",
    author="HiGHS developers",
    author_email="highsopt@gmail.com",
    description = "A thin set of pybind11 wrappers to HiGHS",
    long_description="",
    ext_modules=[CMakeExtension("highspy")],
    cmdclass={"build_ext": CMakeBuild},
    zip_safe=False,
    # extras_require={"test": ["pytest>=6.0"]},
    python_requires=">=3.9",
    install_requires=[
       'numpy',
     ],
)

# class CMakeBuild(build_ext):
#     def run(self):
#         try:
#             out = subprocess.check_output(['cmake', '--version'])
#         except OSError:
#             raise RuntimeError(
#                 "CMake must be installed to build the following extensions: " +
#                 ", ".join(e.name for e in self.extensions))

#         build_directory = os.path.abspath(self.build_temp)

#         cmake_args = [
#             '-DPYTHON=ON'
#             '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + build_directory,
#             # '-DPYTHON_EXECUTABLE=' + sys.executable,
#         ]

#         cfg = 'Debug' if self.debug else 'Release'
#         build_args = ['--config', cfg, '--parallel']

#         # cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]

#         # Assuming Makefiles
#         # build_args += ['--', '-j2']

#         self.build_args = build_args

#         env = os.environ.copy()
#         env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(
#             env.get('CXXFLAGS', ''),
#             self.distribution.get_version())
#         if not os.path.exists(self.build_temp):
#             os.makedirs(self.build_temp)

#         # CMakeLists.txt is in the same directory as this setup.py file
#         cmake_list_dir = os.path.abspath(os.path.dirname(__file__))
#         print('-'*10, 'Running CMake prepare', '-'*40)
#         subprocess.check_call(['cmake', cmake_list_dir] + cmake_args,
#                               cwd=self.build_temp, env=env)

#         print('-'*10, 'Building extensions', '-'*40)
#         cmake_cmd = ['cmake', '--build', '.'] + self.build_args
#         subprocess.check_call(cmake_cmd,
#                               cwd=self.build_temp)

#         # Move from build temp to final position
#         for ext in self.extensions:
#             self.move_output(ext)

#     def move_output(self, ext):
#         build_temp = Path(self.build_temp).resolve()
#         dest_path = Path(self.get_ext_fullpath(ext.name)).resolve()
#         source_path = build_temp / self.get_ext_filename(ext.name)
#         dest_directory = dest_path.parents[0]
#         dest_directory.mkdir(parents=True, exist_ok=True)
#         self.copy_file(source_path, dest_path)


# ext_modules = [
#   CMakeExtension('highspy.highs'),
#   CMakeExtension('highspy.highs_bindings')
# ]

# setup(
#  # ...
#   packages=find_packages(),
#   ext_modules=ext_modules,
#   cmdclass=dict(build_ext=CMakeBuild),
#   zip_safe=False,
# )
