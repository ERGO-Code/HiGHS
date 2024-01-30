from setuptools import setup, find_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext

import os
import sys
import sysconfig

# def path_to_build_folder():
#     """Returns the name of a distutils build directory"""
#     f = "{dirname}.{platform}-{version[0]}.{version[1]}"
#     dir_name = f.format(dirname='lib',
#                         platform=sysconfig.get_platform(),
#                         version=sys.version_info)
#     return os.path.join('build', dir_name, 'highspy')


# def pick_library():
#     my_system = platform.system()
#     if my_system == 'Linux':
#         return "highs_bindings"
#     if my_system == 'Darwin':
#         return "highs_bindings"
#     if my_system == 'Windows':
#         return "highs_bindings"
#     raise ValueError("Unknown platform: " + my_system)


# def get_extra_link_args():
#     if platform.system() == 'Windows':
#         return []
#     else:
#         return ["-Wl,-rpath=$ORIGIN/lib/."]


ext_modules = [
    
    Pybind11Extension(
        "highspy.highs_bindings",
        sources=["highspy/highs_bindings.cpp"],
        language='c++',
        include_dirs=['include/highs'],
        library_dirs=['lib', 'bin'],
        # library_dirs=[os.path.join(path_to_build_folder(), 'lib')],
        libraries=['highs'],
    ),
]


# native_module = Extension(
#     name='highspy.highs_bindinds',
#     sources=["highspy/highs_bindings.cpp"],
#     libraries=['highs_bindings', 'highs'],
#     include_dirs=['highspy/include/highs',
#                   'highspy/include/pybind11',
#                   'highspy/include'],
#     library_dirs=['highspy/lib'],
#     # extra_link_args=get_extra_link_args(),
#     extra_compile_args=['-std=c++11'],
# )

kwargs = {
    'name': 'highspy',
    # 'version': '1.6.0.dev5',
    'packages': find_packages(),
    'include_package_data': True,
    'package_dir': {'highspy': "highspy"},
    'package_data': {'highspy': [
        '*.so',
        '*.pyd',
        'highspy/*.so',
        # '.libs/*.so',
        'lib/*.so',
        'lib/*.dylib',
        'lib/*.lib',
        'bin/*.dll',
        'include/highs/*.h',
        'include/highs/lp_data/*.h',
        'include/highs/util/*.h',
        'include/highs/io/*.h',
        'include/highs/simplex/*.h',
        'include/highs/model/*.h',
        'include/highs/presolve/*.h',
    ]},
    # 'ext_modules':  [native_module],
    'cmdclass': {"build_ext": build_ext},
    'ext_modules': ext_modules,
}

setup(**kwargs)
