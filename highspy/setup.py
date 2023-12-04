import platform
# from setuptools import setup, Extension, find_packages, 
from setuptools import setup, find_packages 
from pybind11.setup_helpers import Pybind11Extension, build_ext

import os
import sys
import sysconfig

def path_to_build_folder():
    """Returns the name of a distutils build directory"""
    f = "{dirname}.{platform}-{version[0]}.{version[1]}"
    dir_name = f.format(dirname='lib',
                        platform=sysconfig.get_platform(),
                        version=sys.version_info)
    return os.path.join('build', dir_name, 'grumbo')

def pick_library():
    my_system = platform.system()
    if my_system == 'Linux':
        return "highs_bindings"
    if my_system == 'Darwin':
        return "highs_bindings"
    if my_system == 'Windows':
        return "highs_bindings"
    raise ValueError("Unknown platform: " + my_system)


# def get_extra_link_args():
#     if platform.system() == 'Windows':
#         return []
#     else:
#         return ["-Wl,-rpath=$ORIGIN/lib/."]
    
ext_modules = [
    Pybind11Extension(
        "highs_bindings",
        ["highspy/highs_bindings.cpp"],
        # include_dirs=['highspy/include']
        include_dirs=['highspy/include/highs'],
        library_dirs=['highspy/lib'],
        libraries=['highs'],
    ),
]


# native_module = Extension(
#     name='highspy.highs_bindinds',
#     sources=["highspy/highs_bindings.cpp"],
#     # include_dirs=[os.path.join(path_to_build_folder(), 'highspy/include')],
#     libraries=[pick_library()],
#     library_dirs=['highspy/lib'],
#     # library_dirs=[os.path.join(path_to_build_folder(), 'highspy/lib')],
#     extra_link_args=get_extra_link_args(),
    
# )

kwargs = {
    'name': 'highspy',
    'version': '1.6.0.dev4',
    # 'ext_modules':  [native_module],
    'packages': find_packages(),
    # 'package_dir': {"": "highspy"},
    # 'package_data': {'highspy': ['highspy/include/*.h', 'highspy/lib/*.so',
    #                              'highspy/lib/*.lib', 'highspy/*.dll',      # for windows
    #                              ]},
    'ext_modules' : ext_modules,
    'cmdclass' : {"build_ext": build_ext},
}

setup(**kwargs)