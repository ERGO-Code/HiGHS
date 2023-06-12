from setuptools import setup, find_packages
import pybind11.setup_helpers
from pybind11.setup_helpers import Pybind11Extension, build_ext
import os
import ctypes.util
import sys

# Simplified version of the pyomo utility
# https://github.com/Pyomo/pyomo/blob/afc9ee5346eabbdfaf8d3802746cfba95682fe91/pyomo/common/fileutils.py#LL309C1-L380C49
def find_library(libname):
    ext = (os.path.splitext(libname)[1]).lower()
    if ext:
        return ctypes.util.find_library(libname)
    if sys.platform.startswith('win'):
        ext = '.dll'
    elif sys.platform.startswith('darwin'):
        ext = '.dylib'
    else:
        ext = '.so'
    full_libname = libname
    if not full_libname.startswith('lib'):
        full_libname = 'lib' + full_libname
    full_libname += ext
    return ctypes.util.find_library(full_libname)

original_pybind11_setup_helpers_macos = pybind11.setup_helpers.MACOS
pybind11.setup_helpers.MACOS = False

try:
    highs_lib = find_library('highs')
    if highs_lib is None:
        raise RuntimeError('Could not find HiGHS library; Please make sure it is in the LD_LIBRARY_PATH environment variable')
    highs_lib_dir = os.path.dirname(highs_lib)
    highs_build_dir = os.path.dirname(highs_lib_dir)
    highs_include_dir = os.path.join(highs_build_dir, 'include', 'highs')
    if not os.path.exists(os.path.join(highs_include_dir, 'Highs.h')):
        raise RuntimeError('Could not find HiGHS include directory')
    
    extensions = list()
    extensions.append(Pybind11Extension('highspy.highs_bindings',
                                        sources=['highspy/highs_bindings.cpp'],
                                        language='c++',
                                        include_dirs=[highs_include_dir],
                                        library_dirs=[highs_lib_dir],
                                        libraries=['highs']))
    
    setup(name='highspy',
          version='1.5.3',
          packages=find_packages(),
          description='Python interface to HiGHS',
          maintainer_email='highsopt@gmail.com',
          license='MIT',
          url='https://github.com/ergo-code/highs',
          include_package_data=True,
          package_data={'highspy': ['highspy/*.so']},
          ext_modules=extensions,
          cmdclass={'build_ext': build_ext},
          python_requires='>=3.7',
          classifiers=["Programming Language :: Python :: 3",
                       "Programming Language :: Python :: 3.7",
                       "Programming Language :: Python :: 3.8",
                       "Programming Language :: Python :: 3.9",
                       "Programming Language :: Python :: 3.10",
                       "Programming Language :: Python :: 3.11",
                       "License :: OSI Approved :: MIT License"]
          )
finally:
    pybind11.setup_helpers.MACOS = original_pybind11_setup_helpers_macos
