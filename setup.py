from setuptools import setup, find_packages, Extension
import pybind11.setup_helpers
import os
import numpy
from pybind11.setup_helpers import Pybind11Extension

# original_pybind11_setup_helpers_macos = pybind11.setup_helpers.MACOS
# pybind11.setup_helpers.MACOS = False

# extensions.append(Pybind11Extension('highspy.highs_bindings',
#                                     sources=['highspy/highs_bindings.cpp'],
#                                     language='c++',
#                                     include_dirs=[highs_include_dir],
#                                     library_dirs=[highs_lib_dir],
#                                     libraries=['highs']))

setup(name='highspy',
      packages=find_packages(),
      include_package_data=True,
    #   package_data={ # 'highspy.highs': ['highs*.so'],
    #     'highspy.hig'highspy.highs': ['highs*.so'],hs_bindings': ['highs_bindings*.so']}
      package_data={ 'highspy.highs': ['build/lib/libhighs.so'],
       'highspy.highs_bindings': ['build/lib/highs_bindings*.so'] }
      )

# finally:
    # pybind11.setup_helpers.MACOS = original_pybind11_setup_helpers_macos