# from glob import glob
# from setuptools import setup, find_packages
# from pybind11.setup_helpers import Pybind11Extension, build_ext

# original_pybind11_setup_helpers_macos = pybind11.setup_helpers.MACOS
# pybind11.setup_helpers.MACOS = False

import os
from setuptools import setup
from setuptools.dist import Distribution

class BinaryDistribution(Distribution):
    def is_pure(self):
        return False


setup(name='highspy', 
    packages=['highspy'],
     package_data={'highspy': ['highs_bindings.so', 'libhighs.dylib']},
    include_package_data=True,
    distclass=BinaryDistribution,
)


# finally:
    # pybind11.setup_helpers.MACOS = original_pybind11_setup_helpers_macos