from glob import glob
from setuptools import setup, find_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext

# original_pybind11_setup_helpers_macos = pybind11.setup_helpers.MACOS
# pybind11.setup_helpers.MACOS = False

ext_modules = [
    Pybind11Extension(
        "highspy.highs_bindings",
         sorted(glob("highspy/*.cpp")),
         ),
]

# how to get setup to recognise highs_bindings.so or highs_bindings.dll
# setup(..., cmdclass={"build_ext": build_ext}, ext_modules=ext_modules)

setup(name='highspy', 
      cmdclass={"build_ext": build_ext}, 
      ext_modules=ext_modules)

      # packages=find_packages(),
      # include_package_data=True,
      # package_data={ # 'highspy.highs': ['highs*.so'],
      #     'highspy.highs_bindinfs': ['highs_bindings*.so']}

# finally:
    # pybind11.setup_helpers.MACOS = original_pybind11_setup_helpers_macos