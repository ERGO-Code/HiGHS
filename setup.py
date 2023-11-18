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
        # language='c++',
        # include_dir='include'
        # library_dirs='lib'
        # libraries=['highs']),
]


# setup(..., cmdclass={"build_ext": build_ext}, ext_modules=ext_modules)



setup(name='highspy', 
      cmdclass={"build_ext": build_ext}, 
      ext_modules=ext_modules)

      # packages=find_packages(),
      # include_package_data=True,
    #   package_data={ # 'highspy.highs': ['highs*.so'],
    #     'highspy.hig'highspy.highs': ['highs*.so'],hs_bindings': ['highs_bindings*.so']}
      # package_data={ 'highspy.highs': ['build/lib/libhighs.so'],
      #  'highspy.highs_bindings': ['build/lib/highs_bindings*.so'] }
      # )

# finally:
    # pybind11.setup_helpers.MACOS = original_pybind11_setup_helpers_macos