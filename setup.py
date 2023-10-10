from setuptools import setup, find_packages
import pybind11.setup_helpers
import os

# original_pybind11_setup_helpers_macos = pybind11.setup_helpers.MACOS
# pybind11.setup_helpers.MACOS = False

# extensions.append(Pybind11Extension('highspy.highs_bindings',
#                                     sources=['highspy/highs_bindings.cpp'],
#                                     language='c++',
#                                     include_dirs=[highs_include_dir],
#                                     library_dirs=[highs_lib_dir],
#                                     libraries=['highs']))

setup(name='highspy',
        version='1.5.3',
        packages=find_packages(),
        description='Python interface to HiGHS',
        maintainer_email='highsopt@gmail.com',
        license='MIT',
        url='https://github.com/ergo-code/highs',
        include_package_data=True,
        package_data={ 'highs': ['build/lib/highs*.so'],
         'highs_bindings': ['build/lib/highs_bindings*.so']
         },
        python_requires='>=3.7',
        classifiers=["Programming Language :: Python :: 3",
                    "Programming Language :: Python :: 3.7",
                    "Programming Language :: Python :: 3.8",
                    "Programming Language :: Python :: 3.9",
                    "Programming Language :: Python :: 3.10",
                    "Programming Language :: Python :: 3.11",
                    "License :: OSI Approved :: MIT License"]
          )
# finally:
    # pybind11.setup_helpers.MACOS = original_pybind11_setup_helpers_macos