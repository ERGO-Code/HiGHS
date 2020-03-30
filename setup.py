'''Create shared library for use within scipy.'''

#from distutils.core import setup
#from distutils.extension import Extension
from Cython.Build import cythonize
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext as _build_ext

import pathlib
from datetime import datetime
import numpy as np

# Define some things for the module
MODULE_NAME = 'pyHiGHS'
VERSION = '0.0.16'

# Dependencies
CYTHON_VERSION = '0.29.16'
NUMPY_VERSION = '1.18.2'
SCIPY_VERSION = '1.4.1'

class build_ext(_build_ext):
    '''Subclass build_ext to bootstrap numpy.'''
    def finalize_options(self):
        _build_ext.finalize_options(self)

        # Prevent numpy from thinking it's still in its setup process
        import numpy as np
        self.include_dirs.append(np.get_include())

# Create HConfig.h: this is usually created by cmake,
# but we just need an empty file and we'll do the
# pound defines here in setup.py
HConfig_h = pathlib.Path('src/HConfig.h')
if not HConfig_h.exists():
    HConfig_h.touch()

def get_sources(CMakeLists, start_token, end_token):
    # Read in sources from CMakeLists.txt
    with open(CMakeLists, 'r') as f:
        s = f.read()

        # Find block where sources are listed
        start_idx = s.find(start_token) + len(start_token)
        end_idx = s[start_idx:].find(end_token) + len(s[:start_idx])
        sources = s[start_idx:end_idx].split('\n')
        sources = [s.strip() for s in sources if s[0] != '#']

    # Make relative to setup.py
    sources = [str(pathlib.Path('src/' + s)) for s in sources]
    return sources

# Preprocess the highs.pyx.in to pull info from the cmake files.
# This step produces src/highs.pyx.
sources = get_sources('src/CMakeLists.txt', 'set(sources\n', ')')
basiclu_sources = get_sources('src/CMakeLists.txt', 'set(basiclu_sources\n', ')')
ipx_sources = get_sources('src/CMakeLists.txt', 'set(ipx_sources\n', ')')

# Grab some more info about HiGHS from root CMakeLists
def get_version(CMakeLists, start_token, end_token=')'):
    with open('CMakeLists.txt', 'r') as f:
        s = f.read()
        start_idx = s.find(start_token) + len(start_token) + 1
        end_idx = s[start_idx:].find(end_token) + len(s[:start_idx])
    return s[start_idx:end_idx].strip()
HIGHS_VERSION_MAJOR = get_version('CMakeLists.txt', 'HIGHS_VERSION_MAJOR')
HIGHS_VERSION_MINOR = get_version('CMakeLists.txt', 'HIGHS_VERSION_MINOR')
HIGHS_VERSION_PATCH = get_version('CMakeLists.txt', 'HIGHS_VERSION_PATCH')

# Get some directory paths
CYTHON_DIR = pathlib.Path(__file__).parent / MODULE_NAME
HIGHS_DIR = str(CYTHON_DIR.parent.resolve())

# Here are the pound defines that HConfig.h would usually provide:
TODAY_DATE = datetime.today().strftime('%Y-%m-%d')
DEFINE_MACROS = [
    ('OPENMP', None),
    ('CMAKE_BUILD_TYPE', '"Release"'),
    ('HiGHSRELEASE', None),
    ('IPX_ON', None),
    ('HIGHS_GITHASH', '"b83a87d9"'),
    ('HIGHS_COMPILATION_DATE', '"' + TODAY_DATE +'"'),
    ('HIGHS_VERSION_MAJOR', HIGHS_VERSION_MAJOR),
    ('HIGHS_VERSION_MINOR', HIGHS_VERSION_MINOR),
    ('HIGHS_VERSION_PATCH', HIGHS_VERSION_PATCH),
    ('HIGHS_DIR', '"' + HIGHS_DIR + '"'),
]
UNDEF_MACROS = [
    'EXT_PRESOLVE',
    'SCIP_DEV',
    'HiGHSDEV',
    'OSI_FOUND',
]

# We use some modern C++, as you should. HiGHS uses C++11, no penalty for going to C++14
EXTRA_COMPILE_ARGS = ['-std=c++14']

extensions = [
    Extension(
        MODULE_NAME + '.linprog',
        [
            str(pathlib.Path(MODULE_NAME + '/src/linprog.pyx'))
        ] + basiclu_sources + ipx_sources + sources,
        include_dirs=[
            str(pathlib.Path(MODULE_NAME+ '/src/')),
            str(pathlib.Path('src/ipm/basiclu/include/')),
            str(pathlib.Path('external/')),
            str(pathlib.Path('src/')),
            str(pathlib.Path('src/ipm/ipx/include/')),
            str(pathlib.Path('src/lp_data/')),
            str(pathlib.Path('src/io/')),
            str(pathlib.Path('src/mip/')),
            str(pathlib.Path('src/interfaces/')),
        ],
        language="c++",
        define_macros=DEFINE_MACROS,
        undef_macros=UNDEF_MACROS,
        extra_compile_args=EXTRA_COMPILE_ARGS + ['-fopenmp'],
        extra_link_args=['-fopenmp'],
    ),
]


setup(
    name='scikit-highs',
    version=VERSION,
    author='Nicholas McKibben',
    author_email='nicholas.bgp@gmail.com',
    packages=find_packages(),
    scripts=[],
    url='https://github.com/mckib2/HiGHS',
    license='MIT',
    description='Cython interface to HiGHS.',
    long_description=open('PYREADME.rst', encoding='utf-8').read(),
    install_requires=[
        "numpy>=" + NUMPY_VERSION,
        "scipy>=" + SCIPY_VERSION,
        "Cython>=" + CYTHON_VERSION,
    ],
    cmdclass={'build_ext': build_ext},
    setup_requires=['numpy', 'Cython'],
    python_requires='>=3',
    include_package_data=True, # include example .mps file

    ext_modules=cythonize(extensions),
)
