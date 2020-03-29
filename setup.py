'''Create shared library for use within scipy.'''

#from setuptools import dist
#dist.Distribution().fetch_build_eggs(['Cython>=0.29.16', 'numpy>=1.18.2'])

from distutils.core import setup
from distutils.extension import Extension
from distutils.ccompiler import new_compiler
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext as _build_ext
from Cython.Build import cythonize

from datetime import datetime
import pathlib
import sysconfig

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
    sources = [str(pathlib.Path('src/' + s).resolve()) for s in sources]
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

# Get path to shared libraries
CYTHON_DIRNAME = 'pyHiGHS'
CYTHON_DIR = pathlib.Path(__file__).resolve().parent / CYTHON_DIRNAME
HIGHS_DIR = str(CYTHON_DIR.parent)
#CYTHON_DIR = str(CYTHON_DIR)
#LIBRARY_DIRS = [CYTHON_DIR]
LIBRARY_DIRS = [str(CYTHON_DIR.parent / 'build/lib.linux-x86_64-3.6/pyHiGHS/')]

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

# Naming conventions of shared libraries differ platform to platform:
SO_PREFIX = str(pathlib.Path(new_compiler().library_filename('', lib_type='shared')).with_suffix(''))
SO_SUFFIX = str(pathlib.Path(sysconfig.get_config_var('EXT_SUFFIX')).with_suffix(''))
if SO_SUFFIX is None:
    # https://bugs.python.org/issue19555
    SO_SUFFIX = str(pathlib.Path(sysconfig.get_config_var('SO')).with_suffix(''))

# We use some modern C++, as you should. HiGHS uses C++11, no penalty for going to C++14
EXTRA_COMPILE_ARGS = ['-std=c++14']

extensions = [
    # BASICLU
    Extension(
        CYTHON_DIRNAME + '.' + SO_PREFIX + 'basiclu',
        basiclu_sources,
        include_dirs=[
            str(pathlib.Path('src/').resolve()),
            str(pathlib.Path('src/ipm/basiclu/include/').resolve()),
        ],
        language="c",
        define_macros=DEFINE_MACROS,
        undef_macros=UNDEF_MACROS,
    ),

    # IPX
    Extension(
        CYTHON_DIRNAME + '.' + SO_PREFIX + 'ipx',
        ipx_sources,
        include_dirs=[
            str(pathlib.Path('src/').resolve()),
            str(pathlib.Path('src/ipm/ipx/include/').resolve()),
            str(pathlib.Path('src/ipm/basiclu/include/').resolve()),
        ],
        language="c++",
        library_dirs=LIBRARY_DIRS,
        libraries=['basiclu' + SO_SUFFIX],
        runtime_library_dirs=LIBRARY_DIRS,
        define_macros=DEFINE_MACROS,
        undef_macros=UNDEF_MACROS,
        extra_compile_args=EXTRA_COMPILE_ARGS,
    ),

    # HiGHS
    Extension(
        CYTHON_DIRNAME + '.libhighs',
        sources,
        include_dirs=[
            str(pathlib.Path(CYTHON_DIRNAME + '/src/').resolve()),
            str(pathlib.Path('src/').resolve()),
            str(pathlib.Path('src/ipm/ipx/include/').resolve()),
            str(pathlib.Path('src/lp_data/').resolve()),
        ],
        language="c++",
        library_dirs=LIBRARY_DIRS,
        libraries=['ipx' + SO_SUFFIX],
        runtime_library_dirs=LIBRARY_DIRS,
        define_macros=DEFINE_MACROS,
        undef_macros=UNDEF_MACROS,

        # Should only be here if using openMP.
        # For Microsoft Visual C++ compiler, use '/openmp' instead of '-fopenmp'.
        extra_compile_args=EXTRA_COMPILE_ARGS + ['-fopenmp'],
        extra_link_args=['-fopenmp'],
    ),

    # Cython wrapper around RunHighs (for solving MPS files)
    Extension(
        CYTHON_DIRNAME + '.linprog_mps',
        [str(pathlib.Path(CYTHON_DIRNAME + '/src/linprog_mps.pyx').resolve())],
        include_dirs=[
            str(pathlib.Path(CYTHON_DIRNAME + '/src/').resolve()),
            str(pathlib.Path('src/').resolve()),
            str(pathlib.Path('src/ipm/ipx/include/').resolve()),
            str(pathlib.Path('src/lp_data/').resolve()),
            str(pathlib.Path('src/io/').resolve()),
            str(pathlib.Path('src/mip/').resolve()),
        ],
        language="c++",
        library_dirs=LIBRARY_DIRS,
        libraries=['highs' + SO_SUFFIX],
        runtime_library_dirs=LIBRARY_DIRS,
        define_macros=DEFINE_MACROS,
        undef_macros=UNDEF_MACROS,
        extra_compile_args=EXTRA_COMPILE_ARGS,
    ),

    # Cython wrapper for Highs_call
    Extension(
        CYTHON_DIRNAME + '.linprog',
        [str(pathlib.Path(CYTHON_DIRNAME + '/src/linprog.pyx').resolve())],
        include_dirs=[
            str(pathlib.Path(CYTHON_DIRNAME + '/src/').resolve()),
            str(pathlib.Path('src/').resolve()),
            str(pathlib.Path('src/interfaces/').resolve()),
            str(pathlib.Path('src/lp_data/').resolve()),
            #np.get_include(),
        ],
        language='c++',
        library_dirs=LIBRARY_DIRS,
        libraries=['highs' + SO_SUFFIX],
        runtime_library_dirs=LIBRARY_DIRS,
        define_macros=DEFINE_MACROS,
        undef_macros=UNDEF_MACROS,
        extra_compile_args=EXTRA_COMPILE_ARGS,
    ),
]

setup(
    name='scikit-highs',
    version='0.0.0',
    author='Nicholas McKibben',
    author_email='nicholas.bgp@gmail.com',
    packages=find_packages(),
    scripts=[],
    url='https://github.com/mckib2/HiGHS',
    description='Cython interface to HiGHS.',
    install_requires=[
        "numpy>=1.18.2",
        "scipy>=1.4.1",
        "Cython>=0.29.16",
    ],
    cmdclass={'build_ext': build_ext},
    setup_requires=['numpy', 'Cython'],
    python_requires='>=3',

    ext_modules=cythonize(extensions),
    #options={'build_ext': {'inplace': True}},
)
