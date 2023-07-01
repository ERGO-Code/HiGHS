project('highs', 'cpp',
  version : '1.5.3',
  default_options : ['warning_level=3', 'cpp_std=c++17'])


# Add C++ compiler options
_args = [] # Extra arguments
_deps = [] # Dependencies
_linkto = [] # All the sub-libraries
_incdirs = [] # All the includes

add_languages('c', required: true)
cc = meson.get_compiler('c')
cppc = meson.get_compiler('cpp')

# Platform detection
host_system = host_machine.system()
is_windows = host_system == 'windows'
is_mingw = is_windows and cc.get_id() == 'gcc'

# Conditional arguments
if host_system == 'linux'
  _args += '-Wno-return-type'
  _args += '-Wno-switch'
  _args += '-Wno-unused-variable'
  _args += '-Wno-unused-const-variable'
endif

# --------------------- Dependencies
# libm for Unix systems
m_dep = cppc.find_library('m', required: false)
_deps += m_dep
# For building with clang
_deps += [declare_dependency(link_args: '-lstdc++')]

# Required
threads_dep = dependency('threads',
                        required: true)
_deps += threads_dep

# Optional
zlib_dep = dependency('zlib',
                      required: get_option('use_zlib'))
if zlib_dep.found()
  _deps += zlib_dep
  _incdirs += include_directories(['extern/zstr'])
endif

# Include Sources
_incdirs += include_directories([
  'extern/',
  'extern/pdqsort/',
  'extern/filereaderlp/',
])

# Optional arguments
if get_option('fast_build')
  add_project_arguments(cc.get_supported_arguments('-fno-omit-frame-pointer'), language : 'c')
  add_project_arguments(cppc.get_supported_arguments('-fno-omit-frame-pointer'), language : 'cpp')
endif

if get_option('debug_sol')
  _args += ['-DHIGHS_DEBUGSOL']
endif

# --------------------- Library

subdir('src') # defines highslib
_linkto += highslib
