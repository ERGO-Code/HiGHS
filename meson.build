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

if cppc.get_id() == 'msvc'
  add_project_arguments(cppc.get_supported_arguments('/W4'), language : 'cpp')
  add_project_arguments('/wd4018', '/wd4061', '/wd4100', '/wd4101', '/wd4127', '/wd4189', '/wd4244', '/wd4245', '/wd4267', '/wd4324', '/wd4365', '/wd4389', '/wd4456', '/wd4457', '/wd4458', '/wd4459', '/wd4514', '/wd4701', '/wd4820', language: 'cpp')
  add_project_arguments('/MP', language: 'cpp')
  _args += '-D_CRT_SECURE_NO_WARNINGS'
endif

cpu_family = host_machine.cpu_family()

if cpu_family in ['x86_64', 'i686']
    add_project_arguments('-mpopcnt', language: 'cpp')
endif
if cpu_family in ['ppc64', 'powerpc64'] and not meson.is_cross_build()
    add_project_arguments('-mpopcntd', language: 'cpp')
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
  add_project_arguments(cc.get_supported_arguments('-fno-omit-frame-pointer'), language : ['c', 'cpp'])
endif

if get_option('debug_sol')
  _args += ['-DHIGHS_DEBUGSOL']
endif

if cppc.get_id() == 'msvc'
  if get_option('with_stdcall')
    add_project_arguments(cppc.get_supported_arguments('/Gz'), language : ['cpp', 'c'])
  endif
endif

# --------------------- Library

subdir('src') # defines highslib
_linkto += highslib
