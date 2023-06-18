project('highspy', 'cpp',
       version: '0.0.13',
  default_options : ['warning_level=2', 'cpp_std=c++17'])

host_system = host_machine.system()

# Add C++ compiler options
cpp_args = ['-Wall', '-Wextra', '-Wpedantic']
_args = [cpp_args] # Extra arguments
_deps = [] # Dependencies
_linkto = [] # All the sub-libraries (potentials)

# Languages
add_languages('c', required: true)
cppc = meson.get_compiler('cpp')

# ---------------------- Library Dependencies
highs_dep = dependency('highs')
_deps += highs_dep
# _linkto += highs_lib

# ---------------------- Bindings
py_mod = import('python')
py = py_mod.find_installation(pure: false)


py = import('python').find_installation()
python_dep = py.dependency()
pyb11f_deps = [ python_dep ]


_deps += [dependency('python3'),
          dependency('pybind11'),
          declare_dependency(link_args: '-lstdc++') ]

# highs_bindings, for the bindings to potlib
py.extension_module(
  'highspy',
  sources : [
    'highspy/highs_bindings.cpp',
  ],
  dependencies: _deps,
  # link_with: _linkto,
  cpp_args : _args,
  install: true,
)
