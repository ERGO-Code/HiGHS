project('highspy', 'cpp',
       version: '0.0.13',
  default_options : ['warning_level=3', 'cpp_std=c++17'])

host_system = host_machine.system()
_deps = []

# Languages
add_languages('c', required: true)
cppc = meson.get_compiler('cpp')

# ---------------------- Library Dependencies
highs_dep = dependency('highs')
_deps += highs_dep

# ---------------------- Bindings
py_mod = import('python')
py = py_mod.find_installation(pure: false)


py = import('python').find_installation()
_deps += [ py.dependency(),
            dependency('pybind11'),
            declare_dependency(link_args: '-lstdc++')
         ]

# highs_bindings, for the bindings to potlib
py.extension_module(
  'highspy',
  sources : [
    'highspy/highs_bindings.cpp',
  ],
  dependencies: _deps,
  install: true,
)
