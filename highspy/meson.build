py_mod = import('python')
py = py_mod.find_installation(pure: false)
pyb11_dep = [
  py.dependency(),
  dependency('pybind11')
]
_deps += pyb11_dep

highs_dep = declare_dependency(link_with: _linkto,
                               dependencies: _deps)
highspyext = py.extension_module(
  'highspy',
  sources : [
    'highs_bindings.cpp',
  ],
  dependencies: highs_dep,
  include_directories: _incdirs,
  cpp_args : _args,
  install: true,
)
