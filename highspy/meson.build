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
  '_highs',
  sources : [
    'highs_bindings.cpp',
  ],
  dependencies: _deps,
  link_with: _linkto,
  include_directories: _incdirs,
  cpp_args : _args,
  install: true,
  subdir: 'highspy/',
)

# highspy, main package
py.install_sources([
    '__init__.py',
  ],
  pure: false, # install next to compiled extension
  subdir: 'highspy'
)
