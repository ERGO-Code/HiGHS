py_mod = import('python')
py = py_mod.find_installation(pure: false)
pyb11_dep = [
  py.dependency(),
  dependency('pybind11')
]

highspyext = py.extension_module(
  '_highs',
  sources : highspy_cpp,
  dependencies: [pyb11_dep, highs_dep],
  cpp_args: _args,
  install: true,
  subdir: 'highspy',
  include_directories: _incdirs,
)

py.install_sources(
  highspy_py,
  subdir: 'highspy',
)
