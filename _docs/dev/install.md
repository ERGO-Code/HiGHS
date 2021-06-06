---
title: Install HiGHS
permalink: /install
---

### Install HiGHS

Running `make` from the `build` folder generates the following by default

* `build/bin/highs` executable 
* `build/lib/libhighs.so` shared library.

#### Linux

To install HiGHS in a directory called `highs_install_dir/` use `-DCMAKE_INSTALL_PREFIX=/path_to/highs_install_dir` when generating the build files. From the `build/` directory, run

``` bash
  cmake -DCMAKE_INSTALL_PREFIX=/path_to/highs_install_dir ..
  make
  make install
```

If no CMake installation prefix was specified HiGHS is installed in the default user application installation directories. The default installation path on Ubuntu 18.04
To install HiGHS at the default location 
`make install` 

HiGHS installation directory structure

* `bin/` containing the executable
* `lib/` containing the shared library
* `include/` header files 

