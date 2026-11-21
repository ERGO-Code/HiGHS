# Third Party Licenses and Acknowledgements

The majority of the HiGHS source code is available under the [MIT license](https://opensource.org/license/MIT).

Third-party code is located in the `/extern` directory and parts of the `/highs`
directory. The corresponding licenses and acknowledgements are listed below.

## amd

The source code in `/extern/amd` is distributed under the [BSD-3 license](https://opensource.org/license/bsd-3-clause)
at `/extern/amd/License.txt`.

It was originally developed by Timothy Davis.

The upstream source code is available at:

 * https://github.com/DrTimothyAldenDavis/SuiteSparse

To avoid compiling this code into HiGHS, use `-DHIPO=OFF` (the default value).

## catch

The source code in `/extern/catch/catch.hpp` is Catch v2.13.8, distributed under
the [Boost Software License 1.0](https://www.boost.org/LICENSE_1_0.txt) at
`/extern/catch/LICENCE_1_0.txt`.

Copyright (c) 2022 Two Blue Cubes Ltd. All rights reserved.

The upstream source code is available at:

 * https://github.com/catchorg/Catch2/tree/v2.13.8

## cli

The source code in `/extern/cli11/CLI11.hpp` is distributed under the
[BSD-3-Clause license](https://opensource.org/license/bsd-3-clause), whose full
text and copyright notice are included in the header.

The upstream source code is available at:

* https://github.com/CLIUtils/CLI11

CLI11 is only used to parse command line input for the HiGHS executable, so does not affect the license status of the HiGHS library or language interfaces to it.

## cuPDLP-C

The source code in `/highs/pdlp/cupdlp` was adapted from cuPDLP-C, whose
[MIT license](https://opensource.org/license/MIT) is reproduced at
`/highs/pdlp/cupdlp/LICENSE`.

Copyright (c) 2023 COPT-Public.

The upstream source code for the revision imported in March 2025 is available at:

 * https://github.com/COPT-Public/cuPDLP-C/tree/7b94c41efca42fa176a129a935b59718e1455fad

## metis

The source code in `/extern/metis` is distributed under the [Apache 2.0 license](https://opensource.org/license/apache-2-0)
at `/extern/metis/LICENSE.txt`.

It was originally developed by George Karypis.

The upstream source code is available at:

 * https://github.com/KarypisLab/METIS
 * https://github.com/KarypisLab/GKlib

To avoid compiling this code into HiGHS, use `-DHIPO=OFF` (the default value).

## OpenBLAS

With `-DBUILD_OPENBLAS=ON`, OpenBLAS is downloaded at build time and linked into
the HiPO library (`highs_extras`), which the `Highs.Native.Extras` NuGet package
distributes on Linux and Windows. It is distributed under the
[BSD-3 license](https://opensource.org/license/bsd-3-clause):

```
Copyright (c) 2011-2014, The OpenBLAS Project
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

   1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the
      distribution.
   3. Neither the name of the OpenBLAS project nor the names of
      its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

The upstream source code is available at:

 * https://github.com/OpenMathLib/OpenBLAS

## pdqsort

The source code in `/extern/pdqsort` is distributed under the [zlib license](https://opensource.org/license/zlib)
at `/extern/pdqsort/license.txt`.

It was originally developed by Orson Peters.

The upstream source code is available at:

 * https://github.com/orlp/pdqsort

## rcm

The source code in `/extern/rcm` is distributed under the [MIT license](https://opensource.org/license/MIT)
at `/extern/rcm/LICENSE`.

It was originally developed by Alan George, Joseph Liu, and John Burkardt.

The upstream source code is available at:

 * https://people.sc.fsu.edu/~jburkardt/cpp_src/rcm/rcm.html

To avoid compiling this code into HiGHS, use `-DHIPO=OFF` (the default value).

## zstr

The source code in `/extern/zstr` is distributed under the [MIT license](https://opensource.org/license/MIT)
at `/extern/zstr/LICENSE`.

It was originally developed by Matei David.

The upstream source code is available at:

 * https://github.com/mateidavid/zstr
