"""Write the vcpkg clang-cl triplet documented in docs/src/installation.md.

The two cmake code blocks in the section "Windows: OpenBLAS from vcpkg, built
with clang-cl" are the triplet and the toolchain file. They are written as
documented, except for the changes that the documentation describes for other
architectures, linkages and scopes:

    python vcpkg_clang_cl_triplet.py <x64|arm64> <dynamic|static> <all-ports|openblas-only> <out dir>

Fails if the documentation no longer contains the text to be changed.
"""

import pathlib
import re
import sys

arch, linkage, scope, out_dir = sys.argv[1:5]

doc = pathlib.Path("docs/src/installation.md").read_text(encoding="utf-8")
section = doc[doc.index("##### Windows: OpenBLAS from vcpkg, built with clang-cl"):]
triplet, toolchain = re.findall(r"```cmake\n(.*?)```", section, re.S)[:2]


def replace(text, old, new):
    if old not in text:
        sys.exit(f"Not found in the documented triplet: {old!r}")
    return text.replace(old, new)


# Architecture
triplet = replace(triplet, "set(VCPKG_TARGET_ARCHITECTURE x64)",
                  f"set(VCPKG_TARGET_ARCHITECTURE {arch})")

# Linkage: the documented triplet builds DLLs
if linkage == "static":
    triplet = replace(triplet, "set(VCPKG_CRT_LINKAGE dynamic)", "set(VCPKG_CRT_LINKAGE static)")
    triplet = replace(triplet, "set(VCPKG_LIBRARY_LINKAGE dynamic)", "set(VCPKG_LIBRARY_LINKAGE static)")

# Scope: the documented triplet builds every port with clang-cl. For OpenBLAS
# only, move the toolchain and vcvars settings (with their comment) into the
# OpenBLAS block.
if scope == "openblas-only":
    chainload = 'set(VCPKG_CHAINLOAD_TOOLCHAIN_FILE "${CMAKE_CURRENT_LIST_DIR}/clang-cl-toolchain.cmake")\n'
    vcvars = re.search(r"(?:#.*\n)*set\(VCPKG_LOAD_VCVARS_ENV ON\)\n", triplet).group(0)
    triplet = replace(triplet, chainload, "")
    triplet = replace(triplet, vcvars, "")
    moved = "".join(f"    {line}\n" for line in (chainload + vcvars).splitlines())
    triplet = replace(triplet, 'if(PORT STREQUAL "openblas")\n',
                      'if(PORT STREQUAL "openblas")\n' + moved)

out = pathlib.Path(out_dir)
out.mkdir(parents=True, exist_ok=True)
files = {f"{arch}-windows-clang.cmake": triplet, "clang-cl-toolchain.cmake": toolchain}
for name, text in files.items():
    (out / name).write_text(text)
    print(f"==> {out / name}\n{text}")
