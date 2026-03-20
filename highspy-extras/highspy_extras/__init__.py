"""External dependencies extension for highspy."""

from __future__ import annotations

import ctypes
import sys
from importlib.metadata import version
from pathlib import Path

__all__ = ["__version__", "get_library_version"]

__version__ = version("highspy-extras")

_PACKAGE_DIR = Path(__file__).resolve().parent
_library_handle: ctypes.CDLL | None = None


def _load_library() -> ctypes.CDLL:
    global _library_handle

    if _library_handle is not None:
        return _library_handle

    if sys.platform == "win32":
        library_name = "highs_extras.dll"
    elif sys.platform == "darwin":
        library_name = "libhighs_extras.dylib"
    else:
        library_name = "libhighs_extras.so"

    library_path = _PACKAGE_DIR / library_name
    if not library_path.is_file():
        raise FileNotFoundError(
            f"Could not find the highs_extras shared library at {library_path}"
        )

    _library_handle = ctypes.CDLL(str(library_path))
    return _library_handle


def get_library_version() -> str:
    """Return the ABI version string exported by the highs_extras library."""

    get_version = _load_library().highs_extras_get_version
    get_version.argtypes = []
    get_version.restype = ctypes.c_char_p

    version_bytes = get_version()
    if version_bytes is None:
        raise RuntimeError("highs_extras_get_version() returned NULL")

    return version_bytes.decode("utf-8")

