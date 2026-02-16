"""HiPO IPM solver extension for highspy."""

import os
import sys
from pathlib import Path

# Version is kept in sync with HiGHS version in Version.txt
# This is updated during the build/release process
__version__ = "1.12.0"


def _get_library_name() -> str:
    """Get the platform-specific library filename."""
    if sys.platform == "win32":
        return "highs_hipo.dll"
    elif sys.platform == "darwin":
        return "libhighs_hipo.dylib"
    else:
        return "libhighs_hipo.so"


def get_library_path() -> Path:
    """
    Get the path to the highs_hipo shared library.
    
    Returns:
        Path to the shared library file.
    
    Raises:
        FileNotFoundError: If the library is not found.
    """
    package_dir = Path(__file__).parent
    lib_name = _get_library_name()
    lib_path = package_dir / lib_name
    
    if not lib_path.exists():
        raise FileNotFoundError(f"HiPO library not found at {lib_path}")
    
    return lib_path


def is_available() -> bool:
    """
    Check if the HiPO library is available.
    
    Returns:
        True if the library file exists.
    """
    try:
        get_library_path()
        return True
    except FileNotFoundError:
        return False


def _setup_environment():
    """
    Set up environment variables for the HiGHS loader to find the library.
    
    This is called when the package is imported to ensure the dynamic
    loader in HiGHS can find the HiPO library.
    """
    try:
        lib_path = get_library_path()
        lib_dir = str(lib_path.parent)
        
        # Set environment variable for HiGHS dynamic loader
        os.environ["HIGHSPY_HIPO_LIBRARY_PATH"] = lib_dir
        
        # Also add to system library path for completeness
        if sys.platform == "win32":
            path = os.environ.get("PATH", "")
            if lib_dir not in path:
                os.environ["PATH"] = lib_dir + os.pathsep + path
        elif sys.platform == "darwin":
            path = os.environ.get("DYLD_LIBRARY_PATH", "")
            if lib_dir not in path:
                os.environ["DYLD_LIBRARY_PATH"] = lib_dir + (os.pathsep + path if path else "")
        else:
            path = os.environ.get("LD_LIBRARY_PATH", "")
            if lib_dir not in path:
                os.environ["LD_LIBRARY_PATH"] = lib_dir + (os.pathsep + path if path else "")
    except FileNotFoundError:
        pass  # Library not found, will be handled later


# Set up environment on import
_setup_environment()
