"""Thin wrapper around scikit-build-core that stages parent-repo files
into the project tree before building the sdist."""

import shutil
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

from scikit_build_core import build as _orig

_HERE = Path(__file__).resolve().parent
_PARENT = _HERE.parent


def _get_sdist_include():
    with open(_HERE / "pyproject.toml", "rb") as f:
        cfg = tomllib.load(f)
    return cfg.get("tool", {}).get("scikit-build", {}).get("sdist", {}).get("include", [])


def build_sdist(sdist_directory, config_settings=None):
    created = []
    try:
        for name in _get_sdist_include():
            dst, src = _HERE / name, _PARENT / name
            print(f"Staging {src} -> {dst}")
            if dst.exists() or not src.exists():
                continue

            # Track any new parent dirs so they can be cleaned up
            new_dirs = []
            for parent in reversed(dst.relative_to(_HERE).parents[:-1]):
                d = _HERE / parent
                if d.exists():
                    break
                new_dirs.append(d)

            if new_dirs:
                new_dirs[0].mkdir(parents=True)
                created.extend(reversed(new_dirs))

            (shutil.copytree if src.is_dir() else shutil.copy2)(src, dst)
            created.append(dst)
        return _orig.build_sdist(sdist_directory, config_settings)
    finally:
        for p in reversed(created):
            shutil.rmtree(p) if p.is_dir() else p.unlink(missing_ok=True)


# Re-export everything else unchanged.
build_wheel = _orig.build_wheel
build_editable = _orig.build_editable
get_requires_for_build_sdist = _orig.get_requires_for_build_sdist
get_requires_for_build_wheel = _orig.get_requires_for_build_wheel
get_requires_for_build_editable = _orig.get_requires_for_build_editable
prepare_metadata_for_build_wheel = _orig.prepare_metadata_for_build_wheel
prepare_metadata_for_build_editable = _orig.prepare_metadata_for_build_editable
