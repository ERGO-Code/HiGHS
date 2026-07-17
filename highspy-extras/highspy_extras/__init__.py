"""External dependencies extension for highspy."""

from __future__ import annotations

import ctypes
import sys
from dataclasses import dataclass
from functools import cached_property
from importlib.metadata import version
from pathlib import Path

__all__ = ["__version__", "HighsExtrasFeatureInfo", "library"]

__version__ = version("highspy-extras")
_PACKAGE_DIR = Path(__file__).resolve().parent


class _HighsExtrasFeatureInfoRaw(ctypes.Structure):
    """Raw feature info layout exposed by the shared library."""

    _fields_ = [
        ("provider", ctypes.c_char_p),
        ("version", ctypes.c_char_p),
        ("license", ctypes.c_char_p),
        ("enabled", ctypes.c_bool),
    ]


@dataclass(frozen=True)
class HighsExtrasFeatureInfo:
    """Feature metadata for an external dependency."""

    provider: str
    version: str
    license: str
    enabled: bool

    @classmethod
    def from_raw(cls, raw: _HighsExtrasFeatureInfoRaw) -> HighsExtrasFeatureInfo:
        return cls(
            provider=raw.provider.decode("utf-8"),
            version=raw.version.decode("utf-8"),
            license=raw.license.decode("utf-8"),
            enabled=bool(raw.enabled),
        )


class HighsExtrasLibrary:
    """Wrapper around the highs_extras shared library."""

    @cached_property
    def handle(self) -> ctypes.CDLL:
        if sys.platform == "win32":
            library_name = "highs_extras.dll"
        elif sys.platform == "darwin":
            library_name = "libhighs_extras.dylib"
        else:
            library_name = "libhighs_extras.so"

        library_path = _PACKAGE_DIR / library_name
        if not library_path.is_file():
            raise FileNotFoundError(f"Could not find the shared library at {library_path}")

        handle = ctypes.CDLL(str(library_path))

        handle.HighsExtras_getVersion.argtypes = []
        handle.HighsExtras_getVersion.restype = ctypes.c_char_p

        handle.HighsExtras_getFeatureCount.argtypes = []
        handle.HighsExtras_getFeatureCount.restype = ctypes.c_size_t

        handle.HighsExtras_getFeatureName.argtypes = [ctypes.c_size_t]
        handle.HighsExtras_getFeatureName.restype = ctypes.c_char_p

        handle.HighsExtras_getFeatureInfo.argtypes = []
        handle.HighsExtras_getFeatureInfo.restype = ctypes.POINTER(_HighsExtrasFeatureInfoRaw)

        return handle

    @property
    def version(self) -> str:
        version_bytes = self.handle.HighsExtras_getVersion()
        if version_bytes is None:
            raise RuntimeError("HighsExtras_getVersion() returned NULL")
        return version_bytes.decode("utf-8")

    def _feature_name(self, index: int) -> str:
        name_bytes = self.handle.HighsExtras_getFeatureName(index)
        if name_bytes is None:
            raise RuntimeError(f"HighsExtras_getFeatureName({index}) returned NULL")
        return name_bytes.decode("utf-8")

    @cached_property
    def features(self) -> dict[str, HighsExtrasFeatureInfo]:
        info_ptr = self.handle.HighsExtras_getFeatureInfo()
        if not info_ptr:
            raise RuntimeError("HighsExtras_getFeatureInfo() returned NULL")

        count = int(self.handle.HighsExtras_getFeatureCount())
        return {self._feature_name(index): HighsExtrasFeatureInfo.from_raw(info_ptr[index]) for index in range(count)}

    def __getitem__(self, name: str) -> HighsExtrasFeatureInfo:
        return self.features[name]

    @cached_property
    def feature_table(self) -> str:
        """Return a human-readable table describing the external dependency features."""

        headers = ("key", "name", "version", "license", "enabled")
        rows = [
            (
                name,
                info.provider,
                info.version,
                info.license,
                "yes" if info.enabled else "no",
            )
            for name, info in self.features.items()
        ]

        widths = [max(len(headers[i]), *(len(row[i]) for row in rows)) if rows else len(headers[i]) for i in range(len(headers))]

        def _fmt(row: tuple[str, ...]) -> str:
            return "  ".join(cell.ljust(widths[i]) for i, cell in enumerate(row))

        separator = "  ".join("-" * w for w in widths)
        lines = [_fmt(headers), separator]
        lines.extend(_fmt(row) for row in rows)
        return "\n".join(lines)

    def __str__(self) -> str:
        return self.feature_table


library = HighsExtrasLibrary()
