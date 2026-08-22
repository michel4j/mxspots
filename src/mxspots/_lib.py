import ctypes
import os
import sys
import sysconfig
from pathlib import Path
from typing import Optional
import numpy as np
from .models import SpotParams


class CMxSpotsParams(ctypes.Structure):
    _fields_ = [
        ("snr_threshold", ctypes.c_float),
        ("min_spot_area", ctypes.c_int),
        ("max_spot_area", ctypes.c_int),
        ("beam_x", ctypes.c_float),
        ("beam_y", ctypes.c_float),
        ("pixel_size_x", ctypes.c_float),
        ("pixel_size_y", ctypes.c_float),
        ("distance", ctypes.c_float),
        ("wavelength", ctypes.c_float),
    ]

    @classmethod
    def from_params(cls, params: SpotParams) -> "CMxSpotsParams":
        return cls(
            snr_threshold=ctypes.c_float(params.snr_threshold),
            min_spot_area=ctypes.c_int(params.min_spot_area),
            max_spot_area=ctypes.c_int(params.max_spot_area),
            beam_x=ctypes.c_float(params.beam_x),
            beam_y=ctypes.c_float(params.beam_y),
            pixel_size_x=ctypes.c_float(params.pixel_size_x),
            pixel_size_y=ctypes.c_float(params.pixel_size_y),
            distance=ctypes.c_float(params.distance),
            wavelength=ctypes.c_float(params.wavelength),
        )


class CMxSpot(ctypes.Structure):
    _fields_ = [
        ("x", ctypes.c_float),
        ("y", ctypes.c_float),
        ("d_spacing", ctypes.c_float),
        ("intensity", ctypes.c_float),
        ("snr", ctypes.c_float),
    ]


class CMxScoreResult(ctypes.Structure):
    _fields_ = [
        ("spot_count", ctypes.c_int),
        ("avg_snr", ctypes.c_float),
        ("d_min", ctypes.c_float),
        ("percentage_indexed", ctypes.c_float),
    ]


def _find_library() -> Path:
    """Locate the compiled libmxspots shared library."""
    current_dir = Path(__file__).resolve().parent
    repo_root = current_dir.parent.parent
    build_dir = repo_root / "build"

    search_dirs = [
        current_dir,
        current_dir / "lib",
        build_dir,
        build_dir / "csrc",
    ]

    # Include virtualenv/site-packages locations (e.g. editable installs)
    for path_name in ("platlib", "purelib"):
        try:
            p = Path(sysconfig.get_path(path_name)) / "mxspots"
            if p not in search_dirs:
                search_dirs.append(p)
        except Exception:
            pass

    lib_names = [
        "libmxspots.so",
        "libmxspots.dylib",
        "mxspots.dll",
        "libmxspots.dll",
        "mxspots.so",
    ]

    for directory in search_dirs:
        for lib_name in lib_names:
            candidate = directory / lib_name
            if candidate.is_file():
                return candidate

    # Search dynamically in system library path
    import ctypes.util
    found = ctypes.util.find_library("mxspots")
    if found:
        return Path(found)

    raise FileNotFoundError(
        "Could not find compiled libmxspots shared library. "
        "Please build the project using 'pip install .' or CMake."
    )


_lib: Optional[ctypes.CDLL] = None


def get_lib() -> ctypes.CDLL:
    """Get or load the libmxspots CDLL handle."""
    global _lib
    if _lib is None:
        lib_path = _find_library()
        _lib = ctypes.CDLL(str(lib_path))
        _configure_signatures(_lib)
    return _lib


def _configure_signatures(lib: ctypes.CDLL) -> None:
    """Configure argtypes and restype for C functions."""
    lib.mxspots_get_version.argtypes = []
    lib.mxspots_get_version.restype = ctypes.c_int

    lib.mxspots_ping.argtypes = [ctypes.POINTER(CMxSpotsParams)]
    lib.mxspots_ping.restype = ctypes.c_int

    # mxspots_find_spots signature
    lib.mxspots_find_spots.argtypes = [
        ctypes.POINTER(ctypes.c_float),          # data
        ctypes.c_int,                            # nx
        ctypes.c_int,                            # ny
        ctypes.POINTER(CMxSpotsParams),          # params
        ctypes.POINTER(CMxSpot),                 # out_spots
        ctypes.c_int,                            # max_spots
    ]
    lib.mxspots_find_spots.restype = ctypes.c_int
