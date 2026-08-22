import ctypes
import os
import sys
import sysconfig
from pathlib import Path
from typing import Optional
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
        ("d_min", ctypes.c_float),
        ("d_max", ctypes.c_float),
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
            d_min=ctypes.c_float(params.d_min),
            d_max=ctypes.c_float(params.d_max),
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


class CMxIndexResult(ctypes.Structure):
    _fields_ = [
        ("unit_cell", ctypes.c_float * 6),
        ("percentage_indexed", ctypes.c_float),
        ("indexed_spot_count", ctypes.c_int),
        ("total_spot_count", ctypes.c_int),
    ]


_LIB_CACHE: Optional[ctypes.CDLL] = None


def _find_library() -> Path:
    """Find the compiled libmxspots shared library."""
    lib_name = "libmxspots.so" if sys.platform != "win32" else "mxspots.dll"
    if sys.platform == "darwin":
        lib_name = "libmxspots.dylib"

    possible_dirs = [
        # In platlib directory where wheel/pip places extension modules
        Path(sysconfig.get_path("platlib")) / "mxspots",
        Path(sysconfig.get_path("purelib")) / "mxspots",
        # In editable install build directory
        Path(__file__).resolve().parent.parent.parent / "build",
        # In package directory directly
        Path(__file__).resolve().parent,
        # System/Current directory
        Path("."),
    ]

    for d in possible_dirs:
        # Check direct path
        candidate = d / lib_name
        if candidate.is_file():
            return candidate
        # Recursively look in build folders
        for match in d.glob(f"**/{lib_name}"):
            if match.is_file():
                return match

    raise FileNotFoundError(
        f"Could not locate {lib_name}. Please build or install the package with `pip install .`"
    )


def get_lib() -> ctypes.CDLL:
    """Load and return the libmxspots ctypes CDLL handle with bound function prototypes."""
    global _LIB_CACHE
    if _LIB_CACHE is not None:
        return _LIB_CACHE

    lib_path = _find_library()
    lib = ctypes.CDLL(str(lib_path))

    # Version check
    lib.mxspots_get_version.argtypes = []
    lib.mxspots_get_version.restype = ctypes.c_int

    # Ping check
    lib.mxspots_ping.argtypes = [ctypes.POINTER(CMxSpotsParams)]
    lib.mxspots_ping.restype = ctypes.c_int

    # Find spots
    lib.mxspots_find_spots.argtypes = [
        ctypes.POINTER(ctypes.c_float),        # data
        ctypes.c_int,                          # nx
        ctypes.c_int,                          # ny
        ctypes.POINTER(CMxSpotsParams),        # params
        ctypes.POINTER(CMxSpot),               # out_spots
        ctypes.c_int,                          # max_spots
    ]
    lib.mxspots_find_spots.restype = ctypes.c_int

    # Score spots
    lib.mxspots_score_spots.argtypes = [
        ctypes.POINTER(CMxSpot),               # spots
        ctypes.c_int,                          # spot_count
        ctypes.POINTER(CMxScoreResult),        # out_score
    ]
    lib.mxspots_score_spots.restype = ctypes.c_int

    # Score frame directly
    lib.mxspots_score_frame.argtypes = [
        ctypes.POINTER(ctypes.c_float),        # data
        ctypes.c_int,                          # nx
        ctypes.c_int,                          # ny
        ctypes.POINTER(CMxSpotsParams),        # params
        ctypes.POINTER(CMxScoreResult),        # out_score
    ]
    lib.mxspots_score_frame.restype = ctypes.c_int

    # Index spots
    lib.mxspots_index_spots.argtypes = [
        ctypes.POINTER(CMxSpot),               # spots
        ctypes.c_int,                          # spot_count
        ctypes.POINTER(CMxSpotsParams),        # params
        ctypes.POINTER(CMxIndexResult),        # out_index
    ]
    lib.mxspots_index_spots.restype = ctypes.c_int

    # Index frame directly
    lib.mxspots_index_frame.argtypes = [
        ctypes.POINTER(ctypes.c_float),        # data
        ctypes.c_int,                          # nx
        ctypes.c_int,                          # ny
        ctypes.POINTER(CMxSpotsParams),        # params
        ctypes.POINTER(CMxIndexResult),        # out_index
    ]
    lib.mxspots_index_frame.restype = ctypes.c_int

    _LIB_CACHE = lib
    return _LIB_CACHE
