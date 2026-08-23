import ctypes
import math
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
        ("num_masked_rings", ctypes.c_int),
        ("masked_rings_r2", (ctypes.c_float * 2) * 16),
    ]

    @classmethod
    def from_params(cls, params: SpotParams) -> "CMxSpotsParams":
        distance = params.distance if params.distance > 0.0 else 200.0
        wavelength = params.wavelength if params.wavelength > 0.0 else 1.0
        masked_r2_arr = ((ctypes.c_float * 2) * 16)()
        num_rings = 0
        if params.masked_rings:
            for d_min_ring, d_max_ring in params.masked_rings[:16]:
                # d_max corresponds to smaller radius r_min, d_min corresponds to larger radius r_max
                r_min = 0.0
                if d_max_ring > 0.0:
                    arg = wavelength / (2.0 * d_max_ring)
                    if arg < 1.0:
                        theta = math.asin(arg)
                        r_min = distance * math.tan(2.0 * theta)
                r_max = 0.0
                if d_min_ring > 0.0:
                    arg = wavelength / (2.0 * d_min_ring)
                    if arg < 1.0:
                        theta = math.asin(arg)
                        r_max = distance * math.tan(2.0 * theta)
                    else:
                        r_max = 1e6
                else:
                    r_max = 1e6

                r2_a = r_min * r_min
                r2_b = r_max * r_max
                masked_r2_arr[num_rings][0] = ctypes.c_float(min(r2_a, r2_b))
                masked_r2_arr[num_rings][1] = ctypes.c_float(max(r2_a, r2_b))
                num_rings += 1

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
            num_masked_rings=ctypes.c_int(num_rings),
            masked_rings_r2=masked_r2_arr,
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
        ("bragg_spots", ctypes.c_int),
        ("bragg_percent", ctypes.c_float),
        ("avg_intensity", ctypes.c_float),
        ("avg_snr", ctypes.c_float),
        ("d_min", ctypes.c_float),
        ("ice_score", ctypes.c_float),
        ("num_ice_rings", ctypes.c_int),
        ("num_lattices", ctypes.c_int),
        ("score", ctypes.c_float),
    ]


class CMxIceRing(ctypes.Structure):
    _fields_ = [
        ("d_spacing", ctypes.c_float),
        ("d_min", ctypes.c_float),
        ("d_max", ctypes.c_float),
        ("score", ctypes.c_float),
    ]


class CMxIceResult(ctypes.Structure):
    _fields_ = [
        ("num_rings", ctypes.c_int),
        ("ice_score", ctypes.c_float),
        ("rings", CMxIceRing * 16),
    ]


class CMxSpotsContext(ctypes.c_void_p):
    pass


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

    # Context management
    lib.mxspots_create_context.argtypes = [ctypes.c_int, ctypes.c_int]
    lib.mxspots_create_context.restype = ctypes.c_void_p

    lib.mxspots_free_context.argtypes = [ctypes.c_void_p]
    lib.mxspots_free_context.restype = None

    # Find spots with context
    lib.mxspots_find_spots_ctx.argtypes = [
        ctypes.c_void_p,                       # ctx
        ctypes.POINTER(ctypes.c_float),        # data
        ctypes.c_int,                          # nx
        ctypes.c_int,                          # ny
        ctypes.POINTER(CMxSpotsParams),        # params
        ctypes.POINTER(CMxSpot),               # out_spots
        ctypes.c_int,                          # max_spots
    ]
    lib.mxspots_find_spots_ctx.restype = ctypes.c_int

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

    # Regularity analysis
    lib.mxspots_analyze_regularity.argtypes = [
        ctypes.POINTER(CMxSpot),               # spots
        ctypes.c_int,                          # spot_count
        ctypes.POINTER(CMxSpotsParams),        # params
        ctypes.POINTER(ctypes.c_float),        # out_bragg_percent
        ctypes.POINTER(ctypes.c_int),          # out_bragg_spots
        ctypes.POINTER(ctypes.c_float),        # out_avg_intensity
        ctypes.POINTER(ctypes.c_int),          # out_num_lattices
        ctypes.POINTER(ctypes.c_float),        # out_d_min
    ]
    lib.mxspots_analyze_regularity.restype = ctypes.c_int

    # Detect ice rings
    lib.mxspots_detect_ice.argtypes = [
        ctypes.POINTER(ctypes.c_float),        # data
        ctypes.c_int,                          # nx
        ctypes.c_int,                          # ny
        ctypes.POINTER(CMxSpotsParams),        # params
        ctypes.POINTER(CMxIceResult),          # out_result
    ]
    lib.mxspots_detect_ice.restype = ctypes.c_int

    _LIB_CACHE = lib
    return _LIB_CACHE
