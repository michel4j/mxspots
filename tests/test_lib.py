import pytest
import ctypes
from mxspots._lib import get_lib, CMxSpotsParams
from mxspots.models import SpotParams


def test_lib_version():
    lib = get_lib()
    assert lib is not None
    version = lib.mxspots_get_version()
    assert version == 100


def test_lib_ping_valid():
    lib = get_lib()
    params = SpotParams(
        snr_threshold=3.5,
        min_spot_area=3,
        max_spot_area=100,
        beam_x=1200.0,
        beam_y=1250.0,
        distance=150.0,
        wavelength=0.979,
    )
    c_params = CMxSpotsParams.from_params(params)
    res = lib.mxspots_ping(ctypes.byref(c_params))
    assert res == 0


def test_lib_ping_invalid():
    lib = get_lib()
    params = SpotParams(snr_threshold=-1.0)
    c_params = CMxSpotsParams.from_params(params)
    res = lib.mxspots_ping(ctypes.byref(c_params))
    assert res == -2
