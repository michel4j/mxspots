import pytest
import ctypes
import numpy as np
from mxspots._lib import get_lib, CMxSpotsParams, CMxSpot, CMxScoreResult
from mxspots.models import SpotParams


def test_c_regularity_clean_frame(clean_frame):
    lib = get_lib()
    params = SpotParams(
        beam_x=clean_frame.cx,
        beam_y=clean_frame.cy,
        pixel_size_x=clean_frame.qx,
        pixel_size_y=clean_frame.qy,
        distance=clean_frame.distance,
        wavelength=clean_frame.wavelength,
    )
    c_params = CMxSpotsParams.from_params(params)
    data = clean_frame.data
    ny, nx = data.shape
    data_ptr = data.ctypes.data_as(ctypes.POINTER(ctypes.c_float))

    max_spots = 5000
    c_spots = (CMxSpot * max_spots)()
    spot_count = lib.mxspots_find_spots(data_ptr, nx, ny, ctypes.byref(c_params), c_spots, max_spots)
    assert spot_count > 20

    bragg_pct = ctypes.c_float(0.0)
    bragg_spots = ctypes.c_int(0)
    avg_intensity = ctypes.c_float(0.0)
    num_lattices = ctypes.c_int(0)

    ret = lib.mxspots_analyze_regularity(
        c_spots,
        spot_count,
        ctypes.byref(c_params),
        ctypes.byref(bragg_pct),
        ctypes.byref(bragg_spots),
        ctypes.byref(avg_intensity),
        ctypes.byref(num_lattices),
    )
    assert ret == 0
    assert bragg_pct.value > 50.0
    assert bragg_spots.value > 10
    assert avg_intensity.value > 0.0
    assert num_lattices.value >= 1


def test_c_regularity_random_noise_spots():
    lib = get_lib()
    params = SpotParams(
        beam_x=1500.0,
        beam_y=1500.0,
        pixel_size_x=0.075,
        pixel_size_y=0.075,
        distance=200.0,
        wavelength=1.0,
    )
    c_params = CMxSpotsParams.from_params(params)

    # Generate 50 random spots with no lattice periodicity
    rng = np.random.default_rng(42)
    n_spots = 50
    c_spots = (CMxSpot * n_spots)()
    for i in range(n_spots):
        c_spots[i].x = float(rng.uniform(200, 2800))
        c_spots[i].y = float(rng.uniform(200, 2800))
        c_spots[i].d_spacing = float(rng.uniform(1.5, 10.0))
        c_spots[i].intensity = float(rng.uniform(500, 2000))
        c_spots[i].snr = 10.0

    bragg_pct = ctypes.c_float(0.0)
    bragg_spots = ctypes.c_int(0)
    avg_intensity = ctypes.c_float(0.0)
    num_lattices = ctypes.c_int(0)

    ret = lib.mxspots_analyze_regularity(
        c_spots,
        n_spots,
        ctypes.byref(c_params),
        ctypes.byref(bragg_pct),
        ctypes.byref(bragg_spots),
        ctypes.byref(avg_intensity),
        ctypes.byref(num_lattices),
    )
    assert ret == 0
    assert bragg_pct.value < 30.0
    assert num_lattices.value <= 1


def test_c_regularity_empty():
    lib = get_lib()
    params = SpotParams()
    c_params = CMxSpotsParams.from_params(params)

    bragg_pct = ctypes.c_float(0.0)
    bragg_spots = ctypes.c_int(0)
    avg_intensity = ctypes.c_float(0.0)
    num_lattices = ctypes.c_int(0)

    ret = lib.mxspots_analyze_regularity(
        None,
        0,
        ctypes.byref(c_params),
        ctypes.byref(bragg_pct),
        ctypes.byref(bragg_spots),
        ctypes.byref(avg_intensity),
        ctypes.byref(num_lattices),
    )
    assert ret == 0
    assert bragg_pct.value == 0.0
    assert bragg_spots.value == 0
    assert avg_intensity.value == 0.0
    assert num_lattices.value == 0
