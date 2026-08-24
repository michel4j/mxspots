import pytest
import ctypes
import math
import numpy as np
from mxspots._lib import get_lib, CMxSpotsParams, CMxSpot, CMxScoreResult
from mxspots.models import SpotParams


def test_c_mxspots_score_spots_98th_percentile():
    lib = get_lib()

    # Create 100 spots with controlled d-spacings:
    # 5 spots at ultra-high resolution (d = 1.0, 1.1, 1.2, 1.3, 1.4 A) - top 5%
    # 95 spots at moderate/low resolution (d = 2.0, 2.1, ..., 6.0 A) - bottom 95%
    n_spots = 100
    c_spots = (CMxSpot * n_spots)()

    # Top 5% highest resolution (indices 0..4)
    for i in range(5):
        c_spots[i].x = float(100 + i)
        c_spots[i].y = float(100 + i)
        c_spots[i].d_spacing = float(1.0 + 0.1 * i)
        c_spots[i].intensity = 1000.0
        c_spots[i].snr = 20.0

    # Remaining 95 spots (indices 5..99) with d in [2.0, 5.0]
    for i in range(5, 100):
        c_spots[i].x = float(100 + i)
        c_spots[i].y = float(100 + i)
        c_spots[i].d_spacing = float(2.0 + (i - 5) * 0.03)
        c_spots[i].intensity = 500.0
        c_spots[i].snr = 15.0

    out_score = CMxScoreResult()
    out_score.bragg_percent = 90.0
    out_score.bragg_spots = 90
    out_score.avg_intensity = 500.0
    out_score.ice_score = 0.0
    out_score.num_ice_rings = 0

    ret = lib.mxspots_score_spots(c_spots, n_spots, ctypes.byref(out_score))
    assert ret == 0

    # Raw minimum is 1.0 A, 98th percentile (where 98% have d >= d_98) should be 1.1 A (at index k = floor(0.02 * 99) = 1)
    assert out_score.d_min == pytest.approx(1.1, rel=1e-3)
    assert out_score.score > 60.0
    assert out_score.score <= 100.0


def test_c_mxspots_score_spots_gating_zero_bragg():
    lib = get_lib()
    n_spots = 50
    c_spots = (CMxSpot * n_spots)()
    for i in range(n_spots):
        c_spots[i].x = float(100 + i)
        c_spots[i].y = float(100 + i)
        c_spots[i].d_spacing = 2.0
        c_spots[i].intensity = 500.0
        c_spots[i].snr = 20.0

    # 0 Bragg spots -> score must be 0.0 regardless of spot count or SNR
    out_score = CMxScoreResult()
    out_score.bragg_spots = 0
    out_score.bragg_percent = 0.0
    out_score.avg_intensity = 0.0

    ret = lib.mxspots_score_spots(c_spots, n_spots, ctypes.byref(out_score))
    assert ret == 0
    assert out_score.score == pytest.approx(0.0)


def test_c_mxspots_score_spots_empty():
    lib = get_lib()
    out_score = CMxScoreResult()
    ret = lib.mxspots_score_spots(None, 0, ctypes.byref(out_score))
    assert ret == 0
    assert out_score.spot_count == 0
    assert out_score.bragg_spots == 0
    assert out_score.d_min >= 999.0
    assert out_score.avg_snr == pytest.approx(0.0)
    assert out_score.score == pytest.approx(0.0)


def test_c_mxspots_score_frame_clean(clean_frame):
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

    out_score = CMxScoreResult()
    ret = lib.mxspots_score_frame(data_ptr, nx, ny, ctypes.byref(c_params), ctypes.byref(out_score))

    assert ret == 0
    assert out_score.spot_count > 0
    assert out_score.bragg_spots > 0
    assert out_score.bragg_percent > 50.0
    assert out_score.avg_snr > 10.0
    assert out_score.d_min < 3.0
    assert out_score.score > 50.0


def test_c_mxspots_score_logistic_monotonicity():
    lib = get_lib()

    def compute_score(n_bragg: int, snr_val: float, d_val: float, num_ice: int) -> float:
        n_spots = 100
        c_spots = (CMxSpot * n_spots)()
        for i in range(n_spots):
            c_spots[i].x = float(100 + i)
            c_spots[i].y = float(100 + i)
            c_spots[i].d_spacing = d_val
            c_spots[i].intensity = 1000.0
            c_spots[i].snr = snr_val

        out_score = CMxScoreResult()
        out_score.bragg_spots = n_bragg
        out_score.bragg_percent = float(n_bragg)
        out_score.avg_intensity = 1000.0
        out_score.num_ice_rings = num_ice
        out_score.ice_score = float(num_ice) * 2.5

        ret = lib.mxspots_score_spots(c_spots, n_spots, ctypes.byref(out_score))
        assert ret == 0
        return out_score.score

    # Increasing Bragg spots increases score
    s10 = compute_score(10, 20.0, 2.0, 0)
    s80 = compute_score(80, 20.0, 2.0, 0)
    assert s80 > s10

    # Better resolution (smaller d) increases score
    s_lowres = compute_score(80, 20.0, 3.5, 0)
    s_hires = compute_score(80, 20.0, 1.5, 0)
    assert s_hires > s_lowres

    # Ice contamination penalizes score
    s_no_ice = compute_score(80, 20.0, 2.0, 0)
    s_with_ice = compute_score(80, 20.0, 2.0, 2)
    assert s_no_ice > s_with_ice


def test_c_mxspots_score_intensity_weight_removed():
    """Verify that Bragg average intensity alone does not affect score when SNR and counts are identical."""
    lib = get_lib()
    n_spots = 50
    c_spots = (CMxSpot * n_spots)()
    for i in range(n_spots):
        c_spots[i].x = float(100 + i)
        c_spots[i].y = float(100 + i)
        c_spots[i].d_spacing = 2.0
        c_spots[i].intensity = 500.0
        c_spots[i].snr = 10.0

    out1 = CMxScoreResult()
    out1.bragg_spots = 40
    out1.bragg_percent = 80.0
    out1.avg_intensity = 100.0
    ret1 = lib.mxspots_score_spots(c_spots, n_spots, ctypes.byref(out1))
    assert ret1 == 0

    out2 = CMxScoreResult()
    out2.bragg_spots = 40
    out2.bragg_percent = 80.0
    out2.avg_intensity = 50000.0
    ret2 = lib.mxspots_score_spots(c_spots, n_spots, ctypes.byref(out2))
    assert ret2 == 0

    assert out1.score == pytest.approx(out2.score, rel=1e-5)
