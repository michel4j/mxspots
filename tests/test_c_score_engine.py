import pytest
import ctypes
import numpy as np
from mxspots._lib import get_lib, CMxSpotsParams, CMxSpot, CMxScoreResult
from mxspots.models import SpotParams


def test_c_mxspots_score_spots_95th_percentile():
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
    out_score.percentage_indexed = 90.0
    out_score.indexed_spot_count = 90
    out_score.ice_score = 0.0
    out_score.num_ice_rings = 0

    ret = lib.mxspots_score_spots(c_spots, n_spots, ctypes.byref(out_score))
    assert ret == 0

    # Raw minimum is 1.0 A, but 95th percentile (where 95% have d >= d_95) should be ~2.0 A (at index k = floor(0.05 * 99) = 4, so d ~ 1.4 to 2.0 A)
    assert out_score.d_min > 1.3
    assert out_score.d_min <= 2.05
    assert out_score.score > 60.0
    assert out_score.score <= 100.0


def test_c_mxspots_score_spots_empty():
    lib = get_lib()
    out_score = CMxScoreResult()
    ret = lib.mxspots_score_spots(None, 0, ctypes.byref(out_score))
    assert ret == 0
    assert out_score.spot_count == 0
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
    assert out_score.avg_snr > 10.0
    assert out_score.d_min < 3.0
    assert out_score.score > 50.0


def test_c_mxspots_score_snr_saturation_and_multi_lattice():
    lib = get_lib()

    # Create identical spots with varying SNR and num_lattices
    def compute_score(snr_val: float, num_lattices: int) -> float:
        n_spots = 100
        c_spots = (CMxSpot * n_spots)()
        for i in range(n_spots):
            c_spots[i].x = float(100 + i)
            c_spots[i].y = float(100 + i)
            c_spots[i].d_spacing = 2.0
            c_spots[i].intensity = 1000.0
            c_spots[i].snr = snr_val

        out_score = CMxScoreResult()
        out_score.percentage_regular = 80.0
        out_score.regular_spot_count = 80
        out_score.num_lattices = num_lattices

        ret = lib.mxspots_score_spots(c_spots, n_spots, ctypes.byref(out_score))
        assert ret == 0
        return out_score.score

    score_snr50 = compute_score(50.0, 1)
    score_snr100 = compute_score(100.0, 1)
    score_snr200 = compute_score(200.0, 1)

    # SNR 100 awards 25 pts (12.5 pts higher than SNR 50 which awards 12.5 pts)
    assert score_snr100 == pytest.approx(score_snr50 + 12.5, abs=0.1)
    # SNR 200 should saturate at same 25 pts as SNR 100
    assert score_snr200 == pytest.approx(score_snr100, abs=0.01)

    # Multi-lattice frame (num_lattices = 3) should NOT deduct any points vs single lattice
    score_multi = compute_score(100.0, 3)
    assert score_multi == pytest.approx(score_snr100, abs=0.01)
