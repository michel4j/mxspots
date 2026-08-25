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


def test_c_mxspots_score_ice_penalty_capped_at_10_points():
    """Verify that ice penalty reduces score by at most 10 points even with extreme ice contamination."""
    lib = get_lib()
    n_spots = 100
    c_spots = (CMxSpot * n_spots)()
    for i in range(n_spots):
        c_spots[i].x = float(100 + i)
        c_spots[i].y = float(100 + i)
        c_spots[i].d_spacing = 1.8
        c_spots[i].intensity = 2000.0
        c_spots[i].snr = 30.0

    # Base score with 0 ice
    out_no_ice = CMxScoreResult()
    out_no_ice.bragg_spots = 90
    out_no_ice.bragg_percent = 90.0
    out_no_ice.avg_intensity = 2000.0
    out_no_ice.num_ice_rings = 0
    out_no_ice.ice_score = 0.0
    ret = lib.mxspots_score_spots(c_spots, n_spots, ctypes.byref(out_no_ice))
    assert ret == 0
    base_score = out_no_ice.score

    # Maximum ice (e.g. 6 ice rings, ice_score = 50.0)
    out_max_ice = CMxScoreResult()
    out_max_ice.bragg_spots = 90
    out_max_ice.bragg_percent = 90.0
    out_max_ice.avg_intensity = 2000.0
    out_max_ice.num_ice_rings = 6
    out_max_ice.ice_score = 50.0
    ret = lib.mxspots_score_spots(c_spots, n_spots, ctypes.byref(out_max_ice))
    assert ret == 0
    max_ice_score = out_max_ice.score

    # Difference must be exactly 10.0 points
    assert base_score - max_ice_score == pytest.approx(10.0, abs=1e-4)

    # Moderate ice (1 ring, ice_score = 3.0) -> p_ice = 2.0*1 + 1.0*(3.0-2.0) = 3.0 points
    out_mod_ice = CMxScoreResult()
    out_mod_ice.bragg_spots = 90
    out_mod_ice.bragg_percent = 90.0
    out_mod_ice.avg_intensity = 2000.0
    out_mod_ice.num_ice_rings = 1
    out_mod_ice.ice_score = 3.0
    ret = lib.mxspots_score_spots(c_spots, n_spots, ctypes.byref(out_mod_ice))
    assert ret == 0

    assert base_score - out_mod_ice.score == pytest.approx(3.0, abs=1e-4)


def test_c_mxspots_score_major_factor_sensitivity():
    """Verify that Bragg spot count and SNR are major co-equal drivers with large dynamic score impact."""
    lib = get_lib()

    def run_eval(n_spots: int, n_bragg: int, snr_val: float, p_bragg: float, d_val: float) -> float:
        c_spots = (CMxSpot * n_spots)()
        for i in range(n_spots):
            c_spots[i].x = float(100 + i)
            c_spots[i].y = float(100 + i)
            c_spots[i].d_spacing = d_val
            c_spots[i].intensity = 1000.0
            c_spots[i].snr = snr_val

        out = CMxScoreResult()
        out.bragg_spots = n_bragg
        out.bragg_percent = p_bragg
        out.d_min = d_val
        out.num_ice_rings = 0
        out.ice_score = 0.0
        ret = lib.mxspots_score_spots(c_spots, n_spots, ctypes.byref(out))
        assert ret == 0
        return out.score

    # 1. SNR sensitivity when Bragg count is fixed (50 spots)
    score_low_snr = run_eval(50, 40, snr_val=2.5, p_bragg=80.0, d_val=2.0)
    score_high_snr = run_eval(50, 40, snr_val=35.0, p_bragg=80.0, d_val=2.0)
    # Delta should be substantial (> 25 points)
    assert score_high_snr - score_low_snr > 25.0

    # 2. Bragg count sensitivity when SNR is fixed (SNR = 15)
    score_few_bragg = run_eval(10, 8, snr_val=15.0, p_bragg=80.0, d_val=2.0)
    score_many_bragg = run_eval(120, 96, snr_val=15.0, p_bragg=80.0, d_val=2.0)
    # Delta should be substantial (> 25 points)
    assert score_many_bragg - score_few_bragg > 25.0


def test_c_mxspots_score_factor_hierarchy_coequal_major_factors():
    """Verify that a high-SNR crystal with moderate spot count outscores a noisy high-count frame with poor SNR."""
    lib = get_lib()

    def run_eval(n_spots: int, n_bragg: int, snr_val: float, p_bragg: float, d_val: float) -> float:
        c_spots = (CMxSpot * n_spots)()
        for i in range(n_spots):
            c_spots[i].x = float(100 + i)
            c_spots[i].y = float(100 + i)
            c_spots[i].d_spacing = d_val
            c_spots[i].intensity = 1000.0
            c_spots[i].snr = snr_val

        out = CMxScoreResult()
        out.bragg_spots = n_bragg
        out.bragg_percent = p_bragg
        out.d_min = d_val
        out.num_ice_rings = 0
        out.ice_score = 0.0
        ret = lib.mxspots_score_spots(c_spots, n_spots, ctypes.byref(out))
        assert ret == 0
        return out.score

    # Frame A: Strong diffraction with moderate spot count (40 Bragg spots, SNR = 30.0, 80% Bragg)
    score_strong = run_eval(50, 40, snr_val=30.0, p_bragg=80.0, d_val=2.0)

    # Frame B: Weak / noisy frame with high spot count (150 Bragg spots, SNR = 3.5, 80% Bragg)
    score_weak_noisy = run_eval(188, 150, snr_val=3.5, p_bragg=80.0, d_val=2.0)

    # Strong high-SNR diffraction frame must decisively outscore the weak noisy high-count frame
    assert score_strong > score_weak_noisy
    assert score_strong - score_weak_noisy > 15.0


def test_c_mxspots_score_secondary_factor_bragg_percent_modulation():
    """Verify that Bragg percentage acts as a secondary factor providing moderate score modulation."""
    lib = get_lib()

    def run_eval(n_spots: int, n_bragg: int, snr_val: float, p_bragg: float, d_val: float) -> float:
        c_spots = (CMxSpot * n_spots)()
        for i in range(n_spots):
            c_spots[i].x = float(100 + i)
            c_spots[i].y = float(100 + i)
            c_spots[i].d_spacing = d_val
            c_spots[i].intensity = 1000.0
            c_spots[i].snr = snr_val

        out = CMxScoreResult()
        out.bragg_spots = n_bragg
        out.bragg_percent = p_bragg
        out.d_min = d_val
        out.num_ice_rings = 0
        out.ice_score = 0.0
        ret = lib.mxspots_score_spots(c_spots, n_spots, ctypes.byref(out))
        assert ret == 0
        return out.score

    # Fixed N_B = 50, SNR = 15.0, d_98 = 2.0 A
    score_low_purity = run_eval(250, 50, snr_val=15.0, p_bragg=20.0, d_val=2.0)
    score_high_purity = run_eval(55, 50, snr_val=15.0, p_bragg=90.9, d_val=2.0)

    # Higher crystalline purity increases score moderately (typically 8-18 points)
    delta = score_high_purity - score_low_purity
    assert 5.0 < delta < 25.0


def test_c_mxspots_score_tertiary_factor_resolution_refinement():
    """Verify that resolution limit acts as a tertiary refinement factor."""
    lib = get_lib()

    def run_eval(n_spots: int, n_bragg: int, snr_val: float, p_bragg: float, d_val: float) -> float:
        c_spots = (CMxSpot * n_spots)()
        for i in range(n_spots):
            c_spots[i].x = float(100 + i)
            c_spots[i].y = float(100 + i)
            c_spots[i].d_spacing = d_val
            c_spots[i].intensity = 1000.0
            c_spots[i].snr = snr_val

        out = CMxScoreResult()
        out.bragg_spots = n_bragg
        out.bragg_percent = p_bragg
        out.d_min = d_val
        out.num_ice_rings = 0
        out.ice_score = 0.0
        ret = lib.mxspots_score_spots(c_spots, n_spots, ctypes.byref(out))
        assert ret == 0
        return out.score

    # Fixed N_B = 50, SNR = 15.0, P_B = 80.0%
    score_low_res = run_eval(62, 50, snr_val=15.0, p_bragg=80.0, d_val=3.5)
    score_high_res = run_eval(62, 50, snr_val=15.0, p_bragg=80.0, d_val=1.4)

    # Resolution refinement yields subtle but positive improvement (typically 3-10 points)
    delta = score_high_res - score_low_res
    assert 2.0 < delta < 15.0


def test_c_mxspots_score_exact_formula_calibration():
    """Verify exact numerical agreement between C engine scoring and ADR 0020 specification."""
    lib = get_lib()

    test_cases = [
        # (n_bragg, total_spots, p_bragg, snr, d_98, n_ice, ice_score)
        (5, 20, 25.0, 3.0, 3.8, 0, 0.0),      # Minimal / weak
        (35, 50, 70.0, 12.0, 2.2, 0, 0.0),    # Moderate
        (120, 140, 85.7, 25.0, 1.6, 0, 0.0),  # Good
        (300, 320, 93.75, 45.0, 1.2, 0, 0.0), # Outstanding
        (80, 100, 80.0, 18.0, 2.0, 2, 4.5),   # Moderate with ice
    ]

    for n_b, n_s, p_b, snr, d, n_ice, ice_score in test_cases:
        c_spots = (CMxSpot * n_s)()
        for i in range(n_s):
            c_spots[i].x = float(100 + i)
            c_spots[i].y = float(100 + i)
            c_spots[i].d_spacing = d
            c_spots[i].intensity = 1000.0
            c_spots[i].snr = snr

        out = CMxScoreResult()
        out.bragg_spots = n_b
        out.bragg_percent = p_b
        out.d_min = d
        out.num_ice_rings = n_ice
        out.ice_score = ice_score
        ret = lib.mxspots_score_spots(c_spots, n_s, ctypes.byref(out))
        assert ret == 0

        # Compute ADR 0020 reference score in Python
        s_res = max(0.0, min(1.0, (4.0 - d) / (4.0 - 1.2))) if d < 4.0 else 0.0
        z = -6.50 + 0.75 * math.log(1.0 + n_b) + 1.20 * math.log(1.0 + snr) + 1.40 * (p_b / 100.0) + 0.50 * s_res
        raw_score = 100.0 / (1.0 + math.exp(-z))
        p_ice = min(10.0, max(0.0, 2.0 * n_ice + max(0.0, ice_score - 2.0)))
        expected_score = max(0.0, min(100.0, raw_score - p_ice))

        assert out.score == pytest.approx(expected_score, abs=1e-3)


def test_c_mxspots_score_bragg_avg_snr_preservation():
    """Verify that when out_score.avg_snr is pre-populated from Bragg regularity,
    injecting low-SNR noise spots does not overwrite or dilute avg_snr or slash the quality score."""
    lib = get_lib()

    n_bragg = 40
    n_noise = 160
    total_spots = n_bragg + n_noise

    c_spots = (CMxSpot * total_spots)()
    # 40 strong Bragg spots with SNR = 35.0
    for i in range(n_bragg):
        c_spots[i].x = float(100 + i)
        c_spots[i].y = float(100 + i)
        c_spots[i].d_spacing = 2.0
        c_spots[i].intensity = 5000.0
        c_spots[i].snr = 35.0

    # 160 low-SNR noise spots with SNR = 3.0
    for i in range(n_bragg, total_spots):
        c_spots[i].x = float(1000 + i)
        c_spots[i].y = float(1000 + i)
        c_spots[i].d_spacing = 4.5
        c_spots[i].intensity = 150.0
        c_spots[i].snr = 3.0

    # Pre-populate out_score with Bragg avg_snr = 35.0
    out = CMxScoreResult()
    out.bragg_spots = n_bragg
    out.bragg_percent = (n_bragg / total_spots) * 100.0
    out.avg_snr = 35.0
    out.d_min = 2.0
    out.num_ice_rings = 0
    out.ice_score = 0.0

    ret = lib.mxspots_score_spots(c_spots, total_spots, ctypes.byref(out))
    assert ret == 0
    # avg_snr must NOT have been overwritten with (40*35 + 160*3)/200 = 9.4
    assert out.avg_snr == pytest.approx(35.0, abs=1e-3)
    # Score should remain high because the genuine Bragg SNR is preserved
    assert out.score > 45.0
