import ctypes
import pytest
import math
from mxspots.models import Spot, SpotParams
from mxspots._lib import get_lib, CMxSpot, CMxScoreResult
from mxspots.scorer import score, score_spots


def test_blank_frame_zero_snr_leakage():
    """Verify that blank frames (N_B == 0) yield score == 0.0 and avg_snr == 0.0,
    even with hot pixel / high-SNR noise spots present."""
    lib = get_lib()

    n_spots = 5
    c_spots = (CMxSpot * n_spots)()
    for i in range(n_spots):
        c_spots[i].x = float(100 + i)
        c_spots[i].y = float(100 + i)
        c_spots[i].d_spacing = 2.0
        c_spots[i].intensity = 150000.0
        c_spots[i].snr = 133000.0  # Simulating detector hot pixel

    out = CMxScoreResult()
    out.bragg_spots = 0
    out.bragg_percent = 0.0
    out.avg_snr = 0.0
    ret = lib.mxspots_score_spots(c_spots, n_spots, None, ctypes.byref(out))
    assert ret == 0
    assert out.score == 0.0
    assert out.d_min == 999.0
    assert out.avg_snr == 0.0


def test_spot_count_dominance_over_extreme_snr():
    """Verify that dense diffraction frames outscore sparse frames with extreme SNR."""
    lib = get_lib()

    def eval_spots(n_bragg: int, snr_val: float) -> float:
        c_spots = (CMxSpot * n_bragg)()
        for i in range(n_bragg):
            c_spots[i].x = float(100 + i)
            c_spots[i].y = float(100 + i)
            c_spots[i].d_spacing = 1.95
            c_spots[i].intensity = 1000.0
            c_spots[i].snr = snr_val

        out = CMxScoreResult()
        out.bragg_spots = n_bragg
        out.bragg_percent = 95.0
        out.d_min = 1.95
        out.avg_snr = snr_val
        ret = lib.mxspots_score_spots(c_spots, n_bragg, None, ctypes.byref(out))
        assert ret == 0
        return out.score

    # Sparse frame with extreme SNR (simulating lyso1 frame 0027: 178 spots, SNR 806.7)
    score_sparse_high_snr = eval_spots(178, 806.7)

    # Dense frame with moderate-to-high SNR (simulating lyso1 frame 0024: 678 spots, SNR 209.0)
    score_dense = eval_spots(678, 209.0)

    # Ultra-dense frame with typical SNR (simulating lyso1 frame 0025: 815 spots, SNR 98.7)
    score_ultra_dense = eval_spots(815, 98.7)

    # Invariant: dense frame must strictly outscore sparse frame
    assert score_dense > score_sparse_high_snr
    assert score_ultra_dense > score_sparse_high_snr


def test_separation_deltas_across_diffraction_regimes():
    """Verify that score deltas between edge, core, and strong diffraction meet minimum acceptance thresholds."""
    lib = get_lib()

    def compute_model_score(n_b: int, snr: float, p_b: float, d: float) -> float:
        c_spots = (CMxSpot * n_b)()
        for i in range(n_b):
            c_spots[i].x = float(100 + i)
            c_spots[i].y = float(100 + i)
            c_spots[i].d_spacing = d
            c_spots[i].intensity = 1000.0
            c_spots[i].snr = snr

        out = CMxScoreResult()
        out.bragg_spots = n_b
        out.bragg_percent = p_b
        out.d_min = d
        out.avg_snr = snr
        ret = lib.mxspots_score_spots(c_spots, n_b, None, ctypes.byref(out))
        assert ret == 0
        return out.score

    # Weak edge diffraction (~75 spots, SNR 25, 80% Bragg, 2.3 A)
    score_weak = compute_model_score(75, 25.0, 80.0, 2.3)
    # Typical core diffraction (~250 spots, SNR 100, 90% Bragg, 2.0 A)
    score_typical = compute_model_score(250, 100.0, 90.0, 2.0)
    # Strong core diffraction (~650 spots, SNR 200, 95% Bragg, 1.9 A)
    score_strong = compute_model_score(650, 200.0, 95.0, 1.9)

    # Acceptance criterion: Delta >= 20% between weak edge and typical core
    delta_edge_core = score_typical - score_weak
    assert delta_edge_core >= 20.0, f"Expected >= 20.0%, got {delta_edge_core:.2f}%"

    # Acceptance criterion: Delta >= 10% between typical core and strong core
    delta_core_strong = score_strong - score_typical
    assert delta_core_strong >= 10.0, f"Expected >= 10.0%, got {delta_core_strong:.2f}%"


def test_strict_monotonicity():
    """Verify that increases in any positive factor strictly increase the score."""
    lib = get_lib()

    def eval_score(nb: int, snr: float, pb: float, d: float) -> float:
        c_spots = (CMxSpot * nb)()
        for i in range(nb):
            c_spots[i].x = float(100 + i)
            c_spots[i].y = float(100 + i)
            c_spots[i].d_spacing = d
            c_spots[i].intensity = 1000.0
            c_spots[i].snr = snr
        out = CMxScoreResult()
        out.bragg_spots = nb
        out.bragg_percent = pb
        out.d_min = d
        out.avg_snr = snr
        lib.mxspots_score_spots(c_spots, nb, None, ctypes.byref(out))
        return out.score

    base = eval_score(100, 50.0, 80.0, 2.0)

    # More Bragg spots strictly increases score
    assert eval_score(150, 50.0, 80.0, 2.0) > base

    # Higher SNR strictly increases score
    assert eval_score(100, 80.0, 80.0, 2.0) > base

    # Higher Bragg purity strictly increases score
    assert eval_score(100, 50.0, 95.0, 2.0) > base

    # Higher resolution (smaller d_min) strictly increases score
    assert eval_score(100, 50.0, 80.0, 1.5) > base
