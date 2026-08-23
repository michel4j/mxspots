import pytest
import time
import numpy as np
from pathlib import Path
from mxspots.models import SpotParams, IceRing
from mxspots.spotfinder import detect_ice_rings, detect_ice_rings_data
from mxspots.synthetic import add_powder_ring, get_cached_synthetic_frame


def test_detect_ice_clean_frame_no_rings(clean_frame):
    # A clean protein frame without ice rings should detect 0 ice rings
    params = SpotParams(
        beam_x=clean_frame.cx,
        beam_y=clean_frame.cy,
        pixel_size_x=clean_frame.qx,
        pixel_size_y=clean_frame.qy,
        distance=clean_frame.distance,
        wavelength=clean_frame.wavelength,
        ice_sensitivity=3.0,
    )
    rings, score = detect_ice_rings(clean_frame.data, params=params)
    assert len(rings) == 0
    assert score < 3.0


def test_detect_ice_synthetic_powder_ring():
    # Construct a flat noisy frame and inject a strong ice ring at 3.897 A
    rng = np.random.default_rng(42)
    ny, nx = 2048, 2048
    cx, cy = 1024.0, 1024.0
    qx, qy = 0.075, 0.075
    distance = 200.0
    wavelength = 1.0
    d_target = 3.897

    data = rng.normal(loc=15.0, scale=3.0, size=(ny, nx)).astype(np.float32)

    # Inject powder ring
    add_powder_ring(
        data,
        cx=cx,
        cy=cy,
        qx=qx,
        qy=qy,
        distance=distance,
        wavelength=wavelength,
        d_spacing=d_target,
        radial_width=2.5,
        peak_intensity=50.0,
    )

    params = SpotParams(
        beam_x=cx,
        beam_y=cy,
        pixel_size_x=qx,
        pixel_size_y=qy,
        distance=distance,
        wavelength=wavelength,
        ice_sensitivity=3.0,
    )

    rings, score = detect_ice_rings(data, params=params)
    assert len(rings) >= 1
    assert score >= 3.0

    # The 3.897 A ring should be detected
    found = any(pytest.approx(r.d_spacing, abs=0.1) == d_target for r in rings)
    assert found
    detected_ring = [r for r in rings if abs(r.d_spacing - d_target) < 0.1][0]
    assert detected_ring.d_min < d_target < detected_ring.d_max
    assert detected_ring.score >= 3.0


def test_detect_ice_multiple_powder_rings():
    # Inject both 3.897 A and 3.669 A rings
    rng = np.random.default_rng(123)
    ny, nx = 2048, 2048
    cx, cy = 1024.0, 1024.0
    qx, qy = 0.075, 0.075
    distance = 200.0
    wavelength = 1.0

    data = rng.normal(loc=12.0, scale=2.5, size=(ny, nx)).astype(np.float32)

    add_powder_ring(
        data,
        cx=cx,
        cy=cy,
        qx=qx,
        qy=qy,
        distance=distance,
        wavelength=wavelength,
        d_spacing=3.897,
        radial_width=2.0,
        peak_intensity=60.0,
    )

    add_powder_ring(
        data,
        cx=cx,
        cy=cy,
        qx=qx,
        qy=qy,
        distance=distance,
        wavelength=wavelength,
        d_spacing=3.669,
        radial_width=2.0,
        peak_intensity=50.0,
    )

    params = SpotParams(
        beam_x=cx,
        beam_y=cy,
        pixel_size_x=qx,
        pixel_size_y=qy,
        distance=distance,
        wavelength=wavelength,
        ice_sensitivity=3.0,
    )

    rings, score = detect_ice_rings(data, params=params)
    assert len(rings) >= 2
    d_spacings = [r.d_spacing for r in rings]
    assert any(pytest.approx(3.897, abs=0.1) == d for d in d_spacings)
    assert any(pytest.approx(3.669, abs=0.1) == d for d in d_spacings)


def test_detect_ice_filepath_input(test_data_dir):
    ice_path = test_data_dir / "ice.yaml"
    rings, score = detect_ice_rings(ice_path, params=SpotParams(ice_sensitivity=2.0))
    assert isinstance(rings, list)
    assert score > 1.5


def test_detect_ice_sensitivity_filtering():
    rng = np.random.default_rng(42)
    ny, nx = 1024, 1024
    data = rng.normal(loc=10.0, scale=2.0, size=(ny, nx)).astype(np.float32)

    # Moderate ring
    add_powder_ring(
        data,
        cx=512.0,
        cy=512.0,
        qx=0.075,
        qy=0.075,
        distance=200.0,
        wavelength=1.0,
        d_spacing=3.897,
        radial_width=2.0,
        peak_intensity=12.0,
    )

    # Very high sensitivity threshold (100.0) should filter it out
    params_high = SpotParams(
        ice_sensitivity=100.0,
        beam_x=512.0,
        beam_y=512.0,
        distance=200.0,
        wavelength=1.0,
        pixel_size_x=0.075,
        pixel_size_y=0.075,
    )
    rings_high, score_high = detect_ice_rings_data(data, params=params_high)
    assert len(rings_high) == 0

    # Low sensitivity threshold (20.0) should include it
    params_low = SpotParams(
        ice_sensitivity=20.0,
        beam_x=512.0,
        beam_y=512.0,
        distance=200.0,
        wavelength=1.0,
        pixel_size_x=0.075,
        pixel_size_y=0.075,
    )
    rings_low, score_low = detect_ice_rings_data(data, params=params_low)
    assert len(rings_low) >= 1
    assert score_low > 20.0


def test_detect_ice_benchmark_scaling(clean_frame):
    # Benchmark radial profile and ice detection on 3072x3072 frame
    params = SpotParams(
        beam_x=clean_frame.cx,
        beam_y=clean_frame.cy,
        pixel_size_x=clean_frame.qx,
        pixel_size_y=clean_frame.qy,
        distance=clean_frame.distance,
        wavelength=clean_frame.wavelength,
    )

    t0 = time.perf_counter()
    n_iters = 5
    for _ in range(n_iters):
        detect_ice_rings(clean_frame.data, params=params)
    elapsed = (time.perf_counter() - t0) / n_iters

    print(f"\n[Ice Detection Benchmark] 3072x3072 processing time: {elapsed * 1000.0:.2f} ms")
    assert elapsed < 0.15  # Sub-150ms for entire 3072x3072 frame OpenMP radial profiling
