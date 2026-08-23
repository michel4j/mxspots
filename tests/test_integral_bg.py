import pytest
import time
import numpy as np
from mxspots.models import SpotParams
from mxspots.spotfinder import findspots_data


def test_integral_bg_synthetic_flat_noise():
    # Frame with uniform Gaussian noise: mean=20.0, std=5.0
    rng = np.random.default_rng(42)
    ny, nx = 1024, 1024
    data = rng.normal(loc=20.0, scale=5.0, size=(ny, nx)).astype(np.float32)

    # Inject 10 distinct bright spots
    for i in range(10):
        x = 100 + i * 80
        y = 100 + i * 80
        data[y-1:y+2, x-1:x+2] += 200.0

    params = SpotParams(snr_threshold=4.0, min_spot_area=2, max_spot_area=100)
    spot_list = findspots_data(data, params=params)

    # All 10 injected spots should be detected
    assert spot_list.count == 10
    for s in spot_list.spots:
        assert s.snr > 10.0


def test_integral_bg_peak_exclusion():
    # Verify that a massive intense spot does not pollute its own local annulus background
    data = np.full((256, 256), fill_value=10.0, dtype=np.float32)
    # Put a huge spot at center
    data[126:130, 126:130] = 50000.0

    params = SpotParams(snr_threshold=3.0, min_spot_area=2, max_spot_area=50)
    spot_list = findspots_data(data, params=params)

    assert spot_list.count == 1
    spot = spot_list.spots[0]
    assert 126 <= spot.x <= 130
    assert 126 <= spot.y <= 130
    assert spot.intensity > 100000.0


def test_integral_bg_benchmark_scaling(clean_frame):
    # Benchmark full integral image background filtering on 3072x3072
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
        findspots_data(clean_frame.data, params=params)
    elapsed = (time.perf_counter() - t0) / n_iters

    print(f"\n[Integral Annulus Benchmark] 3072x3072 processing time: {elapsed * 1000.0:.2f} ms")
    assert elapsed < 0.25  # Sub-250ms for entire 3072x3072 frame including CCL
