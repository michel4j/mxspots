import os
import time
import pytest
from mxspots.models import SpotParams
from mxspots.spotfinder import findspots_data
from mxspots.synthetic import generate_synthetic_frame


def test_multithreading_consistency(clean_frame):
    # Spot finding results should be completely identical regardless of thread count
    params = SpotParams(
        beam_x=clean_frame.cx,
        beam_y=clean_frame.cy,
        pixel_size_x=clean_frame.qx,
        pixel_size_y=clean_frame.qy,
        distance=clean_frame.distance,
        wavelength=clean_frame.wavelength,
    )

    res1 = findspots_data(clean_frame.data, params=params)
    res2 = findspots_data(clean_frame.data, params=params)

    assert res1.count == res2.count
    assert res1.count > 0

    for s1, s2 in zip(res1.spots, res2.spots):
        assert pytest.approx(s1.x, abs=1e-3) == s2.x
        assert pytest.approx(s1.y, abs=1e-3) == s2.y
        assert pytest.approx(s1.intensity, abs=1e-2) == s2.intensity
        assert pytest.approx(s1.snr, abs=1e-3) == s2.snr


def test_multithreading_benchmark(clean_frame):
    params = SpotParams(
        beam_x=clean_frame.cx,
        beam_y=clean_frame.cy,
        pixel_size_x=clean_frame.qx,
        pixel_size_y=clean_frame.qy,
        distance=clean_frame.distance,
        wavelength=clean_frame.wavelength,
    )

    # Warmup
    findspots_data(clean_frame.data, params=params)

    t0 = time.perf_counter()
    n_runs = 5
    for _ in range(n_runs):
        findspots_data(clean_frame.data, params=params)
    avg_time = (time.perf_counter() - t0) / n_runs

    print(f"\n[OpenMP Benchmark] 3072x3072 frame average spot finding time: {avg_time * 1000.0:.2f} ms")
    assert avg_time < 1.0  # Should easily process in < 1 second
