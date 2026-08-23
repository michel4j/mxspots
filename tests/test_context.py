import pytest
import time
import numpy as np
from mxspots.models import SpotParams
from mxspots.spotfinder import findspots_data, SpotFinderContext


def test_context_consistency(clean_frame):
    params = SpotParams(
        beam_x=clean_frame.cx,
        beam_y=clean_frame.cy,
        pixel_size_x=clean_frame.qx,
        pixel_size_y=clean_frame.qy,
        distance=clean_frame.distance,
        wavelength=clean_frame.wavelength,
    )

    # 1. Standard spot finding without preallocated context
    res_no_ctx = findspots_data(clean_frame.data, params=params)

    # 2. Spot finding with reusable SpotFinderContext
    with SpotFinderContext(max_nx=clean_frame.nx, max_ny=clean_frame.ny) as ctx:
        res_with_ctx = findspots_data(clean_frame.data, params=params, context=ctx)

    assert res_no_ctx.count == res_with_ctx.count
    assert res_no_ctx.count > 0

    for s1, s2 in zip(res_no_ctx.spots, res_with_ctx.spots):
        assert pytest.approx(s1.x, abs=1e-3) == s2.x
        assert pytest.approx(s1.y, abs=1e-3) == s2.y
        assert pytest.approx(s1.intensity, abs=1e-2) == s2.intensity
        assert pytest.approx(s1.snr, abs=1e-3) == s2.snr


def test_context_batch_throughput(clean_frame):
    params = SpotParams(
        beam_x=clean_frame.cx,
        beam_y=clean_frame.cy,
        pixel_size_x=clean_frame.qx,
        pixel_size_y=clean_frame.qy,
        distance=clean_frame.distance,
        wavelength=clean_frame.wavelength,
    )

    n_frames = 10
    with SpotFinderContext(max_nx=clean_frame.nx, max_ny=clean_frame.ny) as ctx:
        t0 = time.perf_counter()
        for _ in range(n_frames):
            res = findspots_data(clean_frame.data, params=params, context=ctx)
            assert res.count > 0
        total_time = time.perf_counter() - t0

    avg_ms = (total_time / n_frames) * 1000.0
    print(f"\n[Zero-Allocation Batch Throughput] 3072x3072 avg time per frame: {avg_ms:.2f} ms")
    assert avg_ms < 150.0


def test_bounded_ring_buffer_overflow_protection():
    # Construct an image with a huge saturated patch that would overflow a naive fixed queue
    data = np.zeros((512, 512), dtype=np.float32)
    data[100:400, 100:400] = 5000.0  # 300x300 = 90,000 pixels (>> BFS_QUEUE_CAPACITY 2048)

    # Inject 3 small valid spots
    data[50:53, 50:53] = 500.0
    data[50:53, 450:453] = 500.0
    data[450:453, 50:53] = 500.0

    params = SpotParams(beam_x=256.0, beam_y=256.0, min_spot_area=2, max_spot_area=50)

    with SpotFinderContext(max_nx=512, max_ny=512) as ctx:
        res = findspots_data(data, params=params, context=ctx)

    # The 90,000 pixel blob is correctly rejected (area > max_spot_area) without queue overflow/crash
    # and the 3 valid small spots are detected
    assert res.count == 3
