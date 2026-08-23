import pytest
import numpy as np
from mxspots.models import SpotParams
from mxspots.spotfinder import findspots_data


def test_ccl_ushape_connected_component():
    # Test U-shaped component where two branches merge into a single component
    data = np.zeros((100, 100), dtype=np.float32)

    # U-shape:
    # Left vertical column: y in 20..40, x=30
    # Right vertical column: y in 20..40, x=36
    # Bottom connector: y=40, x in 30..36
    data[20:41, 30:32] = 200.0
    data[20:41, 35:37] = 200.0
    data[39:41, 30:37] = 200.0

    params = SpotParams(
        beam_x=50.0,
        beam_y=50.0,
        snr_threshold=3.0,
        min_spot_area=10,
        max_spot_area=500,
        d_max=0.0,  # Unbounded resolution
    )
    spots = findspots_data(data, params=params)

    # Must resolve into exactly 1 merged spot component
    assert spots.count == 1
    spot = spots.spots[0]
    # Centroid X should be near symmetric center x=33
    assert pytest.approx(33.0, abs=1.0) == spot.x
    # Centroid Y should be below 30 due to bottom heavy bar
    assert 28.0 <= spot.y <= 36.0


def test_ccl_diagonal_8_connectivity():
    # Test 8-connected diagonal stepping
    data = np.zeros((50, 50), dtype=np.float32)
    # Diagonal chain from (10, 10) to (14, 14)
    for i in range(5):
        data[10 + i, 10 + i] = 300.0

    params = SpotParams(
        beam_x=25.0,
        beam_y=25.0,
        snr_threshold=3.0,
        min_spot_area=4,
        max_spot_area=50,
        d_max=0.0,
    )
    spots = findspots_data(data, params=params)

    # Must be 1 connected spot
    assert spots.count == 1
    spot = spots.spots[0]
    assert pytest.approx(12.0, abs=0.5) == spot.x
    assert pytest.approx(12.0, abs=0.5) == spot.y


def test_ccl_multiple_distinct_spots():
    # 4 distinct separate spots
    data = np.zeros((100, 100), dtype=np.float32)
    coords = [(20, 20), (20, 80), (80, 20), (80, 80)]
    for y, x in coords:
        data[y-1:y+2, x-1:x+2] = 400.0

    params = SpotParams(
        beam_x=50.0,
        beam_y=50.0,
        snr_threshold=3.0,
        min_spot_area=4,
        max_spot_area=20,
        d_max=0.0,
    )
    spots = findspots_data(data, params=params)

    assert spots.count == 4
    for spot in spots.spots:
        # Each spot centroid should match one of the coords
        assert any(
            pytest.approx(cx, abs=0.5) == spot.x and pytest.approx(cy, abs=0.5) == spot.y
            for cy, cx in coords
        )
