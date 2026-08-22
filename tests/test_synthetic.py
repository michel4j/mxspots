import pytest
import numpy as np
from pathlib import Path


def test_load_yaml_spots(test_data_dir):
    from mxspots.synthetic import load_synthetic_spec

    yaml_path = test_data_dir / "clean.yaml"
    spec = load_synthetic_spec(yaml_path)

    assert spec.wavelength == pytest.approx(1.7712)
    assert spec.distance == pytest.approx(100.0)
    assert spec.nx == 3072
    assert spec.ny == 3072
    assert spec.qx == pytest.approx(0.07324)
    assert spec.qy == pytest.approx(0.07324)
    assert len(spec.spots) > 0

    first_spot = spec.spots[0]
    # first spot in clean.yaml is [2043.29, 1157.48, 14435.0]
    assert first_spot[0] == pytest.approx(2043.29)
    assert first_spot[1] == pytest.approx(1157.48)
    assert first_spot[2] == pytest.approx(14435.0)


def test_generate_synthetic_frame(test_data_dir):
    from mxspots.synthetic import generate_synthetic_frame

    yaml_path = test_data_dir / "clean.yaml"
    frame = generate_synthetic_frame(yaml_path, max_spots=50, add_noise=False)

    assert isinstance(frame.data, np.ndarray)
    assert frame.data.dtype == np.float32
    assert frame.data.shape == (frame.ny, frame.nx)

    # Spot locations should have high intensity
    for x, y, intensity in frame.spots[:10]:
        ix, iy = int(round(x)), int(round(y))
        if 0 <= iy < frame.ny and 0 <= ix < frame.nx:
            assert frame.data[iy, ix] > 0.0


def test_synthetic_frame_caching(test_data_dir):
    from mxspots.synthetic import get_cached_synthetic_frame

    frame1 = get_cached_synthetic_frame("clean.yaml", max_spots=10)
    frame2 = get_cached_synthetic_frame("clean.yaml", max_spots=10)

    # Should be the exact same object from cache
    assert frame1 is frame2
    assert frame1.data.shape == (3072, 3072)


def test_synthetic_fixture(clean_frame):
    assert clean_frame.data is not None
    assert clean_frame.data.shape == (3072, 3072)
    assert clean_frame.wavelength == pytest.approx(1.7712)
