import pytest
import numpy as np
from mxspots.models import SpotParams, SpotList, Spot
from mxspots.synthetic import get_cached_synthetic_frame, generate_synthetic_frame, load_synthetic_spec


def test_findspots_synthetic_clean(test_data_dir):
    from mxspots.spotfinder import findspots_data

    spec = load_synthetic_spec(test_data_dir / "clean.yaml")
    # Generate frame with top 30 spots
    synth = generate_synthetic_frame(spec, max_spots=30, add_noise=False)

    params = SpotParams(
        snr_threshold=3.0,
        min_spot_area=2,
        max_spot_area=500,
        beam_x=1536.0,  # Center of 3072x3072
        beam_y=1536.0,
        pixel_size_x=spec.qx,
        pixel_size_y=spec.qy,
        distance=spec.distance,
        wavelength=spec.wavelength,
    )

    spot_list = findspots_data(synth.data, params)
    assert isinstance(spot_list, SpotList)
    assert spot_list.count >= 25

    # Check top detected spot matches ground truth
    top_detected = spot_list.spots[0]
    # Ground truth top spot: [2043.29, 1157.48, 14435.0]
    assert top_detected.x == pytest.approx(2043.29, abs=1.5) or top_detected.x == pytest.approx(1825.22, abs=1.5)
    assert top_detected.intensity > 10000.0
    assert top_detected.d_spacing > 0.0


def test_findspots_insulin(insulin_frame):
    from mxspots.spotfinder import findspots

    spot_list = findspots(insulin_frame, max_spots=50)
    assert isinstance(spot_list, SpotList)
    assert spot_list.count > 0
    # Spot centroids should be within image frame bounds
    for s in spot_list.spots:
        assert 0 <= s.x < insulin_frame.nx
        assert 0 <= s.y < insulin_frame.ny
        assert s.d_spacing > 0.0


def test_findspots_empty_image():
    from mxspots.spotfinder import findspots_data

    empty_data = np.zeros((100, 100), dtype=np.float32)
    spot_list = findspots_data(empty_data)
    assert spot_list.count == 0


def test_findspots_cli_json(test_data_dir, capsys, monkeypatch):
    from mxspots.cli import findspots_main
    import json

    yaml_file = str(test_data_dir / "clean.yaml")
    monkeypatch.setattr("sys.argv", ["mxspots.findspots", yaml_file, "--json", "--max-spots", "10"])
    
    findspots_main()
    captured = capsys.readouterr()
    
    data = json.loads(captured.out)
    assert "spot_count" in data
    assert "spots" in data
    assert len(data["spots"]) > 0


def test_findspots_cli_text(test_data_dir, capsys, monkeypatch):
    from mxspots.cli import findspots_main

    yaml_file = str(test_data_dir / "clean.yaml")
    monkeypatch.setattr("sys.argv", ["mxspots.findspots", yaml_file, "--max-spots", "10"])
    
    findspots_main()
    captured = capsys.readouterr()
    assert "Found" in captured.out
    assert "spots" in captured.out
