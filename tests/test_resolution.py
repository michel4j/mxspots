import pytest
from mxspots.models import SpotParams
from mxspots.spotfinder import findspots
from mxspots.scorer import score


def test_findspots_d_max_filter(clean_frame):
    # Test restricting to high resolution only (e.g. d_min=1.0, d_max=3.0)
    params = SpotParams(
        d_min=1.0,
        d_max=3.0,
        beam_x=clean_frame.cx,
        beam_y=clean_frame.cy,
        pixel_size_x=clean_frame.qx,
        pixel_size_y=clean_frame.qy,
        distance=clean_frame.distance,
        wavelength=clean_frame.wavelength,
    )
    spot_list = findspots(clean_frame, params=params)

    assert spot_list.count > 0
    for s in spot_list.spots:
        assert 1.0 <= s.d_spacing <= 3.0 + 0.1  # Allow small margin near border


def test_findspots_d_min_filter(clean_frame):
    # Test restricting to low resolution only (e.g. d_min=4.0, d_max=20.0)
    params = SpotParams(
        d_min=4.0,
        d_max=20.0,
        beam_x=clean_frame.cx,
        beam_y=clean_frame.cy,
        pixel_size_x=clean_frame.qx,
        pixel_size_y=clean_frame.qy,
        distance=clean_frame.distance,
        wavelength=clean_frame.wavelength,
    )
    spot_list = findspots(clean_frame, params=params)

    assert spot_list.count > 0
    for s in spot_list.spots:
        assert s.d_spacing >= 4.0 - 0.1
        assert s.d_spacing <= 20.0 + 0.1


def test_cli_dmin_dmax_flags(test_data_dir, monkeypatch, capsys):
    from mxspots.cli import findspots_main
    import json

    yaml_file = str(test_data_dir / "clean.yaml")
    monkeypatch.setattr(
        "sys.argv",
        ["mxspots.findspots", yaml_file, "--json", "--dmin", "2.0", "--dmax", "4.0", "--max-spots", "20"]
    )

    findspots_main()
    captured = capsys.readouterr()

    data = json.loads(captured.out)
    assert data["spot_count"] > 0
    for s in data["spots"]:
        assert 2.0 - 0.1 <= s["d_spacing"] <= 4.0 + 0.1
