import math
from pathlib import Path
import pytest
import numpy as np
from mxspots.models import SpotParams
from mxspots.spotfinder import findspots, extract_frame_and_params
from mxspots.scorer import score
from mxspots.synthetic import load_synthetic_spec


def _calc_d_spacing(
    x: float,
    y: float,
    cx: float = 1494.0,
    cy: float = 1536.0,
    distance: float = 100.0,
    wavelength: float = 1.7712,
    qx: float = 0.07324,
    qy: float = 0.07324,
) -> float:
    rx = (x - cx) * qx
    ry = (y - cy) * qy
    r = math.sqrt(rx * rx + ry * ry)
    theta = 0.5 * math.atan2(r, distance)
    sin_theta = math.sin(theta)
    return wavelength / (2.0 * sin_theta) if sin_theta > 1e-6 else 999.0


def test_extract_frame_and_params_merges_source_geometry(clean_frame):
    # Pass only resolution filter without detector geometry
    custom_params = SpotParams(d_min=2.0, d_max=5.0, snr_threshold=4.0)
    data, merged_params, detected_index = extract_frame_and_params(clean_frame, custom_params)

    assert data.shape == (clean_frame.ny, clean_frame.nx)
    assert merged_params.d_min == 2.0
    assert merged_params.d_max == 5.0
    assert merged_params.snr_threshold == 4.0
    # Detector geometry should be merged from clean_frame
    assert merged_params.beam_x == pytest.approx(clean_frame.cx)
    assert merged_params.beam_y == pytest.approx(clean_frame.cy)
    assert merged_params.distance == pytest.approx(clean_frame.distance)
    assert merged_params.wavelength == pytest.approx(clean_frame.wavelength)
    assert merged_params.pixel_size_x == pytest.approx(clean_frame.qx)
    assert merged_params.pixel_size_y == pytest.approx(clean_frame.qy)


def test_findspots_with_only_resolution_bounds(test_data_dir):
    yaml_file = test_data_dir / "clean.yaml"
    spec = load_synthetic_spec(yaml_file)

    # Find spots specifying only resolution bounds
    params = SpotParams(d_min=2.0, d_max=4.0)
    spot_list = findspots(yaml_file, params=params)

    assert spot_list.count > 0
    for s in spot_list.spots:
        assert 2.0 - 0.1 <= s.d_spacing <= 4.0 + 0.1
        # Check that d_spacing is computed relative to the true frame beam center and wavelength
        expected_d = _calc_d_spacing(
            s.x,
            s.y,
            cx=spec.cx,
            cy=spec.cy,
            distance=spec.distance,
            wavelength=spec.wavelength,
            qx=spec.qx,
            qy=spec.qy,
        )
        assert s.d_spacing == pytest.approx(expected_d, abs=0.1)


def test_spot_xds_export_resolution_filtering(test_data_dir, tmp_path):
    yaml_file = test_data_dir / "clean.yaml"
    spec = load_synthetic_spec(yaml_file)
    xds_path = tmp_path / "SPOT.XDS"

    params = SpotParams(d_min=2.5, d_max=5.0)
    spot_list = findspots(yaml_file, params=params, xds_output=xds_path, index=3)

    assert xds_path.is_file()
    lines = xds_path.read_text().strip().splitlines()
    assert len(lines) == spot_list.count
    assert len(lines) > 0

    for line in lines:
        parts = line.split()
        assert len(parts) == 4
        x, y, z, intensity = float(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])
        assert z == pytest.approx(2.5)  # 3 - 0.5 = 2.5
        assert intensity > 0.0

        # Verify resolution of coordinates written to SPOT.XDS matches true frame geometry
        d = _calc_d_spacing(
            x,
            y,
            cx=spec.cx,
            cy=spec.cy,
            distance=spec.distance,
            wavelength=spec.wavelength,
            qx=spec.qx,
            qy=spec.qy,
        )
        assert 2.5 - 0.1 <= d <= 5.0 + 0.1


def test_insulin_frame_resolution_filtering(insulin_frame):
    # Test resolution filtering on real detector frame
    params = SpotParams(d_min=3.0, d_max=8.0)
    spot_list = findspots(insulin_frame, params=params, max_spots=100)

    assert spot_list.count > 0
    for s in spot_list.spots:
        assert s.d_spacing >= 3.0 - 0.1
        assert s.d_spacing <= 8.0 + 0.1


def test_score_with_resolution_params(test_data_dir):
    yaml_file = test_data_dir / "clean.yaml"

    params = SpotParams(d_min=1.5, d_max=4.0)
    score_res = score(yaml_file, params=params)
    assert score_res.spot_count > 0
    assert score_res.d_min >= 1.5 - 0.1


def test_cli_dmin_dmax_flags(test_data_dir, monkeypatch, capsys):
    from mxspots.cli import findspots_main
    import json

    yaml_file = str(test_data_dir / "clean.yaml")
    spec = load_synthetic_spec(Path(yaml_file))

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
        expected_d = _calc_d_spacing(
            s["x"],
            s["y"],
            cx=spec.cx,
            cy=spec.cy,
            distance=spec.distance,
            wavelength=spec.wavelength,
            qx=spec.qx,
            qy=spec.qy,
        )
        assert s["d_spacing"] == pytest.approx(expected_d, abs=0.1)
