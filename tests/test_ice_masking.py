import pytest
import json
import numpy as np
from pathlib import Path
from mxspots.models import SpotParams, SpotList, ScoreResult
from mxspots.spotfinder import findspots, findspots_data
from mxspots.scorer import score, score_data
from mxspots.indexer import index, index_data
from mxspots.synthetic import add_powder_ring


def test_ice_masking_filters_powder_ring_spots():
    # Construct a frame with geometry where 3.897 A falls at ~280 px radius
    rng = np.random.default_rng(42)
    ny, nx = 1024, 1024
    cx, cy = 512.0, 512.0
    qx, qy = 0.075, 0.075
    distance = 80.0
    wavelength = 1.0
    d_target = 3.897

    # At distance=80, lambda=1.0, d=3.897:
    # 2*theta = 2 * asin(1.0 / (2 * 3.897)) = 14.74 deg -> r_mm = 80 * tan(14.74) = 21.05 mm -> r_px = 280.6 px
    r_target_px = 280.6

    data = rng.normal(loc=10.0, scale=2.0, size=(ny, nx)).astype(np.float32)

    # Add 4 true spots far outside the ice ring (r_px ~ 450 px, d ~ 2.5 A)
    outside_spots = [
        (int(cy - 450), int(cx)),
        (int(cy + 450), int(cx)),
        (int(cy), int(cx - 450)),
        (int(cy), int(cx + 450)),
    ]
    for y, x in outside_spots:
        data[y-1:y+2, x-1:x+2] += 400.0

    # Add 4 spots located right on the ice ring (r_px ~ 280 px, d ~ 3.897 A)
    on_ring_spots = [
        (int(cy - r_target_px), int(cx)),
        (int(cy + r_target_px), int(cx)),
        (int(cy), int(cx - r_target_px)),
        (int(cy), int(cx + r_target_px)),
    ]
    for y, x in on_ring_spots:
        data[y-1:y+2, x-1:x+2] += 400.0

    # Add powder ring at 3.897 A
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

    # 1. With ice_mask=True (default)
    params_masked = SpotParams(
        snr_threshold=5.0,
        beam_x=cx,
        beam_y=cy,
        pixel_size_x=qx,
        pixel_size_y=qy,
        distance=distance,
        wavelength=wavelength,
        ice_mask=True,
        ice_sensitivity=3.0,
    )
    spots_masked = findspots_data(data, params=params_masked)

    # 2. With ice_mask=False
    params_unmasked = SpotParams(
        snr_threshold=5.0,
        beam_x=cx,
        beam_y=cy,
        pixel_size_x=qx,
        pixel_size_y=qy,
        distance=distance,
        wavelength=wavelength,
        ice_mask=False,
    )
    spots_unmasked = findspots_data(data, params=params_unmasked)

    # Unmasked spot finding detects all 8 spots (outside + on ring)
    assert spots_unmasked.count == 8
    # Masked spot finding masks out the 3.897 A shell, leaving exactly the 4 outside spots
    assert spots_masked.count == 4
    assert len(spots_masked.ice_rings) >= 1
    assert any(pytest.approx(d_target, abs=0.1) == r.d_spacing for r in spots_masked.ice_rings)
    # Check that none of the detected spots in masked mode are near 3.897 A
    for s in spots_masked.spots:
        assert abs(s.d_spacing - d_target) > 0.5


def test_scorer_ice_metrics(ice_frame):
    # Running score() on ice_frame with ice_mask=True
    res = score(ice_frame, params=SpotParams(ice_mask=True, ice_sensitivity=2.5))
    assert isinstance(res, ScoreResult)
    assert res.ice_score > 2.0
    assert len(res.ice_rings_detected) >= 1


def test_findspots_cli_ice_mask_flags(test_data_dir, monkeypatch, capsys):
    from mxspots.cli import findspots_main

    yaml_file = str(test_data_dir / "ice.yaml")

    # 1. Run with default ice masking
    monkeypatch.setattr("sys.argv", ["mxspots.findspots", yaml_file])
    findspots_main()
    captured_masked = capsys.readouterr()
    assert "Found" in captured_masked.out

    # 2. Run with --json output
    monkeypatch.setattr("sys.argv", ["mxspots.findspots", yaml_file, "--json"])
    findspots_main()
    captured_json = capsys.readouterr()
    data = json.loads(captured_json.out)
    assert "spots" in data
    assert "ice_rings" in data


def test_score_cli_ice_mask_flags(test_data_dir, monkeypatch, capsys):
    from mxspots.cli import score_main

    yaml_file = str(test_data_dir / "ice.yaml")

    # Run with default ice masking
    monkeypatch.setattr("sys.argv", ["mxspots.score", yaml_file])
    score_main()
    captured = capsys.readouterr()
    assert "Ice Score:" in captured.out
    assert "Ice Rings Detected:" in captured.out

    # Run with --json
    monkeypatch.setattr("sys.argv", ["mxspots.score", yaml_file, "--json"])
    score_main()
    captured_json = capsys.readouterr()
    data = json.loads(captured_json.out)
    assert "ice_score" in data
    assert "ice_rings_detected" in data
    assert data["ice_score"] > 0.0


def test_index_with_ice_masking(clean_frame):
    # Test index with ice_mask=True on clean frame
    res = index(clean_frame, params=SpotParams(ice_mask=True))
    assert res.total_spot_count > 0
    assert res.percentage_indexed >= 50.0
