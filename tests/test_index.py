import pytest
import json
import numpy as np
from pathlib import Path
from mxspots.models import SpotParams, IndexResult
from mxspots.indexer import index, index_data, index_spots
from mxspots.spotfinder import findspots


def test_index_clean_frame(clean_frame):
    # Perform indexing on clean synthetic diffraction frame
    res = index(clean_frame)

    assert isinstance(res, IndexResult)
    assert len(res.unit_cell) == 6
    a, b, c, alpha, beta, gamma = res.unit_cell

    # Ensure cell lengths and angles are positive and physically reasonable
    assert 10.0 <= a <= 300.0
    assert 10.0 <= b <= 300.0
    assert 10.0 <= c <= 300.0
    assert 30.0 <= alpha <= 150.0
    assert 30.0 <= beta <= 150.0
    assert 30.0 <= gamma <= 150.0

    # Clean frame should have high percentage of indexed spots
    assert res.total_spot_count > 0
    assert res.indexed_spot_count > 0
    assert res.percentage_indexed >= 50.0


def test_index_insulin_frame(insulin_frame):
    res = index(insulin_frame)
    assert isinstance(res, IndexResult)
    assert res.total_spot_count > 0
    assert res.percentage_indexed > 0.0


def test_index_data_numpy():
    data = np.zeros((512, 512), dtype=np.float32)
    # Add a simple periodic 2D grid of spots centered around (256, 256)
    for y in range(106, 410, 30):
        for x in range(106, 410, 30):
            data[y-1:y+2, x-1:x+2] = 500.0

    params = SpotParams(beam_x=256.0, beam_y=256.0, distance=150.0, wavelength=1.0)
    res = index_data(data, params=params)

    assert isinstance(res, IndexResult)
    assert res.total_spot_count > 0
    assert res.percentage_indexed >= 50.0


def test_index_spots_empty():
    params = SpotParams()
    res = index_spots([], params=params)

    assert isinstance(res, IndexResult)
    assert res.total_spot_count == 0
    assert res.indexed_spot_count == 0
    assert res.percentage_indexed == 0.0


def test_index_cli(test_data_dir, monkeypatch, capsys):
    from mxspots.cli import index_main

    yaml_file = str(test_data_dir / "clean.yaml")
    monkeypatch.setattr("sys.argv", ["mxspots.index", yaml_file])

    index_main()
    captured = capsys.readouterr()

    assert "Unit Cell:" in captured.out
    assert "Percentage Indexed:" in captured.out


def test_index_cli_json(test_data_dir, monkeypatch, capsys):
    from mxspots.cli import index_main

    yaml_file = str(test_data_dir / "clean.yaml")
    monkeypatch.setattr("sys.argv", ["mxspots.index", yaml_file, "--json"])

    index_main()
    captured = capsys.readouterr()

    data = json.loads(captured.out)
    assert "unit_cell" in data
    assert len(data["unit_cell"]) == 6
    assert "percentage_indexed" in data
    assert data["percentage_indexed"] > 0
