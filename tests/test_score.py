import pytest
import numpy as np
from mxspots.models import SpotParams, ScoreResult, Spot
from mxspots.synthetic import load_synthetic_spec, generate_synthetic_frame


def test_score_clean_frame(test_data_dir):
    from mxspots.scorer import score

    yaml_path = test_data_dir / "clean.yaml"
    result = score(yaml_path)

    assert isinstance(result, ScoreResult)
    assert result.spot_count > 0
    assert result.avg_snr > 0.0
    assert result.d_min > 0.0
    # On clean diffraction, d_min should be high resolution (e.g. < 5.0 Angstroms)
    assert result.d_min < 5.0
    assert result.score > 50.0
    assert result.score <= 100.0
    assert result.percentage_indexed is not None
    assert result.indexed_spot_count is not None
    assert result.indexed_spot_count > 0


def test_score_data_empty():
    from mxspots.scorer import score_data

    empty_data = np.zeros((200, 200), dtype=np.float32)
    result = score_data(empty_data)

    assert result.spot_count == 0
    assert result.avg_snr == pytest.approx(0.0)
    assert result.d_min >= 999.0
    assert result.score == pytest.approx(0.0)
    assert result.percentage_indexed is None
    assert result.indexed_spot_count is None


def test_score_spots_auto_indexing():
    from mxspots.scorer import score_spots

    spots = [
        Spot(x=1500.0 + i * 20.0, y=1500.0 + i * 20.0, d_spacing=2.5, intensity=1000.0, snr=25.0)
        for i in range(30)
    ]
    params = SpotParams(
        beam_x=1500.0,
        beam_y=1500.0,
        pixel_size_x=0.075,
        pixel_size_y=0.075,
        distance=200.0,
        wavelength=1.0,
    )
    result = score_spots(spots, params=params)

    assert isinstance(result, ScoreResult)
    assert result.spot_count == 30
    assert result.avg_snr == pytest.approx(25.0)
    assert result.d_min == pytest.approx(2.5)
    assert result.score > 0.0


def test_score_cli_json(test_data_dir, capsys, monkeypatch):
    from mxspots.cli import score_main
    import json

    yaml_file = str(test_data_dir / "insulin.yaml")
    monkeypatch.setattr("sys.argv", ["mxspots.score", yaml_file, "--json"])

    score_main()
    captured = capsys.readouterr()

    data = json.loads(captured.out)
    assert "spot_count" in data
    assert "avg_snr" in data
    assert "d_min" in data
    assert "score" in data
    assert data["spot_count"] > 0


def test_score_cli_text(test_data_dir, capsys, monkeypatch):
    from mxspots.cli import score_main

    yaml_file = str(test_data_dir / "clean.yaml")
    monkeypatch.setattr("sys.argv", ["mxspots.score", yaml_file])

    score_main()
    captured = capsys.readouterr()

    assert "Spot Count" in captured.out
    assert "Average SNR" in captured.out
    assert "Resolution Limit" in captured.out
