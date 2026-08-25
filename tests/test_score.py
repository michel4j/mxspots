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
    assert result.bragg_percent > 0.0
    assert result.bragg_spots > 0
    assert result.avg_intensity > 0.0
    assert result.num_lattices >= 1


def test_score_data_empty():
    from mxspots.scorer import score_data

    empty_data = np.zeros((200, 200), dtype=np.float32)
    result = score_data(empty_data)

    assert result.spot_count == 0
    assert result.bragg_spots == 0
    assert result.bragg_percent == pytest.approx(0.0)
    assert result.avg_intensity == pytest.approx(0.0)
    assert result.avg_snr == pytest.approx(0.0)
    assert result.d_min >= 999.0
    assert result.score == pytest.approx(0.0)
    assert result.num_lattices == 0


def test_score_spots_regularity_pipeline():
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
    assert result.bragg_spots >= 0
    assert result.bragg_percent >= 0.0
    assert result.avg_intensity >= 0.0


def test_score_spots_explicit_bragg_values():
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
    result = score_spots(
        spots,
        params=params,
        bragg_spots=25,
        bragg_percent=85.0,
        avg_intensity=950.0,
    )

    assert result.bragg_spots == 25
    assert result.bragg_percent == 85.0
    assert result.avg_intensity == 950.0
    assert result.score > 0.0


def test_score_spots_explicit_avg_snr():
    from mxspots.scorer import score_spots

    spots = [
        Spot(x=1500.0 + i * 20.0, y=1500.0 + i * 20.0, d_spacing=2.5, intensity=1000.0, snr=3.0)
        for i in range(30)
    ]
    # Pass explicit avg_snr=40.0 overriding the spots' nominal SNR of 3.0
    result = score_spots(
        spots,
        bragg_spots=25,
        bragg_percent=83.3,
        avg_snr=40.0,
        avg_intensity=950.0,
    )

    assert result.avg_snr == pytest.approx(40.0)
    assert result.score > 50.0


def test_score_ice_frame_metrics(test_data_dir):
    from mxspots.scorer import score

    yaml_path = test_data_dir / "ice.yaml"
    result = score(yaml_path)

    assert isinstance(result, ScoreResult)
    assert result.ice_score > 0.0
    assert result.ice_rings_detected is not None
    assert len(result.ice_rings_detected) > 0
    assert all(hasattr(r, "d_min") and hasattr(r, "d_max") for r in result.ice_rings_detected)


def test_score_lyso_split_multi_lattice(test_data_dir):
    from mxspots.scorer import score

    yaml_path = test_data_dir / "lyso-split.yaml"
    result = score(yaml_path)

    assert isinstance(result, ScoreResult)
    assert result.spot_count > 100
    assert result.bragg_percent > 60.0
    assert result.bragg_spots > 100
    assert result.avg_intensity > 0.0
    # Multi-lattice frame should detect more than 1 lattice cluster
    assert result.num_lattices > 1


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
    assert "bragg_percent" in data
    assert "bragg_spots" in data
    assert "avg_intensity" in data
    assert "num_lattices" in data
    assert data["spot_count"] > 0


def test_score_cli_text(test_data_dir, capsys, monkeypatch):
    from mxspots.cli import score_main

    yaml_file = str(test_data_dir / "clean.yaml")
    monkeypatch.setattr("sys.argv", ["mxspots.score", yaml_file])

    score_main()
    captured = capsys.readouterr()

    assert "Quality Score for" in captured.out
    assert "Score:" in captured.out
    assert "Spot Count" in captured.out
    assert "Average SNR" in captured.out
    assert "Resolution Limit" in captured.out
    assert "98th percentile" in captured.out
    assert "Bragg Spots:" in captured.out
    assert "Bragg %:" in captured.out
    assert "Avg Intensity:" in captured.out
    assert "Percentage Indexed:" not in captured.out
    assert "Indexed Spots:" not in captured.out


def test_score_cli_split_warning(test_data_dir, capsys, monkeypatch):
    from mxspots.cli import score_main

    yaml_file = str(test_data_dir / "lyso-split.yaml")
    monkeypatch.setattr("sys.argv", ["mxspots.score", yaml_file])

    score_main()
    captured = capsys.readouterr()

    assert "Quality Score for" in captured.out
    assert "Bragg Spots:" in captured.out
    assert "Warning: Multi-lattice / split crystal detected" in captured.out
