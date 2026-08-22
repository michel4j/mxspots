import pytest
from pathlib import Path
from mxspots.models import Spot, SpotList
from mxspots.spotfinder import findspots


def test_spot_list_to_xds_z_coord(tmp_path):
    spots = [
        Spot(x=2043.29, y=1157.48, d_spacing=4.11, intensity=14435.0, snr=67.2),
        Spot(x=1825.22, y=1928.86, d_spacing=5.19, intensity=13466.0, snr=80.1),
    ]
    spot_list = SpotList(spots=spots)

    xds_file = tmp_path / "SPOT.XDS"
    # Testing frame index 5 -> z = 5 - 0.5 = 4.5
    out_path = spot_list.to_xds(xds_file, frame_index=5)

    assert out_path.is_file()
    lines = out_path.read_text().strip().splitlines()
    assert len(lines) == 2

    # Verify line format: 4 columns
    parts1 = lines[0].split()
    assert len(parts1) == 4
    assert float(parts1[0]) == pytest.approx(2043.29)
    assert float(parts1[1]) == pytest.approx(1157.48)
    assert float(parts1[2]) == pytest.approx(4.5)
    assert float(parts1[3]) == pytest.approx(14435.0)


def test_spot_list_to_xds_default_z(tmp_path):
    spots = [
        Spot(x=100.0, y=200.0, d_spacing=3.0, intensity=5000.0, snr=50.0),
    ]
    spot_list = SpotList(spots=spots)

    xds_file = tmp_path / "SPOT.XDS"
    # Default z should be 0.5 (frame 1)
    out_path = spot_list.to_xds(xds_file)

    lines = out_path.read_text().strip().splitlines()
    parts = lines[0].split()
    assert float(parts[2]) == pytest.approx(0.5)


def test_findspots_cli_index_flag(test_data_dir, tmp_path, monkeypatch, capsys):
    from mxspots.cli import findspots_main

    yaml_file = str(test_data_dir / "clean.yaml")
    output_xds = tmp_path / "SPOT.XDS"

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("sys.argv", ["mxspots.findspots", yaml_file, "--xds", "--index", "10", "--max-spots", "5"])

    findspots_main()

    assert output_xds.is_file()
    lines = output_xds.read_text().strip().splitlines()
    assert len(lines) == 5

    parts = lines[0].split()
    # z should be 10 - 0.5 = 9.5
    assert float(parts[2]) == pytest.approx(9.5)
