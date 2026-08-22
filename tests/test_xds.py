import pytest
from pathlib import Path
from mxspots.models import Spot, SpotList
from mxspots.spotfinder import findspots


def test_spot_list_to_xds(tmp_path):
    spots = [
        Spot(x=2043.29, y=1157.48, d_spacing=4.11, intensity=14435.0, snr=67.2),
        Spot(x=1825.22, y=1928.86, d_spacing=5.19, intensity=13466.0, snr=80.1),
    ]
    spot_list = SpotList(spots=spots)

    xds_file = tmp_path / "SPOT.XDS"
    out_path = spot_list.to_xds(xds_file, angle=45.0)

    assert out_path.is_file()
    lines = out_path.read_text().strip().splitlines()
    assert len(lines) == 2

    # Verify line format: 4 columns
    parts1 = lines[0].split()
    assert len(parts1) == 4
    assert float(parts1[0]) == pytest.approx(2043.29)
    assert float(parts1[1]) == pytest.approx(1157.48)
    assert float(parts1[2]) == pytest.approx(45.0)
    assert float(parts1[3]) == pytest.approx(14435.0)


def test_findspots_cli_xds_flag(test_data_dir, tmp_path, monkeypatch, capsys):
    from mxspots.cli import findspots_main

    yaml_file = str(test_data_dir / "clean.yaml")
    output_xds = tmp_path / "SPOT.XDS"

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("sys.argv", ["mxspots.findspots", yaml_file, "--xds", "--max-spots", "20"])

    findspots_main()

    assert output_xds.is_file()
    lines = output_xds.read_text().strip().splitlines()
    assert len(lines) == 20

    captured = capsys.readouterr()
    assert "Found" in captured.out
