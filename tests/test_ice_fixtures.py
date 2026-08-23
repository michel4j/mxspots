import math
import numpy as np
import pytest
from mxspots.models import IceRing, SpotParams, SpotList, Spot, ScoreResult
from mxspots._lib import CMxSpotsParams, get_lib
from mxspots.synthetic import add_powder_ring, SyntheticFrame


def test_ice_ring_dataclass():
    ring = IceRing(d_spacing=3.897, d_min=3.85, d_max=3.95, score=4.5)
    assert ring.d_spacing == pytest.approx(3.897)
    assert ring.d_min == pytest.approx(3.85)
    assert ring.d_max == pytest.approx(3.95)
    assert ring.score == pytest.approx(4.5)
    d = ring.to_dict()
    assert d["d_spacing"] == pytest.approx(3.897)
    assert d["score"] == pytest.approx(4.5)


def test_spot_params_ice_defaults():
    params = SpotParams()
    assert params.ice_mask is True
    assert params.ice_sensitivity == pytest.approx(1.0)
    assert params.masked_rings is None

    custom = SpotParams(ice_mask=False, ice_sensitivity=5.0, masked_rings=[(3.6, 3.7), (3.85, 3.95)])
    assert custom.ice_mask is False
    assert custom.ice_sensitivity == pytest.approx(5.0)
    assert custom.masked_rings == [(3.6, 3.7), (3.85, 3.95)]


def test_spot_list_ice_metadata():
    spots = [Spot(x=100.0, y=100.0, d_spacing=3.0, intensity=500.0, snr=5.0)]
    ice_rings = [IceRing(d_spacing=3.897, d_min=3.85, d_max=3.95, score=4.2)]
    slist = SpotList(spots=spots, ice_rings=ice_rings)
    assert slist.count == 1
    assert slist.ice_rings is not None
    assert len(slist.ice_rings) == 1

    d = slist.to_dict()
    assert "ice_rings" in d
    assert len(d["ice_rings"]) == 1
    assert d["ice_rings"][0]["d_spacing"] == pytest.approx(3.897)

    json_str = slist.to_json()
    assert "3.897" in json_str


def test_score_result_ice_metadata():
    res = ScoreResult(
        spot_count=50,
        avg_snr=6.5,
        d_min=1.8,
        percentage_indexed=90.0,
        ice_score=4.8,
        ice_rings_detected=[3.897, 3.669],
    )
    assert res.ice_score == pytest.approx(4.8)
    assert res.ice_rings_detected == [3.897, 3.669]

    d = res.to_dict()
    assert d["ice_score"] == pytest.approx(4.8)
    assert d["ice_rings_detected"] == [3.897, 3.669]
    json_str = res.to_json()
    assert "ice_score" in json_str


def test_c_mxspots_params_masked_rings_conversion():
    # Distance = 200 mm, Wavelength = 1.0 A
    # For d = 3.897 A:
    # 2*theta = 2 * asin(lambda / (2*d)) = 2 * asin(1.0 / (2 * 3.897)) = 2 * asin(0.1283038) = 0.25732 rad
    # r = distance * tan(2*theta) = 200 * tan(0.25732) = 200 * 0.26315 = 52.63 mm
    # r^2 = 2770 mm^2
    params = SpotParams(
        distance=200.0,
        wavelength=1.0,
        masked_rings=[(3.85, 3.95)],
    )
    c_params = CMxSpotsParams.from_params(params)
    assert c_params.num_masked_rings == 1

    # Check that min_r2 and max_r2 are positive and min_r2 < max_r2
    # Note: smaller d corresponds to larger r, so d_min=3.85 -> max_r, d_max=3.95 -> min_r
    min_r2 = c_params.masked_rings_r2[0][0]
    max_r2 = c_params.masked_rings_r2[0][1]
    assert min_r2 > 0.0
    assert max_r2 > min_r2
    assert 2600.0 < min_r2 < 3000.0
    assert 2600.0 < max_r2 < 3000.0


def test_add_powder_ring_synthetic():
    data = np.zeros((2000, 2000), dtype=np.float32)
    cx, cy = 1000.0, 1000.0
    qx, qy = 0.075, 0.075
    distance = 200.0
    wavelength = 1.0
    d_spacing = 3.897

    # Calculate expected radius in pixels
    theta = math.asin(wavelength / (2.0 * d_spacing))
    r_mm = distance * math.tan(2.0 * theta)
    r_px = r_mm / qx

    out = add_powder_ring(
        data,
        cx=cx,
        cy=cy,
        qx=qx,
        qy=qy,
        distance=distance,
        wavelength=wavelength,
        d_spacing=d_spacing,
        radial_width=2.0,
        peak_intensity=100.0,
    )

    # Pixel exactly on the ring radius should have peak intensity
    sample_x = int(round(cx + r_px))
    sample_y = int(round(cy))
    assert out[sample_y, sample_x] > 50.0

    # Center pixel should have zero added intensity
    assert out[int(cy), int(cx)] == pytest.approx(0.0, abs=1e-3)
