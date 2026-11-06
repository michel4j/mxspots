import pytest
import numpy as np
from mxspots.models import DiffractionFrame, SpotParams, SpotList, ScoreResult
from mxspots.spotfinder import findspots, extract_frame_and_params
from mxspots.scorer import score
from mxspots.synthetic import get_cached_synthetic_frame


def test_diffraction_frame_from_yaml(test_data_dir):
    yaml_path = test_data_dir / "clean.yaml"
    frame = DiffractionFrame.from_source(yaml_path)

    assert isinstance(frame, DiffractionFrame)
    assert isinstance(frame.data, np.ndarray)
    assert frame.data.dtype == np.float32
    assert frame.data.ndim == 2
    assert frame.params.distance > 0.0
    assert frame.params.wavelength > 0.0
    assert frame.params.pixel_size_x > 0.0
    assert frame.frame_index == 1.0


def test_diffraction_frame_from_synthetic():
    synth = get_cached_synthetic_frame("clean.yaml")
    frame = DiffractionFrame.from_source(synth)

    assert isinstance(frame, DiffractionFrame)
    assert frame.data.shape == synth.data.shape
    assert frame.params.beam_x == pytest.approx(synth.cx)
    assert frame.params.beam_y == pytest.approx(synth.cy)
    assert frame.params.distance == pytest.approx(synth.distance)
    assert frame.params.wavelength == pytest.approx(synth.wavelength)


def test_diffraction_frame_from_existing_frame_override_params():
    synth = get_cached_synthetic_frame("clean.yaml")
    base_frame = DiffractionFrame.from_source(synth)

    custom_params = SpotParams(snr_threshold=12.0, d_max=15.0)
    overridden = DiffractionFrame.from_source(base_frame, params=custom_params)

    assert overridden.params.snr_threshold == 12.0
    assert overridden.params.d_max == 15.0
    # Detector geometry preserved from base
    assert overridden.params.beam_x == base_frame.params.beam_x
    assert overridden.params.distance == base_frame.params.distance


def test_findspots_with_diffraction_frame(test_data_dir):
    frame = DiffractionFrame.from_source(test_data_dir / "clean.yaml")
    spots = findspots(frame, max_spots=30)

    assert isinstance(spots, SpotList)
    assert spots.count > 0


def test_score_with_diffraction_frame(test_data_dir):
    frame = DiffractionFrame.from_source(test_data_dir / "clean.yaml")
    res = score(frame)

    assert isinstance(res, ScoreResult)
    assert res.spot_count > 0
    assert res.score > 0.0


def test_extract_frame_and_params_deprecated_compatibility(test_data_dir):
    data, params, idx = extract_frame_and_params(test_data_dir / "clean.yaml")

    assert isinstance(data, np.ndarray)
    assert isinstance(params, SpotParams)
    assert idx == 1.0


def test_ice_rings_metadata_in_score_result(test_data_dir):
    ice_path = test_data_dir / "ice.yaml"
    frame = DiffractionFrame.from_source(ice_path)
    res = score(frame)

    assert isinstance(res, ScoreResult)
    if res.num_ice_rings > 0:
        assert res.ice_rings_detected is not None
        assert len(res.ice_rings_detected) == res.num_ice_rings
        for ring in res.ice_rings_detected:
            assert ring.d_spacing > 0.0
            assert ring.d_min > 0.0
            assert ring.d_max > ring.d_min
