import ctypes
from pathlib import Path
from typing import Optional, Union, Any, List

import numpy as np
from .models import SpotParams, Spot, SpotList, ScoreResult
from ._lib import get_lib, CMxSpotsParams, CMxSpot, CMxScoreResult
from .synthetic import SyntheticFrame
from .spotfinder import extract_frame_and_params, detect_ice_rings_data


def score_spots(
    spots: Union[SpotList, List[Spot]],
) -> ScoreResult:
    """
    Compute quality metrics from a list of detected spots.
    """
    spot_objs = spots.spots if isinstance(spots, SpotList) else list(spots)
    spot_count = len(spot_objs)

    if spot_count == 0:
        return ScoreResult(
            spot_count=0,
            avg_snr=0.0,
            d_min=999.0,
            percentage_indexed=None,
        )

    lib = get_lib()
    c_spots = (CMxSpot * spot_count)()
    for i, s in enumerate(spot_objs):
        c_spots[i].x = ctypes.c_float(s.x)
        c_spots[i].y = ctypes.c_float(s.y)
        c_spots[i].d_spacing = ctypes.c_float(s.d_spacing)
        c_spots[i].intensity = ctypes.c_float(s.intensity)
        c_spots[i].snr = ctypes.c_float(s.snr)

    out_score = CMxScoreResult()
    ret = lib.mxspots_score_spots(
        c_spots,
        spot_count,
        ctypes.byref(out_score),
    )

    if ret != 0:
        raise RuntimeError(f"Error computing spot quality score (code {ret})")

    return ScoreResult(
        spot_count=int(out_score.spot_count),
        avg_snr=float(out_score.avg_snr),
        d_min=float(out_score.d_min),
        percentage_indexed=float(out_score.percentage_indexed) if out_score.percentage_indexed > 0.0 else None,
    )


def score_data(
    data: np.ndarray,
    params: Optional[SpotParams] = None,
) -> ScoreResult:
    """
    Compute quality metrics for a 2D float32 image array.
    """
    if params is None:
        params = SpotParams()

    if not isinstance(data, np.ndarray):
        data = np.asarray(data, dtype=np.float32)
    elif data.dtype != np.float32 or not data.flags.c_contiguous:
        data = np.ascontiguousarray(data, dtype=np.float32)

    if data.ndim != 2:
        raise ValueError(f"Expected 2D image array, got {data.ndim}D shape {data.shape}")

    ice_score: Optional[float] = None
    ice_rings_detected: Optional[List[float]] = None

    if params.ice_mask:
        detected_rings, score_val = detect_ice_rings_data(data, params=params)
        ice_score = score_val
        ice_rings_detected = [r.d_spacing for r in detected_rings]
        if detected_rings:
            active_masked = list(params.masked_rings) if params.masked_rings is not None else []
            for ring in detected_rings:
                active_masked.append((ring.d_min, ring.d_max))
            params = SpotParams(
                snr_threshold=params.snr_threshold,
                min_spot_area=params.min_spot_area,
                max_spot_area=params.max_spot_area,
                beam_x=params.beam_x,
                beam_y=params.beam_y,
                pixel_size_x=params.pixel_size_x,
                pixel_size_y=params.pixel_size_y,
                distance=params.distance,
                wavelength=params.wavelength,
                d_min=params.d_min,
                d_max=params.d_max,
                ice_mask=params.ice_mask,
                ice_sensitivity=params.ice_sensitivity,
                masked_rings=active_masked,
            )

    ny, nx = data.shape
    lib = get_lib()

    c_params = CMxSpotsParams.from_params(params)
    data_ptr = data.ctypes.data_as(ctypes.POINTER(ctypes.c_float))

    out_score = CMxScoreResult()
    ret = lib.mxspots_score_frame(
        data_ptr,
        nx,
        ny,
        ctypes.byref(c_params),
        ctypes.byref(out_score),
    )

    if ret != 0:
        raise RuntimeError(f"Error computing frame quality score (code {ret})")

    return ScoreResult(
        spot_count=int(out_score.spot_count),
        avg_snr=float(out_score.avg_snr),
        d_min=float(out_score.d_min),
        percentage_indexed=float(out_score.percentage_indexed) if out_score.percentage_indexed > 0.0 else None,
        ice_score=ice_score,
        ice_rings_detected=ice_rings_detected,
    )


def score(
    image_source: Union[str, Path, SyntheticFrame, Any],
    params: Optional[SpotParams] = None,
) -> ScoreResult:
    """
    Load an image frame and compute quality score metrics.
    :param image_source: Image source, Supports file paths (.cbf, .h5, .yaml, etc.), mxio
    ImageFrame and DataSet objects, and SyntheticFrame objects.
    :param params: SpotParams object, defaults to None
    """
    data, merged_params, _ = extract_frame_and_params(image_source, params)
    return score_data(data, params=merged_params)
