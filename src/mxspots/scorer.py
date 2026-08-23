import ctypes
from pathlib import Path
from typing import Optional, Union, Any, List

import numpy as np
from .models import SpotParams, Spot, SpotList, ScoreResult
from ._lib import get_lib, CMxSpotsParams, CMxSpot, CMxScoreResult
from .synthetic import SyntheticFrame
from .spotfinder import extract_frame_and_params


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
