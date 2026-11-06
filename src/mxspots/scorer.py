import ctypes
from pathlib import Path
from typing import Optional, Union, Any, List

import numpy as np
from .models import SpotParams, Spot, SpotList, ScoreResult, IceRing, DiffractionFrame
from ._lib import get_lib, CMxSpotsParams, CMxSpot, CMxScoreResult


def score_spots(
    spots: Union[SpotList, List[Spot]],
    params: Optional[SpotParams] = None,
    bragg_spots: Optional[int] = None,
    bragg_percent: Optional[float] = None,
    avg_intensity: Optional[float] = None,
    avg_snr: Optional[float] = None,
    num_lattices: Optional[int] = None,
    ice_score: Optional[float] = None,
    num_ice_rings: Optional[int] = None,
) -> ScoreResult:
    """
    Compute quality metrics from a list of detected spots.
    """
    spot_objs = spots.spots if isinstance(spots, SpotList) else list(spots)
    spot_count = len(spot_objs)

    if spot_count == 0:
        return ScoreResult(
            spot_count=0,
            bragg_spots=0,
            bragg_percent=0.0,
            avg_intensity=0.0,
            avg_snr=0.0,
            d_min=999.0,
            ice_score=ice_score if ice_score is not None else 0.0,
            num_ice_rings=num_ice_rings if num_ice_rings is not None else 0,
            num_lattices=0,
            score=0.0,
            ice_rings_detected=None,
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
    if bragg_spots is not None:
        out_score.bragg_spots = ctypes.c_int(bragg_spots)
    if bragg_percent is not None:
        out_score.bragg_percent = ctypes.c_float(bragg_percent)
    if avg_intensity is not None:
        out_score.avg_intensity = ctypes.c_float(avg_intensity)
    if avg_snr is not None:
        out_score.avg_snr = ctypes.c_float(avg_snr)
    if num_lattices is not None:
        out_score.num_lattices = ctypes.c_int(num_lattices)
    if ice_score is not None:
        out_score.ice_score = ctypes.c_float(ice_score)
    if num_ice_rings is not None:
        out_score.num_ice_rings = ctypes.c_int(num_ice_rings)

    c_params = ctypes.byref(CMxSpotsParams.from_params(params)) if params is not None else None

    ret = lib.mxspots_score_spots(
        c_spots,
        spot_count,
        c_params,
        ctypes.byref(out_score),
    )

    if ret != 0:
        raise RuntimeError(f"Error computing spot quality score (code {ret})")

    return ScoreResult(
        spot_count=int(out_score.spot_count),
        bragg_spots=int(out_score.bragg_spots),
        bragg_percent=float(out_score.bragg_percent),
        avg_intensity=float(out_score.avg_intensity),
        avg_snr=float(out_score.avg_snr),
        d_min=float(out_score.d_min),
        ice_score=ice_score if ice_score is not None else float(out_score.ice_score),
        num_ice_rings=num_ice_rings if num_ice_rings is not None else int(out_score.num_ice_rings),
        num_lattices=int(out_score.num_lattices),
        score=float(out_score.score),
        ice_rings_detected=None,
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

    ice_rings_detected = [
        IceRing(
            d_spacing=float(out_score.rings[i].d_spacing),
            d_min=float(out_score.rings[i].d_min),
            d_max=float(out_score.rings[i].d_max),
            score=float(out_score.rings[i].score),
        )
        for i in range(out_score.num_ice_rings)
    ] if out_score.num_ice_rings > 0 else None

    return ScoreResult(
        spot_count=int(out_score.spot_count),
        bragg_spots=int(out_score.bragg_spots),
        bragg_percent=float(out_score.bragg_percent),
        avg_intensity=float(out_score.avg_intensity),
        avg_snr=float(out_score.avg_snr),
        d_min=float(out_score.d_min),
        ice_score=float(out_score.ice_score),
        num_ice_rings=int(out_score.num_ice_rings),
        num_lattices=int(out_score.num_lattices),
        score=float(out_score.score),
        ice_rings_detected=ice_rings_detected,
    )


def score(
    image_source: Union[str, Path, DiffractionFrame, Any],
    params: Optional[SpotParams] = None,
) -> ScoreResult:
    """
    Load an image frame and compute quality score metrics.
    :param image_source: Image source, Supports file paths (.cbf, .h5, .yaml, etc.), mxio
    ImageFrame and DataSet objects, SyntheticFrame objects, and DiffractionFrame.
    :param params: SpotParams object, defaults to None
    """
    frame = DiffractionFrame.from_source(image_source, params)
    return score_data(frame.data, params=frame.params)
