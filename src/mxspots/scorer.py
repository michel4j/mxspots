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

    reg_d_min = None
    # If regularity analysis not provided but params given, analyze regularity
    if bragg_spots is None and params is not None and spot_count >= 5:
        c_params = CMxSpotsParams.from_params(params)
        c_pct_bragg = ctypes.c_float(0.0)
        c_bragg_cnt = ctypes.c_int(0)
        c_avg_int = ctypes.c_float(0.0)
        c_avg_snr = ctypes.c_float(0.0)
        c_n_lat = ctypes.c_int(0)
        c_d_min = ctypes.c_float(999.0)
        try:
            ret_reg = lib.mxspots_analyze_regularity(
                c_spots,
                spot_count,
                ctypes.byref(c_params),
                ctypes.byref(c_pct_bragg),
                ctypes.byref(c_bragg_cnt),
                ctypes.byref(c_avg_int),
                ctypes.byref(c_avg_snr),
                ctypes.byref(c_n_lat),
                ctypes.byref(c_d_min),
            )
            if ret_reg == 0:
                bragg_percent = float(c_pct_bragg.value)
                bragg_spots = int(c_bragg_cnt.value)
                avg_intensity = float(c_avg_int.value)
                if avg_snr is None:
                    avg_snr = float(c_avg_snr.value)
                num_lattices = int(c_n_lat.value)
                reg_d_min = float(c_d_min.value)
        except Exception:
            pass

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
    if reg_d_min is not None:
        out_score.d_min = ctypes.c_float(reg_d_min)
    if ice_score is not None:
        out_score.ice_score = ctypes.c_float(ice_score)
    if num_ice_rings is not None:
        out_score.num_ice_rings = ctypes.c_int(num_ice_rings)

    ret = lib.mxspots_score_spots(
        c_spots,
        spot_count,
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

    ice_score: Optional[float] = None
    ice_rings_detected: Optional[List[float]] = None

    if params.ice_mask:
        detected_rings, score_val = detect_ice_rings_data(data, params=params)
        ice_score = score_val
        ice_rings_detected = [r.d_spacing for r in detected_rings] if detected_rings else None
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
        bragg_spots=int(out_score.bragg_spots),
        bragg_percent=float(out_score.bragg_percent),
        avg_intensity=float(out_score.avg_intensity),
        avg_snr=float(out_score.avg_snr),
        d_min=float(out_score.d_min),
        ice_score=ice_score if ice_score is not None else float(out_score.ice_score),
        num_ice_rings=int(out_score.num_ice_rings),
        num_lattices=int(out_score.num_lattices),
        score=float(out_score.score),
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
