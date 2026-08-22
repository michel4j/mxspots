import ctypes
from pathlib import Path
from typing import Optional, Union, Any

import mxio
import numpy as np
from .models import SpotParams, ScoreResult
from ._lib import get_lib, CMxSpotsParams, CMxScoreResult
from .synthetic import SyntheticFrame, generate_synthetic_frame


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
    ImageFrame and DataSet, objects, and SyntheticFrame objects.
    :param params: SpotParams object, defaults to None
    """

    if isinstance(image_source, (str, Path)):
        p = Path(image_source)
        if p.is_file() and p.suffix.lower() in (".yaml", ".yml"):
            image_source = generate_synthetic_frame(p)
        else:
            image_source = mxio.DataSet.new_from_file(p)
        if image_source is None:
            raise ValueError(f"Could not load image source: {p}")

    # SyntheticFrame directly
    if isinstance(image_source, SyntheticFrame):
        if params is None:
            params = SpotParams(
                beam_x=image_source.nx / 2.0,
                beam_y=image_source.ny / 2.0,
                pixel_size_x=image_source.qx,
                pixel_size_y=image_source.qy,
                distance=image_source.distance,
                wavelength=image_source.wavelength,
            )
        return score_data(image_source.data, params=params)

    # mxio.DataSet directly
    if isinstance(image_source, mxio.DataSet):
        image_source = image_source.frame

    # mxio.ImageFrame direcly
    if isinstance(image_source, mxio.ImageFrame):
        if params is None:
            params = SpotParams(
                beam_x=image_source.center.x,
                beam_y=image_source.center.y,
                pixel_size_x=image_source.pixel_size.x,
                pixel_size_y=image_source.pixel_size.y,
                distance=image_source.distance,
                wavelength=image_source.wavelength,
            )

        return score_data(image_source.data, params=params)

    raise TypeError(f"Unsupported image source type: {type(image_source)}")
