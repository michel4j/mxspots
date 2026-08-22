import ctypes
from pathlib import Path
from typing import Optional, Union, Any
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
    Supports file paths (.cbf, .h5, .yaml, etc.), mxio ImageFrame objects, and SyntheticFrame objects.
    """
    # 1. SyntheticFrame directly
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

    # 2. String or Path object
    if isinstance(image_source, (str, Path)):
        p = Path(image_source)
        if not p.is_file():
            raise FileNotFoundError(f"Image file not found: {p}")

        if p.suffix.lower() in (".yaml", ".yml"):
            synth = generate_synthetic_frame(p)
            return score(synth, params=params)

        import mxio
        frame = mxio.read_image(str(p))
        if frame is None:
            raise ValueError(f"Could not read image file with mxio: {p}")

        img_data = frame.data
        beam_x = getattr(frame.center, "x", 0.0) if hasattr(frame, "center") else 0.0
        beam_y = getattr(frame.center, "y", 0.0) if hasattr(frame, "center") else 0.0
        pixel_x = getattr(frame.pixel_size, "x", 0.075) if hasattr(frame, "pixel_size") else 0.075
        pixel_y = getattr(frame.pixel_size, "y", 0.075) if hasattr(frame, "pixel_size") else 0.075
        distance = getattr(frame, "distance", 200.0)
        wavelength = getattr(frame, "wavelength", 1.0)

        if params is None:
            params = SpotParams(
                beam_x=beam_x,
                beam_y=beam_y,
                pixel_size_x=pixel_x,
                pixel_size_y=pixel_y,
                distance=distance,
                wavelength=wavelength,
            )

        return score_data(img_data, params=params)

    # 3. Object with .data attribute (e.g. mxio ImageFrame directly)
    if hasattr(image_source, "data") and isinstance(image_source.data, np.ndarray):
        beam_x = getattr(image_source.center, "x", 0.0) if hasattr(image_source, "center") else 0.0
        beam_y = getattr(image_source.center, "y", 0.0) if hasattr(image_source, "center") else 0.0
        pixel_x = getattr(image_source.pixel_size, "x", 0.075) if hasattr(image_source, "pixel_size") else 0.075
        pixel_y = getattr(image_source.pixel_size, "y", 0.075) if hasattr(image_source, "pixel_size") else 0.075
        distance = getattr(image_source, "distance", 200.0)
        wavelength = getattr(image_source, "wavelength", 1.0)

        if params is None:
            params = SpotParams(
                beam_x=beam_x,
                beam_y=beam_y,
                pixel_size_x=pixel_x,
                pixel_size_y=pixel_y,
                distance=distance,
                wavelength=wavelength,
            )

        return score_data(image_source.data, params=params)

    raise TypeError(f"Unsupported image source type: {type(image_source)}")
