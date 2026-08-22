import ctypes
from pathlib import Path
from typing import Optional, Union, Any
import numpy as np
from .models import SpotParams, Spot, SpotList
from ._lib import get_lib, CMxSpotsParams, CMxSpot
from .synthetic import SyntheticFrame, generate_synthetic_frame

import mxio


def findspots_data(
    data: np.ndarray,
    params: Optional[SpotParams] = None,
    max_spots: int = 50000,
) -> SpotList:
    """
    Find spots in a 2D NumPy float32 image array.
    """
    if params is None:
        params = SpotParams()

    # Ensure contiguous 2D float32 array
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

    out_buffer = (CMxSpot * max_spots)()
    num_found = lib.mxspots_find_spots(
        data_ptr,
        nx,
        ny,
        ctypes.byref(c_params),
        out_buffer,
        max_spots,
    )

    actual_count = min(num_found, max_spots)
    spots = [
        Spot(
            x=float(out_buffer[i].x),
            y=float(out_buffer[i].y),
            d_spacing=float(out_buffer[i].d_spacing),
            intensity=float(out_buffer[i].intensity),
            snr=float(out_buffer[i].snr),
        )
        for i in range(actual_count)
    ]

    return SpotList(spots=spots)


def findspots(
    image_source: Union[str, Path, SyntheticFrame, Any],
    params: Optional[SpotParams] = None,
    max_spots: int = 2000,
) -> SpotList:
    """
    Load an image frame and perform spot finding.

    :param image_source: Image source, Supports file paths (.cbf, .h5, .yaml, etc.), mxio
    ImageFrame and DataSet, objects, and SyntheticFrame objects.
    :param params: SpotParams object, defaults to None
    :param max_spots: Maximum number of spots to return, defaults to 2000.
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
        return findspots_data(image_source.data, params, max_spots=max_spots)

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

        return findspots_data(image_source.data, params, max_spots=max_spots)

    raise TypeError(f"Unsupported image source type: {type(image_source)}")
