import ctypes
from pathlib import Path
from typing import Optional, Union, Any, List
import numpy as np
from .models import SpotParams, Spot, SpotList, IndexResult
from ._lib import get_lib, CMxSpotsParams, CMxSpot, CMxIndexResult
from .synthetic import SyntheticFrame, generate_synthetic_frame

import mxio


def index_spots(
    spots: Union[SpotList, List[Spot]],
    params: Optional[SpotParams] = None,
) -> IndexResult:
    """
    Index a list of detected spots using reciprocal lattice FFT search.
    """
    if params is None:
        params = SpotParams()

    spot_objs = spots.spots if isinstance(spots, SpotList) else list(spots)
    spot_count = len(spot_objs)

    if spot_count == 0:
        return IndexResult(
            unit_cell=[50.0, 50.0, 50.0, 90.0, 90.0, 90.0],
            percentage_indexed=0.0,
            indexed_spot_count=0,
            total_spot_count=0,
        )

    lib = get_lib()
    c_params = CMxSpotsParams.from_params(params)

    c_spots = (CMxSpot * spot_count)()
    for i, s in enumerate(spot_objs):
        c_spots[i].x = ctypes.c_float(s.x)
        c_spots[i].y = ctypes.c_float(s.y)
        c_spots[i].d_spacing = ctypes.c_float(s.d_spacing)
        c_spots[i].intensity = ctypes.c_float(s.intensity)
        c_spots[i].snr = ctypes.c_float(s.snr)

    c_result = CMxIndexResult()
    ret = lib.mxspots_index_spots(
        c_spots,
        spot_count,
        ctypes.byref(c_params),
        ctypes.byref(c_result),
    )

    if ret != 0:
        raise RuntimeError(f"mxspots_index_spots failed with error code {ret}")

    cell = [float(c_result.unit_cell[i]) for i in range(6)]
    return IndexResult(
        unit_cell=cell,
        percentage_indexed=float(c_result.percentage_indexed),
        indexed_spot_count=int(c_result.indexed_spot_count),
        total_spot_count=int(c_result.total_spot_count),
    )


def index_data(
    data: np.ndarray,
    params: Optional[SpotParams] = None,
) -> IndexResult:
    """
    Find spots and perform lattice indexing on a 2D float32 NumPy image array.
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

    c_result = CMxIndexResult()
    ret = lib.mxspots_index_frame(
        data_ptr,
        nx,
        ny,
        ctypes.byref(c_params),
        ctypes.byref(c_result),
    )

    if ret != 0:
        raise RuntimeError(f"mxspots_index_frame failed with error code {ret}")

    cell = [float(c_result.unit_cell[i]) for i in range(6)]
    return IndexResult(
        unit_cell=cell,
        percentage_indexed=float(c_result.percentage_indexed),
        indexed_spot_count=int(c_result.indexed_spot_count),
        total_spot_count=int(c_result.total_spot_count),
    )


def index(
    image_source: Union[str, Path, SyntheticFrame, Any],
    params: Optional[SpotParams] = None,
) -> IndexResult:
    """
    Load an image frame, find diffraction spots, and index the reciprocal lattice.

    :param image_source: Image source (.cbf, .h5, .yaml, mxio DataSet/ImageFrame, or SyntheticFrame)
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
                beam_x=image_source.cx,
                beam_y=image_source.cy,
                pixel_size_x=image_source.qx,
                pixel_size_y=image_source.qy,
                distance=image_source.distance,
                wavelength=image_source.wavelength,
            )
        return index_data(image_source.data, params=params)

    # mxio.DataSet directly
    if isinstance(image_source, mxio.DataSet):
        image_source = image_source.frame

    # mxio.ImageFrame directly
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
        return index_data(image_source.data, params=params)

    if hasattr(image_source, "data") and isinstance(image_source.data, np.ndarray):
        return index_data(image_source.data, params=params)

    raise TypeError(f"Unsupported image source type: {type(image_source)}")
