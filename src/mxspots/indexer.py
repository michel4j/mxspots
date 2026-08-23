import ctypes
from pathlib import Path
from typing import Optional, Union, Any, List
import numpy as np
from .models import SpotParams, Spot, SpotList, IndexResult
from ._lib import get_lib, CMxSpotsParams, CMxSpot, CMxIndexResult
from .synthetic import SyntheticFrame
from .spotfinder import extract_frame_and_params


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
    data, merged_params, _ = extract_frame_and_params(image_source, params)
    return index_data(data, params=merged_params)
