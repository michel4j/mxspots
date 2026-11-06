import ctypes
from pathlib import Path
from typing import Optional, Union, Any, Tuple, List
import numpy as np
from .models import SpotParams, Spot, SpotList, IceRing, DiffractionFrame
from ._lib import get_lib, CMxSpotsParams, CMxSpot, CMxIceResult
from .synthetic import SyntheticFrame


class SpotFinderContext:
    """
    Reusable execution scratch context for zero-allocation spot finding loops.
    """

    def __init__(self, max_nx: int = 3100, max_ny: int = 3100):
        self.max_nx = max_nx
        self.max_ny = max_ny
        self._lib = get_lib()
        self._handle = self._lib.mxspots_create_context(max_nx, max_ny)
        if not self._handle:
            raise MemoryError(f"Failed to allocate SpotFinderContext for {max_nx}x{max_ny}")

    @property
    def handle(self):
        return self._handle

    def close(self):
        if self._handle is not None:
            self._lib.mxspots_free_context(self._handle)
            self._handle = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()


def extract_frame_and_params(
    image_source: Union[str, Path, SyntheticFrame, Any],
    params: Optional[SpotParams] = None,
) -> Tuple[np.ndarray, SpotParams, float]:
    """
    Deprecated procedural helper. Use DiffractionFrame.from_source instead.
    """
    frame = DiffractionFrame.from_source(image_source, params)
    return frame.data, frame.params, frame.frame_index


def findspots_data(
    data: np.ndarray,
    params: Optional[SpotParams] = None,
    max_spots: int = 50000,
    context: Optional[SpotFinderContext] = None,
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
    out_ice = CMxIceResult()

    if context is not None and context.handle is not None and nx <= context.max_nx and ny <= context.max_ny:
        num_found = lib.mxspots_find_spots_ctx(
            context.handle,
            data_ptr,
            nx,
            ny,
            ctypes.byref(c_params),
            out_buffer,
            max_spots,
            ctypes.byref(out_ice),
        )
    else:
        num_found = lib.mxspots_find_spots(
            data_ptr,
            nx,
            ny,
            ctypes.byref(c_params),
            out_buffer,
            max_spots,
            ctypes.byref(out_ice),
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

    detected_ice_rings: Optional[List[IceRing]] = None
    if params.ice_mask and out_ice.num_rings > 0:
        detected_ice_rings = [
            IceRing(
                d_spacing=float(out_ice.rings[i].d_spacing),
                d_min=float(out_ice.rings[i].d_min),
                d_max=float(out_ice.rings[i].d_max),
                score=float(out_ice.rings[i].score),
            )
            for i in range(out_ice.num_rings)
        ]

    return SpotList(spots=spots, ice_rings=detected_ice_rings)


def findspots(
    image_source: Union[str, Path, SyntheticFrame, Any],
    params: Optional[SpotParams] = None,
    max_spots: int = 2000,
    xds_output: Optional[Union[str, Path]] = None,
    index: Optional[int] = None,
    z: Optional[float] = None,
    context: Optional[SpotFinderContext] = None,
) -> SpotList:
    """
    Load an image frame and perform spot finding.

    :param image_source: Image source, Supports file paths (.cbf, .h5, .yaml, etc.), mxio
    ImageFrame and DataSet objects, SyntheticFrame objects, and DiffractionFrame.
    :param params: SpotParams object, defaults to None
    :param max_spots: Maximum number of spots to return, defaults to 2000.
    :param xds_output: Optional path to export detected spots to SPOT.XDS format.
    :param index: Frame index (1-based), sets Z coordinate to index - 0.5 in SPOT.XDS.
    :param z: Explicit Z continuous frame coordinate in SPOT.XDS (overrides index).
    :param context: Optional pre-allocated SpotFinderContext for zero-allocation execution.
    """
    frame = DiffractionFrame.from_source(image_source, params)
    spot_list = findspots_data(frame.data, frame.params, max_spots=max_spots, context=context)

    if xds_output is not None:
        if z is not None:
            effective_z = float(z)
        elif index is not None:
            effective_z = float(index) - 0.5
        else:
            effective_z = float(frame.frame_index) - 0.5

        spot_list.to_xds(xds_output, z=effective_z)

    return spot_list


def detect_ice_rings_data(
    data: np.ndarray,
    params: Optional[SpotParams] = None,
) -> Tuple[List[IceRing], float]:
    """
    Detect azimuthal ice powder rings on a 2D NumPy float32 image array.

    :param data: 2D image data array.
    :param params: SpotParams containing detector geometry and ice sensitivity.
    :return: Tuple of (detected IceRing list, overall ice_score).
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
    c_params.snr_threshold = ctypes.c_float(params.ice_sensitivity)

    data_ptr = data.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    out_result = CMxIceResult()

    lib.mxspots_detect_ice(
        data_ptr,
        nx,
        ny,
        ctypes.byref(c_params),
        ctypes.byref(out_result),
    )

    rings = [
        IceRing(
            d_spacing=float(out_result.rings[i].d_spacing),
            d_min=float(out_result.rings[i].d_min),
            d_max=float(out_result.rings[i].d_max),
            score=float(out_result.rings[i].score),
        )
        for i in range(out_result.num_rings)
    ]
    return rings, float(out_result.ice_score)


def detect_ice_rings(
    image_source: Union[str, Path, SyntheticFrame, Any, np.ndarray],
    params: Optional[SpotParams] = None,
) -> Tuple[List[IceRing], float]:
    """
    Detect azimuthal ice powder rings on an image frame or raw numpy data.

    :param image_source: Image file path, mxio DataSet/ImageFrame, SyntheticFrame, DiffractionFrame, or 2D numpy array.
    :param params: SpotParams with detector geometry and ice sensitivity.
    :return: Tuple of (detected IceRing list, overall ice_score).
    """
    if isinstance(image_source, np.ndarray):
        return detect_ice_rings_data(image_source, params=params)

    frame = DiffractionFrame.from_source(image_source, params)
    return detect_ice_rings_data(frame.data, params=frame.params)
