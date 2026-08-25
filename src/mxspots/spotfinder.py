import ctypes
from pathlib import Path
from typing import Optional, Union, Any, Tuple, List
import numpy as np
from .models import SpotParams, Spot, SpotList, IceRing
from ._lib import get_lib, CMxSpotsParams, CMxSpot, CMxIceResult
from .synthetic import SyntheticFrame, generate_synthetic_frame

import mxio


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
    Extract 2D numpy data array, merged SpotParams with metadata from image_source,
    and detected frame continuous index (1-based).
    """
    detected_index = 1.0

    if isinstance(image_source, (str, Path)):
        p = Path(image_source)
        if p.is_file() and p.suffix.lower() in (".yaml", ".yml"):
            image_source = generate_synthetic_frame(p)
        else:
            image_source = mxio.DataSet.new_from_file(p)
        if image_source is None:
            raise ValueError(f"Could not load image source: {p}")

    if isinstance(image_source, mxio.DataSet):
        detected_index = float(image_source.index)
        image_source = image_source.frame

    if isinstance(image_source, SyntheticFrame):
        src_bx = float(image_source.cx)
        src_by = float(image_source.cy)
        src_qx = float(image_source.qx)
        src_qy = float(image_source.qy)
        src_dist = float(image_source.distance)
        src_wvl = float(image_source.wavelength)
        data = image_source.data
    elif isinstance(image_source, mxio.ImageFrame):
        src_bx = float(image_source.center.x)
        src_by = float(image_source.center.y)
        src_qx = float(image_source.pixel_size.x)
        src_qy = float(image_source.pixel_size.y)
        src_dist = float(image_source.distance)
        src_wvl = float(image_source.wavelength)
        data = image_source.data
    elif hasattr(image_source, "data") and isinstance(image_source.data, np.ndarray):
        data = image_source.data
        src_bx = float(getattr(image_source, "cx", 0.0))
        src_by = float(getattr(image_source, "cy", 0.0))
        src_qx = float(getattr(image_source, "qx", 0.075))
        src_qy = float(getattr(image_source, "qy", 0.075))
        src_dist = float(getattr(image_source, "distance", 200.0))
        src_wvl = float(getattr(image_source, "wavelength", 1.0))
    else:
        raise TypeError(f"Unsupported image source type: {type(image_source)}")

    if params is None:
        merged_params = SpotParams(
            beam_x=src_bx,
            beam_y=src_by,
            pixel_size_x=src_qx,
            pixel_size_y=src_qy,
            distance=src_dist,
            wavelength=src_wvl,
        )
    else:
        merged_params = SpotParams(
            snr_threshold=params.snr_threshold,
            min_spot_area=params.min_spot_area,
            max_spot_area=params.max_spot_area,
            beam_x=params.beam_x if params.beam_x != 0.0 else src_bx,
            beam_y=params.beam_y if params.beam_y != 0.0 else src_by,
            pixel_size_x=params.pixel_size_x if params.pixel_size_x > 0.0 else src_qx,
            pixel_size_y=params.pixel_size_y if params.pixel_size_y > 0.0 else src_qy,
            distance=params.distance if params.distance > 0.0 else src_dist,
            wavelength=params.wavelength if params.wavelength > 0.0 else src_wvl,
            d_min=params.d_min,
            d_max=params.d_max,
            ice_mask=params.ice_mask,
            ice_sensitivity=params.ice_sensitivity,
            masked_rings=params.masked_rings,
        )

    return data, merged_params, detected_index


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

    detected_ice_rings: Optional[List[IceRing]] = None
    if params.ice_mask:
        detected_ice_rings, _ = detect_ice_rings_data(data, params=params)
        if detected_ice_rings:
            active_masked = list(params.masked_rings) if params.masked_rings is not None else []
            for ring in detected_ice_rings:
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

    out_buffer = (CMxSpot * max_spots)()

    if context is not None and context.handle is not None and nx <= context.max_nx and ny <= context.max_ny:
        num_found = lib.mxspots_find_spots_ctx(
            context.handle,
            data_ptr,
            nx,
            ny,
            ctypes.byref(c_params),
            out_buffer,
            max_spots,
        )
    else:
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

    # Post-mask validation: ensure no spots survive within any detected or user-specified ice shell
    if (params.masked_rings or detected_ice_rings) and spots:
        active_shells = []
        if params.masked_rings:
            active_shells.extend(params.masked_rings)
        if detected_ice_rings:
            active_shells.extend([(r.d_min, r.d_max) for r in detected_ice_rings])

        validated_spots = []
        for s in spots:
            inside = False
            for d_min_ring, d_max_ring in active_shells:
                low_d = min(d_min_ring, d_max_ring)
                high_d = max(d_min_ring, d_max_ring)
                if low_d <= s.d_spacing <= high_d:
                    inside = True
                    break
            if not inside:
                validated_spots.append(s)
        spots = validated_spots

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
    ImageFrame and DataSet objects, and SyntheticFrame objects.
    :param params: SpotParams object, defaults to None
    :param max_spots: Maximum number of spots to return, defaults to 2000.
    :param xds_output: Optional path to export detected spots to SPOT.XDS format.
    :param index: Frame index (1-based), sets Z coordinate to index - 0.5 in SPOT.XDS.
    :param z: Explicit Z continuous frame coordinate in SPOT.XDS (overrides index).
    :param context: Optional pre-allocated SpotFinderContext for zero-allocation execution.
    """
    data, merged_params, detected_index = extract_frame_and_params(image_source, params)
    spot_list = findspots_data(data, merged_params, max_spots=max_spots, context=context)

    if xds_output is not None:
        if z is not None:
            effective_z = float(z)
        elif index is not None:
            effective_z = float(index) - 0.5
        else:
            effective_z = float(detected_index) - 0.5

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

    :param image_source: Image file path, mxio DataSet/ImageFrame, SyntheticFrame, or 2D numpy array.
    :param params: SpotParams with detector geometry and ice sensitivity.
    :return: Tuple of (detected IceRing list, overall ice_score).
    """
    if isinstance(image_source, np.ndarray):
        return detect_ice_rings_data(image_source, params=params)

    data, merged_params, _ = extract_frame_and_params(image_source, params)
    return detect_ice_rings_data(data, params=merged_params)
