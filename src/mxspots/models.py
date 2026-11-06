from dataclasses import dataclass, asdict
from typing import List, Optional, Union, Tuple, Any
from pathlib import Path
import json
import numpy as np


@dataclass(frozen=True)
class IceRing:
    d_spacing: float           # Nominal d-spacing in Angstroms (e.g. 3.897, 3.669)
    d_min: float               # High-resolution boundary of masked annulus in Angstroms
    d_max: float               # Low-resolution boundary of masked annulus in Angstroms
    score: float = 0.0         # Detection significance / peak-to-background ratio

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class SpotParams:
    snr_threshold: float = 6.0
    min_spot_area: int = 2
    max_spot_area: int = 500
    beam_x: float = 0.0          # px (0.0 for auto/source metadata)
    beam_y: float = 0.0          # px (0.0 for auto/source metadata)
    pixel_size_x: float = 0.0    # mm (0.0 for auto/source metadata, default 0.075)
    pixel_size_y: float = 0.0    # mm (0.0 for auto/source metadata, default 0.075)
    distance: float = 0.0        # mm (0.0 for auto/source metadata, default 200.0)
    wavelength: float = 0.0      # Angstroms (0.0 for auto/source metadata, default 1.0)
    d_min: float = 0.0           # Angstroms (0.0 for unbounded)
    d_max: float = 20.0          # Angstroms
    ice_mask: bool = True        # Automatically detect and mask ice rings
    ice_sensitivity: float = 2.5 # Statistical threshold for ice ring detection
    masked_rings: Optional[List[Tuple[float, float]]] = None  # (d_min, d_max) Angstrom shells to mask


@dataclass(frozen=True)
class DiffractionFrame:
    """
    Diffraction Frame adapter encapsulating 2D pixel array, Detector Geometry,
    and frame metadata.
    """
    data: np.ndarray
    params: SpotParams
    frame_index: float = 1.0

    @classmethod
    def from_source(
        cls,
        image_source: Union[str, Path, "DiffractionFrame", Any],
        params: Optional[SpotParams] = None,
    ) -> "DiffractionFrame":
        if isinstance(image_source, DiffractionFrame):
            if params is None:
                return image_source
            merged = SpotParams(
                snr_threshold=params.snr_threshold,
                min_spot_area=params.min_spot_area,
                max_spot_area=params.max_spot_area,
                beam_x=params.beam_x if params.beam_x != 0.0 else image_source.params.beam_x,
                beam_y=params.beam_y if params.beam_y != 0.0 else image_source.params.beam_y,
                pixel_size_x=params.pixel_size_x if params.pixel_size_x > 0.0 else image_source.params.pixel_size_x,
                pixel_size_y=params.pixel_size_y if params.pixel_size_y > 0.0 else image_source.params.pixel_size_y,
                distance=params.distance if params.distance > 0.0 else image_source.params.distance,
                wavelength=params.wavelength if params.wavelength > 0.0 else image_source.params.wavelength,
                d_min=params.d_min,
                d_max=params.d_max,
                ice_mask=params.ice_mask,
                ice_sensitivity=params.ice_sensitivity,
                masked_rings=params.masked_rings,
            )
            return cls(data=image_source.data, params=merged, frame_index=image_source.frame_index)

        detected_index = 1.0

        if isinstance(image_source, (str, Path)):
            p = Path(image_source)
            if p.is_file() and p.suffix.lower() in (".yaml", ".yml"):
                from .synthetic import generate_synthetic_frame
                image_source = generate_synthetic_frame(p)
            else:
                import mxio
                image_source = mxio.DataSet.new_from_file(p)
            if image_source is None:
                raise ValueError(f"Could not load image source: {p}")

        try:
            import mxio
            is_dataset = isinstance(image_source, mxio.DataSet)
            is_imageframe = isinstance(image_source, mxio.ImageFrame)
        except ImportError:
            is_dataset = False
            is_imageframe = False

        if is_dataset:
            detected_index = float(image_source.index)
            image_source = image_source.frame

        from .synthetic import SyntheticFrame
        if isinstance(image_source, SyntheticFrame):
            src_bx = float(image_source.cx)
            src_by = float(image_source.cy)
            src_qx = float(image_source.qx)
            src_qy = float(image_source.qy)
            src_dist = float(image_source.distance)
            src_wvl = float(image_source.wavelength)
            data = image_source.data
        elif is_imageframe:
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

        if not isinstance(data, np.ndarray):
            data = np.asarray(data, dtype=np.float32)
        elif data.dtype != np.float32 or not data.flags.c_contiguous:
            data = np.ascontiguousarray(data, dtype=np.float32)

        if data.ndim != 2:
            raise ValueError(f"Expected 2D image array, got {data.ndim}D shape {data.shape}")

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

        return cls(data=data, params=merged_params, frame_index=detected_index)


@dataclass(frozen=True)
class Spot:
    x: float
    y: float
    d_spacing: float
    intensity: float
    snr: float = 0.0

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class SpotList:
    spots: List[Spot]
    ice_rings: Optional[List[IceRing]] = None

    @property
    def count(self) -> int:
        return len(self.spots)

    def to_dict(self) -> dict:
        res = {
            "spot_count": self.count,
            "spots": [s.to_dict() for s in self.spots],
        }
        if self.ice_rings is not None:
            res["ice_rings"] = [r.to_dict() for r in self.ice_rings]
        return res

    def to_json(self, indent: Optional[int] = None) -> str:
        return json.dumps(self.to_dict(), indent=indent)

    def to_xds(
        self,
        path: Union[str, Path] = "SPOT.XDS",
        z: float = 0.5,
        frame_index: Optional[int] = None,
    ) -> Path:
        """
        Export spot list to standard XDS ASCII table format (SPOT.XDS).
        Columns: X Y Z INTENSITY
        Where Z represents the continuous frame coordinate (frame_index - 0.5).
        """
        out_path = Path(path)
        effective_z = (float(frame_index) - 0.5) if frame_index is not None else float(z)
        with open(out_path, "w", encoding="utf-8") as f:
            for spot in self.spots:
                f.write(f"{spot.x:10.2f} {spot.y:10.2f} {effective_z:10.2f} {spot.intensity:10.1f}\n")
        return out_path


@dataclass(frozen=True)
class ScoreResult:
    spot_count: int
    bragg_spots: int
    bragg_percent: float
    avg_intensity: float
    avg_snr: float
    d_min: float
    ice_score: float = 0.0
    num_ice_rings: int = 0
    num_lattices: int = 0
    score: float = 0.0
    ice_rings_detected: Optional[List[IceRing]] = None

    def to_dict(self) -> dict:
        return asdict(self)

    def to_json(self, indent: Optional[int] = None) -> str:
        return json.dumps(self.to_dict(), indent=indent)
