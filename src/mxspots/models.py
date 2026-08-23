from dataclasses import dataclass, asdict
from typing import List, Optional, Union, Tuple
from pathlib import Path
import json


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
    snr_threshold: float = 3.0
    min_spot_area: int = 2
    max_spot_area: int = 500
    beam_x: float = 0.0          # px (0.0 for auto/source metadata)
    beam_y: float = 0.0          # px (0.0 for auto/source metadata)
    pixel_size_x: float = 0.0    # mm (0.0 for auto/source metadata, default 0.075)
    pixel_size_y: float = 0.0    # mm (0.0 for auto/source metadata, default 0.075)
    distance: float = 0.0        # mm (0.0 for auto/source metadata, default 200.0)
    wavelength: float = 0.0      # Angstroms (0.0 for auto/source metadata, default 1.0)
    d_min: float = 0.0           # Angstroms (0.0 for unbounded)
    d_max: float = 30.0          # Angstroms
    ice_mask: bool = True        # Automatically detect and mask ice rings
    ice_sensitivity: float = 3.0 # Statistical threshold for ice ring detection
    masked_rings: Optional[List[Tuple[float, float]]] = None  # (d_min, d_max) Angstrom shells to mask


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
    avg_snr: float
    d_min: float
    percentage_indexed: Optional[float] = None
    indexed_spot_count: Optional[int] = None
    ice_score: Optional[float] = None
    ice_rings_detected: Optional[List[float]] = None
    score: float = 0.0

    def to_dict(self) -> dict:
        return asdict(self)

    def to_json(self, indent: Optional[int] = None) -> str:
        return json.dumps(self.to_dict(), indent=indent)


@dataclass(frozen=True)
class IndexResult:
    unit_cell: List[float]
    percentage_indexed: float
    indexed_spot_count: int
    total_spot_count: int

    def to_dict(self) -> dict:
        return asdict(self)

    def to_json(self, indent: Optional[int] = None) -> str:
        return json.dumps(self.to_dict(), indent=indent)
