from dataclasses import dataclass, asdict
from typing import List, Optional
import json


@dataclass(frozen=True)
class SpotParams:
    snr_threshold: float = 3.0
    min_spot_area: int = 2
    max_spot_area: int = 500
    beam_x: float = 0.0
    beam_y: float = 0.0
    pixel_size_x: float = 0.075  # mm (e.g. standard Eiger pixel size)
    pixel_size_y: float = 0.075  # mm
    distance: float = 200.0      # mm
    wavelength: float = 1.0      # Angstroms


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

    @property
    def count(self) -> int:
        return len(self.spots)

    def to_dict(self) -> dict:
        return {
            "spot_count": self.count,
            "spots": [s.to_dict() for s in self.spots],
        }

    def to_json(self, indent: Optional[int] = None) -> str:
        return json.dumps(self.to_dict(), indent=indent)


@dataclass(frozen=True)
class ScoreResult:
    spot_count: int
    avg_snr: float
    d_min: float
    percentage_indexed: Optional[float] = None

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
