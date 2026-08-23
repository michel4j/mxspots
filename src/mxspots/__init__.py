from .models import IceRing, SpotParams, Spot, SpotList, ScoreResult
from .spotfinder import (
    findspots,
    findspots_data,
    detect_ice_rings,
    detect_ice_rings_data,
    SpotFinderContext,
)
from .scorer import score, score_data, score_spots
from .synthetic import (
    SyntheticSpec,
    SyntheticFrame,
    load_synthetic_spec,
    generate_synthetic_frame,
    get_cached_synthetic_frame,
    add_powder_ring,
)

__version__ = "0.1.0"

__all__ = [
    "IceRing",
    "SpotParams",
    "Spot",
    "SpotList",
    "ScoreResult",
    "SpotFinderContext",
    "findspots",
    "findspots_data",
    "detect_ice_rings",
    "detect_ice_rings_data",
    "score",
    "score_data",
    "score_spots",
    "SyntheticSpec",
    "SyntheticFrame",
    "load_synthetic_spec",
    "generate_synthetic_frame",
    "get_cached_synthetic_frame",
    "add_powder_ring",
]
