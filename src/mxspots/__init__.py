"""
mxspots: Spot Analysis for Assessing Quality of MX Diffraction Frames
"""

from .models import (
    SpotParams,
    Spot,
    SpotList,
    ScoreResult,
    IndexResult,
)
from ._lib import get_lib
from .synthetic import (
    SyntheticSpec,
    SyntheticFrame,
    load_synthetic_spec,
    generate_synthetic_frame,
    get_cached_synthetic_frame,
)
from .spotfinder import (
    findspots,
    findspots_data,
)
from .scorer import (
    score,
    score_data,
)

__version__ = "0.1.0"

__all__ = [
    "SpotParams",
    "Spot",
    "SpotList",
    "ScoreResult",
    "IndexResult",
    "get_lib",
    "SyntheticSpec",
    "SyntheticFrame",
    "load_synthetic_spec",
    "generate_synthetic_frame",
    "get_cached_synthetic_frame",
    "findspots",
    "findspots_data",
    "score",
    "score_data",
    "__version__",
]
