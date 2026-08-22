import functools
import math
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Union, Optional
import numpy as np
import yaml


@dataclass(frozen=True)
class SyntheticSpec:
    wavelength: float
    distance: float
    nx: int
    ny: int
    qx: float
    qy: float
    spots: List[Tuple[float, float, float]]  # List of (x, y, intensity)


@dataclass(frozen=True)
class SyntheticFrame:
    wavelength: float
    distance: float
    nx: int
    ny: int
    qx: float
    qy: float
    spots: List[Tuple[float, float, float]]
    data: np.ndarray  # 2D float32 array of shape (ny, nx)


def load_synthetic_spec(yaml_path: Union[str, Path]) -> SyntheticSpec:
    """Load an XDS synthetic diffraction frame YAML spec."""
    path = Path(yaml_path)
    if not path.is_file():
        raise FileNotFoundError(f"Synthetic frame YAML spec not found: {path}")

    with open(path, "r", encoding="utf-8") as f:
        doc = yaml.safe_load(f)

    spots = [
        (float(s[0]), float(s[1]), float(s[2]))
        for s in doc.get("spots", [])
    ]

    return SyntheticSpec(
        wavelength=float(doc["wavelength"]),
        distance=float(doc["distance"]),
        nx=int(doc["nx"]),
        ny=int(doc["ny"]),
        qx=float(doc["qx"]),
        qy=float(doc["qy"]),
        spots=spots,
    )


def generate_synthetic_frame(
    spec_or_path: Union[SyntheticSpec, str, Path],
    max_spots: Optional[int] = None,
    sigma: float = 1.2,
    background: float = 10.0,
    noise_sigma: float = 2.0,
    add_noise: bool = False,
) -> SyntheticFrame:
    """
    Render a 2D synthetic diffraction frame with 2D Gaussian spot profiles.
    
    Returns a SyntheticFrame containing metadata and the float32 2D NumPy array.
    """
    if isinstance(spec_or_path, (str, Path)):
        spec = load_synthetic_spec(spec_or_path)
    else:
        spec = spec_or_path

    # Allocate frame array with shape (ny, nx)
    frame_data = np.zeros((spec.ny, spec.nx), dtype=np.float32)

    spots_to_render = spec.spots if max_spots is None else spec.spots[:max_spots]

    radius = int(math.ceil(3.5 * sigma))
    two_sigma_sq = 2.0 * sigma * sigma
    norm_factor = 1.0 / (2.0 * math.pi * sigma * sigma)

    for x_spot, y_spot, intensity in spots_to_render:
        x_min = max(0, int(math.floor(x_spot - radius)))
        x_max = min(spec.nx, int(math.ceil(x_spot + radius)) + 1)
        y_min = max(0, int(math.floor(y_spot - radius)))
        y_max = min(spec.ny, int(math.ceil(y_spot + radius)) + 1)

        if x_min >= x_max or y_min >= y_max:
            continue

        # Generate grid for local patch
        ys = np.arange(y_min, y_max, dtype=np.float32)
        xs = np.arange(x_min, x_max, dtype=np.float32)
        xx, yy = np.meshgrid(xs, ys)

        dx = xx - x_spot
        dy = yy - y_spot
        dist_sq = dx * dx + dy * dy

        patch = intensity * norm_factor * np.exp(-dist_sq / two_sigma_sq)
        frame_data[y_min:y_max, x_min:x_max] += patch.astype(np.float32)

    if add_noise:
        # Add background and Gaussian noise
        noise = np.random.normal(background, noise_sigma, size=frame_data.shape).astype(np.float32)
        frame_data += noise
        np.clip(frame_data, 0.0, None, out=frame_data)

    return SyntheticFrame(
        wavelength=spec.wavelength,
        distance=spec.distance,
        nx=spec.nx,
        ny=spec.ny,
        qx=spec.qx,
        qy=spec.qy,
        spots=spots_to_render,
        data=frame_data,
    )


@functools.lru_cache(maxsize=16)
def get_cached_synthetic_frame(
    yaml_name_or_path: str,
    max_spots: Optional[int] = None,
    sigma: float = 1.2,
    background: float = 10.0,
    noise_sigma: float = 2.0,
    add_noise: bool = False,
) -> SyntheticFrame:
    """Cached loader for synthetic frames to avoid re-generating large arrays across tests."""
    p = Path(yaml_name_or_path)
    if not p.is_file():
        # Check relative to tests/data
        tests_data_dir = Path(__file__).resolve().parent.parent.parent / "tests" / "data"
        p = tests_data_dir / yaml_name_or_path

    return generate_synthetic_frame(
        p,
        max_spots=max_spots,
        sigma=sigma,
        background=background,
        noise_sigma=noise_sigma,
        add_noise=add_noise,
    )
