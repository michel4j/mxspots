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
    cx: float
    cy: float
    spots: List[Tuple[float, float, float]]  # List of (x, y, intensity)


@dataclass(frozen=True)
class SyntheticFrame:
    wavelength: float
    distance: float
    nx: int
    ny: int
    qx: float
    qy: float
    cx: float
    cy: float
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
        cx=float(doc["cx"]),
        cy=float(doc["cy"]),
        spots=spots,
    )


def generate_2d_gaussian(size: tuple[int, int], center: tuple[float, float] | None = None, sigma: float = 3.0):
    """
    Generates a square 2D array filled with a gaussian distribution centered at center
    :param center: Center pixel coordinates
    :param size: size of the 2D array
    :param sigma: The standard deviation of the gaussian distribution
    """

    nx, ny = size
    cx, cy = (nx//2, ny//2) if center is None else center
    x = np.arange(nx) - cx
    y = np.arange(ny) - cy
    x_grid, y_grid = np.meshgrid(x, y)
    squared_distance = x_grid ** 2 + y_grid ** 2
    data = np.exp(-squared_distance / (2.0 * sigma ** 2))
    return data


def add_powder_ring(
    frame_data: np.ndarray,
    cx: float,
    cy: float,
    qx: float,
    qy: float,
    distance: float,
    wavelength: float,
    d_spacing: float,
    radial_width: float = 2.0,
    peak_intensity: float = 50.0,
) -> np.ndarray:
    """
    Inject a continuous concentric powder diffraction ring (e.g. ice ring) into frame_data.
    :param frame_data: 2D numpy array (ny, nx) to modify in-place (and return).
    :param cx: Beam center X coordinate in pixels.
    :param cy: Beam center Y coordinate in pixels.
    :param qx: Pixel size X in mm.
    :param qy: Pixel size Y in mm.
    :param distance: Detector distance in mm.
    :param wavelength: Incident wavelength in Angstroms.
    :param d_spacing: Ring resolution d-spacing in Angstroms.
    :param radial_width: Standard deviation width of the ring profile in pixels.
    :param peak_intensity: Maximum added intensity at the ring center.
    """
    theta = math.asin(wavelength / (2.0 * d_spacing))
    r_mm = distance * math.tan(2.0 * theta)

    ny, nx = frame_data.shape
    ys = np.arange(ny, dtype=np.float32)
    xs = np.arange(nx, dtype=np.float32)
    xx, yy = np.meshgrid(xs, ys)

    rx_mm = (xx - cx) * qx
    ry_mm = (yy - cy) * qy
    r_grid_mm = np.sqrt(rx_mm * rx_mm + ry_mm * ry_mm)

    sigma_mm = radial_width * (0.5 * (qx + qy))
    two_sigma_sq = 2.0 * sigma_mm * sigma_mm

    diff_r = r_grid_mm - r_mm
    ring_profile = peak_intensity * np.exp(-(diff_r * diff_r) / two_sigma_sq)
    frame_data += ring_profile.astype(np.float32)
    return frame_data


def generate_synthetic_frame(
    spec_or_path: Union[SyntheticSpec, str, Path],
    max_spots: Optional[int] = None,
    sigma: float = 1.2,
    background: float = 10.0,
    noise_sigma: float = 2.0,
    add_noise: bool = True,
) -> SyntheticFrame:
    """
    Render a 2D synthetic diffraction frame with 2D Gaussian spot profiles.
    Returns a SyntheticFrame containing metadata and the float32 2D NumPy array.
    :param spec_or_path: Either a Frame Spec or a path to a yaml file containing the frame spec
    :param max_spots: Maximum number of spots to render, if not provided, use all spots in frame spec
    :param sigma: Standard deviation of the spot profile
    :param background: Background intensity
    :param noise_sigma: Standard deviation of the gaussian noise
    :param add_noise: Whether to add a background noise and scatter to the frame
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
        # Add background noise
        noise = np.random.normal(background, noise_sigma, size=frame_data.shape).astype(np.float32)
        frame_data += noise

        # Add background scatter
        scale_1, scale_2 = 10 * background, background
        scatter = generate_2d_gaussian(size=(spec.nx, spec.ny), center=(spec.cx, spec.cy), sigma=spec.nx) * scale_1
        scatter -= generate_2d_gaussian(size=(spec.nx, spec.ny), center=(spec.cx, spec.cy), sigma=spec.nx//8) * scale_2
        frame_data += scatter

        np.clip(frame_data, 0.0, None, out=frame_data)

    return SyntheticFrame(
        wavelength=spec.wavelength,
        distance=spec.distance,
        nx=spec.nx,
        ny=spec.ny,
        qx=spec.qx,
        qy=spec.qy,
        cx=spec.cx,
        cy=spec.cy,
        spots=spots_to_render,
        data=frame_data,
    )


@functools.lru_cache(maxsize=16)
def get_cached_synthetic_frame(
    yaml_name_or_path: str,
    max_spots: Optional[int] = None,
    sigma: float = 1.5,
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
