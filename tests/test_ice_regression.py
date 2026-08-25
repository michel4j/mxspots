import pytest
import numpy as np
from mxspots.models import SpotParams, SpotList, ScoreResult
from mxspots.spotfinder import findspots_data, detect_ice_rings_data
from mxspots.scorer import score_data
from mxspots.synthetic import add_powder_ring


def test_count_list_parity_under_heavy_ice():
    # Verify spot_list.count == len(spot_list.spots) in all scenarios
    rng = np.random.default_rng(101)
    ny, nx = 1024, 1024
    cx, cy = 512.0, 512.0
    data = rng.normal(loc=10.0, scale=2.0, size=(ny, nx)).astype(np.float32)

    # 1. Zero spots
    params = SpotParams(beam_x=cx, beam_y=cy, distance=200.0, wavelength=1.0, pixel_size_x=0.075, pixel_size_y=0.075)
    slist_empty = findspots_data(data, params=params)
    assert slist_empty.count == 0
    assert len(slist_empty.spots) == 0
    assert slist_empty.count == len(slist_empty.spots)

    # 2. Add ring and 20 spots (10 on ring, 10 outside)
    add_powder_ring(
        data,
        cx=cx,
        cy=cy,
        qx=0.075,
        qy=0.075,
        distance=200.0,
        wavelength=1.0,
        d_spacing=3.897,
        radial_width=3.0,
        peak_intensity=100.0,
    )
    # 10 spots outside ring
    for i in range(10):
        x = 100 + i * 30
        y = 100 + i * 30
        data[y-1:y+2, x-1:x+2] += 400.0

    slist_masked = findspots_data(data, params=params)
    assert slist_masked.count == len(slist_masked.spots)
    assert slist_masked.count == 10


def test_multi_ring_frames_detection_and_masking():
    # Test frame with 4 canonical ice rings simultaneously: 3.897, 3.669, 2.249, 2.071 A
    rng = np.random.default_rng(202)
    ny, nx = 2048, 2048
    cx, cy = 1024.0, 1024.0
    qx, qy = 0.075, 0.075
    distance = 180.0
    wavelength = 1.0

    data = rng.normal(loc=12.0, scale=2.0, size=(ny, nx)).astype(np.float32)
    rings_to_inject = [3.897, 3.669, 2.249, 2.071]

    for d_target in rings_to_inject:
        add_powder_ring(
            data,
            cx=cx,
            cy=cy,
            qx=qx,
            qy=qy,
            distance=distance,
            wavelength=wavelength,
            d_spacing=d_target,
            radial_width=2.5,
            peak_intensity=60.0,
        )

    # Inject 8 genuine spots outside all rings (e.g. at ~5.0 A, ~3.0 A, ~1.5 A)
    # Radii for 5.0 A: 2*asin(1/10) = 11.48 deg -> 180*tan(11.48)/0.075 = 487 px
    # Radii for 1.5 A: 2*asin(1/3) = 38.94 deg -> 180*tan(38.94)/0.075 = 1939 px
    genuine_radii = [300, 487, 800, 1200]
    for r in genuine_radii:
        for angle in [0.3, 1.8]:
            x = int(cx + r * np.cos(angle))
            y = int(cy + r * np.sin(angle))
            if 0 < x < nx - 2 and 0 < y < ny - 2:
                data[y-1:y+2, x-1:x+2] += 500.0

    params = SpotParams(
        beam_x=cx,
        beam_y=cy,
        pixel_size_x=qx,
        pixel_size_y=qy,
        distance=distance,
        wavelength=wavelength,
        ice_mask=True,
        ice_sensitivity=2.5,
    )

    spot_list = findspots_data(data, params=params)
    assert spot_list.ice_rings is not None
    assert len(spot_list.ice_rings) >= 3  # At least 3 of the 4 rings detected
    assert spot_list.count == len(spot_list.spots)

    # Zero spots should survive inside any detected ice ring
    for spot in spot_list.spots:
        for ring in spot_list.ice_rings:
            assert not (ring.d_min <= spot.d_spacing <= ring.d_max), (
                f"Leakage: spot at d={spot.d_spacing:.3f} survived inside ice ring [{ring.d_min:.3f}, {ring.d_max:.3f}]"
            )


def test_borderline_width_rings_containment():
    # Test containment across extreme ring widths: narrow (1.0 px), medium (4.0 px), broad diffuse (10.0 px)
    for width in [1.0, 3.0, 6.0]:
        rng = np.random.default_rng(int(width * 10))
        ny, nx = 1024, 1024
        cx, cy = 512.0, 512.0
        qx, qy = 0.075, 0.075
        distance = 150.0
        wavelength = 1.0
        d_target = 3.897

        data = rng.normal(loc=15.0, scale=2.0, size=(ny, nx)).astype(np.float32)
        add_powder_ring(
            data,
            cx=cx,
            cy=cy,
            qx=qx,
            qy=qy,
            distance=distance,
            wavelength=wavelength,
            d_spacing=d_target,
            radial_width=width,
            peak_intensity=100.0,
        )

        params = SpotParams(
            beam_x=cx,
            beam_y=cy,
            pixel_size_x=qx,
            pixel_size_y=qy,
            distance=distance,
            wavelength=wavelength,
            ice_mask=True,
            ice_sensitivity=2.5,
        )

        detected_rings, score = detect_ice_rings_data(data, params=params)
        assert len(detected_rings) >= 1
        ring = detected_rings[0]

        # Verify envelope: d_min and d_max must strictly bound d_target with safety margin
        assert ring.d_min < d_target
        assert ring.d_max > d_target
        assert ring.d_min <= d_target * 0.98 + 1e-4
        assert ring.d_max >= d_target * 1.02 - 1e-4


def test_residual_spot_leakage_dense_sampling():
    # Densely place 60 spots directly along the circumference of the 3.897 A ice ring
    rng = np.random.default_rng(505)
    ny, nx = 1024, 1024
    cx, cy = 512.0, 512.0
    qx, qy = 0.075, 0.075
    distance = 150.0
    wavelength = 1.0
    d_target = 3.897

    data = rng.normal(loc=10.0, scale=2.0, size=(ny, nx)).astype(np.float32)
    add_powder_ring(
        data,
        cx=cx,
        cy=cy,
        qx=qx,
        qy=qy,
        distance=distance,
        wavelength=wavelength,
        d_spacing=d_target,
        radial_width=3.0,
        peak_intensity=60.0,
    )

    # Ring radius in pixels
    theta = np.arcsin(wavelength / (2.0 * d_target))
    r_px = (distance * np.tan(2.0 * theta)) / qx

    # Inject 60 spots along the ring circumference with slight jitter
    angles = np.linspace(0, 2 * np.pi, 60, endpoint=False)
    for a in angles:
        jitter_r = r_px + rng.uniform(-2.0, 2.0)
        x = int(cx + jitter_r * np.cos(a))
        y = int(cy + jitter_r * np.sin(a))
        if 0 < x < nx - 2 and 0 < y < ny - 2:
            data[y-1:y+2, x-1:x+2] += 400.0

    # Inject 4 control spots far away
    control_spots = [(200, 200), (200, 800), (800, 200), (800, 800)]
    for y, x in control_spots:
        data[y-1:y+2, x-1:x+2] += 500.0

    params = SpotParams(
        beam_x=cx,
        beam_y=cy,
        pixel_size_x=qx,
        pixel_size_y=qy,
        distance=distance,
        wavelength=wavelength,
        ice_mask=True,
        ice_sensitivity=2.5,
    )

    spot_list = findspots_data(data, params=params)
    assert spot_list.count == len(spot_list.spots)
    assert spot_list.ice_rings is not None
    assert len(spot_list.ice_rings) >= 1

    # Ensure zero leakage from the 60 on-ring spots
    for spot in spot_list.spots:
        for ring in spot_list.ice_rings:
            assert not (ring.d_min <= spot.d_spacing <= ring.d_max), (
                f"Leakage: spot at ({spot.x}, {spot.y}), d={spot.d_spacing:.3f} inside ring [{ring.d_min:.3f}, {ring.d_max:.3f}]"
            )

    # Control spots outside the ring should survive
    assert spot_list.count == 4


def test_overlapping_borderline_rings_leakage_prevention():
    # Test overlapping adjacent ice rings (3.897 A and 3.669 A)
    rng = np.random.default_rng(606)
    ny, nx = 1024, 1024
    cx, cy = 512.0, 512.0
    qx, qy = 0.075, 0.075
    distance = 120.0  # smaller distance compresses radial separation
    wavelength = 1.0

    data = rng.normal(loc=10.0, scale=2.0, size=(ny, nx)).astype(np.float32)

    add_powder_ring(
        data,
        cx=cx,
        cy=cy,
        qx=qx,
        qy=qy,
        distance=distance,
        wavelength=wavelength,
        d_spacing=3.897,
        radial_width=4.0,
        peak_intensity=100.0,
    )
    add_powder_ring(
        data,
        cx=cx,
        cy=cy,
        qx=qx,
        qy=qy,
        distance=distance,
        wavelength=wavelength,
        d_spacing=3.669,
        radial_width=4.0,
        peak_intensity=100.0,
    )

    params = SpotParams(
        beam_x=cx,
        beam_y=cy,
        pixel_size_x=qx,
        pixel_size_y=qy,
        distance=distance,
        wavelength=wavelength,
        ice_mask=True,
        ice_sensitivity=2.5,
    )

    spot_list = findspots_data(data, params=params)
    assert spot_list.count == len(spot_list.spots)
    assert spot_list.ice_rings is not None
    assert len(spot_list.ice_rings) >= 2

    for spot in spot_list.spots:
        for ring in spot_list.ice_rings:
            assert not (ring.d_min <= spot.d_spacing <= ring.d_max)
