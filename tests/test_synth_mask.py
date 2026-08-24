import pytest
from mxspots.synthetic import get_cached_synthetic_frame, add_masked_region
from mxspots.spotfinder import findspots_data, SpotParams
import numpy as np

def test_negative_intensities_are_ignored():
    """
    Ensure that a synthetic spot placed entirely inside a masked region
    (where pixel intensities are set to -1.0) is not detected by findspots.
    """
    frame = get_cached_synthetic_frame("clean.yaml", max_spots=50, add_noise=True)
    # create a copy to not corrupt cache
    data = frame.data.copy()

    # Base line: spot 0 is detected
    res_base = findspots_data(data, SpotParams())
    assert len(res_base.spots) > 0

    x_spot, y_spot, _ = frame.spots[0]
    
    # Mask a 40x40 region around the first spot
    x_min = int(max(0, x_spot - 20))
    x_max = int(min(frame.nx, x_spot + 20))
    y_min = int(max(0, y_spot - 20))
    y_max = int(min(frame.ny, y_spot + 20))

    data = add_masked_region(data, x_min, x_max, y_min, y_max, mask_val=-1.0)

    # Now detect spots on masked frame
    res_mask = findspots_data(data, SpotParams())

    # Ensure spot 0 is no longer detected!
    for s in res_mask.spots:
        cx, cy = s.x, s.y
        # It shouldn't be inside the masked bounding box
        assert not (x_min <= cx <= x_max and y_min <= cy <= y_max)
