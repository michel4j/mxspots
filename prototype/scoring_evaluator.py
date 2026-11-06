#!/usr/bin/env python3
"""
Throwaway prototype: Scoring model dynamic range evaluation against real MX sweeps.
Part of Wayfinder ticket michel4j/mxspots#91.
"""
import glob
import math
import os
import sys
from typing import Dict, Tuple
from mxspots.scorer import score

# Calibration parameter sets: (w0, wn, ws, wp, wr)
MODELS: Dict[str, Tuple[float, float, float, float, float]] = {
    "Current (ADR-0020)": (-6.50, 0.75, 1.20, 1.40, 0.50),
    "Calib 1 (w0=-4.40)": (-4.40, 0.52, 0.31, 0.60, 0.35),
    "Calib 2 (w0=-5.00)": (-5.00, 0.58, 0.35, 0.70, 0.40),
    "Calib 3 (w0=-5.60)": (-5.60, 0.64, 0.38, 0.80, 0.45),
}


def eval_logit_score(
    nb: int,
    snr: float,
    pb: float,
    d98: float,
    ice_pen: float,
    w0: float,
    wn: float,
    ws: float,
    wp: float,
    wr: float,
) -> float:
    if nb <= 0:
        return 0.0
    s_res = max(0.0, min(1.0, (4.0 - d98) / (4.0 - 1.2)))
    z = (
        w0
        + wn * math.log(1.0 + nb)
        + ws * math.log(1.0 + snr)
        + wp * (pb / 100.0)
        + wr * s_res
    )
    s_raw = 100.0 / (1.0 + math.exp(-z))
    return max(0.0, min(100.0, s_raw - ice_pen))


def run_sweep(dataset_pattern: str, max_frames: int = 50):
    files = sorted(glob.glob(dataset_pattern))[:max_frames]
    print(f"\nEvaluating dataset pattern: {dataset_pattern} ({len(files)} frames)")
    print(
        f"{'Frame':<12} | {'N_B':<5} | {'SNR':<7} | {'d_min':<5} | {'Current':<7} | {'Calib 1':<7} | {'Calib 2':<7} | {'Calib 3':<7}"
    )
    print("-" * 75)

    for f in files:
        res = score(f)
        bname = os.path.basename(f)
        nb = res.bragg_spots
        snr = res.avg_snr if nb > 0 else 0.0
        pb = res.bragg_percent
        d98 = res.d_min
        ice_pen = res.ice_score
        scores = {
            m: eval_logit_score(nb, snr, pb, d98, ice_pen, *params)
            for m, params in MODELS.items()
        }
        if nb > 0:
            print(
                f"{bname:<12} | {nb:<5} | {snr:<7.1f} | {d98:<5.2f} | {res.score:<7.2f} | {scores['Calib 1 (w0=-4.40)']:<7.2f} | {scores['Calib 2 (w0=-5.00)']:<7.2f} | {scores['Calib 3 (w0=-5.60)']:<7.2f}"
            )
        else:
            print(
                f"{bname:<12} | {'0':<5} | {'0.0':<7} | {'-':<5} | {0.0:<7.2f} | {0.0:<7.2f} | {0.0:<7.2f} | {0.0:<7.2f}"
            )


if __name__ == "__main__":
    run_sweep("/data/xtal/datasets/lyso1/lyso_1-0616-094123_*.cbf", max_frames=30)
