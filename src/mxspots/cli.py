import argparse
import sys
from pathlib import Path
from .spotfinder import findspots
from .scorer import score
from .models import SpotParams


def findspots_main():
    parser = argparse.ArgumentParser(
        prog="mxspots.findspots",
        description="Find diffraction spots in an MX image frame",
    )
    parser.add_argument("image", help="Path to diffraction image file (.cbf, .h5, .yaml, etc.)")
    parser.add_argument("--json", action="store_true", help="Output results formatted as JSON")
    parser.add_argument("--xds", action="store_true", help="Export spots to SPOT.XDS in current directory")
    parser.add_argument("--xds-file", type=str, default="SPOT.XDS", help="Filename for XDS export (default: SPOT.XDS)")
    parser.add_argument("--index", type=int, default=None, help="Frame index for SPOT.XDS export (calculates z = index - 0.5, default: auto or 1)")
    parser.add_argument("--dmin", type=float, default=0.0, help="High-resolution cutoff limit in Angstroms (default: 0.0, unbounded)")
    parser.add_argument("--dmax", type=float, default=30.0, help="Low-resolution cutoff limit in Angstroms (default: 30.0)")
    parser.add_argument("--snr", type=float, default=6.0, help="SNR threshold for spot detection (default: 6.0)")
    parser.add_argument("--min-area", type=int, default=2, help="Minimum connected pixels per spot (default: 2)")
    parser.add_argument("--max-area", type=int, default=500, help="Maximum connected pixels per spot (default: 500)")
    parser.add_argument("--max-spots", type=int, default=5000, help="Maximum spots to return (default: 5000)")
    parser.add_argument("--beam-x", type=float, default=0.0, help="Detector beam center X in pixels (0 for auto)")
    parser.add_argument("--beam-y", type=float, default=0.0, help="Detector beam center Y in pixels (0 for auto)")
    parser.add_argument("--distance", type=float, default=0.0, help="Detector distance in mm (0 for auto)")
    parser.add_argument("--wavelength", type=float, default=0.0, help="Wavelength in Angstroms (0 for auto)")
    parser.add_argument("--no-ice-mask", action="store_true", help="Disable automated ice ring detection and masking")
    parser.add_argument("--ice-sensitivity", type=float, default=1.0, help="Ice ring detection sensitivity threshold (default: 1.0)")

    args = parser.parse_args()

    params = SpotParams(
        snr_threshold=args.snr,
        min_spot_area=args.min_area,
        max_spot_area=args.max_area,
        beam_x=args.beam_x,
        beam_y=args.beam_y,
        distance=args.distance,
        wavelength=args.wavelength,
        d_min=args.dmin,
        d_max=args.dmax,
        ice_mask=not args.no_ice_mask,
        ice_sensitivity=args.ice_sensitivity,
    )

    xds_out = args.xds_file if args.xds else None

    try:
        spot_list = findspots(
            args.image,
            params=params,
            max_spots=args.max_spots,
            xds_output=xds_out,
            index=args.index,
        )
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

    if args.json:
        print(spot_list.to_json(indent=2))
    else:
        print(f"Found {spot_list.count} spots in {args.image}:")
        if spot_list.ice_rings:
            ring_strs = [f"{r.d_spacing:.2f} Å" for r in spot_list.ice_rings]
            print(f"  Detected Ice Rings: {len(spot_list.ice_rings)} ({', '.join(ring_strs)})")
        print(f"{'Index':<6} {'X (px)':<10} {'Y (px)':<10} {'d (Å)':<10} {'Intensity':<12} {'SNR':<8}")
        print("-" * 60)
        display_spots = spot_list.spots[:20]
        for idx, spot in enumerate(display_spots, start=1):
            print(f"{idx:<6} {spot.x:<10.2f} {spot.y:<10.2f} {spot.d_spacing:<10.2f} {spot.intensity:<12.1f} {spot.snr:<8.1f}")
        if spot_list.count > 20:
            print(f"... and {spot_list.count - 20} more spots.")
        if args.xds:
            print(f"Exported spots to {args.xds_file}")


def score_main():
    parser = argparse.ArgumentParser(
        prog="mxspots.score",
        description="Compute quality score metrics for an MX diffraction image",
    )
    parser.add_argument("image", help="Path to diffraction image file (.cbf, .h5, .yaml, etc.)")
    parser.add_argument("--json", action="store_true", help="Output results formatted as JSON")
    parser.add_argument("--dmin", type=float, default=0.0, help="High-resolution cutoff limit in Angstroms (default: 0.0, unbounded)")
    parser.add_argument("--dmax", type=float, default=30.0, help="Low-resolution cutoff limit in Angstroms (default: 30.0)")
    parser.add_argument("--snr", type=float, default=3.0, help="SNR threshold for spot detection (default: 3.0)")
    parser.add_argument("--min-area", type=int, default=2, help="Minimum connected pixels per spot (default: 2)")
    parser.add_argument("--max-area", type=int, default=500, help="Maximum connected pixels per spot (default: 500)")
    parser.add_argument("--beam-x", type=float, default=0.0, help="Detector beam center X in pixels (0 for auto)")
    parser.add_argument("--beam-y", type=float, default=0.0, help="Detector beam center Y in pixels (0 for auto)")
    parser.add_argument("--distance", type=float, default=0.0, help="Detector distance in mm (0 for auto)")
    parser.add_argument("--wavelength", type=float, default=0.0, help="Wavelength in Angstroms (0 for auto)")
    parser.add_argument("--no-ice-mask", action="store_true", help="Disable automated ice ring detection and masking")
    parser.add_argument("--ice-sensitivity", type=float, default=3.0, help="Ice ring detection sensitivity threshold (default: 3.0)")

    args = parser.parse_args()

    params = SpotParams(
        snr_threshold=args.snr,
        min_spot_area=args.min_area,
        max_spot_area=args.max_area,
        beam_x=args.beam_x,
        beam_y=args.beam_y,
        distance=args.distance,
        wavelength=args.wavelength,
        d_min=args.dmin,
        d_max=args.dmax,
        ice_mask=not args.no_ice_mask,
        ice_sensitivity=args.ice_sensitivity,
    )

    try:
        score_res = score(args.image, params=params)
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

    if args.json:
        print(score_res.to_json(indent=2))
    else:
        print(f"Quality Score for {args.image}:")
        print(f"  Score:              {score_res.score:.1f} / 100")
        print(f"  Spot Count:         {score_res.spot_count}")
        print(f"  Bragg Spots:        {score_res.bragg_spots}")
        print(f"  Bragg %:            {score_res.bragg_percent:.1f}%")
        print(f"  Avg Intensity:      {score_res.avg_intensity:.1f}")
        if score_res.num_lattices > 0:
            lattice_note = " (Warning: Multi-lattice / split crystal detected)" if score_res.num_lattices > 1 else ""
            print(f"  Lattices Detected:  {score_res.num_lattices}{lattice_note}")
        print(f"  Average SNR:        {score_res.avg_snr:.2f}")
        d_min_str = f"{score_res.d_min:.2f} Å (95th percentile)" if score_res.d_min < 900.0 else "N/A"
        print(f"  Resolution Limit:   {d_min_str}")
        if score_res.ice_score > 0.0:
            print(f"  Ice Score:          {score_res.ice_score:.2f}")
        if score_res.ice_rings_detected:
            print(f"  Ice Rings Detected: {len(score_res.ice_rings_detected)}")
