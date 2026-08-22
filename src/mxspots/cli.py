import argparse
import sys
from pathlib import Path
from .spotfinder import findspots
from .models import SpotParams


def findspots_main():
    parser = argparse.ArgumentParser(
        prog="mxspots.findspots",
        description="Find diffraction spots in an MX image frame",
    )
    parser.add_argument("image", help="Path to diffraction image file (.cbf, .h5, .yaml, etc.)")
    parser.add_argument("--json", action="store_true", help="Output results formatted as JSON")
    parser.add_argument("--snr", type=float, default=3.0, help="SNR threshold for spot detection (default: 3.0)")
    parser.add_argument("--min-area", type=int, default=2, help="Minimum connected pixels per spot (default: 2)")
    parser.add_argument("--max-area", type=int, default=500, help="Maximum connected pixels per spot (default: 500)")
    parser.add_argument("--max-spots", type=int, default=50000, help="Maximum spots to return (default: 50000)")
    parser.add_argument("--beam-x", type=float, default=0.0, help="Detector beam center X in pixels (0 for auto)")
    parser.add_argument("--beam-y", type=float, default=0.0, help="Detector beam center Y in pixels (0 for auto)")
    parser.add_argument("--distance", type=float, default=0.0, help="Detector distance in mm (0 for auto)")
    parser.add_argument("--wavelength", type=float, default=0.0, help="Wavelength in Angstroms (0 for auto)")

    args = parser.parse_args()

    # If non-default parameters specified, pass them
    params = None
    if any([args.snr != 3.0, args.min_area != 2, args.max_area != 500, args.beam_x != 0.0,
            args.beam_y != 0.0, args.distance != 0.0, args.wavelength != 0.0]):
        params = SpotParams(
            snr_threshold=args.snr,
            min_spot_area=args.min_area,
            max_spot_area=args.max_area,
            beam_x=args.beam_x,
            beam_y=args.beam_y,
            distance=args.distance if args.distance > 0 else 200.0,
            wavelength=args.wavelength if args.wavelength > 0 else 1.0,
        )

    try:
        spot_list = findspots(args.image, params=params, max_spots=args.max_spots)
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

    if args.json:
        print(spot_list.to_json(indent=2))
    else:
        print(f"Found {spot_list.count} spots in {args.image}:")
        print(f"{'Index':<6} {'X (px)':<10} {'Y (px)':<10} {'d (Å)':<10} {'Intensity':<12} {'SNR':<8}")
        print("-" * 60)
        display_spots = spot_list.spots[:20]
        for idx, spot in enumerate(display_spots, start=1):
            print(f"{idx:<6} {spot.x:<10.2f} {spot.y:<10.2f} {spot.d_spacing:<10.2f} {spot.intensity:<12.1f} {spot.snr:<8.1f}")
        if spot_list.count > 20:
            print(f"... and {spot_list.count - 20} more spots.")


def score_main():
    print("mxspots.score - Quality scoring engine")


def index_main():
    print("mxspots.index - FFT lattice indexing engine")
