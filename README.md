# mxspots

**High-Performance Spot Finding, Bragg Lattice Analysis, and Quality Scoring for Macromolecular Crystallography (MX)**

`mxspots` is a fast, multithreaded library and command-line toolkit designed for real-time diffraction frame quality assessment and spot analysis in Macromolecular Crystallography. It combines an optimized C core engine with a Python API, integrating with [`mxio`](https://github.com/michel4j/mxio) to process detector images (`.cbf`, `.h5`, Eiger/Pilatus, etc.) at beamline acquisition rates.

---

## Key Features

- **Fast Multithreaded Spot Finding**: Employs 2D Integral Images (Summed-Area Tables) and OpenMP parallelization for constant-time local background and dispersion estimation, followed by single-pass streaming Connected Component Labeling (CCL).
- **Automated Ice Ring Detection & Masking**: Performs 1D azimuthal radial integration to detect characteristic powder ice rings (e.g., at 3.90 Å, 3.67 Å, 2.25 Å, 2.07 Å, 1.92 Å) and masks candidate spots falling within contaminated resolution shells.
- **Reciprocal Difference-Vector Lattice Analysis**: Uses difference-vector recurrence clustering in reciprocal space to classify regular **Bragg spots** from amorphous scatter or noise without requiring reciprocal lattice FFT indexing.
- **Hybrid Gated-Logistic Quality Scoring**: Generates a unified, normalized quality score ($0 - 100$) combining Bragg spot count, Bragg fraction, average SNR, resolution limit ($d_{98}$), and ice contamination penalties.
- **Zero-Allocation Batch Processing**: Pre-allocates reusable execution contexts (`MxSpotsContext`) for scratch buffers during batch grid scans and mesh screening.
- **XDS Compatibility**: Supports exporting spots directly to `SPOT.XDS` format.

---

## Use Cases

1. **Beamline Rastering & Crystal Screening**: Rapidly evaluate hundreds of grid-scan frames to identify optimal crystal centering and diffraction hotspots.
2. **Real-Time Data Collection Monitoring**: Compute instant quality scores and resolution limits during live rotation data collection.
3. **Automated Ice Contamination Flagging**: Detect crystalline water ice rings early and exclude corrupted resolution shells from downstream processing.
4. **Spot Finding & XDS Export**: Generate filtered spot lists and `SPOT.XDS` coordinate files for downstream data processing pipelines.

---

## Installation

### Binary Wheels (Recommended)

Pre-compiled binary wheels are available for standard Linux systems (`manylinux_2_28` and `musllinux_1_2` across `x86_64` and `aarch64` on Python 3.12+):

```bash
pip install mxspots
```

### Install from Source

When installing from source or building custom binaries:

- **Prerequisites**:
  - Python 3.12 or later
  - C Compiler (GCC, Clang, or MSVC) with OpenMP support
  - CMake 3.18 or later

```bash
# Clone the repository
git clone https://github.com/michel4j/mxspots.git
cd mxspots

# Install in editable mode with development dependencies
pip install -e ".[dev]"
```

---

## Quick Start & CLI

`mxspots` provides command-line tools for rapid diffraction frame analysis:

### 1. Finding Spots (`mxspots.findspots`)

Find Bragg spots in a single frame and optionally export them to `SPOT.XDS`:

```bash
# Basic spot finding
mxspots.findspots /path/to/frame_00001.cbf

# Custom SNR threshold, resolution filtering, and XDS export
mxspots.findspots /path/to/frame_00001.cbf --snr 5.0 --dmin 1.5 --dmax 20.0 --xds
```

#### Example Output
```text
Found 342 spots in /path/to/frame_00001.cbf:
Index  X (px)     Y (px)     d (Å)      Intensity    SNR     
------------------------------------------------------------
1      1248.50    1312.20    2.45       5420.0       42.3    
2      1105.10    1480.00    2.12       4980.0       38.7    
3      1380.25    1150.80    2.80       4320.0       35.1    
... and 322 more spots.
```

### 2. Quality Scoring (`mxspots.score`)

Compute composite quality score, Bragg spot metrics, resolution limit, and ice ring diagnostics:

```bash
# Output formatted summary
mxspots.score /path/to/frame_00001.cbf

# Output detailed JSON for automated beamline pipelines
mxspots.score /path/to/frame_00001.cbf --json
```

#### Example Output
```text
Quality Score for /path/to/frame_00001.cbf:
  Score:              88.4 / 100
  Spot Count:         342
  Bragg Spots:        318
  Bragg %:            93.0%
  Avg Intensity:      1420.5
  Lattices Detected:  1
  Average SNR:        14.22
  Resolution Limit:   1.75 Å (98th percentile)
  Ice Score:          0.00
```

---

## Python API

### Finding Spots

```python
from mxspots import findspots, SpotParams

# Configure spot finding parameters
params = SpotParams(
    snr_threshold=6.0,
    d_min=1.5,
    d_max=20.0,
    ice_mask=True,
)

# Find spots from an image file or NumPy array
spot_list = findspots("frame_00001.cbf", params=params)

print(f"Detected {spot_list.count} spots:")
for spot in spot_list.spots[:10]:
    print(f"Spot at ({spot.x:.1f}, {spot.y:.1f}), d = {spot.d_spacing:.2f} Å, I = {spot.intensity:.1f}")

# Export to SPOT.XDS
spot_list.to_xds("SPOT.XDS", frame_index=1)
```

### Scoring Diffraction Frames

```python
from mxspots import score, SpotParams

# Compute quality metrics for an image frame
result = score("frame_00001.cbf")

print(f"Composite Score:     {result.score:.1f} / 100")
print(f"Bragg Spots:         {result.bragg_spots} / {result.spot_count} ({result.bragg_percent:.1f}%)")
print(f"Resolution Limit:    {result.d_min:.2f} Å")
print(f"Lattices Detected:   {result.num_lattices}")
print(f"Ice Score:           {result.ice_score:.2f}")
```

### Ice Ring Detection

```python
from mxspots import detect_ice_rings, SpotParams

params = SpotParams(ice_sensitivity=1.0)
ice_result = detect_ice_rings("frame_00001.cbf", params=params)

if ice_result.num_rings > 0:
    print(f"Ice contamination detected (score: {ice_result.ice_score:.2f})")
    for ring in ice_result.rings:
        print(f"  Ring at {ring.d_spacing:.2f} Å (SNR: {ring.score:.1f})")
```

### High-Throughput Batch Scoring (`score_data`)

For in-memory batch screening of raw NumPy arrays:

```python
import numpy as np
from mxspots import score_data, SpotParams

# Frame 2D float32 array with detector geometry
params = SpotParams(
    beam_x=1500.0,
    beam_y=1500.0,
    distance=200.0,
    wavelength=1.0,
    pixel_size_x=0.075,
    pixel_size_y=0.075,
)

data = np.load("frame_data.npy").astype(np.float32)
res = score_data(data, params=params)
print(f"Score: {res.score:.1f}")
```

---

## Scoring Model

The Composite Quality Score ($S \in [0, 100]$) uses a **Hybrid Gated-Logistic** model:

$$\text{Score} = \begin{cases} 0.0 & \text{if } N_{\text{bragg}} = 0 \\ \text{clamp}\left(\frac{100}{1 + e^{-z}} - P_{\text{ice}}, 0, 100\right) & \text{if } N_{\text{bragg}} > 0 \end{cases}$$

where the logit $z$ is computed from:
- **Bragg Spot Count ($N_{\text{bragg}}$)**: Logarithmic scaling $\ln(1 + N_{\text{bragg}})$
- **Bragg Spot Fraction ($P_{\text{bragg}}$)**: Linear weighting of lattice conformity
- **Signal-to-Noise Ratio ($\text{SNR}$)**: Logarithmic peak quality $\ln(1 + \text{SNR})$
- **Bragg Resolution Limit ($d_{98}$)**: Linear scaling between $4.0\,\text{Å}$ and $1.2\,\text{Å}$
- **Ice Penalties ($P_{\text{ice}}$)**: Post-logistic penalty capped at at most 10 points ($P_{\text{ice}} \in [0, 10]$) based on detected ice ring count and contamination severity

---

## Testing

Run the test suite with `pytest`:

```bash
pytest
```

---

## License

This project is licensed under the MIT License.
